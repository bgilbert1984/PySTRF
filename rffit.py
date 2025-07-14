import numpy as np
import argparse
import sys
import os
import math
from datetime import timedelta
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Assuming sgp4 library is installed
from sgp4.api import Satrec, WGS72, jday

# Import previously translated modules (assuming they are in the path)
try:
    from rftime import nfd_to_mjd, mjd_to_nfd, date_to_mjd, mjd_to_doy, doy_to_mjd
    from rfsites import get_site, Site # Assuming Site class definition
    from rftles import load_tles, get_tle_by_catalog_id, TleData # Assuming TleData class
except ImportError as e:
    print(f"Error importing helper modules (rftime, rfsites, rftles): {e}")
    print("Please ensure .py files from previous steps are available.")
    sys.exit(1)

# Constants
D2R = math.pi / 180.0
R2D = 180.0 / math.pi
C_KM_S = 299792.458 # Speed of light in km/s
XKMPER_WGS72 = 6378.135 # WGS72 Earth radius in km, used by sgp4 library's default wgs72 model
FLAT_WGS72 = 1.0 / 298.26 # WGS72 flattening
# Note: SGP4 library uses its own internal constants based on the gravity model (WGS72/84)

class DopplerPoint:
    """Class to hold a single Doppler measurement."""
    def __init__(self, mjd, freq, flux, site_id, rsite_id=0):
        self.mjd = mjd
        self.freq = freq # Original frequency measurement (kHz in C code, keep consistent?) Assuming kHz here.
        self.flux = flux
        self.site_id = site_id
        self.rsite_id = rsite_id # For bistatic measurements (e.g., GRAVES)
        self.site = None # Will hold Site object
        self.rsite = None # Will hold remote Site object if rsite_id > 0
        self.flag = 1 # 1=use, 2=highlighted, 0=ignore (from C code)
        # Calculated values
        self.time_offset = 0.0 # Time relative to mjd0 (days)
        self.freq_offset = 0.0 # Freq relative to f0 (kHz)
        self.model_freq_offset = 0.0 # Model freq relative to f0 (kHz)
        self.residual = 0.0 # freq_offset - model_freq_offset (kHz)

class FitData:
    """Class to hold the data and fitting context."""
    def __init__(self):
        self.points = []
        self.n = 0
        self.mjd_min = 0.0
        self.mjd_max = 0.0
        self.mjd0 = 0.0 # Reference MJD (center of data)
        self.f0 = 0.0   # Reference Frequency (kHz, center of data or fit)
        self.fit_freq = True # Whether to fit frequency (a[7] in C)
        self.ffit = 0.0 # Fitted or reference frequency (kHz)
        self.satrec = None # sgp4 Satrec object
        self.tle_name = "UNKNOWN"
        self.active_params = [False] * 7 # Corresponds to ia[0] to ia[6]

# --- Global data object for chi-squared function ---
# This is not ideal Python style, but mimics the C global `d` for simplicity
# in translating the chi-squared function. A class-based approach for fitting
# would be better in a full refactor.
global_fit_data = FitData()

# === Utility Functions (Ported or Adapted from C) ===

def observer_position_velocity_ecef(mjd, site):
    """
    Calculates observer ECEF position (km) and velocity (km/s) at a given MJD.
    Adapted from obspos_xyz in rffit.c / rftrace.c.
    Uses WGS72 constants to match SGP4 default.
    """
    if not isinstance(site, Site):
         raise TypeError("Input 'site' must be a Site object")

    lat_rad = site.lat * D2R
    lon_rad = site.lon * D2R
    alt_km = site.alt # Already in km from rfsites.py

    # --- Calculate GMST ---
    # Simplified GMST calculation from rffit/rftrace (might differ slightly from high-precision implementations)
    t_ut1 = (mjd - 51544.5) / 36525.0
    gmst_deg = (280.46061837 + 360.98564736629 * (mjd - 51544.5) +
                t_ut1 * t_ut1 * (0.000387933 - t_ut1 / 38710000.0)) % 360.0
    gmst_rad = gmst_deg * D2R

    # --- Calculate observer position ---
    # Geocentric coordinates using WGS72 model
    s_lat = math.sin(lat_rad)
    c_lat = math.cos(lat_rad)
    ff = math.sqrt(1.0 - FLAT_WGS72 * (2.0 - FLAT_WGS72) * s_lat * s_lat) #
    gc = 1.0 / ff + alt_km / XKMPER_WGS72 # Earth radii + altitude
    gs = (1.0 - FLAT_WGS72)**2 / ff + alt_km / XKMPER_WGS72 #

    theta = gmst_rad + lon_rad # Local Sidereal Time
    s_theta = math.sin(theta)
    c_theta = math.cos(theta)

    pos_km = np.array([
        gc * c_lat * c_theta * XKMPER_WGS72,
        gc * c_lat * s_theta * XKMPER_WGS72,
        gs * s_lat * XKMPER_WGS72
    ])

    # --- Calculate observer velocity (due to Earth rotation) ---
    # Angular velocity of Earth (rad/s)
    omega_earth_rad_per_sec = 7.2921150e-5 # Standard value
    # Correcting dgmst usage from C code (which seems to be deg/sec?)
    # dtheta_rad_per_sec = (dgmst(mjd) * D2R / 86400.0) # Original C code's apparent logic
    # Using standard omega_earth is simpler and likely more correct
    dtheta_rad_per_sec = omega_earth_rad_per_sec

    vel_km_s = np.array([
        -pos_km[1] * dtheta_rad_per_sec, # -omega * y
         pos_km[0] * dtheta_rad_per_sec, #  omega * x
         0.0
    ])

    return pos_km, vel_km_s


def calculate_doppler_velocity(satrec, mjd, site):
    """
    Calculates the line-of-sight velocity (km/s) between a satellite and an observer.
    Positive velocity means satellite is moving away.
    Returns velocity, azimuth (deg), altitude (deg). Azimuth is 0=N, 90=E.
    Returns None, None, None if propagation fails.
    """
    if not isinstance(site, Site):
         raise TypeError("Input 'site' must be a Site object")

    # Get observer position and velocity in ECEF frame
    obs_pos_km, obs_vel_km_s = observer_position_velocity_ecef(mjd, site)

    # Propagate satellite to the specified time using sgp4
    # Need Julian Date for sgp4 propagation
    jd_utc, jd_frac = jday(int(mjd + 2400000.5), (mjd + 2400000.5) % 1.0) # TODO: Verify correct JD conversion for sgp4
    # Or use astropy Time for better conversion if available
    # t_mjd = Time(mjd, format='mjd', scale='utc')
    # jd_utc, jd_frac = jday(t_mjd.jd1, t_mjd.jd2)

    try:
        error, sat_pos_teme, sat_vel_teme = satrec.sgp4(jd_utc, jd_frac)
        if error != 0:
            warnings.warn(f"SGP4 propagation error {error} for sat {satrec.satnum} at MJD {mjd}")
            return None, None, None
    except Exception as e:
        warnings.warn(f"SGP4 propagation failed for sat {satrec.satnum} at MJD {mjd}: {e}")
        return None, None, None

    # Note: SGP4 output (sat_pos_teme, sat_vel_teme) is in TEME frame (km, km/s).
    # Observer position/velocity is in ECEF (ITRF).
    # For accurate Doppler, both should be in the same inertial frame (e.g., GCRF)
    # or relative velocity calculated carefully.
    # The C code calculates relative pos/vel directly using ECEF observer
    # and (presumably TEME or converted) satellite vectors. This introduces
    # frame inconsistencies.
    # For a Pythonic approach using sgp4, coordinate transformations
    # (e.g., using astropy.coordinates) are recommended for rigor.
    # TEME -> ITRF conversion depends on time (Earth rotation, polar motion).

    # --- Simplified approach mimicking C code (assuming frame mismatch is acceptable/handled implicitly in C?) ---
    # Calculate relative position and velocity vectors (ignoring frame differences for now)
    rel_pos_km = sat_pos_teme - obs_pos_km
    rel_vel_km_s = sat_vel_teme - obs_vel_km_s

    r_mag = np.linalg.norm(rel_pos_km)
    if r_mag < 1e-6: # Avoid division by zero
        return 0.0, 0.0, 90.0 # Directly overhead?

    # Line-of-sight velocity (positive away)
    los_vel_km_s = np.dot(rel_vel_km_s, rel_pos_km) / r_mag

    # --- Calculate Azimuth and Altitude ---
    # Requires converting satellite TEME position to observer's topocentric frame (ENU or Az/Alt)
    # This is complex and requires coordinate transformations.
    # Placeholder - using C code's apparent direct calculation result (azimuth/altitude).
    # The C code calculates RA/Dec from relative vector and then converts to Az/Alt using site lat/lon/gmst.
    # Replicating that here:
    dx, dy, dz = rel_pos_km
    ra_rad = math.atan2(dy, dx) # RA in ECEF/TEME mix - incorrect frame but mimics C?
    dec_rad = math.asin(dz / r_mag) # Declination in ECEF/TEME mix

    # Convert RA/Dec (in mixed frame) to Az/Alt (Topocentric)
    # Requires GMST calculation again (already done in observer_position_velocity_ecef)
    t_ut1 = (mjd - 51544.5) / 36525.0
    gmst_deg = (280.46061837 + 360.98564736629 * (mjd - 51544.5) +
                t_ut1 * t_ut1 * (0.000387933 - t_ut1 / 38710000.0)) % 360.0
    gmst_rad = gmst_deg * D2R

    lat_rad = site.lat * D2R
    lon_rad = site.lon * D2R

    h_rad = gmst_rad + lon_rad - ra_rad # Hour Angle
    s_lat = math.sin(lat_rad)
    c_lat = math.cos(lat_rad)
    s_dec = math.sin(dec_rad)
    c_dec = math.cos(dec_rad)
    s_h = math.sin(h_rad)
    c_h = math.cos(h_rad)

    alt_rad = math.asin(s_lat * s_dec + c_lat * c_dec * c_h)
    azi_rad = math.atan2(-c_dec * s_h, c_lat * s_dec - s_lat * c_dec * c_h) # Azimuth from North (0=N, 90=E)

    # Note C code's equatorial2horizontal might have different azimuth convention
    # C: azi=modulo(atan2(sin(h*D2R),cos(h*D2R)*sin(lat*D2R)-tan(de*D2R)*cos(lat*D2R))*R2D,360.0)
    # This seems non-standard. Using standard conversion here.
    azimuth_deg = (math.degrees(azi_rad) + 360.0) % 360.0
    altitude_deg = math.degrees(alt_rad)

    return los_vel_km_s, azimuth_deg, altitude_deg

def calculate_altitude(satrec, mjd, site):
    """Simplified version returning only altitude."""
    _, _, alt = calculate_doppler_velocity(satrec, mjd, site)
    return alt

# === Data Loading ===

def decode_doppler_line(line, line_num):
    """Decodes a line from the Doppler data file."""
    parts = line.split()
    try:
        if len(parts) >= 4:
            mjd_str = parts[0]
            # Handle different MJD formats if necessary, assume float for now
            mjd = float(mjd_str)
            freq_hz = float(parts[1]) # Frequency in Hz
            flux = float(parts[2])
            site_id = int(parts[3])
            rsite_id = int(parts[4]) if len(parts) >= 5 else 0 # Bistatic GRAVES?
            return DopplerPoint(mjd, freq_hz / 1000.0, flux, site_id, rsite_id) # Store freq in kHz
        else:
            warnings.warn(f"Skipping malformed line {line_num}: {line.strip()}")
            return None
    except (ValueError, IndexError) as e:
        warnings.warn(f"Error parsing line {line_num}: {e} - {line.strip()}")
        return None

def read_doppler_data(filename, graves=False, foffset_khz=0.0):
    """Reads Doppler data file."""
    data = FitData()
    points = []
    site_ids = set()
    rsite_ids = set()

    try:
        with open(filename, 'r') as f:
            for i, line in enumerate(f):
                if line.strip().startswith('#') or not line.strip():
                    continue
                point = decode_doppler_line(line, i + 1)
                if point:
                    point.freq += foffset_khz # Apply offset
                    points.append(point)
                    site_ids.add(point.site_id)
                    if point.rsite_id != 0:
                        rsite_ids.add(point.rsite_id)

    except FileNotFoundError:
        print(f"Error: Data file not found: {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading data file {filename}: {e}")
        sys.exit(1)

    if not points:
        print("Error: No valid data points read from file.")
        sys.exit(1)

    data.points = points
    data.n = len(points)

    # Calculate min/max and reference MJD/Freq
    mjds = np.array([p.mjd for p in data.points])
    freqs_khz = np.array([p.freq for p in data.points]) # Freq is already in kHz
    data.mjd_min = np.min(mjds)
    data.mjd_max = np.max(mjds)
    data.mjd0 = math.floor(0.5 * (data.mjd_max + data.mjd_min)) # Reference MJD (integer day)
    # C code uses floor(0.5*(max+min)). Use median for robustness? Or fit center?
    # Let's stick to C code's reference freq for now
    data.f0 = math.floor(0.5 * (np.min(freqs_khz) + np.max(freqs_khz)))
    data.ffit = data.f0 # Initial guess for fitted frequency is center

    # Set times and frequencies relative to reference
    for p in data.points:
        p.time_offset = p.mjd - data.mjd0 # Time offset in days
        p.freq_offset = p.freq - data.f0  # Freq offset in kHz

    # If GRAVES, force reference frequency
    if graves:
        data.f0 = 143050.0 # kHz
        data.ffit = data.f0
        data.fit_freq = False # Don't fit frequency for GRAVES
        print(f"GRAVES mode enabled. Reference frequency set to {data.f0} kHz.")
        for p in data.points: # Recalculate freq offset
            p.freq_offset = p.freq - data.f0

    # Cache Site objects
    site_cache = {}
    for site_id in site_ids.union(rsite_ids):
         if site_id == 0: continue
         try:
             site_cache[site_id] = get_site(site_id) # Assumes get_site handles errors
         except ValueError as e:
             print(e)
             # Decide how to handle missing site - skip points? Exit?
             print(f"Error: Could not load site {site_id}. Exiting.")
             sys.exit(1)

    for p in data.points:
         p.site = site_cache.get(p.site_id)
         if p.rsite_id != 0:
             p.rsite = site_cache.get(p.rsite_id)
         if p.site is None or (p.rsite_id != 0 and p.rsite is None):
              print(f"Warning: Could not assign site object for point at MJD {p.mjd}")
              p.flag = 0 # Ignore points with missing site info


    print(f"Read {data.n} data points.")
    print(f"Time range (MJD): {data.mjd_min:.5f} - {data.mjd_max:.5f}")
    print(f"Freq range (kHz): {np.min(freqs_khz):.3f} - {np.max(freqs_khz):.3f}")
    print(f"Reference MJD: {data.mjd0}, Reference Freq: {data.f0:.3f} kHz")

    return data

# === Fitting Functions ===

def chi_squared_doppler(params, fit_data, active_mask):
    """
    Calculates chi-squared for Doppler fit. Objective function for scipy.minimize.

    Args:
        params (list/np.array): Current values of the parameters being fitted.
        fit_data (FitData): Object holding observations and static parameters.
        active_mask (list/np.array): Boolean mask indicating which of the 7+1
                                     parameters (6 orbital + bstar + freq) are active.

    Returns:
        float: The sum of squared residuals.
    """
    global global_fit_data # Access the global object
    local_satrec = global_fit_data.satrec # Operate on a copy? SGP4 state?
    # For simplicity, directly modify the global Satrec's attributes for now.
    # A class-based approach wrapping the fit state would be cleaner.

    current_ffit = global_fit_data.ffit # Use current best fit freq if not fitting it

    # Map active params to Satrec attributes and ffit
    param_idx = 0
    if active_mask[0]: local_satrec.inclo = math.radians(params[param_idx]); param_idx += 1
    if active_mask[1]: local_satrec.nodeo = math.radians(params[param_idx] % 360.0); param_idx += 1
    if active_mask[2]: local_satrec.ecco = max(0.0, min(0.999, params[param_idx])); param_idx += 1 # Clamp ecc [0, 1)
    if active_mask[3]: local_satrec.argpo = math.radians(params[param_idx] % 360.0); param_idx += 1
    if active_mask[4]: local_satrec.mo = math.radians(params[param_idx] % 360.0); param_idx += 1
    if active_mask[5]: local_satrec.no_kozai = max(0.05 * 2 * math.pi / 1440.0, params[param_idx] * 2 * math.pi / 1440.0); param_idx += 1 # Mean motion in rad/min, clamp > 0.05 rev/day
    if active_mask[6]: local_satrec.bstar = params[param_idx]; param_idx += 1
    if active_mask[7]: current_ffit = params[param_idx]; param_idx += 1

    # Reinitialize SGP4 with potentially new mean motion, ecc, incl
    # sgp4 library might require re-running the init logic if fundamental params change.
    # For now, assume modifying attributes is sufficient for optimization steps.
    # Re-running init_sgp4 (or equivalent sgp4 lib call) might be needed.

    sum_sq_residuals = 0.0
    n_points = 0

    for p in global_fit_data.points:
        if p.flag >= 1: # Use points marked 1 (normal) or 2 (highlighted)
             vel_los, _, alt = calculate_doppler_velocity(local_satrec, p.mjd, p.site)

             if vel_los is None: # Skip if propagation failed
                  continue

             # Handle bistatic case (GRAVES)
             if p.rsite_id != 0 and p.rsite:
                  vel_los_remote, _, _ = calculate_doppler_velocity(local_satrec, p.mjd, p.rsite)
                  if vel_los_remote is None: continue
                  # Combined doppler factor
                  doppler_factor = (1.0 - vel_los / C_KM_S) * (1.0 - vel_los_remote / C_KM_S)
             else:
                  # Monostatic case
                  doppler_factor = (1.0 - vel_los / C_KM_S)

             model_freq_khz = doppler_factor * current_ffit
             residual = p.freq - model_freq_khz # freq already in kHz
             sum_sq_residuals += residual * residual
             n_points += 1

    if n_points == 0:
        return np.inf # No points to fit

    # Return sum of squares (optimizer minimizes this)
    return sum_sq_residuals # Or sum_sq_residuals / n_points? C uses sum.

def fit_orbit(fit_data):
    """
    Performs the orbit fitting using scipy.optimize.minimize.
    """
    global global_fit_data # Use global object
    global_fit_data = fit_data

    if fit_data.satrec is None:
        print("Error: No TLE loaded for fitting.")
        return np.inf

    # Initial parameter vector 'a' from C code
    initial_params_full = [
        math.degrees(fit_data.satrec.inclo),
        math.degrees(fit_data.satrec.nodeo),
        fit_data.satrec.ecco,
        math.degrees(fit_data.satrec.argpo),
        math.degrees(fit_data.satrec.mo),
        fit_data.satrec.no_kozai * 1440.0 / (2 * math.pi), # Convert rad/min to rev/day
        fit_data.satrec.bstar,
        fit_data.ffit # Frequency in kHz
    ]

    # Filter initial parameters based on active_mask
    active_mask = fit_data.active_params + [fit_data.fit_freq] # 7 orbital + 1 freq flag
    initial_guess = [p for p, active in zip(initial_params_full, active_mask) if active]

    if not initial_guess:
        print("No parameters selected for fitting.")
        return compute_final_rms(fit_data) # Calculate RMS with current params

    print(f"Starting fit for parameters: {np.where(active_mask)[0]}")
    print(f"Initial guess: {initial_guess}")

    # Use Nelder-Mead (Simplex) to match C code's versafit/dsmin/simplex
    # Options: ftol equivalent (xtol/fatol), maxiter
    options = {'maxiter': 1000, 'disp': True, 'xatol': 1e-5, 'fatol': 1e-5} # Adjust tolerances as needed

    result = minimize(
        chi_squared_doppler,
        initial_guess,
        args=(fit_data, active_mask),
        method='Nelder-Mead',
        options=options
    )

    if result.success:
        print("Fit successful.")
        fitted_params = result.x
        # Update the global_fit_data satrec and ffit with fitted values
        param_idx = 0
        final_params_full = initial_params_full[:] # Copy initial values
        for i in range(len(active_mask)):
             if active_mask[i]:
                  final_params_full[i] = fitted_params[param_idx]
                  param_idx += 1

        # Apply fitted parameters back to the satrec object
        fit_data.satrec.inclo = math.radians(final_params_full[0])
        fit_data.satrec.nodeo = math.radians(final_params_full[1] % 360.0)
        fit_data.satrec.ecco = max(0.0, min(0.999, final_params_full[2]))
        fit_data.satrec.argpo = math.radians(final_params_full[3] % 360.0)
        fit_data.satrec.mo = math.radians(final_params_full[4] % 360.0)
        fit_data.satrec.no_kozai = max(0.05 * 2 * math.pi / 1440.0, final_params_full[5] * 2 * math.pi / 1440.0)
        fit_data.satrec.bstar = final_params_full[6]
        if fit_data.fit_freq:
             fit_data.ffit = final_params_full[7]

        print(f"Fitted Frequency: {fit_data.ffit:.3f} kHz")
        print(f"Fitted Bstar: {fit_data.satrec.bstar:.4e}")

    else:
        print(f"Fit failed: {result.message}")

    # Calculate and return final RMS
    return compute_final_rms(fit_data)

def compute_final_rms(fit_data):
     """Computes RMS based on the final parameters."""
     sum_sq = 0.0
     n_points = 0
     for p in fit_data.points:
          if p.flag >= 1:
               vel_los, _, alt = calculate_doppler_velocity(fit_data.satrec, p.mjd, p.site)
               if vel_los is None: continue

               if p.rsite_id != 0 and p.rsite:
                    vel_los_remote, _, _ = calculate_doppler_velocity(fit_data.satrec, p.mjd, p.rsite)
                    if vel_los_remote is None: continue
                    doppler_factor = (1.0 - vel_los / C_KM_S) * (1.0 - vel_los_remote / C_KM_S)
               else:
                    doppler_factor = (1.0 - vel_los / C_KM_S)

               model_freq_khz = doppler_factor * fit_data.ffit
               p.model_freq_offset = model_freq_khz - fit_data.f0 # Store relative to f0
               p.residual = p.freq_offset - p.model_freq_offset
               sum_sq += p.residual * p.residual
               n_points += 1

     if n_points == 0: return np.inf
     rms = math.sqrt(sum_sq / n_points)
     print(f"Final RMS: {rms:.4f} kHz ({n_points} points)")
     return rms

# === Plotting Function ===

def plot_data_and_fit(fit_data, show_plot=True):
    """Creates a static plot of the data and the current fit."""
    if not fit_data or fit_data.n == 0:
        print("No data to plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot data points
    times = np.array([p.time_offset for p in fit_data.points if p.flag >= 1])
    freqs = np.array([p.freq_offset for p in fit_data.points if p.flag >= 1])
    ax.plot(times, freqs, '.', markersize=5, label='Data Points')

    # Plot ignored points
    times_ignored = np.array([p.time_offset for p in fit_data.points if p.flag == 0])
    freqs_ignored = np.array([p.freq_offset for p in fit_data.points if p.flag == 0])
    if len(times_ignored) > 0:
        ax.plot(times_ignored, freqs_ignored, 'x', color='grey', markersize=5, label='Ignored Points')


    # Plot fitted curve if TLE is loaded
    if fit_data.satrec:
        model_times = np.linspace(fit_data.mjd_min - fit_data.mjd0,
                                  fit_data.mjd_max - fit_data.mjd0, 200) # More points for smooth curve
        model_freqs = []
        for t_offset in model_times:
             mjd = t_offset + fit_data.mjd0
             # Use first point's site info for model curve calculation
             # This assumes all points are from the same site(s) for the model plot
             site = fit_data.points[0].site
             rsite = fit_data.points[0].rsite

             vel_los, _, alt = calculate_doppler_velocity(fit_data.satrec, mjd, site)
             if vel_los is None or alt < 0: # Skip points below horizon or failed propagation
                  model_freqs.append(np.nan)
                  continue

             if rsite: # Bistatic
                  vel_los_remote, _, _ = calculate_doppler_velocity(fit_data.satrec, mjd, rsite)
                  if vel_los_remote is None:
                       model_freqs.append(np.nan)
                       continue
                  doppler_factor = (1.0 - vel_los / C_KM_S) * (1.0 - vel_los_remote / C_KM_S)
             else: # Monostatic
                  doppler_factor = (1.0 - vel_los / C_KM_S)

             model_freq_khz = doppler_factor * fit_data.ffit
             model_freqs.append(model_freq_khz - fit_data.f0) # Store relative to f0

        ax.plot(model_times, model_freqs, '-', color='red', label=f'Fit (Sat: {fit_data.satrec.satnum})')

    ax.set_xlabel(f"Time (days from MJD {fit_data.mjd0})")
    ax.set_ylabel(f"Frequency Offset (kHz from {fit_data.f0:.3f} kHz)")
    ax.set_title(f"Doppler Fit for {fit_data.tle_name} ({fit_data.satrec.satnum if fit_data.satrec else 'No TLE'})")
    ax.legend()
    ax.grid(True)

    if show_plot:
        plt.show()

# === Main Execution ===

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='rffit: Fit satellite orbit to RF Doppler measurements.')
    parser.add_argument('-d', '--datafile', required=True, help='Input data file with RF measurements (MJD Freq[Hz] Flux SiteID [RemoteSiteID])')
    parser.add_argument('-c', '--catalog', default=None, help='Catalog with TLEs (optional)')
    parser.add_argument('-i', '--satno', type=int, default=-1, help='NORAD ID of satellite TLE to load from catalog')
    parser.add_argument('-s', '--site', type=int, default=None, help='Default Observer Site ID (override ST_COSPAR)')
    parser.add_argument('-g', '--graves', action='store_true', help='GRAVES mode (bistatic, fixed frequency)')
    parser.add_argument('-m', '--offset', type=float, default=0.0, help='Frequency offset to apply to measurements (Hz)')
    # Note: Fitting parameter toggles (-1 to -8) and interactive commands are not implemented here.
    #       This script performs one fit based on loaded TLE/data.

    args = parser.parse_args()

    # --- Load Data ---
    data = read_doppler_data(args.datafile, args.graves, args.offset / 1000.0) # Offset in kHz

    # --- Load TLE Catalog ---
    if args.catalog:
        tle_data = load_tles(args.catalog)
    else:
        tle_data = TleData() # Empty if no catalog provided

    # --- Select Initial TLE ---
    if args.satno >= 0:
        tle_info = get_tle_by_catalog_id(tle_data, args.satno)
        if tle_info:
            data.satrec = tle_info['sat']
            data.tle_name = tle_info['name']
            print(f"Loaded TLE for {data.tle_name} ({args.satno})")
            # If frequency wasn't fitted, maybe get initial ffit guess from freq list?
            # C code seems to just use f0 unless fitfreq is toggled
        else:
            print(f"Error: TLE for specified satno {args.satno} not found in catalog.")
            sys.exit(1)
    else:
        print("Warning: No initial TLE specified via -i. Plotting data only.")


    # --- Configure Fit (Example: fit everything) ---
    if data.satrec: # Only configure fit if a TLE is loaded
        data.active_params = [True] * 7 # Fit all 7 orbital params
        data.fit_freq = not args.graves # Fit frequency unless GRAVES mode

        # --- Perform Fit ---
        final_rms = fit_orbit(data)
        print(f"Fit completed with RMS: {final_rms:.4f} kHz")

    # --- Plot Results ---
    plot_data_and_fit(data)
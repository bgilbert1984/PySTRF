import numpy as np
import math
import warnings
import os

# Assuming sgp4 library is installed
from sgp4.api import Satrec, WGS72, jday

# Import previously translated modules (assuming they are in the path)
try:
    from rftime import nfd_to_mjd, mjd_to_nfd
    from rfsites import get_site, Site
    from rftles import load_tles, get_tle_by_catalog_id, get_tle_by_index, TleData
    # Import velocity calculation logic (originally derived for rffit.py)
    # Assume it's placed in a shared utility or within rfsites/rftrace
    from rffit import observer_position_velocity_ecef, calculate_doppler_velocity # Or move this logic here/to utils
except ImportError as e:
    print(f"Error importing helper modules (rftime, rfsites, rftles, rffit): {e}")
    print("Please ensure .py files from previous steps are available.")
    sys.exit(1)

# Constants
C_KM_S = 299792.458 # Speed of light in km/s

class Trace:
    """Class to hold computed satellite trace data."""
    def __init__(self, satno, site_id, n_points, freq0=0.0, satname="", classfd=0, graves=False):
        self.satno = satno
        self.site_id = site_id
        self.n = n_points
        self.freq0 = freq0 # Reference frequency (MHz) (read from freqlist) -> Store in Hz
        self.satname = satname # Satellite name from TLE
        self.classfd = classfd # Classified flag (0 or 1)
        self.graves = graves # Is this a GRAVES (bistatic) calculation?
        # Allocate arrays
        self.mjd = np.zeros(n_points, dtype=np.float64)
        self.freq = np.zeros(n_points, dtype=np.float64) # Calculated Doppler-shifted frequency (Hz)
        self.za = np.full(n_points, 180.0, dtype=np.float32) # Zenith angle (degrees), init to below horizon

def read_frequency_list(freqlist_path):
    """Reads the frequency list file."""
    freq_map = {}
    if not os.path.exists(freqlist_path):
        warnings.warn(f"Frequency list file not found: {freqlist_path}")
        return freq_map
    try:
        with open(freqlist_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                try:
                    # Format: satno freq_mhz ... (ignore timestamp/site if present) usage suggests satno freq_mhz
                    satno = int(parts[0])
                    freq_mhz = float(parts[1])
                    freq_map[satno] = freq_mhz * 1e6 # Store frequency in Hz
                except (ValueError, IndexError):
                    warnings.warn(f"Skipping invalid line in frequency list: {line}")
    except IOError as e:
        warnings.warn(f"Error reading frequency list {freqlist_path}: {e}")
    return freq_map

def is_classified(satno, tle_dir=None):
    """
    Checks if a satellite is in the classfd.tle file.
    Mimics logic from rftrace.c.
    """
    if tle_dir is None:
        tle_dir = os.getenv("ST_TLEDIR", ".")
    classfd_path = os.path.join(tle_dir, "classfd.tle")

    if not os.path.exists(classfd_path):
        # warnings.warn(f"Classified TLE file not found: {classfd_path}")
        return False # Default to not classified if file missing

    try:
        with open(classfd_path, 'r') as f:
            for line in f:
                # Check only line 1 for satellite number
                if line.startswith('1 '):
                    try:
                        # Extract satno (columns 3-7)
                        file_satno = int(line[2:7])
                        if file_satno == satno:
                            return True
                    except (ValueError, IndexError):
                        continue # Skip malformed lines
    except IOError as e:
        warnings.warn(f"Error reading classified TLE file {classfd_path}: {e}")
    return False

def compute_trace(tle_data, mjd_times, site_id, center_freq_hz, bandwidth_hz, graves=False, freqlist_path=None, tle_dir=None):
    """
    Computes predicted Doppler traces for multiple satellites over given times.

    Args:
        tle_data (TleData): Loaded TLE data object from rftles.py.
        mjd_times (np.array): Array of MJD times for calculation.
        site_id (int): Observer's site ID.
        center_freq_hz (float): Center frequency of the observation band (Hz).
        bandwidth_hz (float): Bandwidth of the observation band (Hz).
        graves (bool): If True, perform bistatic calculation using GRAVES site (9999).
        freqlist_path (str, optional): Path to the frequency list file.
                                       If None, uses default from env ST_DATADIR.
        tle_dir (str, optional): Path to TLE directory (for classified check).
                                 If None, uses env ST_TLEDIR.

    Returns:
        list[Trace]: A list of Trace objects for satellites potentially visible
                     in the frequency band.
    """
    if freqlist_path is None:
        env_datadir = os.getenv("ST_DATADIR", ".")
        freqlist_path = os.path.join(env_datadir, "data", "frequencies.txt")

    freq_map = read_frequency_list(freqlist_path)
    if not freq_map:
        print("Warning: Frequency list is empty or not found. Cannot compute traces.")
        return []

    try:
        site = get_site(site_id)
        graves_site = get_site(9999) if graves else None
    except ValueError as e:
        print(f"Error getting site info: {e}")
        return []

    valid_mjd_times = mjd_times[mjd_times > 1.0] # Filter out invalid MJD markers (e.g., 0.0)
    n_valid_times = len(valid_mjd_times)
    if n_valid_times == 0:
        print("Warning: No valid MJD times provided.")
        return []

    # Calculate frequency band limits considering max possible Doppler shift
    # Max velocity approx +/- 8 km/s for LEO -> +/- 27 kHz per 100 MHz
    # Use a generous margin, e.g., +/- 15 km/s? C code uses 20 km/s
    max_doppler_shift_abs = (20.0 / C_KM_S) * center_freq_hz * (2 if graves else 1)
    fmin_obs = center_freq_hz - 0.5 * bandwidth_hz
    fmax_obs = center_freq_hz + 0.5 * bandwidth_hz
    fmin_search = fmin_obs - max_doppler_shift_abs
    fmax_search = fmax_obs + max_doppler_shift_abs

    print(f"Searching for satellites with frequencies between {fmin_search/1e6:.3f} and {fmax_search/1e6:.3f} MHz")

    traces = []
    # Iterate through satellites in the frequency map
    for satno, freq0_hz in freq_map.items():
        # Check if satellite frequency is potentially within range
        if not (fmin_search <= freq0_hz <= fmax_search):
             if not (graves and abs(freq0_hz - 143.050e6) < 1e3): # Special case for GRAVES freq
                continue

        # Get TLE info for this satellite
        tle_info = get_tle_by_catalog_id(tle_data, satno)
        if not tle_info:
            # warnings.warn(f"No TLE found for satellite {satno} from frequency list.")
            continue

        satrec = tle_info['sat']
        satname = tle_info['name']
        classified = is_classified(satno, tle_dir)

        # Create trace object
        trace = Trace(satno, site_id, n_valid_times, freq0=freq0_hz,
                      satname=satname, classfd=classified, graves=graves)
        trace.mjd = valid_mjd_times.copy() # Store the valid MJD times

        # Calculate trace points
        visible_points = 0
        for i, mjd in enumerate(valid_mjd_times):
            # Calculate observer-satellite Doppler
            vel_los, _, alt = calculate_doppler_velocity(satrec, mjd, site)

            if vel_los is None: # Propagation failed
                trace.freq[i] = 0.0 # Or np.nan?
                trace.za[i] = 180.0
                continue

            # Calculate Zenith Angle
            trace.za[i] = 90.0 - alt # Zenith Angle = 90 - Altitude

            # Check visibility for GRAVES if needed
            if graves:
                 _, _, alt_g = calculate_doppler_velocity(satrec, mjd, graves_site)
                 if alt_g is None: # Propagation failed for GRAVES site
                      trace.freq[i] = 0.0
                      trace.za[i] = 180.0 # Mark as invalid/below horizon
                      continue
                 # GRAVES visibility criteria from rfpng.c - azimuth check skipped here for simplicity
                 # C code uses complex azi check: !((azi<90.0 || azi>270.0) && alt>15.0 && alt<40.0)
                 # Simplified: check if satellite is above horizon for GRAVES site
                 if alt_g < 0: # If below horizon for GRAVES, treat as invisible
                      trace.za[i] = 180.0 # Use ZA as indicator
                      # This differs slightly from C's azi/alt check, might need refinement

            # Only calculate frequency if potentially visible (ZA <= 90)
            # C code calculates freq anyway, uses ZA > 90 later for plotting.
            # Let's calculate always but store ZA.
            # if trace.za[i] <= 90.0: # Check primary site visibility
            if graves:
                 vel_los_g, _, _ = calculate_doppler_velocity(satrec, mjd, graves_site)
                 if vel_los_g is None: # Should have been caught earlier
                      trace.freq[i] = 0.0
                      trace.za[i] = 180.0
                      continue
                 doppler_factor = (1.0 - vel_los / C_KM_S) * (1.0 - vel_los_g / C_KM_S)
            else:
                 doppler_factor = (1.0 - vel_los / C_KM_S)

            trace.freq[i] = doppler_factor * freq0_hz # Doppler shifted freq in Hz

            # Count points above horizon
            if trace.za[i] <= 90.0:
                visible_points += 1
            # else: # Below horizon
                # trace.freq[i] = 0.0 # Optionally zero out freq if below horizon

        # Only add trace if it has visible points
        if visible_points > 0:
            traces.append(trace)
            # print(f"Computed trace for {satname} ({satno}), {visible_points} points above horizon.")


    print(f"Found {len(traces)} potential satellite traces in the band.")
    return traces


def compute_doppler(tle_data, mjd_times, site_id, satno, graves=False, skip_high_orbits=False, output_filename="out.dat"):
    """
    Computes and saves detailed Doppler info for a specific satellite.
    Mimics rfdop.c.
    """
    try:
        site = get_site(site_id)
        graves_site = get_site(9999) if graves else None
    except ValueError as e:
        print(f"Error getting site info: {e}")
        return

    tle_info = get_tle_by_catalog_id(tle_data, satno)
    if not tle_info:
        print(f"Error: TLE not found for satellite {satno}")
        return

    satrec = tle_info['sat']

    # Skip high orbits if requested (like rfdop -H)
    # SGP4 uses no_kozai (rad/min). Convert to rev/day.
    revs_per_day = satrec.no_kozai * (1440.0 / (2 * math.pi))
    if skip_high_orbits and revs_per_day < 10.0:
        print(f"Skipping high orbit satellite {satno} ({revs_per_day:.2f} revs/day)")
        return

    valid_mjd_times = mjd_times[mjd_times > 1.0]
    n_valid_times = len(valid_mjd_times)
    if n_valid_times == 0:
        print("Warning: No valid MJD times provided.")
        return

    print(f"Calculating Doppler data for satellite {satno} to {output_filename}...")
    try:
        with open(output_filename, "w") as f:
            # Write header
            if graves:
                f.write("# satno mjd r_km v_km_s azi_deg alt_deg rg_km vg_km_s azig_deg altg_deg\n")
            else:
                f.write("# satno mjd r_km v_km_s azi_deg alt_deg\n")

            for mjd in valid_mjd_times:
                # Observer site calculations
                vel_los, azi, alt = calculate_doppler_velocity(satrec, mjd, site)
                if vel_los is None: continue # Skip failed points

                # Calculate range (magnitude of relative position vector)
                # Need pos vectors from inside calculate_doppler_velocity - recalculate or modify func
                obs_pos, _ = observer_position_velocity_ecef(mjd, site)
                jd_utc, jd_frac = jday(int(mjd + 2400000.5), (mjd + 2400000.5) % 1.0)
                error, sat_pos, _ = satrec.sgp4(jd_utc, jd_frac)
                if error != 0: continue
                r_km = np.linalg.norm(sat_pos - obs_pos)

                if graves:
                    # GRAVES site calculations
                    vel_los_g, azi_g, alt_g = calculate_doppler_velocity(satrec, mjd, graves_site)
                    if vel_los_g is None: continue
                    # Calculate range for GRAVES site
                    graves_pos, _ = observer_position_velocity_ecef(mjd, graves_site)
                    rg_km = np.linalg.norm(sat_pos - graves_pos)

                    f.write(f"{satno} {mjd:14.8f} {r_km:.3f} {vel_los:.4f} {azi:.3f} {alt:.3f} "
                            f"{rg_km:.3f} {vel_los_g:.4f} {azi_g:.3f} {alt_g:.3f}\n")
                else:
                    f.write(f"{satno} {mjd:14.8f} {r_km:.3f} {vel_los:.4f} {azi:.3f} {alt:.3f}\n")

        print("Doppler data saved.")

    except IOError as e:
        print(f"Error writing output file {output_filename}: {e}")
    except Exception as e:
         print(f"An unexpected error occurred during Doppler calculation: {e}")

# Note: identify_trace functions are more complex as they involve comparing
# an observed trace (e.g., from fitting points in rfplot/rffit) against
# these computed traces. This requires defining how the observed trace is
# represented and implementing the comparison logic (interpolation + RMS).
# This is omitted here but would be the next step for that functionality.
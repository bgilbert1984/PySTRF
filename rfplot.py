import numpy as np
import argparse
import sys
import os
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.widgets import Cursor # For potential future interactivity
import warnings
import csv
from datetime import datetime

# Import previously translated modules
try:
    from rfio import read_spectrogram, Spectrogram # Assuming class from rfio.py
    from rftime import mjd_to_nfd
    from rftrace import compute_trace, Trace # Assuming Trace class from rftrace.py
    from rftles import load_tles
    from zscale import zscale_samples
    # Fitting/Identify functions might be here or in rffit/rftrace
    # from rffit import fit_orbit # Example if fitting integrated
    # from rftrace import identify_trace # Example if identification integrated
except ImportError as e:
    print(f"Error importing helper modules: {e}")
    print("Please ensure .py files from previous steps are available.")
    sys.exit(1)

# Viridis colormap data (from rfplot.c)
# In Python, we can just use matplotlib.colormaps['viridis']
# viridis_cmap = plt.cm.get_cmap('viridis') # Or plt.colormaps['viridis'] in newer matplotlib

# === Helper Functions ===

def time_axis_formatter(mjd_values, ref_mjd0):
    """
    Creates a formatter function for the time axis (x-axis) to display HH:MM:SS.
    """
    def format_func(x_tick_val, pos):
        # x_tick_val is the subintegration index (float)
        # Need to map index back to MJD
        # This requires interpolation or lookup if x_tick_val isn't integer
        # Simple approach: use integer index, handle edges
        idx = int(round(x_tick_val))
        if 0 <= idx < len(mjd_values):
            mjd = mjd_values[idx]
            if mjd > 1.0: # Check if MJD is valid
                 # Convert MJD fractional day to seconds
                 day_fraction = (mjd - math.floor(mjd))
                 total_seconds = day_fraction * 86400.0
                 hours = int(total_seconds // 3600)
                 minutes = int((total_seconds % 3600) // 60)
                 seconds = int(total_seconds % 60)
                 return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
            else:
                 return "" # No valid time for this index
        else:
            return "" # Tick outside data range
    return format_func

def freq_axis_formatter(center_freq_mhz):
    """
    Creates a formatter function for the frequency axis (y-axis)
    to display offset in kHz relative to a center frequency.
    """
    def format_func(y_tick_val_mhz, pos):
        # y_tick_val_mhz is the absolute frequency in MHz
        offset_khz = (y_tick_val_mhz - center_freq_mhz) * 1000.0
        return f"{offset_khz:.1f}" # Format to 1 decimal place kHz
    return format_func

def plot_traces(ax, traces, center_freq_hz, show_names=False):
    """Plots satellite traces on the given axes."""
    if not traces:
        return
    center_freq_mhz = center_freq_hz / 1e6
    for trace in traces:
        # Convert calculated frequency (Hz) to MHz for plotting against freq axis
        freq_mhz = trace.freq / 1e6
        # Create array of subintegration indices (0, 1, 2, ...)
        time_indices = np.arange(trace.n)

        # Mask points below horizon (ZA > 90) or where frequency is invalid
        mask = (trace.za <= 90.0) & (freq_mhz != 0.0) # Check ZA & non-zero freq

        if not np.any(mask): # Skip if no points are visible
             continue

        # Use different colors/styles based on classified status
        color = 'cyan' if trace.classfd else 'lime' # Example colors
        linestyle = '-'
        linewidth = 1.0

        # Plot segments where trace is above horizon
        # Find contiguous segments of visible points
        masked_time = np.where(mask, time_indices, np.nan)
        masked_freq = np.where(mask, freq_mhz, np.nan)

        ax.plot(masked_time, masked_freq, linestyle=linestyle, color=color, linewidth=linewidth, label=f"{trace.satno}" if not show_names else None)

        # Add satellite label (number or name) near the start of visible segments
        # Find start points of visible segments
        starts = np.where(np.diff(np.concatenate(([False], mask, [False]))))[0]
        visible_starts = starts[::2] # Indices where mask goes False -> True

        label_text = f"{trace.satno}"
        if show_names and trace.satname:
             label_text = f"{trace.satno} - {trace.satname}"


        for start_idx in visible_starts:
             if 0 <= start_idx < trace.n:
                  # Add text slightly offset from the point
                  ax.text(time_indices[start_idx], freq_mhz[start_idx], f" {label_text}",
                          color=color, fontsize=7, ha='left', va='center', clip_on=True)


def filter_data(spectrogram, site_id, sigma_thresh, graves=False, output_filename="filter.dat"):
    """
    Filters spectrogram data based on sigma clipping and writes peaks.
    Based on filter() in rfplot.c / rffind.c.
    """
    if spectrogram.nsub == 0 or spectrogram.nchan == 0:
        print("Warning: Cannot filter empty spectrogram.")
        return

    print(f"Filtering data with sigma > {sigma_thresh}, saving to {output_filename}")
    nsub = spectrogram.nsub
    nchan = spectrogram.nchan
    data = spectrogram.z # Assuming shape (nchan, nsub)

    peaks_found = 0
    try:
        with open(output_filename, "w") as f:
            # Write header based on mode
            if graves:
                 f.write("# MJD Frequency_Hz Sigma SiteID RemoteSiteID\n")
            else:
                 f.write("# MJD Frequency_Hz Sigma SiteID\n")

            for i in range(nsub): # Loop through time (subintegrations)
                subint_data = data[:, i] # Get data for this subintegration
                mask = np.ones(nchan, dtype=bool) # Start with all pixels good

                # --- Iterative Sigma Clipping ---
                # matching loop structure
                for k in range(10): # Max 10 iterations like C code
                    valid_data = subint_data[mask]
                    if len(valid_data) < 2: break # Need >= 2 points for std dev

                    avg = np.nanmean(valid_data)
                    std = np.nanstd(valid_data)

                    if std < 1e-9: break # Avoid division by zero or meaningless std dev

                    # Update mask: remove points further than threshold from mean
                    new_mask = np.abs(subint_data - avg) <= (sigma_thresh * std)
                    # Only consider previously good points for potentially becoming bad
                    mask = mask & new_mask
                    if np.all(~mask == ~new_mask[mask]): # Check if mask changed
                        break # Converged

                # --- Find Peaks Above Threshold ---
                final_valid_data = subint_data[mask]
                if len(final_valid_data) < 2: continue # Need baseline
                final_avg = np.nanmean(final_valid_data)
                final_std = np.nanstd(final_valid_data)
                if final_std < 1e-9: continue

                # Identify potential peaks
                is_peak = (subint_data - final_avg) > (sigma_thresh * final_std)
                peak_indices = np.where(is_peak)[0]

                if len(peak_indices) == 0: continue

                # --- Local Maxima Filter ---
                # logic to keep only local maxima among adjacent peaks
                local_max_mask = np.zeros(nchan, dtype=bool)
                for idx in peak_indices:
                    is_local_max = True
                    # Check left neighbor (if it's also a peak)
                    if idx > 0 and is_peak[idx - 1] and subint_data[idx] < subint_data[idx - 1]:
                        is_local_max = False
                    # Check right neighbor (if it's also a peak)
                    if idx < nchan - 1 and is_peak[idx + 1] and subint_data[idx] < subint_data[idx + 1]:
                        is_local_max = False
                    if is_local_max:
                        local_max_mask[idx] = True

                # --- Write detected peaks ---
                mjd = spectrogram.mjd[i]
                if mjd <= 1.0: continue # Skip invalid MJD

                for j in np.where(local_max_mask)[0]:
                    # Calculate frequency in Hz
                    f_hz = (spectrogram.freq - 0.5 * spectrogram.samp_rate +
                            (j + 0.5) * (spectrogram.samp_rate / nchan))
                    sigma_val = (subint_data[j] - final_avg) / final_std

                    if graves:
                        f.write(f"{mjd:.8f} {f_hz:.3f} {sigma_val:.3f} {site_id} 9999\n")
                    else:
                        f.write(f"{mjd:.8f} {f_hz:.3f} {sigma_val:.3f} {site_id}\n")
                    peaks_found += 1

    except IOError as e:
        print(f"Error writing filter output file {output_filename}: {e}")

    print(f"Found {peaks_found} peaks.")


# === Main Script Logic ===

def main():
    parser = argparse.ArgumentParser(description='rfplot: Plot RF observations interactively (basic).')
    # Arguments similar to rfplot.c
    parser.add_argument('-p', '--path', required=True, help='Input path prefix for /path/prefix_xxxxxx.bin files')
    parser.add_argument('-s', '--start', type=int, default=0, help='Number of starting .bin file index [0]')
    parser.add_argument('-l', '--length', type=int, default=3600, help='Number of subintegrations to read [3600]')
    parser.add_argument('-b', '--nbin', type=int, default=1, help='Number of subintegrations to bin [1]')
    parser.add_argument('-z', '--zmax', type=float, default=None, help='Image scaling upper limit override (contrast factor)')
    parser.add_argument('-f', '--freq', type=float, default=0.0, help='Frequency to zoom into (Hz)')
    parser.add_argument('-w', '--bw', type=float, default=0.0, help='Bandwidth to zoom into (Hz)')
    parser.add_argument('-o', '--offset', type=float, default=0.0, help='Frequency offset to apply to header frequency (Hz)')
    # parser.add_argument('-W', '--width', type=float, default=100.0, help='Track selection width (pixels) [100]') # Interactive
    # parser.add_argument('-S', '--sigma', type=float, default=5.0, help='Track selection/filtering significance [5]') # Interactive/Filter
    parser.add_argument('-C', '--site', type=int, default=None, help='Observer Site ID (default from ST_COSPAR)')
    parser.add_argument('-c', '--catalog', default=None, help='TLE catalog file (default from ST_TLEDIR)')
    parser.add_argument('-F', '--freqlist', default=None, help='List with satellite frequencies (default from ST_DATADIR)')
    parser.add_argument('-g', '--graves', action='store_true', help='GRAVES mode (bistatic)')
    parser.add_argument('-m', '--cmap', type=str, default='viridis', help='Colormap name (e.g., viridis, gray, hot, cool) [viridis]')
    parser.add_argument('--no-overlay', action='store_true', help='Disable satellite trace overlay')
    parser.add_argument('--show-names', action='store_true', help='Show satellite names on overlay')
    parser.add_argument('--filter', action='store_true', help='Run peak filter and save to filter.dat')
    parser.add_argument('--sigma', type=float, default=5.0, help='Sigma threshold for filtering [5.0]')
    # Add ISS integration options
    parser.add_argument('--iss', action='store_true', help='Enable ISS data integration for enhanced RF modeling')
    parser.add_argument('--naval-rf', action='store_true', help='Calculate optimal naval RF positioning based on ISS data')
    parser.add_argument('--naval-coords', type=str, default=None, 
                        help='Naval fleet coordinates in format "lat1,lon1;lat2,lon2;..." (decimal degrees)')
    parser.add_argument('--target-coords', type=str, default=None,
                        help='Target coordinates in format "lat,lon" (decimal degrees)')


    args = parser.parse_args()

    # --- Determine Site ID ---
    if args.site is None:
        env_site = os.getenv("ST_COSPAR")
        if env_site:
            try:
                args.site = int(env_site)
                print(f"Using site ID from ST_COSPAR: {args.site}")
            except ValueError:
                print("Error: ST_COSPAR environment variable is not a valid integer.")
                sys.exit(1)
        else:
            print("Error: Site ID must be provided via -C or ST_COSPAR environment variable.")
            sys.exit(1)

    # --- Load Spectrogram ---
    print(f"Reading spectrogram from prefix: {args.path}")
    spectrogram = read_spectrogram(
        prefix=args.path,
        isub=args.start,
        nsub_to_read=args.length,
        f0=args.freq,
        df0=args.bw,
        nbin=args.nbin,
        foff=args.offset
    )

    if spectrogram.nsub == 0 or spectrogram.nchan == 0:
        print("Failed to read spectrogram data.")
        sys.exit(1)
        
    # --- Apply ISS Data Integration if enabled ---
    if args.iss:
        try:
            # Import the ISS integration module
            sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
            from iss_data_integration import ISSDataClient
            
            print("Integrating ISS data for enhanced RF modeling...")
            iss_client = ISSDataClient()
            
            # Get current ISS position
            iss_position = iss_client.get_current_position()
            if iss_position:
                print(f"ISS Position: Lat {iss_position['latitude']:.4f}, Lon {iss_position['longitude']:.4f}")
                
                # Enhance spectrogram with ISS-derived ionospheric data
                spectrogram = iss_client.integrate_with_rf_scythe(spectrogram, args.site)
                
                # Calculate naval RF positioning if requested
                if args.naval_rf and args.naval_coords and args.target_coords:
                    try:
                        # Parse naval coordinates
                        naval_positions = []
                        for coords in args.naval_coords.split(';'):
                            lat, lon = map(float, coords.split(','))
                            naval_positions.append({"latitude": lat, "longitude": lon})
                            
                        # Parse target coordinates
                        target_lat, target_lon = map(float, args.target_coords.split(','))
                        target_position = {"latitude": target_lat, "longitude": target_lon}
                        
                        # Calculate optimal positions
                        print("\nCalculating optimal naval RF positioning based on ISS ionospheric data...")
                        rf_factors = iss_client.calculate_naval_rf_factors(naval_positions, target_position)
                        
                        print("\nNaval Fleet RF Communication Quality:")
                        for i, vessel in enumerate(rf_factors["fleet_rf_factors"]):
                            print(f"Vessel {i+1}: Overall Quality: {vessel['overall_quality']:.2f}")
                            print(f"  HF: {vessel['rf_quality']['hf']:.2f}, VHF: {vessel['rf_quality']['vhf']:.2f}, " +
                                 f"UHF: {vessel['rf_quality']['uhf']:.2f}, SATCOM: {vessel['rf_quality']['satcom']:.2f}")
                            
                        # Get optimized positions
                        optimized = iss_client.optimize_fleet_rf_positioning(naval_positions, target_position)
                        
                        print("\nOptimized Fleet Positions for RF Communications:")
                        for i, vessel in enumerate(optimized["optimized_fleet"]):
                            improvement = vessel["quality_improvement"] * 100
                            print(f"Vessel {i+1}: Quality Improvement: {improvement:.1f}%")
                            print(f"  Original: {vessel['original']['latitude']:.4f}, {vessel['original']['longitude']:.4f}")
                            print(f"  Optimized: {vessel['optimized']['latitude']:.4f}, {vessel['optimized']['longitude']:.4f}")
                            
                    except Exception as e:
                        print(f"Error calculating naval RF positioning: {e}")
            else:
                print("Warning: Could not fetch ISS position data.")
                
        except Exception as e:
            print(f"Error integrating ISS data: {e}")
            
    # --- Load TLEs and Compute Traces ---
    traces = []
    if not args.no_overlay:
        tle_data = load_tles(args.catalog) # Handles env vars internally
        if tle_data.number_of_elements > 0:
            traces = compute_trace(
                tle_data=tle_data,
                mjd_times=spectrogram.mjd,
                site_id=args.site,
                center_freq_hz=spectrogram.freq, # Use the (potentially zoomed) center freq
                bandwidth_hz=spectrogram.samp_rate, # Use the (potentially zoomed) bw
                graves=args.graves,
                freqlist_path=args.freqlist # Handles env vars internally
            )
        else:
            print("Warning: No TLEs loaded, cannot compute traces.")

    # --- Run Filter if requested ---
    if args.filter:
        filter_data(spectrogram, args.site, args.sigma, args.graves)

    # --- Plotting ---
    print("Generating plot...")
    fig, ax = plt.subplots(figsize=(12, 7))

    # Determine plot limits and aspect ratio
    time_extent = spectrogram.nsub # Number of subintegrations on x-axis
    # Frequency axis in MHz relative to center
    freq_center_mhz = spectrogram.freq / 1e6
    freq_bw_mhz = spectrogram.samp_rate / 1e6
    freq_min_mhz = freq_center_mhz - 0.5 * freq_bw_mhz
    freq_max_mhz = freq_center_mhz + 0.5 * freq_bw_mhz
    freq_extent = [freq_min_mhz, freq_max_mhz]

    # Image extent: [left, right, bottom, top] for imshow
    extent = [0, time_extent, freq_min_mhz, freq_max_mhz]

    # Use zscale results for color limits
    vmin, vmax = spectrogram.zmin, spectrogram.zmax
    if args.zmax is not None: # Allow manual override factor for zmax
        print(f"Applying manual zmax factor: {args.zmax}")
        vmax = spectrogram.zavg.mean() + args.zmax * spectrogram.zstd.mean() # Simple override
        # Or vmax = vmin + (spectrogram.zmax - spectrogram.zmin) * args.zmax # Alternative scaling


    # Display the image (transpose Z if needed: rfio reads as (nchan, nsub))
    # imshow expects (row, col) which translates to (freq, time)
    im = ax.imshow(
        spectrogram.z,
        aspect='auto', # Adjust aspect ratio automatically
        origin='lower', # Place 0 frequency at the bottom
        extent=extent,
        cmap=args.cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation='nearest' # Avoid smoothing pixels
    )
    
    # Add ISS position overlay if enabled
    if args.iss and 'iss_position' in locals():
        # Get ionosphere data if available
        ionosphere_data = None
        if hasattr(spectrogram, 'metadata') and 'ionosphere' in spectrogram.metadata:
            ionosphere_data = spectrogram.metadata['ionosphere']
            
        # Add ISS position to plot
        plot_iss_overlay(ax, iss_position, ionosphere_data, show_label=True)
        
        # Add a note about ISS-enhanced data if integrated
        if hasattr(spectrogram, 'metadata') and 'iss_position' in spectrogram.metadata:
            plt.figtext(0.5, 0.01, "Enhanced with ISS ionospheric data", 
                       ha="center", fontsize=9, 
                       bbox={"facecolor":"orange", "alpha":0.2, "pad":5})

    # --- Configure Axes ---
    ax.set_xlabel(f"Time starting {spectrogram.nfd0} UTC")
    ax.set_ylabel(f"Frequency Offset (kHz) from {freq_center_mhz:.6f} MHz")

    # Custom time formatter for x-axis
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(time_axis_formatter(spectrogram.mjd, spectrogram.nfd0))) # Pass MJD array and start NFD

    # Custom frequency formatter for y-axis (showing offset in kHz)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(freq_axis_formatter(freq_center_mhz)))

    # Set initial view limits (optional, user can zoom/pan)
    ax.set_xlim(0, time_extent)
    ax.set_ylim(freq_min_mhz, freq_max_mhz)

    # --- Plot Traces ---
    if not args.no_overlay:
        plot_traces(ax, traces, spectrogram.freq, args.show_names) # Pass center freq in Hz

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Intensity (Arb. Units / Scaled)')

    ax.set_title(f"Spectrogram starting {spectrogram.nfd0}")
    plt.tight_layout()
    print("Showing plot...")
    plt.show()


def plot_iss_overlay(ax, iss_position, ionosphere_data=None, show_label=True):
    """
    Add ISS position and ionospheric data overlay to the plot
    
    Args:
        ax: Matplotlib axis object
        iss_position: Dictionary with ISS position data
        ionosphere_data: Optional dictionary with ionospheric data
        show_label: Whether to show ISS label
    """
    if not iss_position:
        return
        
    # Extract time and frequency extents from axis
    t_min, t_max = ax.get_xlim()
    f_min, f_max = ax.get_ylim()
    
    # ISS time position (simplified - actual implementation would map MJD properly)
    iss_time = (t_min + t_max) / 2  # Place in center for visibility
    
    # Mark ISS with a star symbol
    ax.plot(iss_time, (f_min + f_max) / 2, 'r*', markersize=12, 
            markeredgecolor='white', markeredgewidth=1, alpha=0.8)
    
    if show_label:
        # Add ISS label with position
        label_text = f"ISS: {iss_position['latitude']:.1f}°, {iss_position['longitude']:.1f}°"
        if ionosphere_data:
            f0F2 = ionosphere_data["critical_frequencies"]["f0F2"]
            label_text += f"\nf0F2: {f0F2:.1f} MHz"
            
        ax.annotate(label_text, 
                   xy=(iss_time, (f_min + f_max) / 2),
                   xytext=(iss_time + (t_max - t_min) * 0.05, (f_min + f_max) / 2 + (f_max - f_min) * 0.1),
                   arrowprops=dict(arrowstyle='->'),
                   bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                   fontsize=10)


if __name__ == "__main__":
    main()
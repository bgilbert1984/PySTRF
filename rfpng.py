import numpy as np
import argparse
import sys
import os
import math
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend for saving files
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import warnings

# Import previously translated modules
try:
    from rfio import read_spectrogram, Spectrogram
    from rftime import mjd_to_nfd
    from rftrace import compute_trace, Trace
    from rftles import load_tles
    from zscale import zscale_samples
    # Import plotting helpers developed for rfplot.py
    # Assume these are refactored into a shared plotting utility module
    # or copy them here if needed. For now, assume they exist.
    from rfplot import time_axis_formatter, freq_axis_formatter, plot_traces
except ImportError as e:
    print(f"Error importing helper modules: {e}")
    print("Please ensure rfio.py, rftime.py, rftrace.py, rftles.py, zscale.py, rfplot.py are available.")
    sys.exit(1)

# Filter function (copied/adapted from rfplot.py / rffind.py if needed for -q option)
# Note: rfpng.c calls filter(), but its usage isn't obvious from just rfpng.c.
# It might be related to the -q option mentioned in rfplot.c's usage (Detect and plot peaks?).
# For now, we'll omit the direct call to filter within rfpng.py unless clarified.
# If needed, the `filter_data` function from rfplot.py translation could be added/imported.

def main():
    parser = argparse.ArgumentParser(description='rfpng: Create PNG plot of RF observations.')
    # Arguments matching rfpng.c (very similar to rfplot.c)
    parser.add_argument('-p', '--path', required=True, help='Input path prefix for /path/prefix_xxxxxx.bin files')
    parser.add_argument('-o', '--output', default=None, help='Output PNG filename override [default: <timestamp>_<freq>.png]')
    parser.add_argument('-O', '--offset', type=float, default=0.0, help='Frequency offset to apply during read (Hz) [0.0]')
    parser.add_argument('-s', '--start', type=int, default=0, help='Number of starting .bin file index [0]')
    parser.add_argument('-l', '--length', type=int, default=3600, help='Number of subintegrations to read [3600]')
    parser.add_argument('-b', '--nbin', type=int, default=1, help='Number of subintegrations to bin [1]')
    parser.add_argument('-z', '--zmax', type=float, default=None, help='Image scaling upper limit override (contrast factor)')
    parser.add_argument('-f', '--freq', type=float, default=0.0, help='Frequency to zoom into (Hz) [0 = full bw]')
    parser.add_argument('-w', '--bw', type=float, default=0.0, help='Bandwidth to zoom into (Hz) [0 = full bw]')
    parser.add_argument('-c', '--catalog', default=None, help='TLE catalog file (default from ST_TLEDIR)')
    parser.add_argument('-F', '--freqlist', default=None, help='List with satellite frequencies (default from ST_DATADIR)')
    parser.add_argument('-C', '--site', type=int, default=None, help='Observer Site ID (default from ST_COSPAR)')
    parser.add_argument('-g', '--graves', action='store_true', help='GRAVES mode (bistatic)')
    # parser.add_argument('-S', '--sigma', type=float, default=5.0, help='Significance for peak detection (for -q) [5.0]') # For filter/peakfind
    # parser.add_argument('-q', action='store_true', help='Detect and plot peaks') # Requires filter/peakfind logic
    parser.add_argument('-W', '--width', type=float, default=14.0, help='Plot width in inches [14.0]')
    parser.add_argument('-A', '--aspect', type=float, default=1.0, help='Plot aspect ratio (height/width) [1.0]')
    parser.add_argument('-m', '--cmap', type=int, default=2, choices=[0,1,2,3], help='Colormap index [0:cool, 1:heat, 2:viridis, 3:gray; default:2(viridis)]')
    parser.add_argument('--no-overlay', action='store_true', help='Disable satellite trace overlay')
    parser.add_argument('--show-names', action='store_true', help='Show satellite names on overlay')


    args = parser.parse_args()

    # --- Determine Site ID ---
    if not args.no_overlay: # Only need site ID if plotting overlays
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
                print("Error: Site ID must be provided via -C or ST_COSPAR for trace overlay.")
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

    # --- Determine Output Filename ---
    output_filename = args.output
    if output_filename is None:
        # Default filename: <timestamp>_<freq_mhz>.png (similar to rfplot default)
        timestamp_str = spectrogram.nfd0.replace(":", "").replace("-", "").replace("T", "_") # Basic formatting
        freq_mhz_str = f"{spectrogram.freq / 1e6:.3f}"
        output_filename = f"{timestamp_str}_{freq_mhz_str}.png"
    print(f"Output filename: {output_filename}")

    # --- Load TLEs and Compute Traces ---
    traces = []
    if not args.no_overlay:
        tle_data = load_tles(args.catalog) # Handles env vars internally
        if tle_data.number_of_elements > 0:
            print("Computing satellite traces...")
            traces = compute_trace(
                tle_data=tle_data,
                mjd_times=spectrogram.mjd,
                site_id=args.site,
                center_freq_hz=spectrogram.freq,
                bandwidth_hz=spectrogram.samp_rate,
                graves=args.graves,
                freqlist_path=args.freqlist
            )
        else:
            print("Warning: No TLEs loaded, cannot compute traces.")

    # --- Setup Plot ---
    plot_height = args.width * args.aspect
    fig, ax = plt.subplots(figsize=(args.width, plot_height))

    # Determine plot limits and aspect ratio
    time_extent = spectrogram.nsub
    freq_center_mhz = spectrogram.freq / 1e6
    freq_bw_mhz = spectrogram.samp_rate / 1e6
    freq_min_mhz = freq_center_mhz - 0.5 * freq_bw_mhz
    freq_max_mhz = freq_center_mhz + 0.5 * freq_bw_mhz

    # Image extent for imshow
    extent = [0, time_extent, freq_min_mhz, freq_max_mhz]

    # Use zscale results for color limits
    vmin, vmax = spectrogram.zmin, spectrogram.zmax
    if args.zmax is not None:
        # Simple override interpretation (same as rfplot.py translation)
        if spectrogram.zstd is not None and np.mean(spectrogram.zstd) > 0:
             mean_avg = np.mean(spectrogram.zavg[np.isfinite(spectrogram.zavg)])
             mean_std = np.mean(spectrogram.zstd[np.isfinite(spectrogram.zstd)])
             vmax = mean_avg + args.zmax * mean_std
             vmin = mean_avg - args.zmax * mean_std / 4.0 # Guessing a sensible vmin offset
        else:
             vmax = args.zmax # If std dev calculation failed, just use value directly? Needs better logic.

    # Select colormap based on index
    cmap_names = ['cool', 'hot', 'viridis', 'gray']
    try:
        cmap = plt.get_cmap(cmap_names[args.cmap])
    except (IndexError, ValueError):
        warnings.warn(f"Invalid cmap index {args.cmap}, defaulting to viridis.")
        cmap = plt.get_cmap('viridis')

    # --- Display Image ---
    im = ax.imshow(
        spectrogram.z,
        aspect='auto',
        origin='lower',
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation='nearest'
    )

    # --- Configure Axes ---
    ax.set_xlabel(f"Time starting {spectrogram.nfd0} UTC")
    ax.set_ylabel(f"Frequency Offset (kHz) from {freq_center_mhz:.6f} MHz")

    # Custom formatters
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(time_axis_formatter(spectrogram.mjd, spectrogram.nfd0)))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(freq_axis_formatter(freq_center_mhz)))

    # Set view limits
    ax.set_xlim(0, time_extent)
    ax.set_ylim(freq_min_mhz, freq_max_mhz)

    # --- Plot Traces ---
    if not args.no_overlay:
        plot_traces(ax, traces, spectrogram.freq, args.show_names)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Intensity (Arb. Units / Scaled)')

    ax.set_title(f"Spectrogram: {os.path.basename(args.path)} [{spectrogram.nfd0}]")
    plt.tight_layout()

    # --- Save PNG ---
    try:
        plt.savefig(output_filename, dpi=150) # Adjust dpi as needed
        print(f"Plot saved to {output_filename}")
    except Exception as e:
        print(f"Error saving PNG file {output_filename}: {e}")

    plt.close(fig) # Close figure to free memory

if __name__ == "__main__":
    main()
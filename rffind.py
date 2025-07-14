import numpy as np
import argparse
import sys
import os
import math
import warnings

# Import previously translated modules
try:
    from rfio import read_spectrogram, Spectrogram # Assuming class from rfio.py
    from rftime import mjd_to_nfd # Needed? C uses it for diagnostics maybe
except ImportError as e:
    print(f"Error importing helper modules (rfio, rftime): {e}")
    print("Please ensure rfio.py and rftime.py are available.")
    sys.exit(1)

def find_peaks_in_spectrogram(spectrogram, site_id, sigma_thresh, graves=False, output_filename="find.dat"):
    """
    Filters spectrogram data based on sigma clipping and writes peaks.
    Based on filter() logic in rffind.c / rfplot.c.

    Args:
        spectrogram (Spectrogram): The loaded spectrogram data object.
        site_id (int): The observer's site ID.
        sigma_thresh (float): The sigma threshold for peak detection.
        graves (bool): If True, output includes GRAVES remote site ID (9999).
        output_filename (str): Name of the file to append peaks to.
    """
    if spectrogram.nsub == 0 or spectrogram.nchan == 0:
        print("Warning: Cannot process empty spectrogram.")
        return

    print(f"Finding peaks with sigma > {sigma_thresh}, appending to {output_filename}")
    nsub = spectrogram.nsub
    nchan = spectrogram.nchan
    # Data shape is assumed (nchan, nsub) from rfio.py translation
    data = spectrogram.z

    peaks_found_total = 0
    try:
        # Append to the output file ('a' mode)
        with open(output_filename, "a") as f:
            # Write header if file is new/empty - check size
            if f.tell() == 0:
                 if graves:
                     f.write("# MJD Frequency_Hz Sigma SiteID RemoteSiteID\n")
                 else:
                     f.write("# MJD Frequency_Hz Sigma SiteID\n")

            for i in range(nsub): # Loop through time (subintegrations)
                subint_data = data[:, i] # Get data for this subintegration (one column)
                mask = np.ones(nchan, dtype=bool) # Start with all pixels good

                # Check for all NaN column
                if np.all(np.isnan(subint_data)):
                    continue

                # --- Iterative Sigma Clipping ---
                # matching loop structure
                last_mask_sum = -1
                for k in range(10): # Max 10 iterations like C code
                    valid_data = subint_data[mask & np.isfinite(subint_data)] # Exclude NaNs too
                    if len(valid_data) < 2: break # Need >= 2 points for std dev

                    avg = np.mean(valid_data) # Use nanmean? No, already filtered NaNs
                    std = np.std(valid_data)

                    if std < 1e-9: break # Avoid division by zero or meaningless std dev

                    # Update mask: remove points further than threshold from mean
                    with np.errstate(invalid='ignore'): # Ignore warnings for NaN comparisons
                        new_mask_cond = np.abs(subint_data - avg) <= (sigma_thresh * std)
                    # Only consider previously good points for potentially becoming bad, and non-NaNs
                    new_mask = mask & new_mask_cond & np.isfinite(subint_data)

                    current_mask_sum = np.sum(new_mask)
                    if current_mask_sum == last_mask_sum: # Check if mask changed
                        break # Converged
                    last_mask_sum = current_mask_sum
                    mask = new_mask # Update mask for next iteration

                # --- Find Peaks Above Threshold in the *original* data ---
                # Use the *final* mask to calculate baseline stats
                final_valid_data = subint_data[mask & np.isfinite(subint_data)]
                if len(final_valid_data) < 2: continue # Need baseline
                final_avg = np.mean(final_valid_data)
                final_std = np.std(final_valid_data)
                if final_std < 1e-9: continue

                # Identify potential peaks in the *original* data using the calculated baseline
                with np.errstate(invalid='ignore'):
                     sigma_values = (subint_data - final_avg) / final_std
                     is_potential_peak = (sigma_values > sigma_thresh) & np.isfinite(subint_data)
                peak_indices = np.where(is_potential_peak)[0]

                if len(peak_indices) == 0: continue

                # --- Local Maxima Filter ---
                # Keep only peaks that are greater than their immediate neighbours (if neighbours are also peaks)
                local_max_mask = np.zeros(nchan, dtype=bool)
                for idx in peak_indices:
                    is_local_max = True
                    # Check left neighbor (if it's also a potential peak)
                    if idx > 0 and is_potential_peak[idx - 1] and subint_data[idx] <= subint_data[idx - 1]:
                        is_local_max = False
                    # Check right neighbor (if it's also a potential peak)
                    if idx < nchan - 1 and is_potential_peak[idx + 1] and subint_data[idx] <= subint_data[idx + 1]:
                        is_local_max = False
                    if is_local_max:
                        local_max_mask[idx] = True

                # --- Write detected local maxima peaks ---
                mjd = spectrogram.mjd[i]
                if mjd <= 1.0: continue # Skip invalid MJD

                peaks_found_this_subint = 0
                for j in np.where(local_max_mask)[0]:
                    # Calculate frequency in Hz for the peak channel j
                    # Use channel center frequency
                    channel_bw_hz = spectrogram.samp_rate / spectrogram.nchan
                    f_hz = (spectrogram.freq - 0.5 * spectrogram.samp_rate +
                            (j + 0.5) * channel_bw_hz)
                    sigma_val = sigma_values[j] # Use pre-calculated sigma

                    if graves:
                        f.write(f"{mjd:.8f} {f_hz:.3f} {sigma_val:.3f} {site_id} 9999\n")
                    else:
                        f.write(f"{mjd:.8f} {f_hz:.3f} {sigma_val:.3f} {site_id}\n")
                    peaks_found_this_subint += 1

                peaks_found_total += peaks_found_this_subint

    except IOError as e:
        print(f"Error writing output file {output_filename}: {e}")

    print(f"Found {peaks_found_total} total peaks.")


def main():
    parser = argparse.ArgumentParser(description='rffind: Find signals in RF observations')
    parser.add_argument('-p', '--path', required=True, help='Input path prefix for /path/prefix_xxxxxx.bin files')
    parser.add_argument('-s', '--start', type=int, default=0, help='Number of starting .bin file index [0]')
    parser.add_argument('-l', '--length', type=int, default=0, help='Number of subintegrations to process (0=all) [0]')
    parser.add_argument('-f', '--freq', type=float, default=0.0, help='Frequency to zoom into (Hz) [0 = full bandwidth]')
    parser.add_argument('-w', '--bw', type=float, default=0.0, help='Bandwidth to zoom into (Hz) [0 = full bandwidth]')
    parser.add_argument('-o', '--output', type=str, default="find.dat", help='Output filename to append peaks to [find.dat]')
    parser.add_argument('-C', '--site', type=int, default=None, help='Observer Site ID (default from ST_COSPAR)')
    parser.add_argument('-g', '--graves', action='store_true', help='GRAVES mode output format')
    parser.add_argument('-S', '--sigma', type=float, default=5.0, help='Sigma threshold for peak detection [5.0]')
    # Note: -b (binning) and -foff (freq offset) from rfplot are not options in rffind.c

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
            # C code doesn't exit if site_id is missing, but uses it in output. Default to 0?
            # Let's require it for clarity.
            print("Error: Site ID must be provided via -C or ST_COSPAR environment variable.")
            sys.exit(1)

    # --- Process Files ---
    # If length is 0, we need to process file by file until none are found
    if args.length == 0:
        current_file_index = args.start
        while True:
            print(f"\nProcessing file index: {current_file_index}")
            # Use nsub_to_read=0 in read_spectrogram to read all subints in this file
            spectrogram = read_spectrogram(
                prefix=args.path,
                isub=current_file_index,
                nsub_to_read=0, # Let read_spectrogram determine from header/EOF
                f0=args.freq,
                df0=args.bw,
                nbin=1, # No binning in rffind
                foff=0.0 # No offset option in rffind
            )
            if spectrogram.nsub == 0:
                print(f"No more files found starting from index {current_file_index} or file empty.")
                break
            find_peaks_in_spectrogram(spectrogram, args.site, args.sigma, args.graves, args.output)
            current_file_index += 1 # Move to next file index implicitly (assuming _000000, _000001...)
    else:
        # Read the specified number of subintegrations starting from isub
        spectrogram = read_spectrogram(
            prefix=args.path,
            isub=args.start,
            nsub_to_read=args.length,
            f0=args.freq,
            df0=args.bw,
            nbin=1, # No binning in rffind
            foff=0.0 # No offset option in rffind
        )
        if spectrogram.nsub > 0:
            find_peaks_in_spectrogram(spectrogram, args.site, args.sigma, args.graves, args.output)
        else:
             print("No data read.")

    print("Processing complete.")

if __name__ == "__main__":
    main()
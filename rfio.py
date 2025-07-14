import os
import re
import numpy as np
import warnings
from datetime import timedelta

# Assuming rftime.py contains the necessary time conversion functions
from rftime import nfd_to_mjd, mjd_to_nfd

# Assuming zscale.py contains the zscale function
# from zscale import zscale_samples # Function name might differ based on translation

class Spectrogram:
    """Class to hold spectrogram data and metadata."""
    def __init__(self):
        self.nsub = 0         # Number of subintegrations in this object
        self.nchan = 0        # Number of channels in this object
        self.msub = 0         # Total number of subintegrations in original file series (from header)
        self.isub = 0         # Starting subintegration index (file number)
        self.mjd = None       # Numpy array of MJDs for each subintegration
        self.freq = 0.0       # Center frequency (Hz)
        self.samp_rate = 0.0  # Bandwidth (Hz)
        self.length = None    # Numpy array of integration lengths (seconds) per subintegration
        self.z = None         # 2D Numpy array (nchan, nsub) of spectrogram data (or nsub, nchan depending on read logic)
        self.zavg = None      # Numpy array of average power per subintegration
        self.zstd = None      # Numpy array of standard deviation per subintegration
        self.zmin = 0.0       # Min value for display scaling (e.g., from zscale)
        self.zmax = 1.0       # Max value for display scaling (e.g., from zscale)
        self.nfd0 = ""        # Timestamp string of the first subintegration read

def _parse_header(header_bytes):
    """Parses the 256-byte header."""
    header_dict = {
        'nfd': "", 'freq': 0.0, 'samp_rate': 0.0, 'length': 0.0,
        'nchan': 0, 'msub': 0, 'nbits': -32, 'zavg': 0.0, 'zstd': 1.0
    }
    try:
        # Decode assuming ASCII, remove null bytes
        header_str = header_bytes.decode('ascii', errors='ignore').rstrip('\x00')
    except UnicodeDecodeError:
        warnings.warn("Header decoding failed, might be corrupted.")
        return None

    # Use regex to capture values robustly, handling potential missing fields
    patterns = {
        'nfd': r"UTC_START\s+(\S+)",
        'freq': r"FREQ\s+([\d\.\-eE]+)\s+Hz",
        'samp_rate': r"BW\s+([\d\.\-eE]+)\s+Hz",
        'length': r"LENGTH\s+([\d\.\-eE]+)\s+s",
        'nchan': r"NCHAN\s+(\d+)",
        'msub': r"NSUB\s+(\d+)",
        'nbits': r"NBITS\s+(\d+)",
        'zavg': r"MEAN\s+([\d\.\-eE]+)",
        'zstd': r"RMS\s+([\d\.\-eE]+)",
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, header_str)
        if match:
            value_str = match.group(1)
            try:
                if key in ['freq', 'samp_rate', 'length', 'zavg', 'zstd']:
                    header_dict[key] = float(value_str)
                elif key in ['nchan', 'msub', 'nbits']:
                    header_dict[key] = int(value_str)
                elif key == 'nfd':
                    header_dict[key] = value_str
            except ValueError:
                warnings.warn(f"Could not parse value for {key}: '{value_str}'")
        # else: # Optional: warn if a field is missing, except for optional ones like NBITS
            # if key not in ['nbits', 'zavg', 'zstd']:
               # warnings.warn(f"Header field '{key}' not found.")

    # Handle case where NBITS is not present (defaults to -32 in C code)
    if 'nbits' not in header_dict or header_dict['nbits'] == 0: # Check if NBITS wasn't found or parsed
         # Check if the string actually contains NBITS 8 to be sure
        if "NBITS         8" not in header_str:
            header_dict['nbits'] = -32
        # else: NBITS was found but parsing failed - already warned


    return header_dict


def read_spectrogram(prefix, isub=0, nsub_to_read=0, f0=0.0, df0=0.0, nbin=1, foff=0.0):
    """
    Reads spectrogram data from a series of .bin files.

    Args:
        prefix (str): The base path and filename prefix (e.g., "/path/to/data/obs").
        isub (int): Starting file index (suffix number). Defaults to 0.
        nsub_to_read (int): Number of subintegrations to read in total across files.
                            If 0, tries to read all available based on the first
                            file's NSUB header. Defaults to 0.
        f0 (float): Center frequency (Hz) to zoom into. If > 0, df0 must also be > 0.
        df0 (float): Bandwidth (Hz) to zoom into. If > 0, f0 must also be > 0.
        nbin (int): Number of subintegrations to bin together. Defaults to 1.
        foff (float): Frequency offset to add to the header frequency (Hz). Defaults to 0.0.

    Returns:
        Spectrogram: A Spectrogram object containing the data, or an empty
                     Spectrogram object with nsub=0 on failure.
    """
    s = Spectrogram()
    s.isub = isub

    first_file_path = f"{prefix}_{isub:06d}.bin"
    if not os.path.exists(first_file_path):
        print(f"Error: First file does not exist: {first_file_path}")
        return s # Return empty spectrogram

    # --- Read header from the first file to get parameters ---
    try:
        with open(first_file_path, "rb") as f:
            header_bytes = f.read(256)
            if len(header_bytes) < 256:
                 print(f"Error: Could not read full header from {first_file_path}")
                 return s
            first_header = _parse_header(header_bytes)
            if not first_header:
                print(f"Error: Could not parse header from {first_file_path}")
                return s

            s.freq = first_header['freq'] + foff # Apply offset immediately
            s.samp_rate = first_header['samp_rate']
            s.nfd0 = first_header['nfd']
            s.msub = first_header['msub'] # Total subints expected in this file series
            nch_total = first_header['nchan'] # Total channels in the file

    except IOError as e:
        print(f"Error reading first file {first_file_path}: {e}")
        return s

    if nch_total <= 0:
        print(f"Error: Invalid number of channels ({nch_total}) in header.")
        return s

    # --- Determine channel range for reading/plotting ---
    j0, j1 = 0, nch_total # Default: read all channels
    if f0 > 0.0 and df0 > 0.0:
        # Calculate channel indices based on desired frequency range
        center_freq_hdr = s.freq # Use the already offset frequency
        bw_hdr = s.samp_rate
        f_start_hdr = center_freq_hdr - 0.5 * bw_hdr
        chan_bw = bw_hdr / nch_total

        f_start_req = f0 - 0.5 * df0
        f_end_req = f0 + 0.5 * df0

        j0 = int(round((f_start_req - f_start_hdr) / chan_bw))
        j1 = int(round((f_end_req - f_start_hdr) / chan_bw))

        # Clamp indices to valid range
        j0 = max(0, j0)
        j1 = min(nch_total, j1)

        if j0 >= j1:
            print(f"Error: Calculated channel range is invalid ({j0} >= {j1}). Check f0/df0.")
            return s
        if j0 < 0 or j1 > nch_total : # Should be caught by clamp, but double check
             print(f"Error: Requested frequency range [{f0-0.5*df0:.3f} - {f0+0.5*df0:.3f}] Hz "
                   f"is outside the file's range [{f_start_hdr:.3f} - {f_start_hdr+bw_hdr:.3f}] Hz.")
             # Reset to full range? Or return error? C code returns error.
             return s

        s.nchan = j1 - j0 # Number of channels to actually store/process
        # Update center frequency and bandwidth to reflect the selected range
        s.freq = f_start_hdr + (j0 + s.nchan / 2.0) * chan_bw
        s.samp_rate = s.nchan * chan_bw
        print(f"Zooming into channels {j0} to {j1-1} ({s.nchan} channels)")
        print(f"Adjusted center frequency: {s.freq*1e-6:.6f} MHz")
        print(f"Adjusted bandwidth: {s.samp_rate*1e-6:.6f} MHz")

    else:
        s.nchan = nch_total # Use all channels


    # --- Determine number of subintegrations to read ---
    if nsub_to_read <= 0 and s.msub > 0:
        nsub_total_expected = s.msub # Read all subints declared in the header
        print(f"Reading all {s.msub} subintegrations specified in header.")
    elif nsub_to_read > 0:
        nsub_total_expected = nsub_to_read
    else:
        print("Warning: nsub_to_read=0 and NSUB not found/invalid in header. Reading limited files.")
        nsub_total_expected = 10000 # Arbitrary limit if NSUB missing

    s.nsub = (nsub_total_expected + nbin -1) // nbin # Ceiling division for binned output


    # --- Allocate Numpy arrays ---
    print(f"Allocating for {s.nsub} binned subintegrations and {s.nchan} channels.")
    try:
        # Using float32 to match C code's float size
        # Array shape: (nchan, nsub) for easier channel access later? Or (nsub, nchan)?
        # C code uses s.z[i+s.nsub*j] -> implies column-major or (nchan, nsub) if j is channel, i is time
        # Let's use (nchan, nsub) for consistency
        s.z = np.zeros((s.nchan, s.nsub), dtype=np.float32)
        s.mjd = np.zeros(s.nsub, dtype=np.float64)
        s.length = np.zeros(s.nsub, dtype=np.float32)
    except MemoryError:
        print(f"Error: Failed to allocate memory for spectrogram ({s.nchan}x{s.nsub})")
        return Spectrogram() # Return empty


    # --- Loop through files and read data ---
    subint_count_total = 0 # Total subints read from files
    output_subint_idx = 0  # Index for the binned output array (s.z, s.mjd, ...)
    current_bin_count = 0  # Subints accumulated in the current bin
    current_file_idx = isub
    keep_reading_files = True

    while keep_reading_files and subint_count_total < nsub_total_expected and output_subint_idx < s.nsub:
        file_path = f"{prefix}_{current_file_idx:06d}.bin"
        if not os.path.exists(file_path):
            print(f"Info: File {file_path} does not exist. Stopping read.")
            keep_reading_files = False
            break

        print(f"Opening {file_path}")
        try:
            with open(file_path, "rb") as f:
                while subint_count_total < nsub_total_expected and output_subint_idx < s.nsub:
                    # Read header for this subint
                    header_bytes = f.read(256)
                    if len(header_bytes) < 256:
                        # End of file reached prematurely
                        keep_reading_files = False
                        break
                    header = _parse_header(header_bytes)
                    if not header:
                        warnings.warn(f"Could not parse header in {file_path} at subint offset, stopping.")
                        keep_reading_files = False
                        break

                    # Determine data type and read data
                    nbits = header['nbits']
                    dtype = np.float32 if nbits == -32 else np.int8
                    itemsize = 4 if nbits == -32 else 1
                    count_to_read = nch_total # Read all channels from file

                    data_bytes = f.read(itemsize * count_to_read)
                    if len(data_bytes) < itemsize * count_to_read:
                        warnings.warn(f"Incomplete data read in {file_path}, stopping.")
                        keep_reading_files = False
                        break

                    raw_data = np.frombuffer(data_bytes, dtype=dtype, count=count_to_read)

                    # Select channels if zooming
                    channel_data = raw_data[j0:j1]

                    # Convert to float32 and scale if needed
                    if nbits == 8:
                        zavg = header['zavg']
                        zstd = header['zstd']
                        # Apply scaling: 6.0/256.0 * val * zstd + zavg (from C)
                        # factor = (6.0 / 256.0) * zstd # Precompute factor
                        # float_data = factor * channel_data.astype(np.float32) + zavg
                        # Alternative: simple cast if scaling is just for display? Check purpose.
                        # Let's assume direct conversion needed for processing
                        float_data = channel_data.astype(np.float32) * (6.0/256.0 * zstd) + zavg
                    else: # Already float32
                        float_data = channel_data # Sliced view, no copy needed yet

                    # Accumulate data for binning
                    # Add to the correct column (output_subint_idx)
                    s.z[:, output_subint_idx] += float_data
                    s.mjd[output_subint_idx] += nfd_to_mjd(header['nfd']) + 0.5 * header['length'] / 86400.0
                    s.length[output_subint_idx] += header['length']
                    current_bin_count += 1
                    subint_count_total += 1

                    # If bin is full, finalize it and move to the next output slot
                    if current_bin_count == nbin:
                        s.z[:, output_subint_idx] /= nbin
                        s.mjd[output_subint_idx] /= nbin
                        # s.length is total length in bin, no division needed? C code sums length.
                        output_subint_idx += 1
                        current_bin_count = 0
                        # Reset accumulation slot if moving to a new one
                        # if output_subint_idx < s.nsub:
                        #    s.z[:, output_subint_idx] = 0.0
                        #    s.mjd[output_subint_idx] = 0.0
                        #    s.length[output_subint_idx] = 0.0

        except IOError as e:
            print(f"Error reading file {file_path}: {e}")
            keep_reading_files = False # Stop reading subsequent files

        current_file_idx += 1 # Move to the next file index


    # --- Finalize last bin if partially filled ---
    if current_bin_count > 0 and output_subint_idx < s.nsub:
        s.z[:, output_subint_idx] /= current_bin_count
        s.mjd[output_subint_idx] /= current_bin_count
        # Adjust nsub if we read fewer than expected due to EOF or file missing
        s.nsub = output_subint_idx + 1
        print(f"Finalizing partially filled bin. Actual number of output subints: {s.nsub}")
    elif output_subint_idx < s.nsub:
         # No partial bin, but fewer output bins than allocated
         s.nsub = output_subint_idx
         print(f"Read fewer subints than expected. Actual number of output subints: {s.nsub}")


    # --- Trim arrays if fewer subints were read/created ---
    if s.nsub < len(s.mjd):
        print(f"Trimming arrays to actual size: {s.nsub}")
        s.z = s.z[:, :s.nsub]
        s.mjd = s.mjd[:s.nsub]
        s.length = s.length[:s.nsub]

    if s.nsub == 0:
        print("Warning: No subintegrations were successfully read.")
        return Spectrogram() # Return empty

    # --- Calculate statistics and scaling ---
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning) # Ignore mean of empty slice etc.

        # Calculate average and std dev across channels for each subintegration
        s.zavg = np.nanmean(s.z, axis=0) # Mean along channel axis
        s.zstd = np.nanstd(s.z, axis=0)  # Std dev along channel axis

        # Calculate overall min/max for scaling (simple approach)
        # valid_z = s.z[np.isfinite(s.z)]
        # if valid_z.size > 0:
        #     avg_all = np.mean(s.zavg[np.isfinite(s.zavg)]) # Avg of averages
        #     std_all = np.mean(s.zstd[np.isfinite(s.zstd)]) # Avg of stddevs
        #     s.zmin = avg_all - 1.0 * std_all # Match C code logic?
        #     s.zmax = avg_all + 1.0 * std_all # Match C code logic?
        # else:
        #     s.zmin, s.zmax = 0.0, 1.0

    # Use zscale for better scaling (requires zscale function)
    # Placeholder: Needs the actual zscale implementation adapted for numpy arrays
    try:
        # Adapt zscale call based on its Python signature
        # Assuming zscale_samples takes a numpy array (e.g., flattened or 2D)
        # The C zscale uses the Spectrogram struct directly.
        # We might need to pass s.z and potentially dimensions.
        # scaled_zmin, scaled_zmax = zscale_samples(s.z, nsamples=s.z.size // 10, contrast=0.25) # Example call
        # s.zmin = scaled_zmin
        # s.zmax = scaled_zmax
        # print(f"Zscale results: zmin={s.zmin}, zmax={s.zmax}")
        # --- Fallback scaling if zscale is not available ---
        if 'scaled_zmin' not in locals():
             valid_z = s.z[np.isfinite(s.z)]
             if valid_z.size > 0:
                  s.zmin = np.percentile(valid_z, 1)  # Robust min
                  s.zmax = np.percentile(valid_z, 99) # Robust max
             else:
                  s.zmin, s.zmax = 0.0, 1.0
             print(f"Using percentile scaling: zmin={s.zmin:.3f}, zmax={s.zmax:.3f}")
    except NameError:
        print("Warning: zscale function not found. Using basic percentile scaling.")
        valid_z = s.z[np.isfinite(s.z)]
        if valid_z.size > 0:
            s.zmin = np.percentile(valid_z, 1)
            s.zmax = np.percentile(valid_z, 99)
        else:
            s.zmin, s.zmax = 0.0, 1.0
        print(f"Using percentile scaling: zmin={s.zmin:.3f}, zmax={s.zmax:.3f}")
    except Exception as e:
        print(f"Error during zscale: {e}. Using basic percentile scaling.")
        # Fallback scaling here too
        valid_z = s.z[np.isfinite(s.z)]
        if valid_z.size > 0:
            s.zmin = np.percentile(valid_z, 1)
            s.zmax = np.percentile(valid_z, 99)
        else:
            s.zmin, s.zmax = 0.0, 1.0
        print(f"Using percentile scaling: zmin={s.zmin:.3f}, zmax={s.zmax:.3f}")


    print(f"Successfully read {s.nsub} binned subintegrations.")
    return s


def write_spectrogram(s, prefix):
    """
    Writes spectrogram data to a .bin file (simplified version).
    NOTE: This only writes ONE file, unlike the C version which might
          implicitly handle splitting based on header info. It also
          only writes float32 data, not the int8 format.
    """
    if s.nsub == 0 or s.nchan == 0:
        print("Error: Spectrogram object is empty, cannot write.")
        return

    filename = f"{prefix}_000000.bin" # Simple filename for the whole dataset
    print(f"Writing spectrogram to {filename}")

    try:
        with open(filename, "wb") as f:
            # Write subintegrations sequentially
            for i in range(s.nsub):
                # --- Create Header ---
                # Get NFD timestamp for the start of the subintegration
                mjd_start = s.mjd[i] - 0.5 * s.length[i] / 86400.0
                nfd = mjd_to_nfd(mjd_start) # Requires rftime module

                # Format header string - Adjust NSUB if needed
                # Using the original full bandwidth/channels before zoom for BW/NCHAN?
                # C code seems to write header per subint using current s state.
                header_str = (
                    f"HEADER\n"
                    f"UTC_START    {nfd}\n"
                    f"FREQ         {s.freq:.6f} Hz\n" # Use potentially adjusted freq
                    f"BW           {s.samp_rate:.6f} Hz\n" # Use potentially adjusted bw
                    f"LENGTH       {s.length[i]:f} s\n"
                    f"NCHAN        {s.nchan}\n" # Use potentially adjusted nchan
                    f"NSUB         {s.nsub}\n" # Total subints being written here
                    # f"NBITS         -32\n" # Implicitly float32
                    f"END\n"
                )
                # Pad header to 256 bytes with null characters
                header_bytes = header_str.encode('ascii').ljust(256, b'\x00')
                f.write(header_bytes)

                # --- Write Data ---
                # Data shape is (nchan, nsub), take the i-th column
                data_slice = s.z[:, i].astype(np.float32)
                f.write(data_slice.tobytes())

    except IOError as e:
        print(f"Error writing spectrogram file {filename}: {e}")
    except Exception as e:
         print(f"An unexpected error occurred during writing: {e}")


# Note: free_spectrogram is not needed as Python handles memory management.
# If using large numpy arrays, ensuring they go out of scope or using `del`
# might help release memory sooner, but it's usually automatic.
import argparse
import sys
import os
import re # Needed for header parsing logic if copied directly
import warnings # Needed for header parsing logic if copied directly

# --- Header Parsing Logic (copied/adapted from rfio.py) ---
# Ideally, this would be imported from rfio.py, but for a standalone
# script, we can include it directly.

def _parse_header(header_bytes):
    """Parses the 256-byte header."""
    header_dict = {
        'nfd': "", 'freq': 0.0, 'samp_rate': 0.0, 'length': 0.0,
        'nchan': 0, 'msub': 0, 'nbits': -32, 'zavg': 0.0, 'zstd': 1.0
    }
    try:
        header_str = header_bytes.decode('ascii', errors='ignore').rstrip('\x00')
    except UnicodeDecodeError:
        warnings.warn("Header decoding failed, might be corrupted.")
        return None

    patterns = {
        'freq': r"FREQ\s+([\d\.\-eE]+)\s+Hz",
        'samp_rate': r"BW\s+([\d\.\-eE]+)\s+Hz",
        'nchan': r"NCHAN\s+(\d+)",
        # Add other patterns if needed, but rfinfo only uses freq, bw, nchan
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, header_str)
        if match:
            value_str = match.group(1)
            try:
                if key in ['freq', 'samp_rate']:
                    header_dict[key] = float(value_str)
                elif key in ['nchan']:
                    header_dict[key] = int(value_str)
            except ValueError:
                 warnings.warn(f"Could not parse value for {key}: '{value_str}'")
        else:
             # Only print warning if essential field missing
             if key in ['freq', 'samp_rate', 'nchan']:
                  warnings.warn(f"Header field '{key}' not found.")
                  return None # Essential field missing

    # Check if essential fields were parsed successfully
    if header_dict['nchan'] == 0 or header_dict['samp_rate'] == 0.0:
         warnings.warn("Essential header fields (NCHAN, BW) could not be parsed.")
         return None

    return header_dict

# --- Main Script Logic ---

def get_spectrogram_info(file_argument):
    """
    Reads header info from the first file and counts subsequent files.

    Args:
        file_argument (str): The path prefix or the full path to the first file
                             (e.g., "path/prefix" or "path/prefix_000000.bin").
    """
    # Determine the file root/prefix
    if file_argument.endswith("_000000.bin"):
        prefix = file_argument[:-11]
    elif re.search(r'_\d{6,}\.bin$', file_argument):
        # If user provided a different numbered file, extract prefix
        prefix = re.sub(r'_\d{6,}\.bin$', '', file_argument)
        print(f"Warning: Input file seems numbered, using prefix: {prefix}")
    else:
        # Assume the input is the prefix directly
        prefix = file_argument

    first_file_path = f"{prefix}_000000.bin"

    if not os.path.exists(first_file_path):
        print(f"Error: First file not found: {first_file_path}")
        return

    # --- Read header from the first file ---
    try:
        with open(first_file_path, "rb") as f:
            header_bytes = f.read(256)
            if len(header_bytes) < 256:
                print(f"Error: Incomplete header read from {first_file_path}")
                return
            header_info = _parse_header(header_bytes)
            if not header_info:
                print(f"Error: Failed to parse header from {first_file_path}")
                return
    except IOError as e:
        print(f"Error reading file {first_file_path}: {e}")
        return

    # --- Count subsequent files ---
    file_count = 0
    current_index = 0
    while True:
        current_file_path = f"{prefix}_{current_index:06d}.bin"
        if os.path.exists(current_file_path):
            file_count += 1
            current_index += 1
        else:
            break # Stop counting when a file is missing

    # --- Print Summary ---
    # Format matches C output: prefix freq_MHz bw_MHz nchan nfiles
    # line 41
    print(f"{prefix} {header_info['freq']/1e6:8.3f} {header_info['samp_rate']/1e6:8.3f} "
          f"{header_info['nchan']:d} {file_count:d}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='rfinfo: Print basic info about a spectrogram file sequence.'
    )
    parser.add_argument(
        'file_prefix',
        metavar='PREFIX',
        type=str,
        help='Path prefix for the spectrogram files (e.g., /path/to/data/obs). '
             'Can also be the full path to the first file (prefix_000000.bin).'
    )

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        print("\nError: Missing file prefix argument.")
        sys.exit(1)

    # Workaround: argparse expects options first. Since rfinfo.c just uses argv[1]
    # directly, we mimic that behaviour slightly differently.
    # args = parser.parse_args() # This won't work as expected without changes
    # Instead, just grab the first argument after the script name
    file_arg = sys.argv[1]

    # Basic help check
    if file_arg in ['-h', '--help']:
         parser.print_help()
         sys.exit(0)

    get_spectrogram_info(file_arg)
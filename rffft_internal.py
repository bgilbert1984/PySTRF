import re
import os
from datetime import datetime
import warnings

# Mapping from month abbreviations to numbers (for SDR Console format)
MONTH_MAP = {
    "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
    "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12
}

# Mapping from file extensions/types to format characters used in rffft
FORMAT_MAP = {
    's8': 'c',  # char
    's16': 'i', # int16
    'f32': 'f', # float32
    'wav': 'w', # wav file (handled by sox library in C)
    'gqrx_fc': 'f' # Gqrx uses float32
}

def rffft_params_from_filename(filename):
    """
    Parses parameters (samplerate, frequency, format, starttime) from
    known SDR filename conventions (SatDump, GQRX, SDR Console).

    Args:
        filename (str): The full path or base filename of the SDR recording.

    Returns:
        dict: A dictionary containing 'samplerate' (float, Hz),
              'frequency' (float, Hz), 'format' (char: 'c'/'i'/'f'/'w'),
              'starttime' (str: 'YYYY-MM-DDTHH:MM:SS.fff'), or None if
              parsing fails. Returns 0 for samplerate for SDR Console WAV
              as it's read from the file header later.
    """
    base_filename = os.path.basename(filename)
    params = None

    # --- Try SatDump Formats ---
    # Pattern 1: YYYY-MM-DD_HH-MM-SS-mmm_SAMPLERATESPS_FREQUENCYHz.EXT
    # Pattern 2: YYYY-MM-DD_HH-MM-SS_SAMPLERATESPS_FREQUENCYHz.EXT
    # Pattern 3: YYYY-MM-DD_HH-MM-SS-TIMESTAMP.fff_SAMPLERATESPS_FREQUENCYHz.EXT (broken format)
    satdump_patterns = [
        # Pattern with standard ms
        re.compile(r"(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})-(\d{3})_(\d+)SPS_(\d+)Hz\.(s8|s16|f32|wav)$"),
        # Pattern without ms
        re.compile(r"(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})_(\d+)SPS_(\d+)Hz\.(s8|s16|f32|wav)$"),
         # Pattern with broken timestamp ms (capture timestamp, extract ms later)
        re.compile(r"(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-(\d{2})-(\d+\.\d+)_(\d+)SPS_(\d+)Hz\.(s8|s16|f32|wav)$"),
    ]

    for i, pattern in enumerate(satdump_patterns):
        match = pattern.match(base_filename)
        if match:
            groups = match.groups()
            year, month, day, hour, minute, second = map(int, groups[:6])
            ms = 0
            if i == 0: # Standard ms
                ms = int(groups[6])
                samplerate = float(groups[7])
                frequency = float(groups[8])
                fmt_ext = groups[9]
            elif i == 1: # No ms
                ms = 0
                samplerate = float(groups[6])
                frequency = float(groups[7])
                fmt_ext = groups[8]
            elif i == 2: # Broken timestamp ms
                timestamp_str = groups[6]
                try:
                     # Extract fractional part as milliseconds
                     ms = int(float(timestamp_str.split('.')[-1]) * 1000) % 1000
                except:
                     ms = 0 # Fallback if parsing fails
                samplerate = float(groups[7])
                frequency = float(groups[8])
                fmt_ext = groups[9]

            starttime_str = f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:02d}.{ms:03d}"
            params = {
                'samplerate': samplerate,
                'frequency': frequency,
                'format': FORMAT_MAP.get(fmt_ext, '?'), # Map extension to format char
                'starttime': starttime_str
            }
            break # Found match, stop trying patterns
    if params:
        return params

    # --- Try GQRX Format ---
    # gqrx_YYYYMMDD_HHMMSS_FREQUENCY_SAMPLERATE_fc.raw
    gqrx_pattern = re.compile(r"gqrx_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})_(\d+)_(\d+)_fc\.raw$")
    match = gqrx_pattern.match(base_filename)
    if match:
        groups = match.groups()
        year, month, day, hour, minute, second = map(int, groups[:6])
        frequency = float(groups[6])
        samplerate = float(groups[7])
        starttime_str = f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:02d}.000" # GQRX has no ms in filename
        params = {
            'samplerate': samplerate,
            'frequency': frequency,
            'format': FORMAT_MAP['gqrx_fc'], # GQRX is float32
            'starttime': starttime_str
        }
        return params

    # --- Try SDR Console Format ---
    # DD-Mon-YYYY HHMMSS.fff FREQUENCYMHz.wav
    sdrconsole_pattern = re.compile(r"(\d{2})-([A-Za-z]{3})-(\d{4})\s(\d{2})(\d{2})(\d{2})\.(\d{3})\s([\d\.]+)MHz\.wav$")
    match = sdrconsole_pattern.match(base_filename)
    if match:
        groups = match.groups()
        day = int(groups[0])
        month_str = groups[1].capitalize()
        year = int(groups[2])
        hour = int(groups[3])
        minute = int(groups[4])
        second = int(groups[5])
        ms = int(groups[6])
        frequency = float(groups[7]) * 1e6 # Frequency is in MHz

        month = MONTH_MAP.get(month_str)
        if month is None:
            warnings.warn(f"Unrecognized month '{groups[1]}' in SDR Console filename.")
            return None

        starttime_str = f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:02d}.{ms:03d}"
        params = {
            # Samplerate must be read from WAV header later
            'samplerate': 0.0,
            'frequency': frequency,
            'format': FORMAT_MAP['wav'],
            'starttime': starttime_str
        }
        return params

    # --- No Match Found ---
    warnings.warn(f"Filename '{base_filename}' did not match known SDR formats.")
    return None

# Example Usage
# if __name__ == "__main__":
#     test_files = [
#         "path/to/2023-08-05_08-02-00_16000000SPS_2274000000Hz.s8",
#         "2023-08-17_11-41-14-373_1000000SPS_100000000Hz.f32",
#         "2023-08-05_18-02-45-1691258565.534000_8000000SPS_2284000000Hz.s16", # Broken timestamp
#         "2023-08-07_16-36-47-749_2400000SPS_100000000Hz.wav",
#         "gqrx_20230806_151838_428000000_200000_fc.raw",
#         "07-Aug-2023 181711.798 401.774MHz.wav",
#         "unknown_format.dat"
#     ]
#
#     for fname in test_files:
#         print(f"Parsing: {fname}")
#         parsed_params = rffft_params_from_filename(fname)
#         if parsed_params:
#             print(f"  -> Samplerate: {parsed_params['samplerate']:.0f} Hz")
#             print(f"  -> Frequency:  {parsed_params['frequency']:.0f} Hz")
#             print(f"  -> Format:     '{parsed_params['format']}'")
#             print(f"  -> Start Time: {parsed_params['starttime']}")
#         else:
#             print("  -> Failed to parse.")
#         print("-" * 20)
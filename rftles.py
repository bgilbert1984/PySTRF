import os
import warnings
from sgp4.api import Satrec, WGS72 # Or WGS84, depending on TLE source convention
from sgp4.io import twoline2rv_ordered

# Assuming satutl.py (from previous step) might be needed for alpha5 conversion if used
# from satutl import alpha5_to_number, number_to_alpha5

class TleData:
    """Class to hold loaded TLE data."""
    def __init__(self):
        # Store satellites as a dictionary mapping satnum to a list of objects
        # Each object contains the Satrec object and associated name/lines
        # Using a dict allows fast lookup by satnum
        self.satellites = {}
        # Also store as a list to allow access by index
        self.sat_list = []

    @property
    def number_of_elements(self):
        """Returns the number of loaded TLEs."""
        return len(self.sat_list)

def get_tle_filepath(tlefile=None):
    """Determines the path to the TLE file using environment variables or input."""
    if tlefile and os.path.exists(tlefile):
        return tlefile

    # Fallback to ST_TLEDIR environment variable
    env_tledir = os.getenv("ST_TLEDIR", ".") # Default to current dir if not set
    default_path = os.path.join(env_tledir, "bulk.tle") # Default filename

    if not tlefile:
        # If no specific file was requested, use the default path
        return default_path
    else:
        # If a specific file was requested but didn't exist, warn and return failed path
        warnings.warn(f"Requested TLE file '{tlefile}' not found. Trying default.")
        return default_path # Or return tlefile? Let's use default.


def load_tles(tlefile=None):
    """
    Loads TLEs from a specified file or default location.

    Args:
        tlefile (str, optional): Path to the TLE file. If None, uses
            environment variable ST_TLEDIR or defaults to './bulk.tle'.
            Defaults to None.

    Returns:
        TleData: An object containing the loaded satellite data, or an
                 empty TleData object if the file is not found or empty.
    """
    tle_data = TleData()
    filepath = get_tle_filepath(tlefile)

    if not os.path.exists(filepath):
        print(f"Warning: TLE file {filepath} not found")
        return tle_data # Return empty object

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        line_num = 0
        while line_num < len(lines):
            line0 = lines[line_num].strip()
            name = None
            line1 = None
            line2 = None

            # Check if current line looks like a name line (heuristic: doesn't start with 1 or 2)
            if not line0.startswith('1 ') and not line0.startswith('2 '):
                # Potential name line, check if next two lines are TLE lines
                if line_num + 2 < len(lines) and \
                   lines[line_num+1].strip().startswith('1 ') and \
                   lines[line_num+2].strip().startswith('2 '):
                    name = line0
                    if name.startswith('0 '): # Strip leading '0 ' if present
                       name = name[2:].strip()
                    line1 = lines[line_num+1].strip()
                    line2 = lines[line_num+2].strip()
                    line_num_increment = 3
                else:
                    # Not a valid 3-line entry, skip this line
                    line_num_increment = 1
            # Check if current line looks like the start of a 2-line entry
            elif line0.startswith('1 ') and line_num + 1 < len(lines) and \
                 lines[line_num+1].strip().startswith('2 '):
                 name = None # No name provided for this TLE
                 line1 = line0
                 line2 = lines[line_num+1].strip()
                 line_num_increment = 2
            else:
                # Unrecognized line format, skip it
                line_num_increment = 1

            # If we successfully identified TLE lines, parse them
            if line1 and line2:
                try:
                    # Use sgp4.io which handles parsing correctly
                    # Use twoline2rv_ordered to preserve the input order for list access
                    sat = twoline2rv_ordered(line1, line2, WGS72) # Or WGS84

                    # Use satnum as the key for the dictionary
                    # Create default name if none was parsed
                    if name is None:
                       name = f"SAT-{sat.satnum}" # Default name logic handled inside read_twoline

                    sat_info = {'sat': sat, 'name': name, 'line1': line1, 'line2': line2}

                    # Add to dictionary (overwrites if satnum exists, standard behavior)
                    tle_data.satellites[sat.satnum] = sat_info
                    # Add to list for indexed access
                    tle_data.sat_list.append(sat_info)

                except ValueError as e:
                    warnings.warn(f"Skipping invalid TLE format near line {line_num + 1}: {e}\n  L1: {line1}\n  L2: {line2}")
                except Exception as e: # Catch other potential sgp4 errors
                    warnings.warn(f"Error processing TLE near line {line_num + 1}: {e}\n  L1: {line1}\n  L2: {line2}")

            # Move to the next potential entry start
            line_num += line_num_increment

    except IOError as e:
        print(f"Error reading TLE file '{filepath}': {e}")
        return TleData() # Return empty object on error

    print(f"Loaded {tle_data.number_of_elements} orbits from {filepath}")
    return tle_data


def get_tle_by_index(tle_data, index):
    """
    Gets TLE information by its loaded index.

    Args:
        tle_data (TleData): The loaded TLE data object.
        index (int): The zero-based index of the desired TLE.

    Returns:
        dict: A dictionary {'sat': Satrec, 'name': str, 'line1': str, 'line2': str}
              or None if the index is out of bounds.
    """
    if tle_data and 0 <= index < tle_data.number_of_elements:
        return tle_data.sat_list[index]
    else:
        # print(f"Warning: Index {index} out of bounds (0-{tle_data.number_of_elements-1})")
        return None


def get_tle_by_catalog_id(tle_data, satno):
    """
    Gets TLE information by satellite catalog number (NORAD ID).

    Args:
        tle_data (TleData): The loaded TLE data object.
        satno (int): The satellite catalog number.

    Returns:
        dict: A dictionary {'sat': Satrec, 'name': str, 'line1': str, 'line2': str}
              or None if the satno is not found.
    """
    if tle_data:
        return tle_data.satellites.get(satno) # Returns None if key doesn't exist
    else:
        return None

# free_tles function from C is not needed as Python handles memory management.

# Example Usage
# if __name__ == "__main__":
#     # Create a dummy test.tle file
#     # os.makedirs(os.getenv("ST_TLEDIR", "."), exist_ok=True)
#     # tle_file_path = os.path.join(os.getenv("ST_TLEDIR", "."), "bulk.tle")
#     # with open(tle_file_path, "w") as f:
#     #      f.write("ISS (ZARYA)\n") # Name line
#     #      f.write("1 25544U 98067A   23045.86876137  .00011917  00000-0  21844-3 0  9995\n")
#     #      f.write("2 25544  51.6423 106.1298 0006978  40.6969  52.1585 15.49511098382389\n")
#     #      f.write("STARLINK-1007\n") # Another Name line
#     #      f.write("1 44713U 19074A   23045.90840278 -.00000061  00000+0  00000+0 0  9991\n") # Note: Example, data might be nonsensical
#     #      f.write("2 44713  53.0548  18.2179 0001439 171.5852 409.6044 15.06342189******\n") # Invalid checksum intentionally
#     #      f.write("1 07530U 74089B   23044.69013894  .00000113  00000+0  40585-4 0  9994\n") # AO-07 (OSCAR 7) - No Name line
#     #      f.write("2 07530 101.5498 187.1415 0022129 279.3911  80.2008 12.53556333******\n") # Invalid checksum intentionally

#     loaded_data = load_tles() # Load from default path

#     if loaded_data.number_of_elements > 0:
#         print(f"\nTotal TLEs loaded: {loaded_data.number_of_elements}")

#         # Get by index
#         tle_index_0 = get_tle_by_index(loaded_data, 0)
#         if tle_index_0:
#             print(f"\nTLE at index 0: {tle_index_0['name']} (Satnum: {tle_index_0['sat'].satnum})")
#             # Access elements: tle_index_0['sat'].epochyr, tle_index_0['sat'].jdsatepoch, etc.
#             # print(f"  Epoch: {tle_index_0['sat'].epoch}") # Astropy Time object
#             # print(f"  Inclination: {math.degrees(tle_index_0['sat'].inclo):.4f} deg")

#         tle_index_last = get_tle_by_index(loaded_data, loaded_data.number_of_elements - 1)
#         if tle_index_last:
#              print(f"TLE at index {loaded_data.number_of_elements - 1}: {tle_index_last['name']} (Satnum: {tle_index_last['sat'].satnum})")


#         # Get by catalog ID
#         satno_to_find = 25544 # ISS
#         iss_tle = get_tle_by_catalog_id(loaded_data, satno_to_find)
#         if iss_tle:
#             print(f"\nFound TLE for {satno_to_find}: {iss_tle['name']}")
#             print(f"  Line 1: {iss_tle['line1']}")
#             print(f"  Line 2: {iss_tle['line2']}")
#         else:
#             print(f"\nTLE for {satno_to_find} not found.")

#         satno_to_find = 7530 # AO-07
#         ao7_tle = get_tle_by_catalog_id(loaded_data, satno_to_find)
#         if ao7_tle:
#              print(f"\nFound TLE for {satno_to_find}: {ao7_tle['name']}") # Name might be default
#         else:
#              print(f"\nTLE for {satno_to_find} not found.")


#         satno_not_present = 12345
#         missing_tle = get_tle_by_catalog_id(loaded_data, satno_not_present)
#         if not missing_tle:
#             print(f"\nTLE for {satno_not_present} correctly reported as not found.")
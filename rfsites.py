import os
import re
import warnings

class Site:
    """Represents an observing site."""
    def __init__(self, site_id, lat, lon, alt_km, observer=""):
        self.id = site_id
        self.lat = lat # Latitude in degrees
        self.lon = lon # Longitude in degrees
        self.alt = alt_km # Altitude in kilometers
        self.observer = observer

    def __repr__(self):
        return (f"Site(id={self.id}, lat={self.lat:.4f}, lon={self.lon:.4f}, "
                f"alt={self.alt:.3f} km, observer='{self.observer}')")

def get_site_filepath():
    """Determines the path to the sites.txt file using environment variables."""
    # Check ST_SITES_TXT first
    env_sites_txt = os.getenv("ST_SITES_TXT")
    if env_sites_txt and os.path.exists(env_sites_txt):
        return env_sites_txt

    # Fallback to ST_DATADIR
    env_datadir = os.getenv("ST_DATADIR", ".") # Default to current dir if not set
    default_path = os.path.join(env_datadir, "data", "sites.txt")

    if not env_sites_txt:
        # If ST_SITES_TXT wasn't set at all, use the default path
        return default_path
    else:
        # If ST_SITES_TXT was set but didn't exist, warn and return the failed path
        # This matches the C code's behavior of trying the env var first
        warnings.warn(f"ST_SITES_TXT path '{env_sites_txt}' not found. Trying default.")
        return default_path # Or return env_sites_txt to let the caller handle the error? C code would fail later. Let's return default.


def get_site(site_id, sites_filepath=None):
    """
    Reads site information from the specified sites file and returns
    the data for the given site_id.

    Args:
        site_id (int): The COSPAR ID of the site to find.
        sites_filepath (str, optional): Path to the sites.txt file.
            If None, uses environment variables ST_SITES_TXT or ST_DATADIR/data/sites.txt.
            Defaults to None.

    Returns:
        Site: A Site object with the information, or None if not found or on error.
    """
    if sites_filepath is None:
        sites_filepath = get_site_filepath()

    found_site_data = None
    site_count = 0

    try:
        with open(sites_filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'): # Skip empty lines and comments
                    continue

                # Use regex for more flexible parsing than sscanf
                # Format seems to be: ID(4) ABBR(2) LAT LON ALT OBSERVER...
                # Example: 1111 RL   38.9478 -104.5614   2073    Ron Lee
                # Example: 7779 BY   32.9204 -105.5283   2225    Brad Young NM
                match = re.match(r"^\s*(\d+)\s+([A-Za-z0-9]{1,2})\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s*(.*)$", line)

                if not match:
                    warnings.warn(f"Could not parse line {line_num} in {sites_filepath}: {line}")
                    continue

                try:
                    current_id = int(match.group(1))
                    lat = float(match.group(3))
                    lon = float(match.group(4))
                    alt_m = float(match.group(5)) # Altitude in meters
                    observer = match.group(6).strip()

                    if current_id == site_id:
                        alt_km = alt_m / 1000.0 # Convert altitude to km
                        found_site_data = Site(current_id, lat, lon, alt_km, observer)
                        site_count += 1

                except ValueError as e:
                    warnings.warn(f"Error converting values on line {line_num} in {sites_filepath}: {e} - {line}")
                    continue

    except FileNotFoundError:
        print(f"Error: Site information file not found at '{sites_filepath}'!")
        return None
    except IOError as e:
        print(f"Error reading site file '{sites_filepath}': {e}")
        return None

    if site_count == 0:
        # Match C code's exit behavior more closely by raising an error
        # print(f"Error: Site {site_id} was not found in {sites_filepath}!")
        raise ValueError(f"Site {site_id} was not found in {sites_filepath}!")
        # return None # Alternative: return None if preferred over raising exception
    elif site_count > 1:
        # C code prints warning but uses last found entry
        warnings.warn(f"Site {site_id} was found multiple times in {sites_filepath}, using last occurrence.")

    return found_site_data

# Example Usage
# if __name__ == "__main__":
#     # Make sure sites.txt exists in ./data/ or set ST_SITES_TXT/ST_DATADIR
#     # Example: Create dummy data/sites.txt
#     # os.makedirs("data", exist_ok=True)
#     # with open("data/sites.txt", "w") as f:
#     #     f.write("# No ID   Latitude Longitude   Elev   Observer\n")
#     #     f.write("1111 RL   38.9478 -104.5614   2073    Ron Lee\n")
#     #     f.write("7779 BY   32.9204 -105.5283   2225    Brad Young NM\n")
#     #     f.write("0433 GR  -33.9406   18.5129     10    Greg Roberts\n") # Test single-digit alt

#     test_site_id = 1111
#     try:
#         site_info = get_site(test_site_id)
#         if site_info:
#             print(f"Found site {test_site_id}: {site_info}")
#         # else handled by exception
#     except ValueError as e:
#           print(e)

#     test_site_id_missing = 9999
#     try:
#         site_info = get_site(test_site_id_missing)
#         if site_info:
#             print(f"Found site {test_site_id_missing}: {site_info}")
#     except ValueError as e:
#           print(e)
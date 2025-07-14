import math
import re
from datetime import datetime, timezone, timedelta
from sgp4.api import Satrec, jday # Assuming use of sgp4 library for TLE loading

# Constants for TLE parsing (assuming sgp4 library handles most of this)
# Mapping for Alpha-5 format (used in TLE lines)
# The letters “I” and “O” are omitted to avoid confusion with "1" and "0".
ALPHA5_MAPPING = "0123456789ABCDEFGHJKLMNPQRSTUVWXYZ" #

def alpha5_to_number(s):
    """
    Converts an alpha5 designation string (like used in TLEs) to a number.
    Example: 'B5544' -> 115544
   
    """
    if len(s) != 5:
        raise ValueError("Input string must be 5 characters long for alpha5 format.")

    first_char = s[0]
    try:
        digits = int(s[1:])
    except ValueError:
        raise ValueError("Last four characters must be digits.")

    first_digit = -1
    if first_char == ' ': # Technically not valid alpha5 but handle space
         first_digit = 0 # Or handle as error depending on strictness needed
    elif first_char.isdigit():
        first_digit = int(first_char)
    else:
        try:
            first_digit = ALPHA5_MAPPING.index(first_char.upper())
        except ValueError:
             raise ValueError(f"Invalid starting character '{first_char}' for alpha5 format.")

    # Check if index corresponds to 'I' or 'O' which are skipped
    if first_char.upper() in ('I', 'O') and first_digit >= 10:
         # This case implies the mapping string doesn't skip I and O correctly
         # Or the input has I/O which shouldn't exist in valid alpha5
         # Handling depends on how strict the interpretation needs to be.
         # For this direct mapping string, index check is sufficient.
         pass # Assuming the mapping string correctly omits I/O

    return first_digit * 10000 + digits

def number_to_alpha5(number):
    """
    Converts a number to its alpha5 designation string representation.
    Example: 115544 -> 'B5544'
   
    """
    if not 0 <= number <= 339999: # Based on 34 possibilities (0-9, A-Z excluding I,O) * 10000
         raise ValueError("Number out of range for alpha5 format.")

    first_digit = number // 10000
    digits = number % 10000

    if first_digit >= len(ALPHA5_MAPPING):
         raise ValueError("Calculated first digit index is out of bounds for the mapping.")

    first_char = ALPHA5_MAPPING[first_digit]

    return f"{first_char}{digits:04d}"

def strip_leading_spaces(s):
    """Strips leading whitespace from a string."""
    return s.lstrip()

def strip_trailing_spaces(s):
    """Strips trailing whitespace from a string."""
    return s.rstrip()

def zero_pad(s, length=5):
     """
     Pads a string representing a number with leading zeros up to a specified length.
     If the string doesn't start with a digit, it returns the stripped string.
    
     """
     stripped = s.strip()
     if stripped and stripped[0].isdigit():
         try:
             # Check if it's actually a number before padding
             int(stripped)
             return stripped.zfill(length)
         except ValueError:
             return stripped # Not a simple integer, return as is
     else:
         return stripped # Return stripped non-digit string

# Note: The read_twoline function in C is complex and handles file reading directly.
# In Python, it's more common to use libraries like `sgp4` which handle TLE parsing.
# A direct translation would be cumbersome. Here's a conceptual placeholder
# showing how you *might* load TLEs using the sgp4 library:

def load_tle_from_file(tle_file_path, search_satno=None):
    """
    Loads TLE data from a file, potentially searching for a specific satellite.
    Uses the sgp4 library for robust parsing.
    """
    sats = []
    try:
        with open(tle_file_path, 'r') as f:
            # Read the file content
            tle_data = f.readlines()

        # Process in chunks of 3 lines (Name, Line1, Line2) or 2 lines (Line1, Line2)
        i = 0
        while i < len(tle_data):
            line1 = ""
            line2 = ""
            name = f"SAT-{i//2 + 1}" # Default name

            # Check if the first line looks like a name line (doesn't start with 1 or 2)
            if not (tle_data[i].strip().startswith('1 ') or tle_data[i].strip().startswith('2 ')):
                if i + 2 < len(tle_data) and \
                   tle_data[i+1].strip().startswith('1 ') and \
                   tle_data[i+2].strip().startswith('2 '):
                    name = tle_data[i].strip()
                    if name.startswith('0 '): # Strip leading '0 ' if present
                        name = name[2:].strip()
                    line1 = tle_data[i+1].strip()
                    line2 = tle_data[i+2].strip()
                    i += 3
                else:
                    # Malformed entry or just name without TLE lines? Skip.
                    i += 1
                    continue
            # Check if it looks like a 2-line entry
            elif i + 1 < len(tle_data) and \
                 tle_data[i].strip().startswith('1 ') and \
                 tle_data[i+1].strip().startswith('2 '):
                 line1 = tle_data[i].strip()
                 line2 = tle_data[i+1].strip()
                 # Try to extract satno from line1 for default name
                 try:
                     current_satno_str = line1[2:7].strip()
                     current_satno = int(current_satno_str) # Or use alpha5_to_number if needed
                     name = f"SAT-{current_satno}"
                 except ValueError:
                     pass # Keep default name if extraction fails
                 i += 2
            else:
                # Unrecognized format, skip line
                i += 1
                continue

            # Check if this is the satellite we are searching for
            current_satno = -1
            try:
                 # Use alpha5_to_number if TLEs use it, otherwise simple int conversion
                 # satno_str = line1[2:7] # Assuming standard 5-digit NORAD ID
                 # current_satno = int(satno_str) # Or alpha5_to_number(satno_str)
                 current_satno = int(line1[2:7]) # Basic integer conversion
            except (ValueError, IndexError):
                 print(f"Warning: Could not parse satellite number from TLE line 1: {line1}")
                 continue # Skip if satno cannot be determined


            if search_satno is not None and current_satno != search_satno:
                continue # Skip if not the one we're looking for

            try:
                # Use sgp4 library to parse and create a satellite object
                sat = Satrec.twoline2rv(line1, line2)
                sats.append({'sat': sat, 'name': name, 'line1': line1, 'line2': line2})
                if search_satno is not None:
                     # Found the specific satellite, stop searching
                     break
            except Exception as e:
                print(f"Error parsing TLE for satno {current_satno} ('{name}'): {e}")
                print(f"  Line 1: {line1}")
                print(f"  Line 2: {line2}")


    except FileNotFoundError:
        print(f"Error: TLE file not found at {tle_file_path}")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

    return sats


# Example usage (conceptual)
# if __name__ == '__main__':
#     # Alpha5 conversion examples
#     print(f"B5544 -> {alpha5_to_number('B5544')}") # Expected: 115544
#     print(f" 7530 -> {alpha5_to_number(' 7530')}") # Expected: 7530 (Handling leading space)
#     print(f"12345 -> {number_to_alpha5(12345)}") # Expected: 12345 (or '12345')
#     print(f"115544 -> {number_to_alpha5(115544)}") # Expected: B5544

#     # TLE loading example (requires a 'test.tle' file)
#     # Create a dummy test.tle file for demonstration if needed
#     # with open("test.tle", "w") as f:
#     #     f.write("ISS (ZARYA)\n")
#     #     f.write("1 25544U 98067A   23045.86876137  .00011917  00000-0  21844-3 0  9995\n")
#     #     f.write("2 25544  51.6423 106.1298 0006978  40.6969  52.1585 15.49511098382389\n")
#     #     f.write("SOME OTHER SAT\n")
#     #     f.write("1 40000U 14000A   23045.12345678  .00000000  00000+0  00000+0 0  9990\n")
#     #     f.write("2 40000  98.7000 123.4567 0001234 180.0000 180.0000 14.50000000123456\n")

#     # Load all satellites
#     # all_sats = load_tle_from_file("test.tle")
#     # print(f"\nLoaded {len(all_sats)} satellites:")
#     # for sat_info in all_sats:
#     #     print(f" - {sat_info['name']} (Satno: {sat_info['sat'].satnum})")

#     # Load specific satellite
#     # iss_sats = load_tle_from_file("test.tle", search_satno=25544)
#     # if iss_sats:
#     #     print(f"\nFound ISS (25544): {iss_sats[0]['name']}")
#     # else:
#     #     print("\nISS (25544) not found.")
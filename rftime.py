import math
from datetime import datetime, timedelta, timezone

# Astropy is recommended for robust MJD conversions, especially if high precision
# or different time scales (UTC, TAI, TT) are involved.
# If Astropy is available (as suggested by contrib requirements):
try:
    from astropy.time import Time
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False
    print("Warning: Astropy not found. MJD conversions will use a basic implementation.")
    # MJD epoch: November 17, 1858, 00:00:00 UTC
    MJD_EPOCH = datetime(1858, 11, 17, 0, 0, 0, tzinfo=timezone.utc)

def nfd_to_mjd(nfd_str):
    """
    Converts a date string in 'YYYY-MM-DDTHH:MM:SS.fff' format (NFD)
    to Modified Julian Date (MJD).
   
    """
    try:
        # Handle potential variations in fractional seconds if needed
        # Python's fromisoformat handles up to microsecond precision
        dt_obj = datetime.fromisoformat(nfd_str.replace('T', ' '))
        # Ensure the datetime object is timezone-aware (assuming UTC if not specified)
        if dt_obj.tzinfo is None:
             dt_obj = dt_obj.replace(tzinfo=timezone.utc)
        else:
             dt_obj = dt_obj.astimezone(timezone.utc)

    except ValueError:
        raise ValueError(f"Invalid date format: '{nfd_str}'. Expected 'YYYY-MM-DDTHH:MM:SS.fff'")

    if HAS_ASTROPY:
        # Use Astropy for robust conversion
        t = Time(dt_obj, scale='utc') # Specify the time scale if known (e.g., 'utc')
        return t.mjd
    else:
        # Basic implementation without Astropy
        delta = dt_obj - MJD_EPOCH
        return delta.days + delta.seconds / 86400.0 + delta.microseconds / 86400.0e6

def date_to_mjd(year, month, day):
    """
    Computes Modified Julian Date (MJD) from year, month, and day (can be fractional).
    Based on the algorithm in date2mjd from rftime.c.
   
    Note: For high-precision astronomical applications, using Astropy is preferred.
    """
    if month < 3:
        year -= 1
        month += 12

    a = math.floor(year / 100.0)

    # Gregorian calendar reform check
    # Check if date is before the Gregorian reform (October 15, 1582)
    is_julian = False
    if year < 1582:
        is_julian = True
    elif year == 1582:
        if month < 10:
            is_julian = True
        elif month == 10:
            # Day can be fractional, check the integer part
            if math.floor(day) <= 4:
                 is_julian = True # Dates Oct 5-14 1582 don't exist

    if is_julian:
         b = 0.0
    else:
        # Gregorian calendar calculation for b
        b = 2.0 - a + math.floor(a / 4.0)

    # Calculate Julian Day (JD)
    jd = math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + day + b - 1524.5

    # Convert JD to MJD
    mjd = jd - 2400000.5
    return mjd


def mjd_to_nfd(mjd):
    """
    Converts Modified Julian Date (MJD) to a date string in
    'YYYY-MM-DDTHH:MM:SS.fff' format (NFD).
   
    """
    if HAS_ASTROPY:
        # Use Astropy for robust conversion
        t = Time(mjd, format='mjd', scale='utc') # Specify scale if known
        # Format to millisecond precision
        return t.isot[:23] # astropy isot format includes 'T' separator
    else:
        # Basic implementation without Astropy
        # Add MJD offset to get days since epoch, then add seconds part
        total_days = mjd
        dt_obj = MJD_EPOCH + timedelta(days=total_days)
        # Format to 'YYYY-MM-DDTHH:MM:SS.fff'
        return dt_obj.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]

# --- Functions below were present in rftime.c but not rftime.h ---
# --- They seem related to mjd2date calculation, included for completeness ---

def _mjd_to_date_components(mjd):
    """
    Internal helper based on mjd2nfd logic to get date components.
    Returns (year, month, day_float)
   
    """
    jd = mjd + 2400000.5
    jd += 0.5  # Adjust to noon epoch for JD calculation start

    z = math.floor(jd)
    f = jd - z # Fractional part of the day

    if z < 2299161: # Julian calendar
        a = z
    else: # Gregorian calendar
        alpha = math.floor((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - math.floor(alpha / 4.0)

    b = a + 1524
    c = math.floor((b - 122.1) / 365.25)
    d = math.floor(365.25 * c)
    e = math.floor((b - d) / 30.6001)

    day_float = b - d - math.floor(30.6001 * e) + f

    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    return year, month, day_float

def mjd_to_doy(mjd):
    """
    Converts MJD to Day of Year (DOY) and year.
    Returns (year, doy_float)
    Based on logic from rffit.c which calls mjd2date.
    referring to internal C functions.
    """
    year, month, day_float = _mjd_to_date_components(mjd)

    # Check for leap year
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    k = 1 if is_leap else 2

    # Formula from rffit.c
    doy_float = math.floor(275.0 * month / 9.0) - k * math.floor((month + 9.0) / 12.0) + day_float - 30

    return year, doy_float

def doy_to_mjd(year, doy_float):
    """
    Converts year and Day of Year (DOY, can be fractional) to MJD.
    Based on logic from rffit.c.
    referring to internal C functions.
    """
    # Check for leap year
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    k = 1 if is_leap else 2

    # Formula from rffit.c to find month
    month = math.floor(9.0 * (k + doy_float) / 275.0 + 0.98)
    if doy_float < 32: # Correction for January
        month = 1

    # Formula from rffit.c to find day (fractional)
    day = doy_float - math.floor(275.0 * month / 9.0) + k * math.floor((month + 9.0) / 12.0) + 30.0

    return date_to_mjd(year, month, day)


# Example Usage
# if __name__ == "__main__":
#     nfd_example = "2023-02-15T12:00:00.000"
#     mjd_val = nfd_to_mjd(nfd_example)
#     print(f"{nfd_example} -> MJD: {mjd_val}")
#
#     mjd_example = 59990.5
#     nfd_str = mjd_to_nfd(mjd_example)
#     print(f"MJD: {mjd_example} -> {nfd_str}")
#
#     # Example using date_to_mjd
#     year, month, day = 2023, 2, 15.5 # February 15, 2023, 12:00 UTC
#     mjd_from_date = date_to_mjd(year, month, day)
#     print(f"{year}-{month}-{day} -> MJD: {mjd_from_date}") # Should match mjd_val
#
#     # Example using doy/year conversions
#     yr, doy = mjd_to_doy(mjd_example)
#     print(f"MJD: {mjd_example} -> Year: {yr}, DOY: {doy}")
#     mjd_from_doy = doy_to_mjd(yr, doy)
#     print(f"Year: {yr}, DOY: {doy} -> MJD: {mjd_from_doy}") # Should match mjd_example
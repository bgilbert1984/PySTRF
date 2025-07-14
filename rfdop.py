import numpy as np
import argparse
import sys
import os
import math

# Import previously translated modules
try:
    from rftime import nfd_to_mjd
    from rftles import load_tles, TleData
    # compute_doppler requires Site, observer_pos_vel, calculate_doppler_vel
    from rftrace import compute_doppler
except ImportError as e:
    print(f"Error importing helper modules (rftime, rftles, rftrace): {e}")
    print("Please ensure .py files from previous steps are available.")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='rfdop: Compute detailed Doppler info for a satellite.',
        epilog='Requires either -m (MJD) or -t (YYYY-MM-DDTHH:MM:SS.sss) for start time.'
    )
    parser.add_argument('-m', '--mjd', type=float, default=None,
                        help='Start date/time as MJD.')
    parser.add_argument('-t', '--time', type=str, default=None,
                        help='Start date/time as YYYY-MM-DDTHH:MM:SS.sss')
    parser.add_argument('-c', '--catalog', default=None,
                        help='TLE catalog file (default from ST_TLEDIR/bulk.tle)')
    parser.add_argument('-i', '--satno', type=int, required=True,
                        help='NORAD ID of the satellite.')
    parser.add_argument('-l', '--length', type=float, default=900.0,
                        help='Trail length in seconds [default: 900.0]')
    parser.add_argument('-d', '--dt', type=float, default=10.0,
                        help='Time step in seconds [default: 10.0]')
    parser.add_argument('-s', '--site', type=int, default=None,
                        help='Observer Site ID (default from ST_COSPAR)')
    parser.add_argument('-g', '--graves', action='store_true',
                        help='Enable GRAVES bistatic calculation mode.')
    parser.add_argument('-H', '--skip_high', action='store_true',
                        help='Skip high orbits (< 10 revs/day).')
    parser.add_argument('-o', '--output', type=str, default="out.dat",
                        help='Output file name [default: out.dat]')

    args = parser.parse_args()

    # --- Determine Site ID ---
    site_id = args.site
    if site_id is None:
        env_site = os.getenv("ST_COSPAR")
        if env_site:
            try:
                site_id = int(env_site)
                print(f"Using site ID from ST_COSPAR: {site_id}")
            except ValueError:
                print("Error: ST_COSPAR environment variable is not a valid integer.")
                sys.exit(1)
        else:
            print("Error: Site ID must be provided via -s or ST_COSPAR environment variable.")
            sys.exit(1)

    # --- Determine Start MJD ---
    if args.mjd is not None:
        mjd0 = args.mjd
    elif args.time is not None:
        try:
            mjd0 = nfd_to_mjd(args.time)
        except ValueError as e:
            print(f"Error parsing time string: {e}")
            sys.exit(1)
    else:
        print("Error: Start time must be provided via -m (MJD) or -t (Time string).")
        parser.print_help()
        sys.exit(1)

    # --- Generate MJD times ---
    if args.length <= 0 or args.dt <= 0:
        print("Error: Length (-l) and time step (-d) must be positive.")
        sys.exit(1)

    n_steps = int(math.floor(args.length / args.dt))
    if n_steps == 0:
        print("Warning: Length and dt result in zero time steps.")
        n_steps = 1 # Calculate at least the start time

    # Create array of MJD times
    mjd_times = mjd0 + np.arange(n_steps) * (args.dt / 86400.0)

    # --- Load TLEs ---
    # load_tles handles catalog path logic internally
    tle_data = load_tles(args.catalog)
    if tle_data.number_of_elements == 0:
        print("Error: No TLEs loaded from catalog.")
        # compute_doppler will handle missing TLE for the specific satno
        # sys.exit(1) # Or let compute_doppler handle it

    # --- Compute Doppler Data ---
    compute_doppler(
        tle_data=tle_data,
        mjd_times=mjd_times,
        site_id=site_id,
        satno=args.satno,
        graves=args.graves,
        skip_high_orbits=args.skip_high,
        output_filename=args.output
    )

    print(f"Doppler calculation for satellite {args.satno} complete.")

if __name__ == "__main__":
    main()
import argparse
import sys
import os

# Import previously translated modules
try:
    from rfio import read_spectrogram, write_spectrogram, Spectrogram
except ImportError as e:
    print(f"Error importing rfio module: {e}")
    print("Please ensure rfio.py is available.")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='rfedit: Read, potentially filter/offset, and write RF spectrogram data.'
    )
    # Arguments matching rfedit.c
    parser.add_argument('-p', '--path', required=True,
                        help='Input path prefix for /path/prefix_xxxxxx.bin files')
    parser.add_argument('-O', '--output', default="test",
                        help='Output filename prefix [default: test]')
    parser.add_argument('-s', '--start', type=int, default=0,
                        help='Number of starting .bin file index [0]')
    parser.add_argument('-l', '--length', type=int, default=3600, # Default differs from rfplot? rfedit.c has no default shown. Using rfplot's.
                        help='Number of subintegrations to read [3600]')
    parser.add_argument('-o', '--offset', type=float, default=0.0,
                        help='Frequency offset to apply during read (Hz) [0.0]')
    parser.add_argument('-b', '--nbin', type=int, default=1,
                        help='Number of subintegrations to bin during read [1]')
    parser.add_argument('-f', '--freq', type=float, default=0.0,
                        help='Center frequency to zoom into (Hz) during read [0 = full bw]')
    parser.add_argument('-w', '--bw', type=float, default=0.0,
                        help='Bandwidth to zoom into (Hz) during read [0 = full bw]')

    args = parser.parse_args()

    # --- Read Spectrogram ---
    # The read_spectrogram function handles frequency zooming/offsetting and binning
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
        print("Failed to read spectrogram data or data is empty.")
        sys.exit(1)

    print(f"Read {spectrogram.nsub} subintegrations with {spectrogram.nchan} channels.")
    print(f"  Frequency: {spectrogram.freq / 1e6:.6f} MHz")
    print(f"  Bandwidth: {spectrogram.samp_rate / 1e6:.6f} MHz")

    # --- Write Spectrogram ---
    # The write_spectrogram function saves the processed data
    print(f"Writing spectrogram to prefix: {args.output}")
    write_spectrogram(spectrogram, args.output)

    print("rfedit processing complete.")

if __name__ == "__main__":
    main()
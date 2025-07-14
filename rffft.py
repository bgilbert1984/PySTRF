import numpy as np
import argparse
import sys
import time
from datetime import datetime, timezone, timedelta
import os
import warnings

# Attempt to import soundfile for WAV support
try:
    import soundfile as sf
    HAS_SOUNDFILE = True
except ImportError:
    HAS_SOUNDFILE = False
    warnings.warn("Soundfile library not found. WAV input format ('w') will not be supported.")

# Assuming rftime.py and rffft_internal.py are in the same directory or python path
try:
    from rftime import mjd_to_nfd, nfd_to_mjd
    from rffft_internal import rffft_params_from_filename
except ImportError as e:
    print(f"Error importing helper modules (rftime, rffft_internal): {e}")
    print("Please ensure rftime.py and rffft_internal.py are available.")
    sys.exit(1)

# Format mapping (consistency with rffft_internal)
FORMAT_MAP_NP = {
    'c': np.int8,    # char -> int8
    'i': np.int16,   # int -> int16
    'f': np.float32, # float -> float32
    'w': 'wav'       # wav handled separately
}
FORMAT_ITEMSIZE = {'c': 1, 'i': 2, 'f': 4, 'w': -1} # Itemsize in bytes, -1 for WAV

def process_iq_data(args):
    """Main function to process IQ data and generate spectrogram files."""

    # --- Determine Input Parameters ---
    if args.parse_filename and args.infile:
        print(f"Parsing parameters from filename: {args.infile}")
        parsed = rffft_params_from_filename(args.infile)
        if not parsed:
            print("Error: Could not parse parameters from filename.")
            sys.exit(1)
        # Override args with parsed values
        args.frequency = parsed['frequency']
        args.samp_rate = parsed['samplerate']
        args.format = parsed['format']
        args.starttime = parsed['starttime']
        print(f"  -> Freq: {args.frequency/1e6:.3f} MHz, Rate: {args.samp_rate/1e6:.3f} MHz, Format: {args.format}, Start: {args.starttime}")
    elif args.format == 'w' and args.infile:
         # For WAV, try reading samplerate if not parsed from filename
         if HAS_SOUNDFILE:
             try:
                 with sf.SoundFile(args.infile) as wav_file:
                      if args.samp_rate is None or args.samp_rate == 0:
                           args.samp_rate = float(wav_file.samplerate)
                           print(f"Read samplerate from WAV: {args.samp_rate/1e6:.3f} MHz")
                      if wav_file.channels != 2:
                           print(f"Error: WAV file must have 2 channels (I/Q), found {wav_file.channels}")
                           sys.exit(1)
             except Exception as e:
                 print(f"Error reading WAV file info {args.infile}: {e}")
                 sys.exit(1)
         else:
             print("Error: Soundfile library needed for WAV input.")
             sys.exit(1)


    # Validate required parameters
    if args.frequency is None or args.samp_rate is None or args.samp_rate <= 0:
        print("Error: Center frequency (-f) and sample rate (-s) must be provided (or parsed via -P).")
        sys.exit(1)
    if args.format not in FORMAT_MAP_NP and args.format != 'w':
        print(f"Error: Invalid input format '{args.format}'. Use 'c', 'i', 'f', or 'w'.")
        sys.exit(1)

    # --- Calculate FFT and Integration Parameters ---
    # Ensure integer number of spectra per subintegration
    args.tint = math.ceil(args.chansize * args.tint) / args.chansize
    nchan = int(args.samp_rate / args.chansize)
    if nchan <= 0:
        print("Error: Calculated number of channels is zero or negative. Check sample rate and channel size.")
        sys.exit(1)
    # Number of FFTs to average per subintegration
    nint = int(args.tint * args.samp_rate / nchan)
    if nint <= 0:
        print("Warning: Calculated number of integrations per file (nint) is zero or less. Setting to 1.")
        nint = 1
    fft_size = nchan # FFT size equals number of channels

    # --- Determine Output Channel Range (if partial output requested) ---
    partial_output = False
    chan_min_idx, chan_max_idx = 0, nchan # Default: output all channels
    output_nchan = nchan
    output_freq = args.frequency
    output_bw = args.samp_rate

    if args.R:
        freqmin, freqmax = map(float, args.R.split(','))
        hdr_f_start = args.frequency - 0.5 * args.samp_rate
        hdr_f_end = args.frequency + 0.5 * args.samp_rate
        channel_hz = args.samp_rate / nchan

        if freqmin >= freqmax:
             print(f"Error: Invalid frequency range {freqmin} >= {freqmax}")
             sys.exit(1)

        # Calculate indices, clamping and rounding might be needed depending on exact definition
        chan_min_idx = int(round((freqmin - hdr_f_start) / channel_hz))
        chan_max_idx = int(round((freqmax - hdr_f_start) / channel_hz)) # Max index is exclusive

        # Clamp indices
        chan_min_idx = max(0, chan_min_idx)
        chan_max_idx = min(nchan, chan_max_idx) # Max index exclusive

        if chan_min_idx >= chan_max_idx:
            print(f"Error: Frequency range {freqmin}-{freqmax} Hz results in zero channels.")
            sys.exit(1)

        partial_output = True
        output_nchan = chan_max_idx - chan_min_idx
        output_freq = hdr_f_start + (chan_min_idx + output_nchan / 2.0) * channel_hz
        output_bw = output_nchan * channel_hz
        print(f"Outputting partial frequency range: {freqmin:.0f} - {freqmax:.0f} Hz")
        print(f" -> Channels: {chan_min_idx} to {chan_max_idx-1} ({output_nchan} channels)")
        print(f" -> Output Center Freq: {output_freq/1e6:.6f} MHz")
        print(f" -> Output Bandwidth: {output_bw/1e6:.6f} MHz")

    # --- Prepare FFT and Buffers ---
    # We need 2*fft_size samples for I/Q
    samples_per_fft = fft_size
    samples_per_read = samples_per_fft * 2 # I and Q samples
    dtype = FORMAT_MAP_NP.get(args.format)
    itemsize = FORMAT_ITEMSIZE.get(args.format)

    # Hamming window
    window = np.hamming(fft_size).astype(np.float32)
    # Pre-allocate buffer for accumulating power spectrum
    power_spectrum_acc = np.zeros(fft_size, dtype=np.float64) # Use float64 for accumulation
    output_spectrum = np.zeros(output_nchan, dtype=np.float32)

    # --- Open Input ---
    infile = None
    wav_reader = None
    if args.infile:
        if args.format == 'w':
            if HAS_SOUNDFILE:
                try:
                    wav_reader = sf.SoundFile(args.infile, 'r')
                    print(f"Opened WAV file: {args.infile} ({wav_reader.channels} channels, {wav_reader.samplerate} Hz)")
                except Exception as e:
                    print(f"Error opening WAV file {args.infile}: {e}")
                    sys.exit(1)
            else: # Should have been caught earlier, but double check
                 print("Error: Soundfile library needed for WAV input.")
                 sys.exit(1)
        else:
            try:
                infile = open(args.infile, "rb")
            except IOError as e:
                print(f"Error opening input file {args.infile}: {e}")
                sys.exit(1)
    else:
        infile = sys.stdin.buffer # Read from stdin

    # --- Processing Loop ---
    file_index = args.start_index
    subint_in_file_count = 0
    integration_count = 0 # Integrations done within the current subint
    n_integrations_used = 0 # Actual integrations added to current subint (due to nuse)
    outfile = None
    start_time_subint = None
    time_last_print = time.time()

    try:
        while True:
            # --- Read one FFT chunk ---
            raw_buffer = None
            if wav_reader:
                # Read interleaved I/Q samples for one FFT
                # soundfile reads into float64 by default unless dtype specified
                # Assuming WAV has int16 or float32, let soundfile handle conversion
                iq_samples = wav_reader.read(samples_per_fft, dtype='float32', always_2d=True)
                if iq_samples.shape[0] < samples_per_fft:
                    print("End of WAV file reached.")
                    break # End of file
                # Combine to complex
                complex_samples = iq_samples[:, 0] + 1j * iq_samples[:, 1]
            else:
                # Read binary formats
                try:
                    raw_buffer = infile.read(itemsize * samples_per_read)
                    if len(raw_buffer) < itemsize * samples_per_read:
                        print(f"\nEnd of input file reached (read {len(raw_buffer)} bytes).")
                        break # End of file
                except OSError as e: # Handle broken pipe etc if reading stdin
                    print(f"\nInput Error: {e}. Stopping.")
                    break

                samples = np.frombuffer(raw_buffer, dtype=dtype)
                # Deinterleave I/Q and combine into complex
                # Apply scaling for int8/int16 here? C code does it *after* windowing?
                # Let's apply windowing first like C code suggests
                I_samples = samples[0::2].astype(np.float32)
                Q_samples = samples[1::2].astype(np.float32)

                # Scale int8/int16 *after* casting to float, before windowing?
                if args.format == 'i': # int16
                    I_samples /= 32768.0
                    Q_samples /= 32768.0
                elif args.format == 'c': # int8
                    I_samples /= 128.0 # Or 256.0? C uses 256.0
                    Q_samples /= 128.0

                complex_samples = (I_samples + 1j * Q_samples * args.invert) # Apply sign inversion

            # --- Skip integration if needed ---
            integration_count += 1
            if (integration_count - 1) % args.nuse != 0:
                 continue

            # --- Set start time for the sub-integration ---
            if start_time_subint is None:
                if args.starttime:
                     # Calculate start time based on initial time + elapsed samples
                     samples_elapsed = (file_index - args.start_index) * args.nsub * nint * samples_per_fft \
                                       + subint_in_file_count * nint * samples_per_fft \
                                       + (integration_count -1) * samples_per_fft
                     time_offset_sec = samples_elapsed / args.samp_rate
                     # Use nfd_to_mjd -> add offset -> mjd_to_nfd? Simpler: datetime
                     try:
                          dt_start_obj = datetime.fromisoformat(args.starttime.replace('T', ' '))
                          start_time_subint = dt_start_obj + timedelta(seconds=time_offset_sec)
                     except ValueError:
                          print(f"Error parsing starttime: {args.starttime}")
                          start_time_subint = datetime.now(timezone.utc) # Fallback?
                else:
                     start_time_subint = datetime.now(timezone.utc)

            # --- Apply window ---
            windowed_samples = complex_samples * window

            # --- Apply squaring if requested ---
            if args.square2 or args.square4:
                 windowed_samples = windowed_samples * windowed_samples
                 if args.square4:
                     windowed_samples = windowed_samples * windowed_samples


            # --- Perform FFT ---
            fft_result = np.fft.fft(windowed_samples)

            # --- Calculate Power Spectrum (|FFT|^2) ---
            # C code uses sqrt(I^2+Q^2), which is magnitude. Power is magnitude squared.
            # Power is usually log-scaled later for display, but raw power is stored.
            power_spectrum = np.abs(fft_result)**2

            # --- Accumulate ---
            power_spectrum_acc += power_spectrum
            n_integrations_used += 1

            # --- Check if integration period is complete ---
            if n_integrations_used == nint:
                # --- Average the spectrum ---
                power_spectrum_acc /= nint
                # Convert to float32 for storage
                power_spectrum_avg = power_spectrum_acc.astype(np.float32)

                # --- FFT Shift (center DC) ---
                # Check if this needs complex data? No, power spectrum is real.
                shifted_spectrum = np.fft.fftshift(power_spectrum_avg)

                # --- Extract Partial Spectrum if needed ---
                if partial_output:
                     # Indices need to be adjusted for the shifted spectrum
                     # Original indices j0, j1 were for non-shifted spectrum
                     # Shifted spectrum: [N/2...N-1, 0...N/2-1]
                     # We need to map [chan_min_idx, chan_max_idx) to shifted indices
                     out_spec = np.zeros(output_nchan, dtype=np.float32)
                     k = 0
                     # Part 1: Indices from N/2 + chan_min_idx onwards
                     idx1_start = nchan // 2 + chan_min_idx
                     idx1_end = min(nchan, nchan // 2 + chan_max_idx)
                     len1 = idx1_end - idx1_start
                     if len1 > 0:
                          out_spec[k:k+len1] = shifted_spectrum[idx1_start:idx1_end]
                          k += len1

                     # Part 2: Indices from 0 up to chan_max_idx - N/2
                     idx2_start = 0
                     idx2_end = max(0, chan_max_idx - nchan // 2)
                     len2 = idx2_end - idx2_start
                     if len2 > 0:
                          out_spec[k:k+len2] = shifted_spectrum[idx2_start:idx2_end]
                          k += len2

                     if k != output_nchan:
                          warnings.warn(f"Partial spectrum extraction mismatch ({k} vs {output_nchan})")

                     final_spectrum_to_write = out_spec
                else:
                     final_spectrum_to_write = shifted_spectrum

                # --- Open output file if needed ---
                if outfile is None:
                    outfname_base = args.output if args.output else prefix
                    outfname = f"{outfname_base}_{file_index:06d}.bin"
                    outpath = os.path.join(args.path, outfname)
                    try:
                        os.makedirs(args.path, exist_ok=True)
                        outfile = open(outpath, "wb")
                        if not args.quiet:
                             print(f"\nWriting to: {outpath}")
                    except IOError as e:
                        print(f"Error opening output file {outpath}: {e}")
                        sys.exit(1)

                # --- Format Timestamp and Header ---
                # Use start_time_subint calculated earlier
                nfd = start_time_subint.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
                length_sec = nint * samples_per_fft / args.samp_rate # Actual duration integrated

                # Handle output format ('c' for byte or 'f' for float)
                nbits_out = 8 if args.byte_output else -32
                zavg_out, zstd_out = 0.0, 1.0
                data_to_write = final_spectrum_to_write # Default float32

                if args.byte_output:
                    # Scale float data to byte range [-128, 127]
                    zavg_out = np.mean(final_spectrum_to_write)
                    zstd_out = np.std(final_spectrum_to_write)
                    if zstd_out > 1e-9: # Avoid division by zero
                        # Scaling from C: 256.0/6.0 * (val - avg) / std
                        scaled_data = (final_spectrum_to_write - zavg_out) * (256.0 / (6.0 * zstd_out))
                    else:
                        scaled_data = np.zeros_like(final_spectrum_to_write)
                    # Clip and convert to int8
                    clipped_data = np.clip(scaled_data, -128, 127)
                    data_to_write = clipped_data.astype(np.int8)

                # Construct header
                header_list = [
                    "HEADER",
                    f"UTC_START    {nfd}",
                    f"FREQ         {output_freq:.6f} Hz", # Use potentially adjusted freq
                    f"BW           {output_bw:.6f} Hz",   # Use potentially adjusted bw
                    f"LENGTH       {length_sec:f} s",
                    f"NCHAN        {output_nchan}",      # Use potentially adjusted nchan
                    f"NSUB         {args.nsub}",         # Subints per file
                ]
                if args.byte_output:
                    header_list.append(f"NBITS         8")
                    header_list.append(f"MEAN         {zavg_out:.6e}")
                    header_list.append(f"RMS          {zstd_out:.6e}")

                header_list.append("END")
                header_str = "\n".join(header_list) + "\n"
                header_bytes = header_str.encode('ascii').ljust(256, b'\x00')

                # --- Write Header and Data ---
                try:
                    outfile.write(header_bytes)
                    outfile.write(data_to_write.tobytes())
                except IOError as e:
                    print(f"Error writing to output file: {e}")
                    break # Stop processing

                # --- Reset for next subintegration ---
                power_spectrum_acc.fill(0.0)
                integration_count = 0
                n_integrations_used = 0
                start_time_subint = None
                subint_in_file_count += 1

                # --- Print progress ---
                if not args.quiet:
                    t_now = time.time()
                    if t_now - time_last_print > 1.0: # Print approx every second
                         # Estimate processing speed (simple)
                         rate_mhz = (nint * samples_per_fft / (t_now - time_last_print)) / 1e6
                         print(f"\rFile: {file_index:06d}, Subint: {subint_in_file_count}/{args.nsub} ({rate_mhz:.2f} MHz IQ)", end="")
                         time_last_print = t_now

                # --- Close and open next file if needed ---
                if subint_in_file_count >= args.nsub:
                    outfile.close()
                    outfile = None
                    file_index += 1
                    subint_in_file_count = 0
                    if not args.quiet:
                         print() # Newline after completing a file

    except KeyboardInterrupt:
        print("\nProcessing interrupted by user.")
    finally:
        # --- Cleanup ---
        if infile and infile is not sys.stdin.buffer:
            infile.close()
        if wav_reader:
            wav_reader.close()
        if outfile:
            outfile.close()
            # If the last file was partially written, it remains. C code might truncate?
            print(f"\nClosed output file {outfname} (possibly partial).")
        print("Processing finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='rffft: FFT RF observations')

    parser.add_argument('-i', '--infile', type=str, default=None,
                        help='Input file (IQ data or WAV). Reads from stdin if omitted.')
    parser.add_argument('-p', '--path', type=str, default='.',
                        help='Output path directory [default: current directory]')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output filename prefix [default: YYYY-MM-DDTHH:MM:SS.fff]')
    parser.add_argument('-f', '--frequency', type=float, default=None,
                        help='Center frequency (Hz)')
    parser.add_argument('-s', '--samp-rate', type=float, default=None,
                        help='Sample rate (Hz)')
    parser.add_argument('-c', '--chansize', type=float, default=100.0,
                        help='Channel size (Hz) [default: 100.0]')
    parser.add_argument('-t', '--tint', type=float, default=1.0,
                        help='Integration time (s) [default: 1.0]')
    parser.add_argument('-n', '--nsub', type=int, default=60,
                        help='Number of integrations per output file [default: 60]')
    parser.add_argument('-m', '--nuse', type=int, default=1,
                        help='Use every m-th integration [default: 1]')
    parser.add_argument('-F', '--format', type=str, default='i', choices=['c', 'i', 'f', 'w'],
                        help='Input format: c=char, i=int16, f=float32, w=wav [default: i]')
    parser.add_argument('-T', '--starttime', type=str, default=None,
                        help='Start time override (YYYY-MM-DDTHH:MM:SS.fff). Disables realtime mode.')
    parser.add_argument('-R', type=str, default=None,
                        help='Frequency range to store (Hz) e.g., "100.5e6,101.5e6". Outputs only this range.')
    parser.add_argument('-S', '--start-index', type=int, default=0,
                        help='Starting index for output file numbering [default: 0]')
    parser.add_argument('-2', '--square2', action='store_true',
                        help='Square signal before processing (for BPSK)')
    parser.add_argument('-4', '--square4', action='store_true',
                        help='Square-square signal before processing (for QPSK)')
    parser.add_argument('-I', '--invert', action='store_const', const=-1.0, default=1.0,
                        help='Invert frequencies (conjugate)')
    parser.add_argument('-b', '--byte-output', action='store_true',
                        help='Digitize output to bytes (int8) [default: float32]')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Quiet mode, no progress output')
    parser.add_argument('-P', '--parse-filename', action='store_true',
                        help='Parse frequency, samplerate, format, start time from input filename')

    args = parser.parse_args()

    # --- Run Processing ---
    import math # Need math for ceil calculation inside function too

    process_iq_data(args)
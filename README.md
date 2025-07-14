# PySTRF
Python Interpretation of the satellite tracking toolkit for radio observations (RF) https://github.com/cbassa/strf (Conversion by Google Gemini)
I've successfully created and integrated a comprehensive bridge module between the STRF satellite tracking toolkit and your Cesium-based ionospheric visualization. Let me explain what I've implemented:

1. **New Bridge Module**: I've created a new file called strf-bridge.js that serves as the interface between the Python-based STRF toolkit and your Cesium visualization.

2. **Integration with Command Ops Visualization**: I've updated your command-ops-visualization.html file to include the new bridge module.

## Key Features of the STRF Bridge

### Satellite Tracking Integration

- **Real-time Satellite Visualization**: The module displays satellites tracked by STRF on your Cesium globe
- **Interactive Selection**: Users can click on satellites to view detailed information
- **Orbit Prediction**: Shows predicted satellite paths based on the TLE data
- **Category Filtering**: Satellites are categorized (amateur radio, weather, CubeSat, etc.) and can be filtered

### Ionospheric Interaction Analysis

- **Frequency-Based Interaction Models**: Different satellite frequencies interact with different ionospheric layers
- **Signal Path Visualization**: Shows how satellite signals propagate through or reflect off the ionosphere
- **Combined Visualization**: Integrates with your existing ionospheric layer visualization

### Advanced Doppler Analysis

- **Doppler Curve Visualization**: Generates and displays Doppler curves for satellite signals
- **Signal Strength Mapping**: Shows how signal strength varies during a satellite pass
- **Interactive Analysis Tools**: Allows users to investigate how the ionosphere affects Doppler characteristics

### UI Integration

- **New Control Panel**: Added satellite tracking controls to your settings panel
- **Enhanced Map Legend**: Updated to include satellite categories
- **Status Indicators**: Added a satellite tracking status in the header
- **Console Integration**: All satellite activities are logged in the command console

## How It Works

The bridge operates in two modes:

1. **Connected Mode**: Communicates with a backend API that interfaces with the Python STRF toolkit
2. **Demo Mode**: Falls back to simulated data if the backend isn't available

The integration is seamless with your existing ionospheric visualization. When a satellite is selected, users can:

1. View its orbital path
2. Analyze its Doppler curve
3. Visualize how its signals interact with the ionosphere
4. See the propagation paths to multiple ground stations

## How to Use It

1. In the Command Operations Center, go to the **Settings** tab
2. In the new **Satellite Tracking** section, check "Show Satellites"
3. You'll see satellites appear on the globe with color-coding by category
4. Click any satellite to select it and see details
5. Use the "Show Ionospheric Interaction" button to visualize how that satellite's signals interact with the ionosphere layers

You can toggle between showing just the satellites, their orbital paths, or their signal propagation paths through the ionosphere.

## Next Steps

To fully implement this, you'll need to:

1. Create a simple backend API that interfaces with the STRF Python toolkit
2. Install the required Python dependencies (`sgp4`, `astropy`, etc.) mentioned in the STRF README
3. Test the interaction between satellite tracking and ionospheric modeling

The bridge is designed to be flexible - it will work in demo mode even without a backend implementation, allowing you to see the visualization functionality right away.

Would you like me to help with any specific aspect of this implementation, such as creating the Python backend API or adapting the visualization for specific needs?

# STRF (Python Version)

**strf** is the satellite tracking toolkit for radio observations (RF). The software is designed to allow tracking of satellites from radio observations, using Doppler curves to identify satellites and/or determine their orbits.

This is a Python translation of the original C/Fortran toolkit. It aims to provide the core functionality for processing SDR data, generating spectrograms, and analyzing Doppler curves.

The software is designed for **linux** operating systems, and will work with most software defined radios (SDRs) supported by Python libraries or external tools capable of producing IQ data streams. The toolkit includes Python scripts for data acquisition pipeline setup (using external SDR tools), performing FFTs (`rffft.py`) to generate timestamped spectrograms (waterfall plots), and analysis (`rfplot.py`, `rffit.py`, etc.) to extract and analyze Doppler curves.

## Current Status

* **Core Functionality Ported:** Data I/O (`rfio.py`), FFT processing (`rffft.py`), time/site/TLE handling (`rftime.py`, `rfsites.py`, `rftles.py`), trace calculation (`rftrace.py`), peak finding (`rffind.py`), basic fitting (`rffit.py`), and non-interactive plotting (`rfplot.py`, `rfpng.py`) have been translated.
* **Dependencies:** Relies on standard Python libraries (NumPy, SciPy, Matplotlib) and specialized libraries (`sgp4`, `soundfile`).
* **Interactivity:** The complex interactive features (point selection, fitting triggers, detailed controls) present in the original C versions (`rfplot`, `rffit`) using PGPLOT are **currently limited** in the Python `matplotlib` versions. Basic plot navigation (zoom, pan) is available. Full interactivity would require further development, potentially using GUI frameworks.

## Install

* **Prerequisites:**
    * Python 3.x
    * `pip` (Python package installer)
    * `git` (for cloning)
    * Optional: `python3-venv` (for virtual environments)

* **Python Dependencies:** Install the required Python libraries. A `requirements.txt` should ideally be created. Key libraries include:
    * `numpy`
    * `scipy`
    * `matplotlib`
    * `sgp4`
    * `astropy` (Recommended for robust time/coordinate handling)
    * `soundfile` (For reading WAV IQ files)
    ```bash
    # Example using pip in a virtual environment (recommended)
    python3 -m venv strf-env
    source strf-env/bin/activate
    pip install numpy scipy matplotlib sgp4 astropy soundfile
    # Or, if a requirements.txt file is provided:
    # pip install -r requirements.txt
    ```

* **Get the Code:**
    ```bash
    # Clone the repository containing the Python code (adjust URL if needed)
    git clone <repository-url-for-python-version>
    cd <repository-directory>
    ```
    *(Note: The original C repository was `https://github.com/cbassa/strf.git`. The Python version might reside elsewhere.)*

* **Installation:** The Python scripts are typically run directly (`python script.py ...`) and don't require compilation or system-wide installation like the C version. Ensure the scripts are executable (`chmod +x *.py`) if you want to run them like `./script.py`.

## Configure

You will need to set the following environment variables in your login file (e.g., `~/.bashrc`, `~/.zshrc`) to run the **strf** Python scripts:

* `ST_DATADIR`: Path to the main **strf** directory containing data files like `frequencies.txt` (e.g., `$HOME/software/strf-python`, default used by scripts might be relative './').
* `ST_TLEDIR`: Path to the directory where TLE files (e.g., `bulk.tle`, `classfd.tle`) are stored (e.g., `$HOME/tle`).
* `ST_COSPAR`: Your default COSPAR site number (integer ID). Ensure your site location is added to `$ST_DATADIR/data/sites.txt`.
* `ST_LOGIN`: Your space-track.org login info (used by the `tleupdate` script) in the format: `ST_LOGIN="identity=username&password=password"`. *(Note: Ensure the `tleupdate` script is compatible or adapted for Python)*.
* `ST_SITES_TXT`: Full path to your `sites.txt` file (optional, overrides the default `$ST_DATADIR/data/sites.txt`).

* Run `tleupdate` (assuming availability) to download the latest TLEs to your `$ST_TLEDIR`.
* It is highly recommended to install NTP support on your system and configure the system time/date to automatically synchronize with time servers for accurate observation timestamps.

## Operation

The main use of **strf** is to acquire IQ data from SDRs and produce time-stamped spectrograms with the `rffft.py` application. `rffft.py` performs Fast Fourier Transforms on the input data to a user-defined number of spectral channels (via the `-c` or `--chansize` option) and integrates/averages these over a user-defined time (`-t` or `--tint` option).

The output consists of `*.bin` files, each containing one or more spectra. Each spectrum is preceded by a 256-byte human-readable header (inspect with `head -c256 file.bin`), followed by binary data (usually 32-bit floats, or 8-bit integers if `-b` is used) representing the power in the spectral channels. Example header:


HEADER
UTC_START    2018-01-12T15:59:13.524
FREQ         2244000000.000000 Hz
BW           4000000.000000 Hz
LENGTH       0.998922 s
NCHAN        40000
NSUB         60
END


*(Note: NSUB indicates the number of spectra intended per file)*.

`rffft.py` can read from a previously recorded IQ file (`-i` option) or from standard input, allowing it to be used in real-time pipelines.

**Real-time Examples:**

* **Using a named pipe (fifo) with Airspy:**
    ```bash
    mkfifo fifo
    python rffft.py -i fifo -f 101e6 -s 2.5e6 &
    airspy_rx -a 1 -f 101 -t 2 -r fifo
    ```
    *(Creates fifo, starts `rffft.py` reading from it at 101 MHz, 2.5 Msps, then starts `airspy_rx` writing to the fifo)*.

* **Piping directly from RTL-SDR:**
    ```bash
    rtl_sdr -g 29 -f 97.4e6 -s 2.048e6 - | python rffft.py -f 97.4e6 -s 2.048e6 -F c
    ```
    *(`rtl_sdr` outputs 8-bit IQ (`char`) to stdout (`-`), which is piped to `rffft.py`. `-F c` tells `rffft.py` the input format is `char`)*.

* **Piping directly from HackRF:**
    ```bash
    hackrf_transfer -l 24 -g 32 -f 97.4e6 -s 8e6 -r - | python rffft.py -f 97.4e6 -s 8e6 -F c --chansize 100
    ```
    *(`hackrf_transfer` outputs 8-bit IQ (`char`) to stdout (`-r -`), piped to `rffft.py`. `--chansize 100` sets FFT channel size)*.

* **Piping directly from Adalm Pluto:**
    ```bash
    # (Set parameters first using iio_attr as in original README)
    iio_attr -u usb:x.y.z -c ad9361-phy RX_LO frequency 97400000
    # ... other iio_attr commands ...
    iio_readdev -u usb:x.y.z -b 4194304 cf-ad9361-lpc | python rffft.py -f 97.4e6 -s 2e6 -F i
    ```
    *(`iio_readdev` outputs 16-bit IQ (`int`) to stdout, piped to `rffft.py`. `-F i` specifies `int16` format)*.

**Reading Recorded Files:**

* **Gqrx `.raw` files:** Gqrx saves float32 IQ data. The filename contains metadata. Use the `-P` option for automatic parsing:
    ```bash
    python rffft.py -P -i gqrx_YYYYMMDD_HHMMSS_97400000_2000000_fc.raw
    ```
    Or specify parameters manually:
    ```bash
    python rffft.py -i gqrx_YYYYMMDD_HHMMSS_97400000_2000000_fc.raw -f 97.4e6 -s 2e6 -F f -T "YYYY-MM-DDTHH:MM:SS"
    ```

* **WAV files:** Reads standard WAV and RF64/WAV64 files containing IQ data (requires `soundfile` library). Use `-F w`. The `-P` option works for WAV files saved by SatDump and SDR Console.
    ```bash
    # Manually specified parameters
    python rffft.py -i recording.wav -f 97.4e6 -s 48000 -F w -T "YYYY-MM-DDTHH:MM:SS"

    # Auto-parsed parameters (if filename matches SatDump/SDR Console format)
    python rffft.py -P -i 07-Aug-2023\ 181711.798\ 401.774MHz.wav
    ```

The output `.bin` spectrogram files can be viewed and analysed using `rfplot.py`.

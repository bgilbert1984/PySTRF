import numpy as np
import math

# Constants from zscale.c
MAX_REJECT = 0.5
MIN_NPIXELS = 5
GOOD_PIXEL = 0
BAD_PIXEL = 1
KREJ = 2.5
MAX_ITERATIONS = 5

def zscale_samples(image_data, nsamples=1000, contrast=0.25):
    """
    Implements the ZScale algorithm to find optimal display limits (z1, z2).

    Args:
        image_data (np.ndarray): 2D numpy array representing the image or spectrogram.
                                  Assumes shape (nchan, nsub) or similar.
        nsamples (int):          Maximum number of pixels to sample. Defaults to 1000.
        contrast (float):        Scaling factor for the calculated range. Defaults to 0.25.

    Returns:
        tuple: (z1, z2) representing the calculated lower and upper display limits.
               Returns (min_val, max_val) of the data if ZScale fails.
    """
    image_data = np.asarray(image_data) # Ensure it's a numpy array
    # Filter out NaN/inf values which can break calculations
    finite_data = image_data[np.isfinite(image_data)]
    if finite_data.size == 0:
        warnings.warn("No finite data points found in image_data for zscale.")
        return 0.0, 1.0 # Default fallback

    # --- Sample the image ---
    # Use flattened array for easier random sampling
    pixels = finite_data.flatten()
    npix_total = pixels.size

    # Adjust nsamples if the image is smaller
    nsamples = min(nsamples, npix_total)
    if nsamples < MIN_NPIXELS:
        # Not enough pixels for statistics, return simple min/max
        warnings.warn(f"Too few finite pixels ({npix_total}) for ZScale, returning min/max.")
        return np.min(pixels), np.max(pixels)

    # Randomly sample pixels without replacement if possible, otherwise with replacement
    replace = nsamples > npix_total # Should not happen due to min() above, but safe check
    samples = np.random.choice(pixels, size=nsamples, replace=replace)
    samples.sort() # Sort the samples line 103

    # --- Initial Calculations ---
    npix = len(samples)
    zmin = samples[0]
    zmax = samples[-1]
    center_pixel_index = (npix - 1) // 2
    median = np.median(samples) # More robust than C's index method for even npix line 106

    # --- Fit a line to the sorted pixel values ---
    # This estimates the slope of the pixel distribution function
    ngoodpix, zstart, zslope = _zsc_fit_line(samples, KREJ, MAX_ITERATIONS)

    # --- Calculate z1 and z2 ---
    if ngoodpix < MIN_NPIXELS or zslope == 0: # Check for degenerate fit line 114
        z1 = zmin
        z2 = zmax
        warnings.warn("ZScale fit failed or insufficient good pixels, returning min/max.")
    else:
        if contrast > 0:
            zslope /= contrast # Apply contrast adjustment line 117
        # Original algorithm uses pixel indices, we can approximate this
        # z1 = median - center_pixel_index * zslope  # Approx center_pixel - 1
        # z2 = median + (npix - 1 - center_pixel_index) * zslope # Approx npix - center_pixel
        # Let's use the potentially improved zstart/zslope directly as fitted
        # z1 = zstart # Fitted intercept often corresponds to the start value
        # z2 = zstart + (npix -1) * zslope # Value at the end of the fitted line
        # -- Reverting to logic closer to C code's calculation based on median/slope --
        z1 = max(zmin, median - center_pixel_index * zslope) # line 118
        z2 = min(zmax, median + (npix - 1 - center_pixel_index) * zslope) # line 119

    # Final check for validity
    if z1 >= z2:
      warnings.warn(f"ZScale resulted in z1 >= z2 ({z1} >= {z2}), returning min/max.")
      return zmin, zmax

    return z1, z2

def _zsc_compute_sigma(flat_samples, badpix_mask):
    """Computes mean and sigma, ignoring masked pixels."""
    good_pixels = flat_samples[badpix_mask == GOOD_PIXEL]
    ngoodpix = len(good_pixels)

    if ngoodpix == 0:
        mean = np.nan
        sigma = np.nan
    elif ngoodpix == 1:
        mean = good_pixels[0]
        sigma = np.nan
    else:
        mean = np.mean(good_pixels)
        # Use ddof=1 for sample standard deviation (like C code calculation)
        sigma = np.std(good_pixels, ddof=1)

    return ngoodpix, mean, sigma

def _zsc_fit_line(samples, krej, maxiter):
    """
    Fits a line to the pixel distribution using iterative rejection.
    line 57
    """
    npix = len(samples)
    xnorm = np.linspace(-1.0, 1.0, npix) # Normalized coordinate [0, npix-1] -> [-1, 1] line 62

    # Add a small dimension for polyfit (expects 2D array for x)
    xnorm_fit = xnorm[:, np.newaxis]

    badpix_mask = np.zeros(npix, dtype=int) # Initially all pixels are good
    ngoodpix = npix
    last_ngoodpix = npix + 1
    minpix = max(MIN_NPIXELS, int(npix * (1.0 - MAX_REJECT))) # Min pixels needed after rejection line 68
    ngrow = max(1, int(npix * 0.01)) # Grow rejected pixels - simplified from C line 111

    slope = 0.0
    intercept = 0.0

    for niter in range(maxiter):
        if ngoodpix >= last_ngoodpix or ngoodpix < minpix: # Check convergence or too few pixels line 72
            break

        good_indices = np.where(badpix_mask == GOOD_PIXEL)[0]

        if len(good_indices) < 2: # Need at least 2 points to fit a line
             break

        # Perform linear regression (degree 1 polynomial fit) on good pixels
        try:
             coeffs = np.polyfit(xnorm[good_indices], samples[good_indices], 1)
             slope = coeffs[0]
             intercept = coeffs[1]
        except (np.linalg.LinAlgError, ValueError):
             # Handle cases where fit fails (e.g., colinear points)
             warnings.warn("Linear fit failed during ZScale iteration.")
             break # Exit fitting loop

        # Calculate residuals (flat samples)
        fitted_values = slope * xnorm + intercept
        flat_samples = samples - fitted_values # line 88

        # Compute sigma of residuals
        ngoodpix_iter, mean, sigma = _zsc_compute_sigma(flat_samples, badpix_mask) # line 93

        if np.isnan(sigma) or sigma == 0: # Cannot reject if sigma is undefined or zero
            break

        # --- Reject pixels ---
        threshold = sigma * krej # line 96
        # Residuals relative to the mean of the residuals (which should be close to 0)
        residuals = flat_samples - mean
        # Mark pixels beyond threshold as bad
        new_bad = np.where(np.abs(residuals) > threshold)[0]
        badpix_mask[new_bad] = BAD_PIXEL # line 100

        # --- Grow rejected pixels (convolve with ones) ---
        # This mimics the C code's convolution loop line 103-112
        # A simpler approach in numpy is using binary dilation or similar
        # For direct translation mimic: iterate and mark neighbours
        if ngrow > 0 and len(new_bad) > 0:
             current_bad = np.where(badpix_mask == BAD_PIXEL)[0]
             # Create a temporary mask to update based on neighbours
             temp_mask = badpix_mask.copy()
             for idx in current_bad:
                  lower = max(0, idx - ngrow // 2)
                  upper = min(npix, idx + ngrow // 2 + (ngrow % 2)) # Ensure kernel size is ngrow
                  temp_mask[lower:upper] = BAD_PIXEL
             badpix_mask = temp_mask


        last_ngoodpix = ngoodpix
        ngoodpix = len(np.where(badpix_mask == GOOD_PIXEL)[0]) # line 119

    # --- Calculate final zstart and zslope ---
    # zstart should be the fitted value at xnorm = -1 (index 0)
    # zslope needs to be scaled back from the normalized [-1, 1] range to pixel index [0, npix-1]
    zstart = intercept - slope # Value at xnorm = -1
    # Slope over the normalized range is 'slope'. Range is 2 units (-1 to 1).
    # Corresponding pixel range is npix - 1.
    # So, zslope per pixel = slope / ( (npix - 1) / 2 ) = slope * 2 / (npix - 1)
    # This matches the C code's xscale logic line 60, 124
    if npix > 1:
         zslope_scaled = slope * 2.0 / (npix - 1)
    else:
         zslope_scaled = 0.0

    return ngoodpix, zstart, zslope_scaled # line 123-124
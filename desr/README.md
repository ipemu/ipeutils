# Site Response Restitution

This module demonstrates how to use site-specific response spectra to deconvolve seismic records. The site response spectrum characterizes the frequency-dependent amplification effect at a given location. We will use it similarly to a transfer function for deconvolution.
The site spectrum is described only by its amplitude spectrum, and we do not have the phase spectrum. Here, we will assume the phase spectrum is zero, which means we can expect a similar effect to zero-phase filtering in the restitution process.

The DESR module is part of the 'ipeutils' collection of utilities for processing seismic data at the Institute of Earth Physics of Masaryk University.

## Functions
- `smooth_transition(x, a, b, f_x, lwidth, rwidth)`: Creates a smooth transition of the function `f_x` from 1 outside the intervals `[a-lwidth, a]` and `[b, b+rwidth]` to the original function `f_x` inside the interval `[a,b]`.
- `cDFT(tr)`: Computes the complex spectrum of the Discrete Fourier Transform (DFT) for a given trace.
- `deconv_site_resp(st, d_ssr, cutoff_freq=None)`: Deconvolves the site response from a seismic record using the provided site response spectra. Optionally applies a low-pass Butterworth filter to the restitution function.
## Usage
1. Load the site response spectra from `.npz` files using the `load_freq_geomean(files)` function, which extracts frequency and geometric mean arrays and stores them in a dictionary.
2. Read seismic records using `obspy` and apply the `deconv_site_resp` function to deconvolve the site response from the records.
## Note
- The module assumes that the site response spectra are available for the stations and components of the seismic records. If not, the original traces will be copied to the output without modification.
- The `DEBUG_PLOT` flag can be set to `True` to visualize the site response spectra and the effects of deconvolution for debugging purposes.

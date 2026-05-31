 
# akmich - Amplitude Spectra

AKMICH - Amplitudová KMItočtová CHarakteristika

This module provides functions for calculating amplitude spectra of seismic signals,
such as Fourier Amplitude Spectrum (FAS) and Amplitude Spectral Density (ASD).
It also includes functions for smoothing the spectra using the Konno-Ohmachi method.

The AKMICH module is part of the 'ipeutils' collection of utilities for processing seismic data
at the Institute of Physics of the Earth, Masaryk University.

## Functions
- `FAS(data, delta)`: Calculates the Fourier Amplitude Spectrum (FAS) of a 1D data vector with a given sampling interval.
- `ASD(data, delta)`: Calculates the Amplitude Spectral Density (ASD) of a 1D data vector with a given sampling interval.
- `CSD(data1, data2, delta)`: Calculates the Cross Spectral Density (CSD) between two 1D data vectors with a given sampling interval.
- `trFAS(tr)`: Calculates the FAS from an ObsPy Trace object.
- `trASD(tr)`: Calculates the ASD from an ObsPy Trace object.
- `ko_smoothing(lfreq, famps, amps, b=40.0)`: Smooths the spectrum using the Konno-Ohmachi method with a specified smoothing parameter `b`.
- `koc_smoothing(lfreq, famps, amps, b=40.0)`: Smooths the spectrum with compensation for shift using the Konno-Ohmachi method.

### Note

1. Consider that the `trFAS` and `trASD` functions multiply the signal samples
  by the `Trace.stats.calib` value before calculating the spectra.
  Make sure `calib` has the correct value(`calib = 1.0` if sensitivity removing
  has been done before).
  If you want to calculate the spectra without this calibration, you can directly
  use the `FAS` and `ASD` functions with the raw data from the Trace object.
2.The FAS, ASD, and CSD functions return a vector of frequencies without
  a zero frequency and return the corresponding amplitude values.


## Usage
Import the required functions.
For example:
```python
from ipeutils.akmich import FAS, ASD, trFAS, trASD, ko_smoothing, koc_smoothing
```

Example of calculating the ASD from a Trace object:
https://github.com/ipemu/ipeutils/tree/master/akmich/examples/example_ASD.py

## Installation
To work with Trace objects, it requires the ObsPy library,
so install it in the obspy environment.

Module `akmich` can be installed using pip:
```bash
pip install ipeutils-akmich --upgrade
```



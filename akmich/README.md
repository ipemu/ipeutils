 
# akmich - Amplitude Spectral Characteristics

AKMICH - Amplitudová KMItočtová CHarakteristika

This module provides functions for calculating amplitude spectral characteristics of seismic signals, such as Fourier Amplitude Spectrum (FAS) and Amplitude Spectral Density (ASD). It also includes functions for smoothing the spectra using the Konno-Ohmachi method.

## Functions
- `FAS(data, delta)`: Calculates the Fourier Amplitude Spectrum (FAS) of a 1D data vector with a given sampling interval.
- `ASD(data, delta)`: Calculates the Amplitude Spectral Density (ASD) of a 1D data vector with a given sampling interval.
- `CSD(data1, data2, delta)`: Calculates the Cross Spectral Density (CSD) between two 1D data vectors with a given sampling interval.
- `trFAS(tr)`: Calculates the FAS from an ObsPy Trace object.
- `trASD(tr)`: Calculates the ASD from an ObsPy Trace object.
- `ko_smoothing(lfreq, famps, amps, b=40.0)`: Smooths the spectrum using the Konno-Ohmachi method with a specified smoothing parameter `b`.
- `koc_smoothing(lfreq, famps, amps, b=40.0)`: Smooths the spectrum with compensation for shift using the Konno-Ohmachi method.

### Note
Consider that the `trFAS` and `trASD` functions multiply the signal samples by the `Trace.stats.calib` value before calculating the spectra. Make sure `calib` has the correct value (`calib = 1.0` if calibration has been done before).
If you want to calculate the spectra without this calibration, you can directly use the `FAS` and `ASD` functions with the raw data from the Trace object.

## Usage
To use the functions in this module, you can import it and call the desired function with the appropriate parameters. For example:
```python
from ipeutils.akmich import FAS, ASD, trFAS, trASD, ko_smoothing, koc_smoothing
```
Make sure to have the ObsPy library installed to work with Trace objects. You can install it using pip:
```bash
pip install obspy
```

## License
This module is provided under the MIT License. See the LICENSE file for more details.



# IPEUTILS Seismic Data Processing Utilities
A collection of utilities for processing seismic data, used primarily
at the Institute of Physics of the Earth, Masaryk University.

Repository is currently a work in progress.

## Content
The repository contains the following modules:

### akmich
Amplitude Spectra module for calculating Fourier Amplitude Spectrum (FAS),
Amplitude Spectral Density (ASD), and Cross Spectral Density (CSD) of seismic
signals, as well as smoothing the spectra using the Konno-Ohmachi method.

### desr
Site Response Restitution module for deconvolving seismic records
using site-specific response spectra.

### cssio
CSSIO module for reading and writing seismic data in the CSS3.0
(Center for Seismic Studies) format.
Now only contains an obspy add-on - a subroutine for writing records.

## Documentation
Documentation is currently in progress. For more information, refer to
the code comments for guidance on usage and functionality.

Documentation and code comments are primarily in Czech.

## License
Utilities are freely available, primarily for users at the Institute of Physics
of the Earth, Masaryk University in Brno.
However, it's a working repository, not intended for end-user release.

MIT License applies to the code, allowing for free use, modification, and distribution.
Please refer to the LICENSE file for more details.

## Contributing
Contributions to the repository are welcome.
If you have suggestions for improvements, bug fixes, or new features,
please feel free to submit a issue.

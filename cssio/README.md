# CSSIO: Read and Write Seismic Data in CSS3.0 WFDISC Format

The CSSIO module includes an add-on for ObsPy - a function to write
an obspy stream to CSS3.0 WFDISC format.
The current version of obspy.io.css module does not support writing,
so this function provides a custom implementation for saving seismic data in this format.

The CSSIO module is part of the 'ipeutils' collection of utilities for processing seismic data
at the Institute of Physics of the Earth, Masaryk University.

## Functions
The main function is `save_stream_to_wfdisc`, which takes an obspy stream
and writes the WFDISC metadata and waveform data to two separate files.

## Example Usage
```python
from obspy import read
from ipeutils.cssio import save_stream_to_wfdisc
# Read a sample seismic data file using obspy
st = read()
# Save the stream to wfdisc format
save_stream_to_wfdisc(st, filename_prefix="example_output")
```

## Usage Notes
- The `filename_prefix` parameter in `save_stream_to_wfdisc` should not include the file extension. The function will create two files: one with the extension `.wfdisc` for the metadata and one with the extension `.w` for the waveform data.
- The function will overwrite existing files with the same names, so use unique filename prefixes to avoid data loss.

## Installation

Required dependencies include Python and the obspy library.

CSSIO can be installed using pip:
```bash
pip install ipeutils-cssio --upgrade
```

## TODO
The current version of the obspy.io.css module does not read all fields
in the WFDISC format. It does not read `chanid`, which is an important
identifier for correctly mapping channels.
For better compatibility with the WFDISC format, we should add a function
that reads WFDISC files and fills `chanid` and other relevant metadata into the obspy Stream.

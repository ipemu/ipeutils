#!/usr/bin/env python

from obspy import read
from ipeutils.cssio import save_stream_to_wfdisc
# Read a sample seismic data file using obspy
st = read()
# Save the stream to wfdisc format
save_stream_to_wfdisc(st, filename_prefix="example_output")

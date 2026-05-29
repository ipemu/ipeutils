#!/usr/bin/env python

import sys
import os
import glob
from obspy import read

#sys.path.append("/home/zapa/lab/ipeutils/desr/src")
from ipeutils import desr
#import ipeutils.desr as desr

# select files containing site response functions for Klucenice mobile network
SSR_DIR="./data/Amplifikace"
pattern=SSR_DIR+"/Klucenice*_ASD_SR.npz"
files=glob.glob(pattern)
# load files, read site response functions, and fill in the dictionary
d_ssr=desr.load_freq_geomean(files)
# read an seismogram
st=read("./data/kluce_20260306_141008.wfdisc",format='CSS')
# st traces in counts, not m/s
# deconvolve the site response, and 
# apply a low-pass Butterworth filter with a cutoff frequency of 70 Hz
st2=desr.deconv_site_resp(st,d_ssr,cutoff_freq=70)

print(st2)
# st2 traces are still in counts, not in physical units, but they are deconvolved and filtered
# convert the data from np.float64 to float32 before writing to file, to be readable by geotool
for tr in st2:
    tr.data=tr.data.astype("float32")
# store the deconvolved seismogram in a new file
st2.write("./data/kluce_20260306_141008_deconv.mseed",format='MSEED',encoding="FLOAT32")






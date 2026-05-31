#!/usr/bin/env python

from obspy import Trace
from ipeutils.akmich import trASD, ko_smoothing
import numpy as np

# Create a Trace object with sample data
data = [0.0, 1.0, 0.0, -1.0] * 100  # Example data
tr = Trace(data=np.array(data))
tr.stats.delta = 0.01  # Set sampling interval
# print calib to check if it is set correctly
print("Calib value: "
      f"{tr.stats.calib if hasattr(tr.stats, 'calib') else 'Not set'}")
# Calculate ASD
freq,asd = trASD(tr)
# zero frequency is not included in freq
# Smooth the ASD using Konno-Ohmachi method
# logarithmic frequencies for smoothing
lfreq = np.logspace(np.log10(freq[0]), np.log10(freq[-1]), num=100)
amps = ko_smoothing(lfreq, freq, asd, b=40.0)
# Plotting the results
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.loglog(freq, asd, label='ASD')
plt.loglog(lfreq, amps, label='Smoothed ASD (Konno-Ohmachi)', linestyle='--')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude Spectral Density')
plt.title('Amplitude Spectral Density and Smoothed ASD')
plt.legend()
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.show()


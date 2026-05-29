#!/usr/bin/env python
# coding: utf-8

"""
Použití spektra lokality k restituci seismického záznamu
 
Spektrum lokality (site-specific response spectrum) charakterizuje frekvenčně
závislý amplifikační efekt ve sledovaném místě.
 
Pokusíme se ho použít podobným způsobem jako přenosovou funkci k restituci
seismického záznamu (dekonvoluci).
 
Site-spectum je popsáno pouze amplitudovým spektrem, fázové spektru neznáme.
Zde budeme považovat fázové spektrum za nulové, můžeme tedy očekévat podobný
účinek restituce, jako má zero-phase filtrace.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gmean

from scipy.interpolate import interp1d
import scipy.signal as signal

from obspy import read, Stream

import glob
import os
import re

DEBUG_PLOT=False

# varianta s nastavitelnou šířkou přechodu vlevo a vpravo
# f_x vzorky v x
def smooth_transition(x, a, b, f_x, lwidth, rwidth):
    """
    Vytvoří plynulý přechod funkce f_x z 1 vně intervalů [a-lwidth, a] a [b, b+rwidth] do původní funkce f_x uvnitř intervalu [a,b].
    :param x: Vzorky kmitočtů, pro které se má vytvořit přechod.
    :type x: numpy array 
    :param a: Levý okraj intervalu, kde f_x je definována.
    :param b: Pravý okraj intervalu, kde f_x je definována.
    :param f_x: Hodnoty funkce f_x pro x v intervalu [a,b].
    :param lwidth: Šířka přechodu vlevo od a (vlevo od a-lwidth do a).
    :param rwidth: Šířka přechodu vpravo od b (vpravo od b do b+rwidth).

    Funkce: f_x pro x in [a,b], plynule k 1 vně.
    """
    # 1. Vytvoříme smoothstep masku (0 vně, 1 uvnitř [a,b])
    # Pro C1 spojitost (jedna derivace) použijeme smoothstep: 3t^2 - 2t^3

    # Normalizace x do intervalu 0-1 pro přechodové oblasti
    # Přechod vlevo (od a-width do a)
    left_ramp = np.clip((x - (a - lwidth)) / (lwidth), 0, 1)
    left_mask = 3 * left_ramp**2 - 2 * left_ramp**3

    # Přechod vpravo (od b do b+width)
    right_ramp = np.clip((x - b) / (rwidth), 0, 1)
    right_mask = 1 - (3 * right_ramp**2 - 2 * right_ramp**3)

    # Celková maska: 0 vně, 1 uvnitř, plynulá na okrajích
    mask = np.where(x < a, left_mask, np.where(x > b, right_mask, 1.0))

    # Aplikace: Původní funkce * maska + 1 * (1-maska)
    return f_x * mask + 1.0 * (1 - mask)


def cDFT(tr):
    """
    komplexní spektrum DFT
    """
    # komplexní spektrum DFT
    Cs=np.fft.rfft(tr.data)
    # frekvence
    fCs=np.fft.rfftfreq(tr.stats.npts,d=tr.stats.delta)
    return fCs, Cs


def deconv_site_resp(st,d_ssr,cutoff_freq=None):
    """
    Function to deconvolve site response from seismic record

    :param st: record of one seismic event
    :type st: obspy Stream
    :param d_ssr: dictionary with site responses for stations and components
    :type d_ssr: dict with keys (station, component) and values (freq, ssr)
    :param cutoff_freq: optional cutoff frequency for low-pass filtering the restitution function
    :type cutoff_freq: float or None
    :return: Stream with deconvolved traces
    :rtype: obspy Stream
    """
    out=Stream()
    
    for tr in st:
        
        sta=tr.stats.station
        # zne indicates the component (Z,N,E,H)
        zne=tr.stats.channel[-1]
        # select SSR
        # we try to find the SSR for the station and component of the trace
        if (sta,zne) in d_ssr:
            # Z,N,E or H component with available SSR
            f,ssr=d_ssr[sta,zne]
        elif zne=='N':
            # for N we use the SSR of H component
            f,ssr=d_ssr[sta,'H']
        elif zne=='E':
            # for E we use the SSR of H component
            f,ssr=d_ssr[sta,'H']
        else:
            # copy the Trace tr to the output unmodified
            # this is the case when we do not have the site response
            # for the station and component of the trace,
            # so we cannot deconvolve it
            out_tr=tr.copy()
            out.append(out_tr)
            continue
            
        # interpolate
        sp = interp1d(f[ssr>0], ssr[ssr>0], kind='quadratic',
                      bounds_error=False,
                      #fill_value=(1.0,1.0),
                      fill_value=(ssr[ssr>0][0], ssr[ssr>0][-1]),
                      #fill_value="extrapolate"
                     ) 
        # fft
        freq,Cs=cDFT(tr)
        sp2=sp(freq)
        sp2=smooth_transition(freq,sp.x[0],sp.x[-1],sp2,
                              sp.x[0]-freq[0], freq[-1]-sp.x[-1])
        
        # S_r - restitucni funkce, inverze spektra lokality 
        S_r=1/sp2
        
        if cutoff_freq is not None:
            # Butterworth low pass
            b,a=signal.butter(3,2*np.pi*cutoff_freq,'low',analog=True)
            w, h = signal.freqs(b, a, worN=freq*2*np.pi)
            S_r=S_r*h   # S_r in complex now
            
        if DEBUG_PLOT:
            plt.plot(f,sp(f), label='SSR log')
            plt.plot(freq,sp2, label='SSR lin')
            plt.plot(freq,sp2, label='SSR taper')
            if cutoff_freq is not None:
                plt.plot(freq,abs(h), label=f"LP {cutoff_freq:.1f}")
            plt.plot(freq,abs(S_r), label="restit. fce")
            plt.xscale('log')
            #plt.yscale('log')
            #plt.ylim(0.25,None)
            plt.grid(which='both')
            plt.legend()
            plt.xlabel("f [Hz]")
            plt.ylabel("spectral ratio")
            plt.title(f"{sta} {zne}")
            plt.show()

        # modify complex spectrum
        # inverse fft
        y_r=np.fft.irfft(Cs*S_r)
        # create output trace
        out_tr=tr.copy()
        out_tr.data=y_r
        # output trace is added to the output stream
        out.append(out_tr)
        
        if DEBUG_PLOT:
            t=tr.stats.starttime
            tr.plot(starttime=t+DEBUG_tA,endtime=t+DEBUG_tB)
            out_tr.plot(starttime=t+DEBUG_tA,endtime=t+DEBUG_tB)
        
    return out

# from Klucenice_OR01A_OR01A-ref_b30_Z_ASD_SR.npz
# extract the station name which is OR01A
def extract_station_name(file_name):
    # use regular expression to extract the station name
    pattern = r'_(\w+)_'
    match = re.search(pattern, file_name)
    if match:
        station_name = match.group(1)
        return station_name
    else:
        return None


# from Klucenice_OR01A_OR01A-ref_b30_Z_ASD_SR.npz
# find the orientation component
def extract_orientation(file_name):
    # use regex to extract the orientation
    match = re.search(r'_(Z|N|E|H)_', file_name)
    if match:
        return match.group(1)
    else:
        return None


def load_freq_geomean(files):
    """
    Extract frequency and geometric mean arrays from npz files
    and store them in a dictionary with keys (station, component).
    :param files: list of npz file paths
    :type files: list of str
    :return: dictionary with keys (station, component) and values (freq, geomean)
    :rtype: dict
    """
    d_ssr=dict()

    for file in files:
        print(os.path.basename(file))
        sta=extract_station_name(file)
        if sta is None: continue
        zneh=extract_orientation(file)
        if zneh is None: continue
        print(sta,zneh)
        npz=np.load(file)
        freq=npz["freq"]
        geomean=npz['geomean']
        d_ssr[sta,zneh]=(freq,geomean)
    return d_ssr


# test main code

if __name__ == "__main__":
    
    SSR_DIR="/home/zapa/dat/zemetreseni_KLUCENICE/Amplifikace"
    
    pattern=SSR_DIR+"/Klucenice*_ASD_SR.npz"
    files=glob.glob(pattern)

    d_ssr=load_freq_geomean(files)
    
    st=read("/home/zapa/dat/klucenice/css26sel/kluce_20260306_141008.wfdisc",format='CSS')
    # in counts, not in m/s

    # filrace Butterworth LP, 3. řádu, 70 Hz
    st2=deconv_site_resp(st,d_ssr,cutoff_freq=70)

    print(st2)






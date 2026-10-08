#!/usr/bin/python3
# coding: utf-8

# akmich - Aplitudová KMItočtová CHarakteristika

import numpy as np
from obspy.core.trace import Trace

def FAS(data,delta):
    '''
    Výpočet Fourierova amplitudového spektra FAS
    FAS = sqrt(ESD), jednotky [U/Hz], např. [m/s/Hz]
    Vstup:
        data - 1D vektor dat, např. kalibrovaný seis. signál,
               tr.data*tr.stats.calib
        delta - vzorkovací interval [s],
                např. tr.stats.delta
    Výstup:
        freqs - řada kmitočtů
        amps  - hodnoty amplitudového spektra FAS
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # FAS = sqrt(2)*dt*|DFT|
    # backward normalization: DFT = sum(x[n]*exp(-2*pi*i*k*n/N))
    # rfft - DFT posloupnosti reálných čísel
    amps = np.sqrt(2)*delta*np.abs(np.fft.rfft(data))[1:]
    freqs = np.fft.rfftfreq(len(data), d=delta)[1:]
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz,
    # protože nelze zobrazit na logaritmické ose f
    return freqs, amps

def ASD(data,delta):
    '''
    Výpočet Amplitudové spektrální hustoty ASD
    ASD = sqrt(PSD), jednotky [U/sqrt(Hz)], např. [m/s/sqrt(Hz)]
    Vstup:
        data - 1D vektor dat, např. kalibrovaný seis. signál,
               tr.data*tr.stats.calib
        delta - vzorkovací interval [s],
                např. tr.stats.delta
    Výstup:
        freqs - řada kmitočtů
        amps  - hodnoty amplitudového spektra ASD
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # PSD = 2*dt/N*|DFT|^2
    # ASD = sqrt(2*dt/N)*|DFT|
    # rfft - DFT posloupnosti reálných čísel
    N = len(data)
    amps = np.sqrt(2*delta/N)*np.abs(np.fft.rfft(data))[1:]
    freqs = np.fft.rfftfreq(N, d=delta)[1:]
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz,
    # protože nelze zobrazit na logaritmické ose f
    return freqs, amps

def CSD(data1,data2,delta):
    '''
    Výpočet vzájemné výkonové spektrální hustoty CSD
    Vstup:
        data1, data2 - 1D vektor dat, např. kalibrovaný seis. signál,
        delta - vzorkovací interval [s],
    Výstup:
        freqs - řada kmitočtů
        pamps  - hodnoty PSD
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # PSD = 2*dt/N*|ZDFT|^2
    # rfft - DFT posloupnosti reálných čísel
    N1 = len(data1)
    N2 = len(data2)
    N=N1    #N1=N2
    S1 = np.fft.rfft(data1)
    S2 = np.fft.rfft(data2)
    S1S2 = S1*np.conj(S2)
    pamps = 2*delta/N*np.abs(S1S2)
    freqs = np.fft.rfftfreq(N, d=delta)
    return freqs[1:], pamps[1:]

def trFAS(tr:Trace):
    '''
    Výpočet FAS záznamu seis. signálu Trace
    Vstup:
        tr - obspy.core.trace.Trace
    Výstup:
        freqs - řada kmitočtů
        amps  - hodnoty amplitudového spektra FAS
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # vzorkovací interval delta
    delta=tr.stats.delta
    # kalibrace
    calib=tr.stats.calib
    # FAS = sqrt(2)*dt*|DFT|
    freqs,amps=FAS(tr.data*calib,delta)
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz
    return freqs, amps

def trASD(tr:Trace):
    '''
    Výpočet ASD záznamu seis. signálu Trace
    Vstup:
        tr - obspy.core.trace.Trace
    Výstup:
        freqs - řada kmitočtů
        amps  - hodnoty amplitudového spektra ASD
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # vzorkovací interval delta
    delta=tr.stats.delta
    # kalibrace
    calib=tr.stats.calib
    # ASD = sqrt(2*dt/N)*|DFT|
    freqs,amps=ASD(tr.data*calib,delta)
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz
    return freqs, amps

def ko_smoothing(ofreq,freqs,amps,b=40.0):
    '''
    Zhlazení spektra
    Konno-Ohmachi váhová funkce, limitace šířky filtru
    Vstup:
        ofreq - řada kmitočtů zhlazeného spektra
        freqs - lineární řada kmitočtů amplitudového spektra
        amps  - amplitudové spektrum odpovídající kmitočtové řadě freqs
        b     - parametr zhlazení
    Návratová hodnota:
        oamps - zhlazené amplitudy v řadě ofreq
    '''
    
    nN=len(freqs)
    Deltaf=(freqs[-1]-freqs[0])/(nN-1) # vzorkovací interval ve spektru

    w_f=2*np.pi/b        # šířka hlavního laloku [zlomek dekády]
    #w_f=0.7*w_f          # zúžení filtru
    alpha=10**(w_f/2)    # polovina kmitočtového intervalu
    
    c=b/np.pi
    oamps=np.empty_like(ofreq)
    if freqs[0] <= 1e-20: freqs[0]=1e-20
    lfreqs=np.log10(freqs)

    for i,f_c in enumerate(ofreq):
        if f_c <= 1e-20: f_c=1e-20
        lf_c=np.log10(f_c)
        # vstup freqs, amps může být výřezem z původního DFT spektrálního rozsahu
        # indexy pro výřez spektra mezi fA a fB, kde fA=f_c/alpha a fB=f_c*alpha
        fA=f_c/alpha                # šířka pásma od f_A
        fB=f_c*alpha                # šířka pásma do f_B
        #idx=np.where(np.logical_and(freqs>=fA-Deltaf/2,freqs<=fB+Deltaf/2))[0]
        idx=np.where(np.logical_and(freqs>=fA,freqs<=fB))[0]
        if len(idx) == 0:
            oamps[i]=np.nan
            continue
        w_lfreqs=lfreqs[idx]                # výřez kmitočtů
        w_amps=amps[idx]                    # výřez amplitud
        w_kos = np.sinc(c*(w_lfreqs-lf_c))**4 # váhová funkce
        w_kos /= w_kos.sum()                  # normování vah
        oamps[i]=w_kos.dot(w_amps)

    return oamps

def koc_smoothing(ofreq,freqs,amps,b=40.0):
    '''
    Zhlazení spektra s kompenzací posunu
    Konno-Ohmachi váhová funkce, limitace šířky filtru
    Místo váženého průměru se počítá integrál spektra
        podle logaritmu f 
    Vstup:
        ofreq - řada kmitočtů zhlazeného spektra
        freqs - lineární řada kmitočtů amplitudového spektra
        amps  - amplitudové spektrum odpovídající kmitočtové řadě freqs
        b     - parametr zhlazení
    Návratová hodnota:
        oamps - zhlazené amplitudy v řadě ofreq
    '''
    
    nN=len(freqs)
    Deltaf=(freqs[-1]-freqs[0])/(nN-1) # vzorkovací interval ve spektru

    w_f=2*np.pi/b        # šířka hlavního laloku [zlomek dekády]
    #w_f=0.7*w_f          # zúžení filtru
    alpha=10**(w_f/2)    # polovina kmitočtového intervalu
    
    c=b/np.pi
    oamps=np.empty_like(ofreq)
    if freqs[0] <= 1e-20: freqs[0]=1e-20
    lfreqs=np.log10(freqs)
    lfdiff=np.diff(lfreqs,append=lfreqs[-1])

    for i,f_c in enumerate(ofreq):
        if f_c <= 1e-20: f_c=1e-20
        lf_c=np.log10(f_c)
        # vstup freqs, amps může být výřezem z původního DFT spektrálního rozsahu
        # indexy pro výřez spektra mezi fA a fB, kde fA=f_c/alpha a fB=f_c*alpha
        fA=f_c/alpha                # šířka pásma od f_A
        fB=f_c*alpha                # šířka pásma do f_B
        #idx=np.where(np.logical_and(freqs>=fA-Deltaf/2,freqs<=fB+Deltaf/2))[0]
        idx=np.where(np.logical_and(freqs>=fA,freqs<=fB))[0]
        if len(idx) == 0:
            oamps[i]=np.nan
            continue
        w_lfd=lfdiff[idx]                   # váhy bez normování
        w_lfreqs=lfreqs[idx]                # výřez kmitočtů
        w_amps=amps[idx]                    # výřez amplitud
        w_kos = w_lfd*np.sinc(c*(w_lfreqs-lf_c))**4 # váhová funkce
        w_kos /= w_kos.sum()                  # normování vah
        oamps[i]=w_kos.dot(w_amps)

    return oamps


def main():
    help(FAS)
    help(ASD)
    help(trFAS)
    help(trASD)
    help(ko_smoothing)
    help(koc_smoothing)

if __name__ == "__main__":
    main()

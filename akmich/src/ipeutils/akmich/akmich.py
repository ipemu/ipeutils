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
        famps - řada kmitočtů
        amps  - hodnoty amplitudového spektra FAS
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # FAS = sqrt(2)*dt*|DFT|
    # backward normalization: DFT = sum(x[n]*exp(-2*pi*i*k*n/N))
    # rfft - DFT posloupnosti reálných čísel
    amps = np.sqrt(2)*delta*np.abs(np.fft.rfft(data))[1:]
    famps = np.fft.rfftfreq(len(data), d=delta)[1:]
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz,
    # protože nelze zobrazit na logaritmické ose f
    return famps, amps

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
        famps - řada kmitočtů
        amps  - hodnoty amplitudového spektra ASD
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # PSD = 2*dt/N*|DFT|^2
    # ASD = sqrt(2*dt/N)*|DFT|
    # rfft - DFT posloupnosti reálných čísel
    N = len(data)
    amps = np.sqrt(2*delta/N)*np.abs(np.fft.rfft(data))[1:]
    famps = np.fft.rfftfreq(N, d=delta)[1:]
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz,
    # protože nelze zobrazit na logaritmické ose f
    return famps, amps

def CSD(data1,data2,delta):
    '''
    Výpočet vzájemné výkonové spektrální hustoty CSD
    Vstup:
        data1, data2 - 1D vektor dat, např. kalibrovaný seis. signál,
        delta - vzorkovací interval [s],
    Výstup:
        famps - řada kmitočtů
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
    famps = np.fft.rfftfreq(N, d=delta)
    return famps[1:], pamps[1:]

def trFAS(tr:Trace):
    '''
    Výpočet FAS záznamu seis. signálu Trace
    Vstup:
        tr - obspy.core.trace.Trace
    Výstup:
        famps - řada kmitočtů
        amps  - hodnoty amplitudového spektra FAS
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # vzorkovací interval delta
    delta=tr.stats.delta
    # kalibrace
    calib=tr.stats.calib
    # FAS = sqrt(2)*dt*|DFT|
    famps,amps=FAS(tr.data*calib,delta)
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz
    return famps, amps

def trASD(tr:Trace):
    '''
    Výpočet ASD záznamu seis. signálu Trace
    Vstup:
        tr - obspy.core.trace.Trace
    Výstup:
        famps - řada kmitočtů
        amps  - hodnoty amplitudového spektra ASD
        Výstup bez stejnosměrné složky 0 Hz
    '''
    # vzorkovací interval delta
    delta=tr.stats.delta
    # kalibrace
    calib=tr.stats.calib
    # ASD = sqrt(2*dt/N)*|DFT|
    famps,amps=ASD(tr.data*calib,delta)
    # Vypustili jsme první člen amps[0], který odpovídá kmitočtu 0 Hz
    return famps, amps

def ko_smoothing(lfreq,famps,amps,b=40.0):
    '''
    Zhlazení spektra
    Konno-Ohmachi váhová funkce, limitace šířky filtru
    Vstup:
        lfreq - řada kmitočtů zhlazeného spektra
        famps - lineární řada kmitočtů amplitudového spektra
        amps  - amplitudové spektrum
        b     - parametr zhlazení
    Návratová hodnota:
        lamps - zhlazené amplitudy v řadě lfreq
    '''
    
    nN=len(famps)
    Deltaf=(famps[-1]-famps[0])/(nN-1) # vzorkovací interval ve spektru

    w_f=2*np.pi/b        # šířka hlavního laloku [zlomek dekády]
    w_f=0.7*w_f          # zúžení filtru
    alpha=10**(w_f/2)    # polovina kmitočtového intervalu
    
    c=b/np.pi
    lamps=np.empty_like(lfreq)
    lfamps=np.log10(famps)

    for i,f_c in enumerate(lfreq):
        lf_c=np.log10(f_c)
        
        fA=f_c/alpha                # šířka pásma od f_A
        iA=int(np.floor(fA/Deltaf)) # index spektra od
        if iA < 0: iA=0
        fB=f_c*alpha                # šířka pásma do f_B
        iB=int(np.ceil(fB/Deltaf))  # index spektra do
        if iB > nN: iB=nN  # max index iN-1, rozsah [:iN]

        w_lfamps=lfamps[iA:iB]                # výřez kmitočtů
        w_amps=amps[iA:iB]                    # výřez amplitud
        w_kos = np.sinc(c*(w_lfamps-lf_c))**4 # váhová funkce
        w_kos /= w_kos.sum()                  # normování vah
        lamps[i]=w_kos.dot(w_amps)

    return lamps

def koc_smoothing(lfreq,famps,amps,b=40.0):
    '''
    Zhlazení spektra s kompenzací posunu
    Konno-Ohmachi váhová funkce, limitace šířky filtru
    Místo váženého průměru se počítá integrál spektra
        podle logaritmu f 
    Vstup:
        lfreq - řada kmitočtů zhlazeného spektra
        famps - lineární řada kmitočtů amplitudového spektra
        amps  - amplitudové spektrum
        b     - parametr zhlazení
    Návratová hodnota:
        lamps - zhlazené amplitudy v řadě lfreq
    '''
    
    nN=len(famps)
    Deltaf=(famps[-1]-famps[0])/(nN-1) # vzorkovací interval ve spektru

    w_f=2*np.pi/b        # šířka hlavního laloku [zlomek dekády]
    w_f=0.7*w_f          # zúžení filtru
    alpha=10**(w_f/2)    # polovina kmitočtového intervalu
    
    c=b/np.pi
    lamps=np.empty_like(lfreq)
    lfamps=np.log10(famps)
    lfdiff=np.diff(lfamps)

    for i,f_c in enumerate(lfreq):
        lf_c=np.log10(f_c)
        
        fA=f_c/alpha                # šířka pásma od f_A
        iA=int(np.floor(fA/Deltaf)) # index spektra od
        if iA < 0: iA=0
        fB=f_c*alpha                # šířka pásma do f_B
        iB=int(np.ceil(fB/Deltaf))  # index spektra do
        #if iB > nN: iB=nN  # max index iN-1, rozsah [:iN]
        if iB > nN-1: iB=nN-1

        w_lfd=lfdiff[iA:iB]         # váhy bez normování

        w_lfamps=lfamps[iA:iB]                # výřez kmitočtů
        #print(len(w_lfamps),len(w_lfd))
        w_amps=amps[iA:iB]                    # výřez amplitud
        w_kos = w_lfd*np.sinc(c*(w_lfamps-lf_c))**4 # váhová funkce
        w_kos /= w_kos.sum()                  # normování vah
        #w_kos=np.flip(w_kos)
        lamps[i]=w_kos.dot(w_amps)

    return lamps


def main():
    help(FAS)
    help(ASD)
    help(trFAS)
    help(trASD)
    help(ko_smoothing)
    help(koc_smoothing)

if __name__ == "__main__":
    main()

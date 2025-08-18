#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import re

from scipy.signal import welch 
from scipy.integrate import trapezoid
from pathlib import Path

dir = Path()
if not dir.joinpath('corotatingData_0.0in.pkl').exists():
    dir = dir.joinpath('validate_data', 'corotatingBlades')
    
#%% Parameters 
Pref = 2e-5

# Blades
gaps = [0.0, 0.25, 0.5, 1.2, 2.0]
phases = [0, 45, 90]
RPM = [4000, 5000, 6000]
B = 2

# FFT Parameters
df = 3.125 # Hz
window = 'hann'
overlap = 0.75
freqCut = 1e4
around = 5 # range to integrate spectrum and get the tonals

pattern1 = r'(\d+(?:\.\d+)?)\s*RPM'
pattern2 = r'(\d+(?:\.\d+)?)\s*deg'



#%% Functions
def get_tonals(freq:np.ndarray, spectrum:np.ndarray, bpf:float) -> np.ndarray:
    freq_ref = freq/bpf
    bpf_pos = np.argmin(np.abs(freq_ref - 1))
    points = np.arange(bpf_pos-around, bpf_pos+around+1)
    tonals  = trapezoid(spectrum[points], freq[points], axis=1)
    return tonals

def mag2dB(mag):
    dB = 10*np.log10(mag/(Pref**2))
    return dB
#%% load data and reduce data
reduceData = np.zeros((len(gaps)*len(phases)*len(RPM), 28)) # gap | phase | RPM | BPF | Mic1 | Mic2 ... 

for k in range(len(gaps)):

    case = dir.joinpath(f'corotatingData_{gaps[k]}in.pkl')
    with open(case, 'rb') as file:
        cluster = pickle.load(file)

    # compute fft
    for data in cluster:
        t0 = data['time'][0]
        dt = data['time'][1]
        n_values = data['time'][2]
        time = np.arange(t0, dt*n_values, dt)
        nfft = int(1/(dt*df))
        noverlap = overlap*nfft
        
        freq, spectrum = welch(
            data['signal'], 
            fs = 1/dt, 
            window=window, 
            nperseg=nfft,
            noverlap=noverlap,
            axis=1
            )
        
        pos = freq <= freqCut
        
        freq = freq[pos]
        spectrum = spectrum[:, pos]
    
        if gaps[k] == 0:
            data['case'] = data['case'].replace('12x8', '0in')
            phase = 0
        else:
            phase = int(re.search(pattern2, data['case'].split('_')[0]).group(1))  # type: ignore
        
        rpm = int(re.search(pattern1, data['case'].split('_')[-1]).group(1)) # type: ignore
        BPF = rpm/60 * B
        tonals = get_tonals(freq, spectrum, BPF)
        
        ## Update matrix
        reduceData[0] = gaps[k]
        reduceData[1] = phase
        reduceData[2] = rpm
        reduceData[3] = BPF
        reduceData[4:] = mag2dB(tonals)
        
        
        # PLOT
        plt.figure(figsize=(11, 5))
        plt.loglog(np.repeat(cluster[0]['freq'].reshape(1,-1), 24, 0).T, cluster[0]['spectrum'].T,)
        plt.xlabel('f [Hz]')
        plt.ylabel('PSD [Pa²/Hz]')
        plt.legend([f'mic{i}' for i in range(1, 25)], ncols = 8, bbox_to_anchor =(0.5, 1.25), loc = 'upper center')

        plt.title('Individual' if gaps[k] != 0 else f'Gap = {gaps[k]} in' )

        plt.grid()
        plt.tight_layout()

        plt.show(block = False)
# %%

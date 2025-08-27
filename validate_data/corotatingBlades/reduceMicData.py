#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle

from matplotlib.colors import TABLEAU_COLORS
from matplotlib.gridspec import GridSpec
from scipy.signal import welch 
from scipy.integrate import trapezoid
from pathlib import Path


colors = list(TABLEAU_COLORS.keys())

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
around = 10 # range to integrate spectrum and get the tonals

#%% Functions
def get_tonals(freq:np.ndarray, spectrum:np.ndarray, bpf:float):
    freq_ref = freq/bpf
    bpf_pos = np.argmin(np.abs(freq_ref - 1))
    points = np.arange(bpf_pos-around, bpf_pos+around+1)
    tonals  = trapezoid(spectrum[:, points], freq[points], axis=1)
    return tonals, bpf_pos

def mag2dB(mag):
    dB = 10*np.log10(mag/(Pref**2))
    return dB
#%% load data and reduce data

reduceData = np.zeros((len(gaps)*len(phases)*len(RPM), 28)) # gap | phase | RPM | BPF | Mic1 | Mic2 ... 

# fig1 = plt.figure(1)
# fig2 = plt.figure(2)


count = 0
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
    
      
        rpm = data['RPM']
        phase = data['phase']
        
        BPF = rpm/60 * B
        tonals, bpf_pos = get_tonals(freq, spectrum, BPF)
        
        ## Update matrix
        reduceData[count, 0] = gaps[k]
        reduceData[count, 1] = phase
        reduceData[count, 2] = rpm
        reduceData[count, 3] = BPF
        reduceData[count, 4:] = mag2dB(tonals)
        count+=1
        
        # PLOT
        dB = mag2dB(spectrum[[0, 11, -1], :]*df)
        
        # plt.figure( figsize=(11, 5))
        # plt.semilogx(freq, dB[0] ,colors[0],  alpha = 0.7,label = 'Mic 1')
        # plt.semilogx(freq, dB[1],colors[1], alpha = 0.7,label = 'Mic 12')
        # plt.semilogx(freq, dB[2],colors[2], alpha = 0.7,label = 'Mic 24')
        
        # plt.semilogx(freq[bpf_pos], mag2dB(tonals[0]),  colors[0] , marker = 'o',)
        # plt.semilogx(freq[bpf_pos], mag2dB(tonals[11]), colors[1], marker = 'o',)
        # plt.semilogx(freq[bpf_pos], mag2dB(tonals[-1]), colors[2], marker = 'o',)
        
        # plt.xlabel('f [Hz]')
        # plt.ylabel('SPL [dB]')
        # plt.legend()

        # title = 'Individual' if k == 0 else f'Gap = {gaps[k]} in - Phase = {int(phase)}'
        # title = title + f' - {int(rpm)}RPM' 
        # plt.title(title)

        # plt.grid()
        # plt.tight_layout()

        # plt.show(block = False)
# %% save data in sheet
sheet =  pd.read_excel('corotating_data.xlsx')

for row in range(len(sheet.index)):
    gap = sheet.iloc[row, 1]
    phase = sheet.iloc[row, 2]
    rpm = sheet.iloc[row, 3]
    
    pos = (reduceData[:, 0] == gap) * (reduceData[:, 1] == phase) * (reduceData[:, 2] == rpm)   
    
    sheet.iloc[row, 7:] = reduceData[pos, 4:]
    # print(f' ------------- row = {row} -----------')
    # print(f' Reduce Data: {reduceData[pos,0:3]}')
    # print(f'Sheet {sheet.iloc[row, 1:4]}')

sheet.to_excel('corotating_data2.xlsx', index=False)

# %%

#%%
import numpy as np
import pickle

# from scipy.signal import welch 
# from scipy.integrate import trapezoid
from scipy.io import loadmat
from pathlib import Path

#%% Paths
dir = Path('D:','Enrico','20241206')
paths = [i.joinpath('Ruído') for i in list(dir.glob('12x8*'))]

#%% get gaps
gap = [i.name.split(' ')[-1] for i in list(dir.glob('12x8*'))]
gap = ['0' if i=='Individual' else i for i in gap]
gap = [0.25 if i ==  'stacked' else float(i.split('in')[0]) for i in gap] 
# %% load data
# dataBase  = {g: None for g in gap}

for i, path in enumerate(paths):
    print(path.parent   .name)
    files =list(path.glob('*.mat'))
    cluster = np.empty(len(files), dtype=dict)
    for k, file in enumerate(files):
        
        data = loadmat(file, simplify_cells = True)
        t0 = data['Signal_00']['x_values']['start_value']
        dt = data['Signal_00']['x_values']['increment']
        n_values = data['Signal_00']['x_values']['number_of_values']
        # time = np.arange(t0, (n_values)*dt, dt)
        time = np.array([t0, dt, n_values])
        signal = np.array([data[key]['y_values']['values'] for key in list(data.keys())[3:]])        
        cluster[k] = dict(time = time,signal = signal,case = file.stem.split(' ')[-1])
        
    with open(f'corotatingData_{gap[i]}in.pkl', 'wb') as file:
        pickle.dump(cluster, file)
    # dataBase[gap[i]] = cluster # type: ignore
         
#%% write

        
        
# %%

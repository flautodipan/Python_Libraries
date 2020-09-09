#%%

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt

##%
# GYRATIO
WTC_identifier = ('wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6','wtc7_16')
gyrs = {}
norm_gyrs = {}


for key in WTC_identifier:
    gyrs[key] = BA.Parse_xvg_skip('gyration_'+key+'_BS_RNA.xvg', path = '../GROMACS/'+key.upper()+'/' if key != 'wtc1_h' else '../GROMACS/WTC1_h/', skip_lines= 28, )[1]
    norm_gyrs[key] = (gyrs[key] - np.min(gyrs[key]))/(np.max(gyrs[key]) - np.min(gyrs[key]))
#%%
f, ax = plt.subplots(1,1)
for key in WTC_identifier:
    sampling = np.arange(0, len(gyrs[key]), 100)
    ax.plot(norm_gyrs[key][sampling], '--', label = key if key == 'wtc7_16' else None, color = 'firebrick' if key == 'wtc7_16' else 'grey', alpha = .7 if key == 'wtc7_16' else 0.4)
    ax.set_title('Gyration Radii comparison')
    ax.set_xlabel('Time (a.u.)')
    ax.set_ylabel('Gyration radium normalized')
    plt.legend()
    plt.xlim(0,100)


# %%
f, ax = plt.subplots(1,1)
for key in WTC_identifier:
    sampling = np.arange(0, len(gyrs[key]), 100)
    ax.plot(gyrs[key][sampling], '--', label = key if key == 'wtc7_16' else None, color = 'firebrick' if key == 'wtc7_16' else 'grey', alpha = .7 if key == 'wtc7_16' else 0.4)
    ax.set_title('Gyration Radii comparison')
    ax.set_xlabel('Time (a.u.)')
    ax.set_ylabel('Gyration radium (nm)')
    plt.legend()
    #plt.xlim(0,100)

# %%

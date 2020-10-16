#%%
import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.spatial.distance import euclidean
from scipy.stats import pearsonr

path = '../GROMACS/'
wtc_keys = ['wtc1_h', 'wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
wtc_keys_red = [wtc_keys[1]] + wtc_keys[3:]
cov_keys = ['cov2', 'cov3', 'cov4']
mtc_keys = ['mtc1', 'mtc2']

exp_df = pd.read_excel(path+'MD_experimental_data.xlsx')

eq = 'eq1'

#%%
# inputs

correlators = []
x_vals = []
x_errs = []

now_keys = wtc_keys_red+mtc_keys+cov_keys

for now_name in now_keys:
    print(now_name)
    now_exp_df = exp_df[exp_df.identifier == now_name]
    now_path = path + now_name.upper()  + '/'
    if eq=='eq':
        correlators.append(np.mean(np.load(now_path+now_name+'_min_CAPdist_'+eq+'.npy')))
    else:        
        correlators.append(np.mean(np.load(now_path+now_name+'_min_CAPdist_'+eq+'.npy')) if now_name == 'wtc1_h_new' else np.mean(np.load(now_path+now_name+'_min_CAPdist_eq.npy')))


    x_vals.append(now_exp_df.Kd[now_exp_df.identifier == now_name])
    x_errs.append(now_exp_df.Kd_err[now_exp_df.identifier == now_name])

y_err = (correlators[4]-correlators[1])/np.sqrt(2)

#%%
fig, ax = plt.subplots()

#for x, x_err, y, key in zip(x_vals[1:], x_errs[1:], correlators[1:], now_keys[1:]):
for x, x_err, y, key in zip(x_vals, x_errs, correlators, now_keys):

    ax.errorbar(x, y, xerr = x_err, yerr = y_err , fmt = 'o', ecolor= 'maroon', color = 'green', label = 'mtc series' if key in wtc_keys[-1] else None)
    plt.annotate('NMR' if key == 'wtc1' else 'wtc1' if key == 'wtc1_h' else 'wtc1_new' if key == 'wtc1_h_new' else key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_title('CA-P minimum mean distance vs Kd')
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('minimum distance (ang)')
plt.legend()
plt.show()

#%%
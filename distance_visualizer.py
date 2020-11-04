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
wtc_keys = ['wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
cov_keys = ['cov2', 'cov3', 'cov4']
mtc_keys = ['MTC1', 'mtc2', 'mtc3', 'mtc4']

exp_df = pd.read_excel(path+'MD_experimental_data.xlsx')


#%%
# inputs

correlators = []
x_vals = []
x_errs = []

now_keys = wtc_keys + mtc_keys
now_eqs1 = ['wtc1_h_new', 'MTC1',  'mtc3', 'mtc4', 'wtc7_16']
now_eqs2 = [ 'mtc2']


for key in now_keys:
    print(key)
    now_exp_df = exp_df[exp_df.identifier == key]
    now_path = path + key.upper()  + '/'

    eq = '_eq2' if key in now_eqs2 else '_eq1' if key in now_eqs1 else '_eq' 
    if ((key in now_eqs1) & (key in now_eqs2)):
        raise ValueError('{} in both eq1 and eq2'.format(key))

    correlators.append(np.mean(np.load(now_path+key+'_min_CAPdist'+eq+'.npy')))
    print('Acquisisco {} con eq = {}'.format(key, eq, ))


    x_vals.append(now_exp_df.Kd[now_exp_df.identifier == key])
    x_errs.append(now_exp_df.Kd_err[now_exp_df.identifier == key])

y_err = (correlators[4]-correlators[1])/np.sqrt(2)

#%%
fig, ax = plt.subplots()

#for x, x_err, y, key in zip(x_vals[1:], x_errs[1:], correlators[1:], now_keys[1:]):
for x, x_err, y, key in zip(x_vals, x_errs, correlators, now_keys):

    ax.errorbar(x, y, xerr = x_err, fmt = 'o', ecolor= 'maroon', color = 'green' if key in mtc_keys else 'firebrick' if key in wtc_keys else 'k', label = 'wtc_series' if key == wtc_keys[-1] else 'cov_series' if key == cov_keys[-1] else 'mtc series' if key == mtc_keys[-1] else None)
    plt.annotate('NMR' if key == 'wtc1' else 'wtc1' if key == 'wtc1_h' else 'wtc1_new' if key == 'wtc1_h_new' else key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_title('CA-P minimum mean distance vs Kd')
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('minimum distance (ang)')
plt.legend()
plt.show()

#%%
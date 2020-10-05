#%%
#0) preambolo

import  BioAlessandria as BA
from    Alessandria import Find_Nearest
import  numpy as np
import  warnings
import  pandas as pd
import  matplotlib.pyplot as plt
import  os
warnings.filterwarnings("ignore")


path = '../GROMACS/'
wtc_keys = ['wtc1_h', 'wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
cov_keys = ['cov2', 'cov3', 'cov4']
mtc_keys = ['mtc1']

#%%
# FORMO IL DATAFRAME con i DATI CHE VOGLIO STAMPARE ORA
now_keys = wtc_keys

dfs = {}
for now_name in now_keys:

    now_path    = path + now_name.upper() +'/'
    dfs[now_name] = pd.read_json(now_path+now_name+'_df.json')

DF = pd.DataFrame()
DF = pd.concat([dfs[key] for key in dfs.keys()], ignore_index=True)

#%%
# riferimenti e z_score

ref_mean = np.mean(DF.Covariance_Mean_Prot_BS[DF.identifier == 'wtc7_16'][DF.Is_BS == True])
ref_std  = np.std(DF.Covariance_Mean_Prot_BS[DF.identifier == 'wtc7_16'][DF.Is_BS == True])

DF['z_Covariance_Mean_Prot_BS_12'] = (DF.Covariance_Mean_Prot_BS - ref_mean)/ref_std
z_covave_rmsf_BS_BS_12 = {key : DF.z_Covariance_Mean_Prot_BS_12[DF.identifier == key][DF.Is_BS == True][DF.Is_Prot == True] for key in now_keys }

# %%
#PLOT di Covarianza
# preparo le variabili

x_vals = [DF.Kd[DF.identifier == key].values[0] for key in now_keys]
x_errs = [DF.Kd[DF.identifier == key].values[0] for key in now_keys]

correlators = [np.mean(z_covave_rmsf_BS_BS_12[key]) for key in now_keys]
normalizers = [np.mean(DF.RMSF[DF.identifier == key][DF.Is_Prot == True][DF.Is_BS==True])  for key in now_keys]

#stampo figura
# non normalizzata

fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd')

for x,x_err, y, key in zip(x_vals, x_errs, correlators, now_keys):
    ax.errorbar(x, y, fmt = 'o', color = 'k', ecolor = 'orange', mew = 0.1)
    plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
plt.show()

#normalizzata
fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd')

for x,x_err, y, key, norm in zip(x_vals, x_errs, correlators, now_keys, normalizers):
    ax.errorbar(x, y/norm, fmt = 'o', color = 'k', ecolor = 'orange', mew = 0.1)
    plt.annotate(key, (x,y/norm), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS / <RMSF> ')
plt.show()

#%%



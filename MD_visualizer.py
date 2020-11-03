#%%
#0) preambolo

import  BioAlessandria as BA
from    Alessandria import Find_Nearest
import  numpy as np
import  warnings
import  pandas as pd
import  matplotlib.pyplot as plt
import  os
from os.path import join
warnings.filterwarnings("ignore")


path = '../GROMACS/'
wtc_keys = ['wtc1_h', 'wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5',  'wtc6', 'wtc7_16']
wtc_keys_red = [wtc_keys[1]]+wtc_keys[3:]
cov_keys = ['cov2', 'cov3', 'cov4']
mtc_keys = ['MTC1', 'mtc2', 'mtc3', 'mtc4']

all_keys = wtc_keys+cov_keys+mtc_keys
gian_keys = [wtc_keys[1]] + wtc_keys[3:]

#%%
# FORMO IL DATAFRAME con i DATI CHE VOGLIO STAMPARE ORA
now_keys = wtc_keys_red +mtc_keys
now_eqs1 = ['wtc1_h_new', 'MTC1', 'mtc3']
now_eqs2 = ['mtc2',]

dfs = {}
for now_name in now_keys:  

    now_path    = join(path, now_name.upper())

    #gestione eq
    eq = '_eq2' if now_name in now_eqs2 else '_eq1' if now_name in now_eqs1 else '_eq' 
    if ((now_name in now_eqs1) & (now_name in now_eqs2)):
        raise ValueError('{} in both eq1 and eq2'.format(now_name))
    dfs[now_name] = pd.read_json(join(now_path, now_name+'_df'+eq+'.json'))
    print('Acquisisco {} con eq = {}'.format(now_name, eq, ))

DF = pd.DataFrame()
DF = pd.concat([dfs[key] for key in dfs.keys()], ignore_index=True)

#DF = DF.rename(columns={'Is_BS_05' : 'Is_BS_O5'})

#%%
# riferimenti e z_score
# O5
ref_mean_O5 = np.mean(DF.Covariance_Mean_Prot_BS_O5[DF.identifier == 'wtc1_h_new'][DF.Is_BS_O5 == True][DF.Is_Prot == True])
ref_std_O5  = np.std(DF.Covariance_Mean_Prot_BS_O5[DF.identifier == 'wtc1_h_new'][DF.Is_BS_O5 == True][DF.Is_Prot == True])
DF['z_Covariance_Mean_Prot_BS_12_O5'] = (DF.Covariance_Mean_Prot_BS_O5 - ref_mean_O5)/ref_std_O5
z_covave_rmsf_BS_BS_12_O5 = {key : DF.z_Covariance_Mean_Prot_BS_12_O5[DF.identifier == key][DF.Is_BS_O5 == True][DF.Is_Prot == True] for key in now_keys }
pearson_BS_BS_O5 = {key : DF.Pearson_Mean_Prot_BS_O5[DF.identifier == key][DF.Is_BS_O5 == True][DF.Is_Prot == True] for key in now_keys }

#normalizzati da rmsf
norm_covave_BS_BS_12_O5 = {key : (DF.Covariance_Mean_Prot_BS_O5[DF.identifier == key][DF.Is_BS_O5 == True][DF.Is_Prot == True])/np.mean(DF.RMSF[DF.identifier == key][DF.Is_Prot == True][DF.Is_BS_O5==True]) for key in now_keys }


# P
ref_mean_P = np.mean(DF.Covariance_Mean_Prot_BS_P[DF.identifier == 'wtc1_h_new'][DF.Is_BS_P == True][DF.Is_Prot == True])
ref_std_P = np.std(DF.Covariance_Mean_Prot_BS_P[DF.identifier == 'wtc1_h_new'][DF.Is_BS_P == True][DF.Is_Prot == True])
DF['z_Covariance_Mean_Prot_BS_12_P'] = (DF.Covariance_Mean_Prot_BS_P - ref_mean_P)/ref_std_P

z_covave_rmsf_BS_BS_12_P = {key : DF.z_Covariance_Mean_Prot_BS_12_P[DF.identifier == key][DF.Is_BS_P == True][DF.Is_Prot == True] for key in now_keys }
pearson_BS_BS_P = {key : DF.Pearson_Mean_Prot_BS_P[DF.identifier == key][DF.Is_BS_P == True][DF.Is_Prot == True] for key in now_keys }
#normalizzati da rmsf
norm_covave_BS_BS_12_P = {key : (DF.Covariance_Mean_Prot_BS_P[DF.identifier == key][DF.Is_BS_P == True][DF.Is_Prot == True])/np.mean(DF.RMSF[DF.identifier == key][DF.Is_Prot == True][DF.Is_BS_P==True]) for key in now_keys }

# %%
"""
#PLOT di Covarianza per O5

# preparo le variabili

x_vals = [DF.Kd[DF.identifier == key].values[0] for key in now_keys]
x_errs = [DF.Kd_err[DF.identifier == key].values[0] for key in now_keys]

correlators_O5 = [np.mean(z_covave_rmsf_BS_BS_12_O5[key]) for key in now_keys]
normaletors_O5 = [np.mean((norm_covave_BS_BS_12_O5[key] - ref_mean_O5)/ref_std_O5)  for key in now_keys]
#stampo figura
# non normalizzata

fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd\nRNA by O5 atoms')


for x,x_err, y, key in zip(x_vals, x_errs, correlators_O5, now_keys):
    ax.errorbar(x, y, fmt = 'o', ecolor = 'orange', mew = 0.1, color = 'black' if key in wtc_keys else 'green' if key in cov_keys else 'goldenrod', label = 'wtc_series' if key == wtc_keys[-1] else 'cov_series' if key == cov_keys[-1] else 'mtc series' if key == mtc_keys[-1] else None)
    plt.annotate('NMR' if key == 'wtc1' else 'wtc1' if key == 'wtc1_h' else 'wtc1_new' if key == 'wtc1_h_new' else key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
plt.legend()
plt.show()

#normalizzata
fig, ax = plt.subplots()
ax.set_title('Protein BS Normalized Covariance with Protein BS  vs Kd\nRNA by O5 atoms')

for x,x_err, y, key in zip(x_vals, x_errs, normaletors_O5, now_keys,):
    ax.errorbar(x, y, fmt = 'o', ecolor = 'orange', mew = 0.1, color = 'black' if key in wtc_keys else 'green' if key in cov_keys else 'goldenrod', label = 'wtc_series' if key == wtc_keys[-1] else 'cov_series' if key == cov_keys[-1] else 'mtc series' if key == mtc_keys[-1] else None)
    plt.annotate('NMR' if key == 'wtc1' else 'wtc1' if key == 'wtc1_h' else 'wtc1_new' if key == 'wtc1_h_new' else key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS / <RMSF> ')
plt.legend()
plt.show()
"""
#%%


#PLOT di Covarianza per P

# preparo le variabili

x_vals = [DF.Kd[DF.identifier == key].values[0] for key in now_keys]
x_errs = [DF.Kd_err[DF.identifier == key].values[0] for key in now_keys]


#pearson

correlators_P = [np.mean(pearson_BS_BS_P[key]) for key in now_keys]

#stampo figura
# non normalizzata

fig, ax = plt.subplots()
ax.set_title('Protein BS Pearson with Protein BS  vs Kd\nRNA by P atoms')

for x,x_err, y, key in zip(x_vals, x_errs, correlators_P, now_keys):
#for x,x_err, y, key in zip(x_vals[1:], x_errs[1:], correlators_P[1:], now_keys[1:]):
    ax.errorbar(x, y, xerr = x_err, fmt = 'o', ecolor = 'orange', mew = 0.1, color = 'black' if key in wtc_keys else 'green' if key in cov_keys else 'goldenrod', label = 'wtc_series' if key == wtc_keys[-1] else 'cov_series' if key == cov_keys[-1] else 'mtc series' if key == mtc_keys[-1] else None)
    plt.annotate('NMR' if key == 'wtc1' else 'wtc1' if key == 'wtc1_h' else 'wtc1_new' if key == 'wtc1_h_new' else key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS Pearson correlation with Prot BS ')
plt.legend()
plt.show()

# covariance
correlators_P = [np.mean(z_covave_rmsf_BS_BS_12_P[key]) for key in now_keys]

"""
same_phys = [('wtc1_h_new', 'wtc7_16',), ('wtc2','wtc5'), ('wtc3', 'wtc4'), ]
errors = []
for which in same_phys:
    now_delta = np.abs(np.mean(z_covave_rmsf_BS_BS_12_P[which[0]])- np.mean(z_covave_rmsf_BS_BS_12_P[which[1]]))
    print(now_delta)
    errors.append(now_delta)
print(errors)
"""
y_err_P = 0.1#np.mean(errors)/np.sqrt(2)
print('Media dello scarto in y su punti con stessa fisica (stessa x ) =\n{}\n'.format(y_err_P))

fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd\nRNA by P atoms')

for x,x_err, y, key in zip(x_vals, x_errs, correlators_P, now_keys):
    ax.errorbar(x, y, xerr= x_err, fmt = 'o', ecolor = 'orange', mew = 0.1, color = 'black' if key in wtc_keys else 'green' if key in cov_keys else 'goldenrod', label = 'wtc_series' if key == wtc_keys[-1] else 'cov_series' if key == cov_keys[-1] else 'mtc series' if key == mtc_keys[-1] else None)
    plt.annotate('NMR' if key == 'wtc1' else 'wtc1' if key == 'wtc1_h' else 'wtc1_new' if key == 'wtc1_h_new' else key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
plt.legend()
plt.show()



# %%


# %%

#%%
# 0) librerie e definizioni

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr



now_path = '../GROMACS/'
save_path = now_path+'COMPLEX/'

cov_identifier = ('cov1', 'cov2', 'cov3', 'cov4')
Kds = [5.4, 1500, 115, 6.9]

dodici = 12

#%%

DF = pd.read_json(now_path+'cov_data_frame_{}.json'.format(str(dodici)))


z_covave_rmsf_BS_BS = {key : DF.z_Covariance_Mean_Prot_BS[DF.cov_identifier == key][DF.Is_BS == True][DF.Is_Prot == True] for key in cov_identifier }
z_covave_rmsf_BS_BS_correlators = [np.mean(z_covave_rmsf_BS_BS[key]) for key in cov_identifier]
z_covave_rmsf_BS_BS_correlations = pearsonr(Kds, z_covave_rmsf_BS_BS_correlators)


z_covave_rmsf_BS_BS_12 = {key : DF.z_Covariance_Mean_Prot_BS_12[DF.cov_identifier == key][DF.Is_BS == True][DF.Is_Prot == True] for key in cov_identifier }
z_covave_rmsf_BS_BS_correlators_12 = [np.mean(z_covave_rmsf_BS_BS_12[key]) for key in cov_identifier]
z_covave_rmsf_BS_BS_correlations_12 = pearsonr(Kds, z_covave_rmsf_BS_BS_correlators)


#%%

correlators = z_covave_rmsf_BS_BS_correlators_12
correlations = z_covave_rmsf_BS_BS_correlations_12
normalizers = [np.mean(DF.RMSF[DF.cov_identifier == key][DF.Is_Prot == True][DF.Is_BS==True])  for key in cov_identifier]



fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd')
ax.errorbar(Kds , correlators/np.array(normalizers), fmt = 'o', color = 'k', ecolor = 'orange', label = 'simulated data', mew = 0.1)
#ax.errorbar(Kds[0] , correlators[0], fmt = 'o', xerr = Kds_errs[0], color = 'orange', ecolor = 'k', label = 'NMR data', mew = 0.1)

for x, key, ll in zip(Kds, cov_identifier, range(len(cov_identifier))):
        y = correlators[ll]/float(normalizers[ll])
        if '2' in key:
                plt.annotate(key, (x,y), xytext = (-10, -10), textcoords = 'offset points')
        else:
                plt.annotate(key, (x,y), xytext = (10, 0), textcoords = 'offset points')


ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
#ax.vlines(np.mean([Kds[kk] for kk in to_exclude]), min_err, max_err, linestyles = 'dashed', color = 'firebrick', label = 'max error bar', linewidths = 2.)
ax.legend()
fig.savefig(save_path+'corr_plot_cov.pdf', format = 'pdf')
plt.show()



# %%
#insieme alla tesi

DF_1 = pd.read_json(now_path+'WTC_data_frame_9ang_16_all.json')

Kds_new =  [4, 4, 1360, 650, 750, 1320, 1800, 70]
Kds_errs = [0.9, 0.9, 600, 135, 150, 350, 400, 10]

WTC_identifier = ('wtc1_h', 'wtc1', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6', 'wtc7')


z_covave_rmsf_BS_BS_1 = {key : DF_1.z_Covariance_Mean_Prot_BS[DF_1.WTC_identifier == key][DF_1.Is_BS == True][DF_1.Is_Prot == True] for key in WTC_identifier }
z_covave_rmsf_BS_BS_correlators_1 = [np.mean(z_covave_rmsf_BS_BS_1[key]) for key in WTC_identifier]
z_covave_rmsf_BS_BS_correlations_1 = pearsonr(Kds_new, z_covave_rmsf_BS_BS_correlators_1)


z_covave_rmsf_BS_BS_12_1 = {key : DF_1.z_Covariance_Mean_Prot_BS_12[DF_1.WTC_identifier == key][DF_1.Is_BS == True][DF_1.Is_Prot == True] for key in WTC_identifier }
z_covave_rmsf_BS_BS_correlators_12_1 = [np.mean(z_covave_rmsf_BS_BS_12_1[key]) for key in WTC_identifier]
z_covave_rmsf_BS_BS_correlations_12_1 = pearsonr(Kds_new, z_covave_rmsf_BS_BS_correlators_1)



correlators_1 = z_covave_rmsf_BS_BS_correlators_12_1
correlations_1 = z_covave_rmsf_BS_BS_correlations_12_1

normalizers_1 = [np.mean(DF_1.RMSF[DF_1.WTC_identifier == key][DF_1.Is_Prot == True][DF_1.Is_BS==True])  for key in WTC_identifier]



#STAMPO I PUNTI normalizzati con rmsf medio della BS di ogni dinamica
all_correlators = correlators_1+correlators
fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd')
ax.errorbar(Kds_new[1:] , correlators_1[1:]/np.array(normalizers_1[1:]), fmt = 'o', color = 'k', ecolor = 'orange', label = 'thesis data', mew = 0.1)
ax.errorbar(Kds_new[0] , correlators_1[0]/np.array(normalizers_1[0]), fmt = 'o', color = 'orange', ecolor = 'k', label = 'NMR data', mew = 0.1)
ax.errorbar(Kds , correlators/np.array(normalizers), fmt = 'o', color = 'g', ecolor = 'orange', label = 'new data', mew = 0.1)
for x, key, ll, norm in zip(Kds_new + Kds, WTC_identifier+cov_identifier, range(len(WTC_identifier+cov_identifier)), normalizers_1 + normalizers):
        y = all_correlators[ll]/float(norm)
        if key in ['wtc1_h_new']:
            key = 'wtc1'
            plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )
        elif key == 'wtc1': 
                key = 'NMR'
                plt.annotate(key, (x,y), xytext = (-10, 10), textcoords = 'offset points', )
        elif '5' in key:
                plt.annotate(key, (x,y), xytext = (-10, -10), textcoords = 'offset points', )

        elif key in ['wtc3', 'wtc6']:plt.annotate(key, (x,y), xytext = (-25, -15), textcoords = 'offset points', )
        elif key == 'wtc7': 
                plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )
        elif key in ['cov1', 'cov4']:
                plt.annotate(key, (x,y), xytext = (10, 0), textcoords = 'offset points', )




        else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS / < RMSF > ')
#ax.vlines(np.mean([Kds[kk] for kk in to_exclude]), min_err, max_err, linestyles = 'dashed', color = 'firebrick', label = 'max error bar', linewidths = 2.)
ax.legend()
plt.show()


#STAMPO I PUNTI NON normalizzati 
all_correlators = correlators_1+correlators
fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd')
ax.errorbar(Kds_new[1:] , correlators_1[1:], fmt = 'o', color = 'k', ecolor = 'orange', label = 'thesis data', mew = 0.1)
ax.errorbar(Kds_new[0] , correlators_1[0], fmt = 'o', color = 'orange', ecolor = 'k', label = 'NMR data', mew = 0.1)
ax.errorbar(Kds , correlators, fmt = 'o', color = 'g', ecolor = 'orange', label = 'new data', mew = 0.1)
for x, key, ll, norm in zip(Kds_new + Kds, WTC_identifier+cov_identifier, range(len(WTC_identifier+cov_identifier)), normalizers_1 + normalizers):
        y = all_correlators[ll]
        if key in ['wtc1_h_new']:
            key = 'wtc1'
            plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )
        elif key == 'wtc1': 
                key = 'NMR'
                plt.annotate(key, (x,y), xytext = (-10, 10), textcoords = 'offset points', )
        elif key in ['wtc3', 'wtc6']:plt.annotate(key, (x,y), xytext = (-25, -15), textcoords = 'offset points', )
        elif key == 'wtc7': 
                plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )
        elif key in ['cov1', 'cov4']:
                plt.annotate(key, (x,y), xytext = (10, 0), textcoords = 'offset points', )




        else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')

ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
#ax.vlines(np.mean([Kds[kk] for kk in to_exclude]), min_err, max_err, linestyles = 'dashed', color = 'firebrick', label = 'max error bar', linewidths = 2.)
ax.legend()
plt.show()

# %%

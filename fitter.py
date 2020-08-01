# Versione del programma correlator_new che include il settimo punto e 
# considera già di per sé 
# - che escludo il wtc1_home
# - che mi concentro sui 9 angstrom di definizione per il BS
# - che uso z_covave_BS_BS (v. programma correlator_new per significato)


#%%

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import scipy.odr as odr


def f(A, x):
    return A[0]*x+A[1]

now_path = '../GROMACS/'
save_path = now_path+'CORRELATIONS/'

treshold = 9
print('\n\n\nTreshold = 9 Angstrom, chosen as the best one\n\n\n')

DF = pd.read_json(now_path+'WTC_data_frame_{}ang.json'.format(str(treshold)))
DF['z_Covariance_Mean'] = (DF.Covariance_Mean - np.mean(DF.Covariance_Mean))/(np.std(DF.Covariance_Mean))
DF['z_Covariance_Mean_Prot_BS'] = (DF.Covariance_Mean_Prot_BS - np.mean(DF.Covariance_Mean_Prot_BS))/(np.std(DF.Covariance_Mean_Prot_BS))
DF['z_Covariance_Mean_RNA_BS'] = (DF.Covariance_Mean_RNA_BS - np.mean(DF.Covariance_Mean_RNA_BS))/(np.std(DF.Covariance_Mean_RNA_BS))

#DF['z_RMSF'] = (DF.RMSF.values - np.mean(DF.RMSF.values))/np.std(DF.RMSF.values)

#WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6')
WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6', 'wtc7')

Kds = [DF.Kd[DF.WTC_identifier == key].values[0] for key in WTC_identifier]

Kds = [4, 4, 1360, 650, 750, 1320, 1800, 70]
Kds_errs = [0.9, 0.9, 600, 135, 150, 350, 400, 10]

colori = ['royalblue', 'cornflowerblue', 'forestgreen', 'goldenrod', 'orange', 'darkorchid', 'firebrick', 'black']
colors = {wtc : color for (wtc, color) in zip(WTC_identifier, colori)}

#%%

z_covave_rmsf_BS_BS = {key : DF.z_Covariance_Mean_Prot_BS[DF.WTC_identifier == key][DF.Is_BS == True][DF.Is_Prot == True] for key in WTC_identifier }
z_covave_rmsf_BS_BS_correlators = [np.mean(z_covave_rmsf_BS_BS[key]) for key in WTC_identifier]
z_covave_rmsf_BS_BS_correlations = pearsonr(Kds[1:], z_covave_rmsf_BS_BS_correlators[1:])

#%%

# FIT

linear = odr.Model(f)

# ERRORI presi sui 2 punti che dovrebbero essere uguali,
#ciclo sugli step per calcolarli

to_exclude = [3,4]

max_err = np.max([z_covave_rmsf_BS_BS_correlators[jj] for jj in to_exclude])
min_err = np.min([z_covave_rmsf_BS_BS_correlators[jj] for jj in to_exclude])
err_max = (max_err - min_err)
err = err_max/3 #procedura standard
print('Errore stimato = 0', err/2)


correlators = z_covave_rmsf_BS_BS_correlators[1:]
correlations = z_covave_rmsf_BS_BS_correlations[1:]



#STAMPO I PUNTI

fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd\ncorrelation = {:3.2f} p-value = {:3.2f}'.format(*pearsonr(Kds[1:], correlators)))
ax.errorbar(Kds[1:] , correlators, fmt = 'o', xerr = Kds_errs[1:], color = 'k', ecolor = 'orange', label = 'simulated data', mew = 0.1)
ax.errorbar(Kds[0] , z_covave_rmsf_BS_BS_correlators[0], fmt = 'o', xerr = Kds_errs[0], color = 'orange', ecolor = 'k', label = 'NMR data', mew = 0.1)

for x,y, key in zip(Kds , z_covave_rmsf_BS_BS_correlators, WTC_identifier):

    if key in ['wtc1_h']:
        key = 'wtc1'
        plt.annotate(key, (x,y), xytext = (5, -12), textcoords = 'offset points', )
    elif key == 'wtc3': plt.annotate(key, (x,y), xytext = (-22, -15), textcoords = 'offset points', )
    elif key == 'wtc1': continue
    else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')
    
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
ax.vlines(np.mean([Kds[kk] for kk in to_exclude]), min_err, max_err, linestyles = 'dashed', color = 'firebrick', label = 'max error bar', linewidths = 2.)
ax.legend()
plt.show()

#FACCIO E STAMPO IL FIT

mydata = odr.RealData(Kds[1:], correlators, sy = err/2, sx = Kds_errs[1:])
myodr = odr.ODR(mydata, linear, beta0=[1.,2.])
myoutput = myodr.run()

fig, ax = plt.subplots()
x = np.linspace(np.min(Kds), np.max(np.array(Kds)+np.array(Kds_errs)), 1000)
ax.set_title(r'$\bf{y = mx + q}$ fit'+'\nProtein BS Covariance with Protein BS vs Kd'+'\ncorrelation = {:3.2f} p-value = {:3.2f}'.format(*pearsonr(Kds[1:], correlators)))
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
ax.errorbar(Kds[1:] , correlators, fmt = 'o', xerr = Kds_errs[1:], yerr = err, color = 'k', ecolor = 'orange', label = 'correlation data', mew = 0.1)
ax.plot(x,f(myoutput.beta, x), color = 'yellowgreen', label = 'linear fit')
ax.legend(title = 'm = {:3.2e} $\pm$ {:3.2e}\nq = {:3.2e} $\pm$ {:3.2e}'.format(myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1]))
plt.show()




# %%

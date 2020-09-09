# Versione del programma correlator_new che include il settimo punto e 
# considera già di per sé 
# - che escludo il wtc1_home
# - che mi concentro sui 9 angstrom di definizione per il BS
# - che uso z_covave_BS_BS (v. programma correlator_new per significato)


#%%
# 0) librerie e definizioni

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import scipy.odr as odr


def f(A, x):
    return A[0]*x+A[1]

now_path = '../GROMACS/'
save_path = now_path+'CORRELATIONS_postLM/'
treshold = 9
print('\n\n\nTreshold = 9 Angstrom, chosen as the best one\n\n\n')

WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6','wtc7')
Kds = [4, 4, 1360, 650, 750, 1320, 1800, 70]
Kds_errs = [0.9, 0.9, 600, 135, 150, 350, 400, 50]

Kds_new =  [4, 4, 1360, 650, 750, 1320, 1800, 70]

colori = ['royalblue', 'cornflowerblue', 'forestgreen', 'goldenrod', 'orange', 'darkorchid', 'firebrick',  'black']
colors = {wtc : color for (wtc, color) in zip(WTC_identifier, colori)}

print('\n\n\nStiamo procedendo con {} dinamiche\n\n\n'.format(len(Kds)))

red = False
wtc7_7 = False
wtc7_8 = False
wtc7_16 = True
wtc7_24 = False
wtc7_24_1 = False

#%%
# 1)            CORRELAZIONE

#1) importo dati dal matricione delle dinamiche e costruisco variabili correlazione

if red:
    print('Sto aprendo il red')
    DF = pd.read_json(now_path+'WTC_data_frame_{}ang_red.json'.format(str(treshold)))
elif wtc7_7:
    wtc7=7
    print('Sto aprendo il 7')

    DF = pd.read_json(now_path+'WTC_data_frame_{}ang_7.json'.format(str(treshold)))
elif wtc7_8:
    wtc7=8
    print('Sto aprendo il 8')
    DF = pd.read_json(now_path+'WTC_data_frame_{}ang_8.json'.format(str(treshold)))
elif wtc7_16:
    wtc7=16
    print('Sto aprendo il 16')
    DF = pd.read_json(now_path+'WTC_data_frame_{}ang_16_all.json'.format(str(treshold)))
elif wtc7_24:
    wtc7=24
    print('Sto aprendo il 24')
    DF = pd.read_json(now_path+'WTC_data_frame_{}ang_24_all.json'.format(str(treshold)))
elif wtc7_24_1:
    wtc7=25
    print('Sto aprendo il 24_1')
    DF = pd.read_json(now_path+'WTC_data_frame_{}ang_24_1_all.json'.format(str(treshold)))
else: 
    wtc7='tot'
    print('Sto aprendo il tot')
    DF = pd.read_json(now_path+'WTC_data_frame_{}ang.json'.format(str(treshold)))


z_covave_rmsf_BS_BS = {key : DF.z_Covariance_Mean_Prot_BS[DF.WTC_identifier == key][DF.Is_BS == True][DF.Is_Prot == True] for key in WTC_identifier }
z_covave_rmsf_BS_BS_correlators = [np.mean(z_covave_rmsf_BS_BS[key]) for key in WTC_identifier]
z_covave_rmsf_BS_BS_correlations = pearsonr(Kds[1:], z_covave_rmsf_BS_BS_correlators[1:])


z_covave_rmsf_BS_BS_12 = {key : DF.z_Covariance_Mean_Prot_BS_12[DF.WTC_identifier == key][DF.Is_BS == True][DF.Is_Prot == True] for key in WTC_identifier }
z_covave_rmsf_BS_BS_correlators_12 = [np.mean(z_covave_rmsf_BS_BS_12[key]) for key in WTC_identifier]
z_covave_rmsf_BS_BS_correlations_12 = pearsonr(Kds[1:], z_covave_rmsf_BS_BS_correlators[1:])

#%%
#FIT

linear = odr.Model(f)

# ERRORI presi sui 2 punti che dovrebbero essere uguali,
#ciclo sugli step per calcolarli

to_exclude = [3,4]

max_err = np.max([z_covave_rmsf_BS_BS_correlators[jj] for jj in to_exclude])
min_err = np.min([z_covave_rmsf_BS_BS_correlators[jj] for jj in to_exclude])
err_max = (max_err - min_err)
err = err_max/3 #procedura standard
print('Errore stimato = 0', err/2)


correlators = z_covave_rmsf_BS_BS_correlators_12
correlations = z_covave_rmsf_BS_BS_correlations_12



#STAMPO I PUNTI

fig, ax = plt.subplots()
ax.set_title('Protein BS Covariance with Protein BS  vs Kd\ncorrelation = {:3.2f} p-value = {:3.2f}'.format(*pearsonr(Kds[1:], correlators[1:])))
ax.errorbar(Kds[1:] , correlators[1:], fmt = 'o', xerr = Kds_errs[1:], yerr= err_max, color = 'k', ecolor = 'orange', label = 'simulated data', mew = 0.1)
#ax.errorbar(Kds[0] , correlators[0], fmt = 'o', xerr = Kds_errs[0], color = 'orange', ecolor = 'k', label = 'NMR data', mew = 0.1)

for x,y, key in zip(Kds[1:] , correlators[1:], ('wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6','wtc7')):

    if key in ['wtc1_h']:
        key = 'wtc1'
        plt.annotate(key, (x,y), xytext = (5, -12), textcoords = 'offset points', )
    elif key == 'wtc3': plt.annotate(key, (x,y), xytext = (-22, -15), textcoords = 'offset points', )
    elif key == 'wtc1': continue
    elif key == 'wtc7': 
        key = 'wtc7_'+str(wtc7)
        plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points', )

    else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')
    
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
#ax.vlines(np.mean([Kds[kk] for kk in to_exclude]), min_err, max_err, linestyles = 'dashed', color = 'firebrick', label = 'max error bar', linewidths = 2.)
ax.legend()
fig.savefig(save_path+'corr_plot'+str(wtc7)+'_all.pdf', format = 'pdf')
plt.show()


#FACCIO E STAMPO IL FIT

mydata = odr.RealData(Kds[1:], correlators[1:], sy = err/2, sx = Kds_errs[1:])
myodr = odr.ODR(mydata, linear, beta0=[1.,2.])
myoutput = myodr.run()

fig, ax = plt.subplots()
x = np.linspace(np.min(Kds), np.max(np.array(Kds)+np.array(Kds_errs)), 1000)
ax.set_title(r'$\bf{y = mx + q}$ fit'+'\nProtein BS Covariance with Protein BS vs Kd'+'\ncorrelation = {:3.2f} p-value = {:3.2f}'.format(*pearsonr(Kds[1:], correlators[1:])))
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('Prot BS z Covariance with Prot BS ')
ax.errorbar(Kds[1:] , correlators[1:], fmt = 'o', xerr = Kds_errs[1:], yerr = err, color = 'k', ecolor = 'orange', label = 'correlation data', mew = 0.1)
ax.plot(x,f(myoutput.beta, x), color = 'yellowgreen', label = 'linear fit')
ax.legend(title = 'm = {:3.2e} $\pm$ {:3.2e}\nq = {:3.2e} $\pm$ {:3.2e}'.format(myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1]))
ax.set_ylim(0, 7)
plt.show()
fig.savefig(save_path+'corr_fit'+str(wtc7)+'_all.pdf', format = 'pdf')




# %%
# 2) CAP_DISTANCES

import scipy.odr as odr
def f_lin(A, x):
    return A[0]*x+A[1]

linear = odr.Model(f_lin)



if red:
    df = pd.read_json('../GROMACS/df_CAPdist_new_red.json')
elif wtc7_7:
    wtc7=7
    df = pd.read_json('../GROMACS/df_CAPdist_new_7.json')
elif wtc7_8:
    wtc7=8
    df = pd.read_json('../GROMACS/df_CAPdist_new_8.json')
elif wtc7_16:
    wtc7=16
    df = pd.read_json('../GROMACS/df_CAPdist_new_16_all.json')
elif wtc7_24:
    wtc7=24
    df = pd.read_json('../GROMACS/df_CAPdist_new_24_all.json')
elif wtc7_24_1:
    wtc7=25
    df = pd.read_json('../GROMACS/df_CAPdist_new_24_1_all.json')
else:
    wtc7='tot'
    df = pd.read_json('../GROMACS/df_CAPdist_new.json')


df = df.rename(index = { old : new for old, new in zip(range(8), WTC_identifier)})

treshold = '9 ang'
for to_exclude in [[3,4]]:#, [2,5]]:

    print(" Treshold = {}  \nPunti con stessa fisica = {}".format(treshold, to_exclude))

    max_err = np.max([df[treshold][jj] for jj in to_exclude])
    min_err = np.min([df[treshold][jj] for jj in to_exclude])
    err_max = (max_err - min_err)
    err = err_max/np.sqrt(2) #procedura standard
    print(err/2)

    f, ax = plt.subplots()
    ax.set_title(r'$C_\alpha$ - $P$ average min distance vs Kd'+'\nBS treshold = {} $\AA$\npearson = {:3.2f} p-value = {:3.2f}'.format(treshold[:2], *pearsonr(Kds_new[1:], [df[treshold][key]for key in WTC_identifier][1:] )))
    ax.errorbar(Kds_new[1:], [df[treshold][key] for key in WTC_identifier][1:], yerr = err/2, xerr = Kds_errs[1:], fmt = 'o',  color = 'green', ecolor = 'firebrick',mew = 0.1, label = 'simulated data')
    #ax.errorbar(Kds_new[0], [df[treshold][key] for key in WTC_identifier][0], xerr = Kds_errs[0], fmt = 'o',  color = 'magenta', ecolor = 'green', label = 'NMR', mew = 0.1)
    #ax.vlines(np.mean([Kds_new[kk] for kk in to_exclude]), min_err, max_err, linestyles = 'dashed', color = 'yellowgreen', label = 'max error bar', linewidths = 2.)
    for x, key in zip(Kds_new, WTC_identifier):
        y = df[treshold][key]
        if key in ['wtc1_h']:
            key = 'wtc1'
            plt.annotate(key, (x,y), xytext = (5, 12), textcoords = 'offset points', )
        elif key == 'wtc1': continue
        elif key in ['wtc3', 'wtc6']:plt.annotate(key, (x,y), xytext = (-25, -15), textcoords = 'offset points', )
        elif key == 'wtc7': 
                key = 'wtc7_'+str(wtc7)
                plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )


        else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')
    ax.legend()
    ax.set_xlabel('Kd (nM)')
    ax.set_ylabel('CA-P Mean Dist (Ang)')
    f.savefig(save_path+'CAP_dist_plot_'+str(wtc7)+'_all.pdf', format = 'pdf')


    mydata = odr.RealData(Kds_new[1:], [df[treshold][key] for key in WTC_identifier][1:], sy = err/2, sx = Kds_errs[1:])
    myodr = odr.ODR(mydata, linear, beta0=[1.,2.])
    myoutput = myodr.run()  

    fig, ax = plt.subplots()
    x = np.linspace(np.min(Kds_new), np.max(np.array(Kds_new)+np.array(Kds_errs)), 1000)
    ax.set_title(r'$\bf{y = mx + q}$ fit '+r' $C_\alpha$ - $P$ average min distance vs Kd'+'\nBS treshold = {} $\AA$\npearson = {:3.2f} p-value = {:3.2f}'.format(treshold[:2], *pearsonr(Kds_new[1:],[df[treshold][key] for key in WTC_identifier][1:])))
    ax.set_xlabel('Kd (nM)')
    ax.set_ylabel('Prot BS z Covariance with Prot BS ')
    ax.errorbar(Kds_new[1:], [df[treshold][key] for key in WTC_identifier][1:], xerr = Kds_errs[1:], yerr = err/2, fmt = 'o',  color = 'green', ecolor = 'magenta',mew = 0.1, label = 'simulated data')

    ax.plot(x,f_lin(myoutput.beta, x), color = 'yellowgreen', label = 'linear fit')
    ax.legend(title = 'm = {:3.2e} $\pm$ {:3.2e}\nq = {:3.2e} $\pm$ {:3.2e}'.format(myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1]))
    plt.show()

# %%
# SALVO DATI    



covars = correlators[1:]
dists = [df[treshold][key] for key in WTC_identifier][1:]

df_cov = pd.DataFrame(index=WTC_identifier[1:])
df_cov['Kd'] = [Kd for Kd in Kds_new[1:]]
df_cov['Kd_errs'] = [Kd_err for Kd_err in Kds_errs[1:]]
df_cov['Cov'] = [cov for cov in covars]
df_cov['cov_errs'] = [err_max]*len(covars)
df_cov['dist'] = [dist for dist in dists]
df_cov['dist_errs'] = [err/2]*len(covars)

if red:
    df_cov.to_csv('../GROMACS/df_final.csv')
elif wtc7_7:
    wtc7=7
    df_cov.to_csv('../GROMACS/df_final_7.csv')
elif wtc7_8:
    wtc7=8
    df_cov.to_csv('../GROMACS/df_final_8.csv')
elif wtc7_16:
    wtc7=16
    df_cov.to_csv('../GROMACS/df_final_16_all.csv')
elif wtc7_24:
    wtc7=24
    df_cov.to_csv('../GROMACS/df_final_24_all.csv')
elif wtc7_24_1:
    wtc7=25
    df_cov.to_csv('../GROMACS/df_final_24_1_all.csv')
else:
    wtc7='tot'
    df_cov.to_csv('../GROMACS/df_final.csv')

#%%
# FACCIO UN GRAFICO UNICO a distribuzione normalizzate


z_covars = (covars - np.mean(covars))/np.std(covars)
f, ax = plt.subplots()
ax.set_title('Normalized Cov vs Kd'+'\nBS treshold = {} $\AA$\npearson = {:3.2f} p-value = {:3.2f}'.format(treshold[:2], *pearsonr(Kds_new[1:], z_covars )))
ax.errorbar(Kds_new[1:], z_covars, fmt = 'o',  color = 'black', ecolor = 'firebrick',mew = 0.1, label = 'simulated data')

for x, key, ii in zip(Kds_new[1:], WTC_identifier[1:], range(len(Kds_new[1:]))):
    y = z_covars[ii]
    if key in ['wtc1_h']:
        key = 'wtc1'
        plt.annotate(key, (x,y), xytext = (5, 12), textcoords = 'offset points', )
    elif key == 'wtc1': continue
    elif key in ['wtc3', 'wtc6']:plt.annotate(key, (x,y), xytext = (-25, -15), textcoords = 'offset points', )
    elif key == 'wtc7': 
            key = 'wtc7_'+str(wtc7)
            plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )


    else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')
ax.legend()
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('z_score covariance')
plt.tight_layout()

# %%

z_dists = (dists - np.mean(dists))/np.std(dists)

f, ax = plt.subplots()
ax.set_title('Normalized Ca-P dist vs Kd'+'\nBS treshold = {} $\AA$\npearson = {:3.2f} p-value = {:3.2f}'.format(treshold[:2], *pearsonr(Kds_new[1:], z_dists)))
ax.errorbar(Kds_new[1:], z_dists, fmt = 'o',  color = 'green', ecolor = 'firebrick',mew = 0.1, label = 'simulated data')

for x, key, ii in zip(Kds_new[1:], WTC_identifier[1:], range(len(Kds_new[1:]))):
    y = z_dists[ii]
    if key in ['wtc1_h']:
        key = 'wtc1'
        plt.annotate(key, (x,y), xytext = (5, 12), textcoords = 'offset points', )
    elif key == 'wtc1': continue
    elif key in ['wtc3', 'wtc6']:plt.annotate(key, (x,y), xytext = (-25, -15), textcoords = 'offset points', )
    elif key == 'wtc7': 
            key = 'wtc7_'+str(wtc7)
            plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )


    else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')
ax.legend()
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('z_score dist')
plt.tight_layout()


#%%
mean_overall = (z_covars + z_dists)/2

f, ax = plt.subplots()
ax.set_title('Covariance and Ca-P Dist Normalized Mean vs Kd'+'\nBS treshold = {} $\AA$\npearson = {:3.2f} p-value = {:3.2f}'.format(treshold[:2], *pearsonr(Kds_new[1:], mean_overall)))
ax.errorbar(Kds_new[1:], mean_overall, fmt = 'o',  color = 'firebrick', ecolor = 'firebrick',mew = 0.1, label = 'simulated data')

for x, key, ii in zip(Kds_new[1:], WTC_identifier[1:], range(len(Kds_new[1:]))):
    y = mean_overall[ii]
    if key in ['wtc1_h']:
        key = 'wtc1'
        plt.annotate(key, (x,y), xytext = (5, 12), textcoords = 'offset points', )
    elif key == 'wtc1': continue
    elif key in ['wtc3', 'wtc6']:plt.annotate(key, (x,y), xytext = (-25, -15), textcoords = 'offset points', )
    elif key == 'wtc7': 
            key = 'wtc7_'+str(wtc7)
            plt.annotate(key, (x,y), xytext = (5, -10), textcoords = 'offset points', )


    else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')
ax.legend()
ax.set_xlabel('Kd (nM)')
ax.set_ylabel('z_dist and z_cov mean')
plt.tight_layout()



# %%

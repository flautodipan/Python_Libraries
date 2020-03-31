#%%

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

now_path = '../GROMACS/'
DF = pd.read_json(now_path+'WTC_data_frame.json')
DF['z_RMSF'] = (DF.RMSF.values - np.mean(DF.RMSF.values))/np.std(DF.RMSF.values)

Kds = np.array([4., 4., 700., 1320., 1350., 1360., 1800.])
WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6')
colori = ['royalblue', 'cornflowerblue', 'forestgreen', 'goldenrod', 'orange', 'darkorchid', 'firebrick' ]
colors = {wtc : color for (wtc, color) in zip(WTC_identifier, colori)}

# %%
factor = 3.
#cutoffs

#z_RMSF
z_rmsf_cutoffs, z_rmsf_step = np.linspace(-1.5,3., 25, retstep = True )
#plt.hist(DF.z_RMSF[DF.Is_BS == True], bins = 50, rwidth=.8)
z_rmsf_BS_cutoffs = {}
z_rmsf_BS_step = {}

for key in WTC_identifier:
    a = DF.z_RMSF[DF.WTC_identifier == key][DF.Is_BS == True]
    z_rmsf_BS_cutoffs[key], z_rmsf_BS_step[key] = np.linspace(np.min(a), np.max(a), 25, retstep=True)

#z_rmsf_BS_cutoffs, z_rmsf_BS_step = np.linspace(-1.5, 3., 25, retstep=True)



#variabili di correlazione

z_rmsf_rmsf_correlators = []
z_rmsf_rmsf_correlations = []
z_rmsf_rmsf_population = []

z_rmsf_rmsf_BS_correlators = []
z_rmsf_rmsf_BS_correlations = []
z_rmsf_rmsf_BS_population = []


n_step = 0

for (z_rmsf_cutoff, z_rmsf_BS_cutoff) in zip(z_rmsf_cutoffs, pd.DataFrame(z_rmsf_BS_cutoffs).T.items()):
    print(n_step)
    #variabili interesse

    z_rmsf_rmsf = {}
    z_rmsf_rmsf_howmany = {}

    z_rmsf_rmsf_BS = {}
    z_rmsf_rmsf_BS_howmany = {}


    #¢alcolo
    for key in WTC_identifier:

        z_rmsf_rmsf[key]          = DF.z_RMSF[DF.WTC_identifier == key].values[DF.z_RMSF[DF.WTC_identifier == key].values >= z_rmsf_cutoff]
        z_rmsf_rmsf_howmany[key] = len(z_rmsf_rmsf[key])

        z_rmsf_rmsf_BS[key]          = DF.z_RMSF[DF.WTC_identifier == key][DF.Is_BS == True].values[DF.z_RMSF[DF.WTC_identifier == key][DF.Is_BS == True].values >= z_rmsf_BS_cutoff[1][key]]
        z_rmsf_rmsf_BS_howmany[key] = len(z_rmsf_rmsf_BS[key])

    #salvo 

    #zRMSF-RMSF
    z_rmsf_rmsf_correlators.append([np.mean(z_rmsf_rmsf[key]) for key in WTC_identifier])
    z_rmsf_rmsf_correlations.append(pearsonr(Kds, z_rmsf_rmsf_correlators[n_step]))
    z_rmsf_rmsf_population.append({ key : z_rmsf_rmsf_howmany[key] for key in WTC_identifier})
    #RMSF-RMSF-BS
    z_rmsf_rmsf_BS_correlators.append([np.mean(z_rmsf_rmsf_BS[key]) for key in WTC_identifier])
    z_rmsf_rmsf_BS_correlations.append(pearsonr(Kds, z_rmsf_rmsf_BS_correlators[n_step]))
    z_rmsf_rmsf_BS_population.append({ key : z_rmsf_rmsf_BS_howmany[key] for key in WTC_identifier})


    n_step+=1

#%%
#Dataframe popolazioni

z_rmsf_rmsf_pop = pd.DataFrame(z_rmsf_rmsf_population)
z_rmsf_rmsf_BS_pop = pd.DataFrame(z_rmsf_rmsf_BS_population)
z_rmsf_rmsf_BS_pop
#%%
steps = np.arange(0,n_step, 1)

#GRAFICO DATI

f, ax = plt.subplots()
ax.plot(steps, [corr[0] for corr in z_rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'z_RMSF_RMSF', color = 'b')
ax.plot(steps, [corr[0] for corr in z_rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'z_RMSF_BS_RMSF', color = 'r')


ax.legend()
plt.tight_layout(rect= (0,0,2,2))

#%%
#GRAFICO POPOLAZIONI

f1,ax1 = plt.subplots()
ax1.set_facecolor('aquamarine')
for n in range(n_step):
    [ax1.plot(n+1, z_rmsf_rmsf_population[n][key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]

ax1.legend([key for key in WTC_identifier])
ax1.set_xlabel('N steps', fontsize = 18)
ax1.set_ylabel('Population', fontsize = 18)

ax1.set_title('z_rmsf_rmsf population', fontsize= 25)
plt.tight_layout(rect= (0,0,2,2))

#%%
f2,ax2 = plt.subplots()
ax2.set_facecolor('indianred')
for n in range(n_step):
    [ax2.plot(n+1, z_rmsf_rmsf_BS_population[n][key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]

ax2.legend([key for key in WTC_identifier])
ax2.set_xlabel('N steps', fontsize = 18)
ax2.set_ylabel('Population', fontsize = 18)

ax2.set_title('z_rmsf_rmsf_BS population', fontsize= 25)
plt.tight_layout(rect= (0,0,2,2))


# %%
# ISTOGRAMMA
var = DF.z_RMSF

plt.figure()
plt.hist(var, bins = 50, rwidth=0.8, color = 'maroon')
plt.title('{}\nMax = {:4f}, Min = {:4f}'.format(var.name, np.max(var), np.min(var)))

#

#spunti
rmsf_cutoffs, rmsf_step = np.linspace(np.mean(DF.RMSF.values)-(factor*np.std(DF.RMSF.values)), np.mean(DF.RMSF.values) +(factor*np.std(DF.RMSF.values)), 25, retstep = True )
rmsf_BS_cutoffs, rmsf_BS_step= np.linspace(np.mean(DF.RMSF[DF.Is_BS == True].values)-(factor*np.std(DF.RMSF[DF.Is_BS == True].values)), np.mean(DF.RMSF[DF.Is_BS == True].values) +(factor*np.std(DF.RMSF[DF.Is_BS == True].values)), 25, retstep = True )


rmsf_rmsf_correlators = []
rmsf_rmsf_correlations = []


rmsf_BS_rmsf_correlators = []
rmsf_BS_rmsf_correlations = []

for key in WTC_identifier:

        #rmsf_rmsf[key]      = DF.RMSF[DF.WTC_identifier == key].values[DF.RMSF[DF.WTC_identifier == key].values >= rmsf_cutoff]
        #rmsf_BS_rmsf[key]   = DF.RMSF[DF.WTC_identifier == key][DF.Is_BS == True].values[DF.RMSF[DF.WTC_identifier == key][DF.Is_BS == True].values >= rmsf_BS_cutoff]
    #salvo 
    """
    #RMSF-RMSF
    rmsf_rmsf_correlators.append([np.mean(rmsf_rmsf[key]) for key in WTC_identifier])
    rmsf_rmsf_correlations.append(pearsonr(Kds, rmsf_rmsf_correlators[n_step]))
    #RMSF-RMSF-BS
    rmsf_BS_rmsf_correlators.append([np.mean(rmsf_BS_rmsf[key]) for key in WTC_identifier])
    rmsf_BS_rmsf_correlations.append(pearsonr(Kds, rmsf_BS_rmsf_correlators[n_step]))
    """

# %%

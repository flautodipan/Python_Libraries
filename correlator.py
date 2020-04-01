#%%

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

now_path = '../GROMACS/'
save_path = now_path+'CORRELATIONS/'
DF = pd.read_json(now_path+'WTC_data_frame.json')
#DF['z_RMSF'] = (DF.RMSF.values - np.mean(DF.RMSF.values))/np.std(DF.RMSF.values)

Kds = np.array([4., 4., 700., 1320., 1350., 1360., 1800.])
WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6')
colori = ['royalblue', 'cornflowerblue', 'forestgreen', 'goldenrod', 'orange', 'darkorchid', 'firebrick' ]
colors = {wtc : color for (wtc, color) in zip(WTC_identifier, colori)}
populations = ['rmsf', 'rmsf_BS', ]
#%%
#1) definisco cutoffs
numerosity = 20

######  RMSF
#totale
rmsf_min = np.min(DF.RMSF)
rmsf_max = 0.38#np.min([np.max(DF.RMSF[DF.WTC_identifier == key]) for key in WTC_identifier]) #il minore tra i massimi delle distribuzioni delle singole dinamiche
rmsf_cutoffs = np.linspace(rmsf_min, rmsf_max, numerosity)
#BS
rmsf_BS_min = np.min(DF.RMSF[DF.Is_BS == True])
rmsf_BS_max = 0.15#np.min([np.max(DF.RMSF[DF.Is_BS == True][DF.WTC_identifier == key]) for key in WTC_identifier])
rmsf_BS_cutoffs = np.linspace(rmsf_BS_min, rmsf_BS_max, numerosity)

#%%
#quantità
#   NOMENCLATURA
#   prima var = media seconda = condizione
#   rmsf_rmsf = media su rmsf condizione su rmsf
#   

### variabili correlazione
rmsf_rmsf_correlators = []
rmsf_rmsf_correlations = []

rmsf_rmsf_BS_correlators = []
rmsf_rmsf_BS_correlations = []

#popolazioni
population_min = 3
rmsf_populations = []
rmsf_BS_populations = []

#da svuotare
populations_step = []
rmsf = {}
rmsf_BS = {}

n_step = 0

for rmsf_cutoff, rmsf_BS_cutoff in zip(rmsf_cutoffs, rmsf_BS_cutoffs):

    stop = False
    n_step += 1


    rmsf = { key :DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    rmsf_BS = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_BS_cutoff][DF.Is_BS == True].values for key in WTC_identifier}
    
    #popolazioni
    rmsf_populations.append({ key : len(rmsf[key]) for key in WTC_identifier})
    rmsf_BS_populations.append({ key : len(rmsf_BS[key]) for key in WTC_identifier})

    count=0
    for pop, pop_name in zip([rmsf_populations[n_step-1], rmsf_BS_populations[n_step-1],], populations):

        populations_step = [pop[key] for key in WTC_identifier]
        check_idx = np.argmin(populations_step)

        if populations_step[check_idx] < population_min:
            print("{} Popolazione di {} con {} ha superato popolazione minima pari a {}\nInterrompo a step {} su {}".format(pop_name, WTC_identifier[check_idx], populations_step[check_idx], population_min, n_step-1, numerosity ))
            print("Valore dei cutoff attuali\n{} per rmsf \n{} per rmsf_BS\n".format(rmsf_cutoff, rmsf_BS_cutoff))
            stop = True
            break
        count += 1
    if stop: 
        for pop in [rmsf_populations, rmsf_BS_populations,]:
            pop.remove(pop[-1])
        break

    else:

        rmsf_rmsf_correlators.append([np.mean(rmsf[key]) for key in WTC_identifier])
        rmsf_rmsf_correlations.append(pearsonr(Kds, rmsf_rmsf_correlators[n_step-1]))

        rmsf_rmsf_BS_correlators.append([np.mean(rmsf_BS[key]) for key in WTC_identifier])
        rmsf_rmsf_BS_correlations.append(pearsonr(Kds, rmsf_rmsf_BS_correlators[n_step-1]))
        

#%%

rmsf_populations = pd.DataFrame(rmsf_populations)
rmsf_BS_populations = pd.DataFrame(rmsf_BS_populations)


#%%
#PLOTS

steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
f, ax = plt.subplots()
#plt.tight_layout(rect= (0,0,2,2))

ax.plot(steps, [corr[0] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf', color = 'b')
ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', color = 'r')

f.set_size_inches(11, 8)
ax.legend()
plt.savefig(save_path+'BS_correlations.pdf', format = 'pdf', bbox_to_inches = 'tight')
plt.show()
plt.close()

#%%
#GRAFICO POPOLAZIONI
till = len(steps)
for pop, pop_name, cutoffs in zip([rmsf_populations, rmsf_BS_populations,], populations, (rmsf_cutoffs[:till], rmsf_BS_cutoffs[:till]), ):
    
    f1,ax1 = plt.subplots()
    f1.set_size_inches(11, 8)
    [ax1.plot(steps, pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]
    ax1.legend([key for key in WTC_identifier])

    ax2 = ax1.twiny()
    [ax2.plot(cutoffs, pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]


    ax2.set_xlabel('{} cutoff (nm)'.format(pop_name), fontsize = 16, labelpad = 20)
    ax2.xaxis.set_label_position('top')    
    ax1.set_xlabel('N steps', fontsize = 16, labelpad = 20)
    ax1.xaxis.set_label_position('bottom')

    ax1.set_ylabel('Population', fontsize = 18)
    ax1.set_title(pop_name+' Population', fontsize= 25, pad = 60)

    plt.tight_layout()
    plt.savefig(save_path+pop_name+'_population.pdf', format = 'pdf', bbox_to_inches = 'tight')
    plt.show()
    plt.close()
    




##################  ##  ######  ######################
##############  ##  ##  ##  ##  ##  ###### #######  ##
#########   ##  ##  ##  ##  ##  ##  ##  ## ##   ##  ##
####    ##  ##  ##  ##  ##  ##  ##  ##  ## ##   ##  ##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## ##   ##  ##
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## #### ##  ##
####    ##  ##  ##  ##  ##  ##  ##  ##  ##      ##  ##
########    ##  ##  ##  ##  ##  ##  ##  ##########  ##
#############   ##  ##  ##  ##  ##  ##              ##
################    ######  ######  ##################



#%%
#VARIE ed EVENTUALI

#istogramma completo

_ =[plt.hist(DF.RMSF[DF.WTC_identifier == key].values, bins = 6, label = key, rwidth=.8, histtype = 'stepfilled', alpha = 0.5, color = colors[key])for key in WTC_identifier]
plt.legend()


_ =[plt.hist(DF.RMSF[DF.WTC_identifier == key][DF.Is_BS == True].values, bins = 6, label = key, rwidth=.8, histtype = 'stepfilled', alpha = 0.5, color = colors[key])for key in WTC_identifier]
plt.legend()
#%%
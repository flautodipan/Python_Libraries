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
#RNA
rmsf_RNA_min = np.min(DF.RMSF[DF.Is_Prot == False])
rmsf_RNA_max = 0.2
rmsf_RNA_cutoffs = np.linspace(rmsf_RNA_min, rmsf_RNA_max, numerosity)


######  covave = mean (cov)
#totale
covave_min = np.min(DF.Covariance_Mean)
covave_max = np.min([np.max(DF.Covariance_Mean[DF.WTC_identifier == key]) for key in WTC_identifier])
covave_cutoffs = np.linspace(covave_min, covave_max, numerosity)
#BS
covave_BS_min = np.min(DF.Covariance_Mean[DF.Is_BS == True])
covave_BS_max = np.min([np.max(DF.Covariance_Mean[DF.WTC_identifier == key][DF.Is_BS == True]) for key in WTC_identifier])
covave_BS_cutoffs = np.linspace(covave_BS_min, covave_BS_max, numerosity)
#RNA
covave_RNA_min = np.min(DF.Covariance_Mean[DF.Is_Prot == False])
covave_RNA_max = np.min([np.max(DF.Covariance_Mean[DF.WTC_identifier == key][DF.Is_Prot == False]) for key in WTC_identifier])
covave_RNA_cutoffs = np.linspace(covave_min, covave_max, numerosity)

######  covstd = std(cov)
#totale
covave_min = np.min(DF.Covariance_Mean)
covave_max = np.min([np.max(DF.Covariance_Mean[DF.WTC_identifier == key]) for key in WTC_identifier])
covave_cutoffs = np.linspace(covave_min, covave_max, numerosity)
#BS
covave_BS_min = np.min(DF.Covariance_Mean[DF.Is_BS == True])
covave_BS_max = np.min([np.max(DF.Covariance_Mean[DF.WTC_identifier == key][DF.Is_BS == True]) for key in WTC_identifier])
covave_BS_cutoffs = np.linspace(covave_BS_min, covave_BS_max, numerosity)
#RNA
covave_RNA_min = np.min(DF.Covariance_Mean[DF.Is_Prot == False])
covave_RNA_max = np.min([np.max(DF.Covariance_Mean[DF.WTC_identifier == key][DF.Is_Prot == False]) for key in WTC_identifier])
covave_RNA_cutoffs = np.linspace(covave_min, covave_max, numerosity)

#%%
#quantità
#   NOMENCLATURA
#   prima var = media (da che descrittore viene) seconda = condizione
#   rmsf_rmsf = media su rmsf (quindi tutti i valori di RMSF tali che vale )condizione su rmsf
#   

### variabili correlazione
#rmsf-rmsf
rmsf_rmsf_correlators = []
rmsf_rmsf_correlations = []

rmsf_rmsf_BS_correlators = []
rmsf_rmsf_BS_correlations = []

rmsf_rmsf_RNA_correlators = []
rmsf_rmsf_RNA_correlations = []

#rmsf-covave
rmsf_covave_correlators = []
rmsf_covave_correlations = []

rmsf_covave_BS_correlators = []
rmsf_covave_BS_correlations = []

rmsf_covave_RNA_correlators = []
rmsf_covave_RNA_correlations = []


#popolazioni
population_min = 3
populations = ['rmsf_RNA_covave', 'rmsf_BS_covave', 'rmsf_covave', 'rmsf_rmsf', 'rmsf_rmsf_BS', 'rmsf_rmsf_RNA' ]

#rmsf-rmsf
rmsf_rmsf_populations = []
rmsf_rmsf_BS_populations = []
rmsf_rmsf_RNA_populations = []

#rmsf-mean(cov)
rmsf_covave_populations = []
rmsf_covave_BS_populations = []
rmsf_covave_RNA_populations = []


#initialize
populations_step = []
n_step = 0

for covave_RNA_cutoff, covave_BS_cutoff, rmsf_cutoff, rmsf_BS_cutoff, rmsf_RNA_cutoff, covave_cutoff in zip(covave_RNA_cutoffs, covave_BS_cutoffs, rmsf_cutoffs, rmsf_BS_cutoffs, rmsf_RNA_cutoffs, covave_cutoffs):

    stop = False
    n_step += 1

    #variabili da RMSF-rmsf
    rmsf_rmsf = { key :DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    rmsf_rmsf_BS = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_BS_cutoff][DF.Is_BS == True].values for key in WTC_identifier}
    rmsf_rmsf_RNA = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_RNA_cutoff][DF.Is_Prot == False].values for key in WTC_identifier}
    #variabili da RMSF-covave
    rmsf_covave = { key :DF.RMSF[DF.WTC_identifier == key][DF.Covariance_Mean > covave_cutoff].values for key in WTC_identifier}
    rmsf_covave_BS = { key : DF.RMSF[DF.WTC_identifier == key][DF.Covariance_Mean > covave_BS_cutoff][DF.Is_BS == True].values for key in WTC_identifier}
    rmsf_covave_RNA = { key : DF.RMSF[DF.WTC_identifier == key][DF.Covariance_Mean > covave_RNA_cutoff][DF.Is_Prot == False].values for key in WTC_identifier}
   
    
    #popolazioni
    rmsf_rmsf_populations.append({ key : len(rmsf_rmsf[key]) for key in WTC_identifier})
    rmsf_rmsf_BS_populations.append({ key : len(rmsf_rmsf_BS[key]) for key in WTC_identifier})
    rmsf_rmsf_RNA_populations.append({ key : len(rmsf_rmsf_RNA[key]) for key in WTC_identifier})
    rmsf_covave_populations.append({ key : len(rmsf_covave[key]) for key in WTC_identifier})
    rmsf_covave_BS_populations.append({ key : len(rmsf_covave_BS[key]) for key in WTC_identifier})
    rmsf_covave_RNA_populations.append({ key : len(rmsf_covave_RNA[key]) for key in WTC_identifier})

    count=0
    for pop, pop_name in zip([rmsf_covave_RNA_populations[n_step-1],rmsf_covave_BS_populations[n_step-1], rmsf_covave_populations[n_step-1], rmsf_rmsf_populations[n_step-1], rmsf_rmsf_BS_populations[n_step-1], rmsf_rmsf_RNA_populations[n_step-1], ], populations):

        populations_step = [pop[key] for key in WTC_identifier]
        check_idx = np.argmin(populations_step)

        if populations_step[check_idx] < population_min:
            print("{} Popolazione di {} con {} ha superato popolazione minima pari a {}\nInterrompo a step {} su {}".format(pop_name, WTC_identifier[check_idx], populations_step[check_idx], population_min, n_step-1, numerosity ))
            print("Valore dei cutoff attuali per rmsf\n{} per tutto \n{} per BS\n{} per RNA\n\n".format(rmsf_cutoff, rmsf_BS_cutoff, rmsf_RNA_cutoff))
            print("Valore dei cutoff attualiper covave\n{} per tutto \n{} per BS\n{} per RNA\n\n".format(covave_cutoff, covave_BS_cutoff, covave_RNA_cutoff))

            stop = True
            break
        count += 1
    if stop: 
        for pop in [rmsf_covave_RNA_populations,rmsf_covave_BS_populations,rmsf_covave_populations, rmsf_rmsf_populations, rmsf_rmsf_BS_populations,rmsf_rmsf_RNA_populations,] :
            pop.remove(pop[-1])
        break

    else:


        #rmsf-rmsf
        rmsf_rmsf_correlators.append([np.mean(rmsf_rmsf[key]) for key in WTC_identifier])
        rmsf_rmsf_correlations.append(pearsonr(Kds, rmsf_rmsf_correlators[n_step-1]))

        rmsf_rmsf_BS_correlators.append([np.mean(rmsf_rmsf_BS[key]) for key in WTC_identifier])
        rmsf_rmsf_BS_correlations.append(pearsonr(Kds, rmsf_rmsf_BS_correlators[n_step-1]))

        rmsf_rmsf_RNA_correlators.append([np.mean(rmsf_rmsf_RNA[key]) for key in WTC_identifier])
        rmsf_rmsf_RNA_correlations.append(pearsonr(Kds, rmsf_rmsf_RNA_correlators[n_step-1]))
        
        #rmsf-covave

        rmsf_covave_correlators.append([np.mean(rmsf_covave[key]) for key in WTC_identifier])
        rmsf_covave_correlations.append(pearsonr(Kds, rmsf_covave_correlators[n_step-1]))

        rmsf_covave_BS_correlators.append([np.mean(rmsf_covave_BS[key]) for key in WTC_identifier])
        rmsf_covave_BS_correlations.append(pearsonr(Kds, rmsf_covave_BS_correlators[n_step-1]))

        rmsf_covave_RNA_correlators.append([np.mean(rmsf_covave_RNA[key]) for key in WTC_identifier])
        rmsf_covave_RNA_correlations.append(pearsonr(Kds, rmsf_covave_RNA_correlators[n_step-1]))
        

#%%

rmsf_rmsf_populations = pd.DataFrame(rmsf_rmsf_populations)
rmsf_rmsf_BS_populations = pd.DataFrame(rmsf_rmsf_BS_populations)
rmsf_rmsf_RNA_populations = pd.DataFrame(rmsf_rmsf_RNA_populations)

rmsf_covave_populations = pd.DataFrame(rmsf_covave_populations)
rmsf_covave_BS_populations = pd.DataFrame(rmsf_covave_BS_populations)
rmsf_covave_RNA_populations = pd.DataFrame(rmsf_covave_RNA_populations)

#%%
#PLOTS

#plot BS - RNA confronto
steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
f, ax = plt.subplots()
#plt.tight_layout(rect= (0,0,2,2))

ax.plot(steps, [corr[0] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf',)
ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', )
ax.plot(steps, [corr[0] for corr in rmsf_rmsf_RNA_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_RNA',)
ax.plot(steps, [corr[0] for corr in rmsf_covave_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covave',)
ax.plot(steps, [corr[0] for corr in rmsf_covave_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covave_BS', )
ax.plot(steps, [corr[0] for corr in rmsf_covave_RNA_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covave_RNA', )

f.set_size_inches(11, 8)
ax.legend()
plt.savefig(save_path+'BS_RNA_correlations.pdf', format = 'pdf', bbox_to_inches = 'tight')
plt.show()
plt.close()

#%%
#plot BS da solo

#plot BS - RNA confronto
steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
f, ax = plt.subplots()
#plt.tight_layout(rect= (0,0,2,2))

ax.plot(steps, [corr[0] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf',)
ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', )
ax.plot(steps, [corr[0] for corr in rmsf_covave_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covave',)
ax.plot(steps, [corr[0] for corr in rmsf_covave_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covave_BS', )


f.set_size_inches(11, 8)
ax.legend()
plt.savefig(save_path+'BS_correlations.pdf', format = 'pdf', bbox_to_inches = 'tight')
plt.show()
plt.close()


#%%
#GRAFICO POPOLAZIONI
till = len(steps)
for pop, pop_name, cutoffs in zip([rmsf_covave_RNA_populations, rmsf_covave_BS_populations, rmsf_covave_populations, rmsf_rmsf_populations, rmsf_rmsf_BS_populations, rmsf_rmsf_RNA_populations], populations, (covave_RNA_cutoffs, covave_BS_cutoffs, covave_cutoffs, rmsf_cutoffs, rmsf_BS_cutoffs, rmsf_RNA_cutoffs)):
    
    f1,ax1 = plt.subplots()
    f1.set_size_inches(11, 8)
    [ax1.plot(steps, pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]
    ax1.legend([key for key in WTC_identifier])

    ax2 = ax1.twiny()
    [ax2.plot(cutoffs[:till], pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]


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
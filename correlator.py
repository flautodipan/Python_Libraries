#%%

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

now_path = '../GROMACS/'
save_path = now_path+'CORRELATIONS/'
treshold = 10
DF = pd.read_json(now_path+'WTC_data_frame_{}ang.json'.format(str(treshold)))
DF_noBS_RNA = pd.concat([DF[DF.Is_BS == False], DF[DF.Is_Prot == False]], ignore_index= True)
DF_BS_RNA = pd.concat([DF[DF.Is_BS == True], DF[DF.Is_Prot == False]], ignore_index = True)

#DF['z_RMSF'] = (DF.RMSF.values - np.mean(DF.RMSF.values))/np.std(DF.RMSF.values)

Kds = np.array([4., 4., 700., 1320., 1350., 1360., 1800.])
Kds_sort = (np.array([4.,4., 1360, 1350, 1320, 700, 1800]))
WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6')
colori = ['royalblue', 'cornflowerblue', 'forestgreen', 'goldenrod', 'orange', 'darkorchid', 'firebrick' ]
colors = {wtc : color for (wtc, color) in zip(WTC_identifier, colori)}
#%%
#1) definisco cutoffs
numerosity = 20

######  RMSF
#totale
rmsf_min = np.min(DF.RMSF)
rmsf_max = 0.23#np.min([np.max(DF.RMSF[DF.WTC_identifier == key]) for key in WTC_identifier]) #il minore tra i massimi delle distribuzioni delle singole dinamiche
rmsf_cutoffs = np.linspace(rmsf_min, rmsf_max, numerosity)
#BS
rmsf_BS_min = np.min(DF.RMSF[DF.Is_BS == True])
rmsf_BS_max = 0.25#np.min([np.max(DF.RMSF[DF.Is_BS == True][DF.WTC_identifier == key]) for key in WTC_identifier])
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
covstd_min = np.min(DF.Covariance_Std)
covstd_max = 0.0018559#np.min([np.max(DF.Covariance_Std[DF.WTC_identifier == key]) for key in WTC_identifier])
covstd_cutoffs = np.linspace(covstd_min, covstd_max, numerosity)
#BS
covstd_BS_min = np.min(DF.Covariance_Std[DF.Is_BS == True])
covstd_BS_max = 0.00054#np.min([np.max(DF.Covariance_Std[DF.WTC_identifier == key][DF.Is_BS == True]) for key in WTC_identifier])
covstd_BS_cutoffs = np.linspace(covstd_BS_min, covstd_BS_max, numerosity)
#RNA
covstd_RNA_min = np.min(DF.Covariance_Std[DF.Is_Prot == False])
covstd_RNA_max = 0.00114#np.min([np.max(DF.Covariance_Std[DF.WTC_identifier == key][DF.Is_Prot == False]) for key in WTC_identifier])
covstd_RNA_cutoffs = np.linspace(covstd_RNA_min, covstd_RNA_max, numerosity)

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

rmsf_rmsf_BS_RNA_correlators = []
rmsf_rmsf_BS_RNA_correlations = []

#rmsf-covstd
rmsf_covstd_correlators = []
rmsf_covstd_correlations = []

rmsf_covave_correlators = []
rmsf_covave_correlations = []

#covave-rmsf

covave_rmsf_correlators = []
covave_rmsf_correlations = []

covave_rmsf_noBS_RNA_correlators = []
covave_rmsf_noBS_RNA_correlations = []

covave_rmsf_BS_RNA_correlators = []
covave_rmsf_BS_RNA_correlations = []

#covave-covstd

covave_covstd_correlators = []
covave_covstd_correlations = []


covave_covstd_BS_correlators = []
covave_covstd_BS_correlations = []


#covstd-covave
covstd_covave_correlators = []
covstd_covave_correlations = []


#covstd_rmsf

covstd_rmsf_correlators = []
covstd_rmsf_correlations = []

#popolazioni
population_min = 3
populations = ['rmsf_rmsf_BS_RNA', 'covave_covstd_BS','covstd_rmsf', 'covstd_covave','rmsf_covave', 'covave_covstd', 'covave_rmsf','covave_rmsf_noBS_RNA', 'covave_rmsf_BS_RNA', 'rmsf_covstd', 'rmsf_rmsf', 'rmsf_rmsf_BS', 'rmsf_rmsf_RNA' ]

#rmsf-rmsf
rmsf_rmsf_populations = []
rmsf_rmsf_BS_populations = []
rmsf_rmsf_RNA_populations = []
rmsf_rmsf_BS_RNA_populations = []


#rmsf-std(cov)
rmsf_covstd_populations = []
rmsf_covave_populations = []

#covave-rmsf
covave_rmsf_populations = []
covave_rmsf_noBS_RNA_populations = []
covave_rmsf_BS_RNA_populations = []
#covave-covstd
covave_covstd_populations = []
covave_covstd_BS_populations = []

#covstd - covave
covstd_covave_populations = []

#covstd-rmsf
covstd_rmsf_populations = []




#initialize
populations_step = []
n_step = 0

for rmsf_cutoff, rmsf_BS_cutoff, rmsf_RNA_cutoff, covstd_cutoff, covstd_BS_cutoff, covave_cutoff in zip(rmsf_cutoffs, rmsf_BS_cutoffs, rmsf_RNA_cutoffs, covstd_cutoffs,covstd_BS_cutoffs, covave_cutoffs):

    stop = False
    n_step += 1

    #variabili da RMSF-rmsf
    rmsf_rmsf = { key :DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    rmsf_rmsf_BS = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_BS_cutoff][DF.Is_BS == True].values for key in WTC_identifier}
    rmsf_rmsf_RNA = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_Prot == False].values for key in WTC_identifier}
    rmsf_rmsf_BS_RNA = {key : DF_BS_RNA.RMSF[DF_BS_RNA.WTC_identifier == key][DF_BS_RNA.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    #variabili da RMSF-covstd
    rmsf_covstd = { key :DF.RMSF[DF.WTC_identifier == key][DF.Covariance_Std > covstd_cutoff].values for key in WTC_identifier}
    #variabili da RMSF-covave
    rmsf_covave = { key :DF.RMSF[DF.WTC_identifier == key][DF.Covariance_Mean > covave_cutoff].values for key in WTC_identifier}
    #variavili da covave_rmsf
    covave_rmsf = {key : DF.Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    covave_rmsf_noBS_RNA = {key : DF_noBS_RNA.Covariance_Mean[DF_noBS_RNA.WTC_identifier == key][DF_noBS_RNA.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    covave_rmsf_BS_RNA = {key : DF_BS_RNA.Covariance_Mean[DF_BS_RNA.WTC_identifier == key][DF_BS_RNA.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    #variabili da covave-cvostd
    covave_covstd = {key : DF.Covariance_Mean[DF.WTC_identifier == key][DF.Covariance_Std > covstd_cutoff].values for key in WTC_identifier}
    covave_covstd_BS = {key : DF.Covariance_Mean[DF.WTC_identifier == key][DF.Is_BS == True][DF.Covariance_Std > covstd_BS_cutoff].values for key in WTC_identifier}
    #covstd-covave
    covstd_covave = {key : DF.Covariance_Std[DF.WTC_identifier == key][DF.Covariance_Mean > covave_cutoff].values for key in WTC_identifier}
    #covstd-rmsf
    covstd_rmsf = {key : DF.Covariance_Std[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}

    
    
    
    #popolazioni
    rmsf_rmsf_populations.append({ key : len(rmsf_rmsf[key]) for key in WTC_identifier})
    rmsf_rmsf_BS_populations.append({ key : len(rmsf_rmsf_BS[key]) for key in WTC_identifier})
    rmsf_rmsf_RNA_populations.append({ key : len(rmsf_rmsf_RNA[key]) for key in WTC_identifier})
    rmsf_rmsf_BS_RNA_populations.append({ key : len(rmsf_rmsf_BS_RNA[key]) for key in WTC_identifier})
    rmsf_covstd_populations.append({ key : len(rmsf_covstd[key]) for key in WTC_identifier})

    rmsf_covave_populations.append({ key : len(rmsf_covave[key]) for key in WTC_identifier})
    
    covave_rmsf_populations.append({ key : len(covave_rmsf[key]) for key in WTC_identifier})
    covave_rmsf_noBS_RNA_populations.append({ key : len(covave_rmsf_noBS_RNA[key]) for key in WTC_identifier})
    covave_rmsf_BS_RNA_populations.append({ key : len(covave_rmsf_BS_RNA[key]) for key in WTC_identifier})

    covave_covstd_populations.append({ key : len(covave_covstd[key]) for key in WTC_identifier})
    covave_covstd_BS_populations.append({ key : len(covave_covstd_BS[key]) for key in WTC_identifier})

    covstd_covave_populations.append({ key : len(covstd_covave[key]) for key in WTC_identifier})

    covstd_rmsf_populations.append({key : len(covstd_rmsf[key]) for key in WTC_identifier})


    count=0
    for pop, pop_name in zip([rmsf_rmsf_BS_RNA_populations[n_step -1], covave_covstd_BS_populations[n_step -1],covstd_rmsf_populations[n_step -1], covstd_covave_populations[n_step -1], rmsf_covave_populations[n_step -1], covave_covstd_populations[n_step-1], covave_rmsf_populations[n_step-1], covave_rmsf_noBS_RNA_populations[n_step-1],covave_rmsf_BS_RNA_populations[n_step-1], rmsf_covstd_populations[n_step-1], rmsf_rmsf_populations[n_step-1], rmsf_rmsf_BS_populations[n_step-1], rmsf_rmsf_RNA_populations[n_step-1], ], populations):

        populations_step = [pop[key] for key in WTC_identifier]
        check_idx = np.argmin(populations_step)

        if populations_step[check_idx] < population_min:
            print("{} Popolazione di {} con {} ha superato popolazione minima pari a {}\nInterrompo a step {} su {}".format(pop_name, WTC_identifier[check_idx], populations_step[check_idx], population_min, n_step-1, numerosity ))
            print("Valore dei cutoff attuali per rmsf\n{} per tutto \n{} per BS\n{} per RNA\n\n".format(rmsf_cutoff, rmsf_BS_cutoff, rmsf_RNA_cutoff))
            print("Valore dei cutoff attualiper covstd\n{} per tutto \n\n".format(covstd_cutoff,))
            print("Valore dei cutoff attualiper covave\n{} per tutto \n\n".format(covave_cutoff,))

            stop = True
            break
        count += 1
    if stop: 
        for pop in [rmsf_rmsf_BS_RNA_populations, covave_covstd_BS_populations, covstd_rmsf_populations, covstd_covave_populations, rmsf_covave_populations, covave_covstd_populations, covave_rmsf_populations, covave_rmsf_noBS_RNA_populations,covave_rmsf_BS_RNA_populations,rmsf_covstd_populations, rmsf_rmsf_populations, rmsf_rmsf_BS_populations,rmsf_rmsf_RNA_populations,] :
            pop.remove(pop[-1])
        break

    else:


        #rmsf-rmsf
        rmsf_rmsf_correlators.append([np.mean(rmsf_rmsf[key]) for key in WTC_identifier])
        rmsf_rmsf_correlations.append(pearsonr(Kds_sort, rmsf_rmsf_correlators[n_step-1]))

        rmsf_rmsf_BS_correlators.append([np.mean(rmsf_rmsf_BS[key]) for key in WTC_identifier])
        rmsf_rmsf_BS_correlations.append(pearsonr(Kds_sort, rmsf_rmsf_BS_correlators[n_step-1]))

        rmsf_rmsf_RNA_correlators.append([np.mean(rmsf_rmsf_RNA[key]) for key in WTC_identifier])
        rmsf_rmsf_RNA_correlations.append(pearsonr(Kds_sort, rmsf_rmsf_RNA_correlators[n_step-1]))

        rmsf_rmsf_BS_RNA_correlators.append([np.mean(rmsf_rmsf_BS_RNA[key]) for key in WTC_identifier])
        rmsf_rmsf_BS_RNA_correlations.append(pearsonr(Kds_sort, rmsf_rmsf_BS_RNA_correlators[n_step-1]))
        
        #rmsf-covstd

        rmsf_covstd_correlators.append([np.mean(rmsf_covstd[key]) for key in WTC_identifier])
        rmsf_covstd_correlations.append(pearsonr(Kds_sort, rmsf_covstd_correlators[n_step-1]))

        #rmsf_covave

        rmsf_covave_correlators.append([np.mean(rmsf_covave[key]) for key in WTC_identifier])
        rmsf_covave_correlations.append(pearsonr(Kds_sort, rmsf_covave_correlators[n_step-1]))

        #covave-rmsf

        covave_rmsf_correlators.append([np.mean(covave_rmsf[key]) for key in WTC_identifier])
        covave_rmsf_correlations.append(pearsonr(Kds_sort, covave_rmsf_correlators[n_step-1]))

        covave_rmsf_noBS_RNA_correlators.append([np.mean(covave_rmsf_noBS_RNA[key]) for key in WTC_identifier])
        covave_rmsf_noBS_RNA_correlations.append(pearsonr(Kds_sort, covave_rmsf_noBS_RNA_correlators[n_step-1]))
        
        covave_rmsf_BS_RNA_correlators.append([np.mean(covave_rmsf_BS_RNA[key]) for key in WTC_identifier])
        covave_rmsf_BS_RNA_correlations.append(pearsonr(Kds_sort, covave_rmsf_BS_RNA_correlators[n_step-1]))


        #covave-covstd
        covave_covstd_correlators.append([np.mean(covave_covstd[key]) for key in WTC_identifier])
        covave_covstd_correlations.append(pearsonr(Kds_sort, covave_covstd_correlators[n_step-1]))
        
        covave_covstd_BS_correlators.append([np.mean(covave_covstd_BS[key]) for key in WTC_identifier])
        covave_covstd_BS_correlations.append(pearsonr(Kds_sort, covave_covstd_BS_correlators[n_step-1]))

        #covstd covave
        covstd_covave_correlators.append([np.mean(covstd_covave[key]) for key in WTC_identifier])
        covstd_covave_correlations.append(pearsonr(Kds_sort, covstd_covave_correlators[n_step-1]))

        covstd_rmsf_correlators.append([np.mean(covstd_rmsf[key]) for key in WTC_identifier])
        covstd_rmsf_correlations.append(pearsonr(Kds_sort, covstd_rmsf_correlators[n_step-1]))


#%%

rmsf_rmsf_populations = pd.DataFrame(rmsf_rmsf_populations)
rmsf_rmsf_BS_populations = pd.DataFrame(rmsf_rmsf_BS_populations)
rmsf_rmsf_RNA_populations = pd.DataFrame(rmsf_rmsf_RNA_populations)
rmsf_rmsf_BS_RNA_populations = pd.DataFrame(rmsf_rmsf_BS_RNA_populations)

rmsf_covstd_populations = pd.DataFrame(rmsf_covstd_populations)
rmsf_covave_populations = pd.DataFrame(rmsf_covave_populations)


covave_rmsf_populations = pd.DataFrame(covave_rmsf_populations)
covave_rmsf_noBS_RNA_populations = pd.DataFrame(covave_rmsf_noBS_RNA_populations)
covave_rmsf_BS_RNA_populations = pd.DataFrame(covave_rmsf_BS_RNA_populations)

covave_covstd_populations = pd.DataFrame(covave_covstd_populations)

covave_covstd_BS_populations = pd.DataFrame(covave_covstd_BS_populations)

covstd_covave_populations = pd.DataFrame(covstd_covave_populations)

covstd_rmsf_populations = pd.DataFrame(covstd_rmsf_populations)


#%%
#PLOTS

#plot completo
steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
f, ax = plt.subplots()
#plt.tight_layout(rect= (0,0,2,2))

ax.plot(steps, [corr[0] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf',)
ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', )
#ax.plot(steps, [corr[0] for corr in rmsf_covstd_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covstd',)
#ax.plot(steps, [corr[0] for corr in rmsf_covave_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covave',)
#ax.plot(steps, [corr[0] for corr in covave_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf',)
#ax.plot(steps, [corr[0] for corr in covave_covstd_correlations], marker = 'o', ls = 'dashed', label = 'covave_covstd',)
#ax.plot(steps, [corr[0] for corr in covave_covstd_BS_correlations], marker = 'o', ls = 'dashed', label = 'covave_covstd_BS',)
ax.plot(steps, [corr[0] for corr in covstd_covave_correlations], marker = 'o', ls = 'dashed', label = 'covstd_covave',)
#ax.plot(steps, [corr[0] for corr in covstd_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'covstd_rmsf',)
ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_RNA_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS_RNA',)

ax.set_title('Pearson correlation with Kd increasing cutoff\nBS treshold = {} Ang'.format(treshold), fontsize = 20, pad = 18)
ax.set_xlabel('N steps', fontsize = 16)
ax.set_ylabel('Correlation', fontsize = 16)
ax.set_ylim(-1,1)

f.set_size_inches(11, 8)
ax.legend(loc = 'upper left', title = '<population>_<cutoff>')
plt.savefig(save_path+'complete_correlations.pdf', format = 'pdf', bbox_to_inches = 'tight')
plt.show()
plt.close()

#%%
# Per le variabili di intereses allungo i apssi
######  RMSF
#totale
rmsf_min = np.min(DF.RMSF)
rmsf_max = 0.239#np.min([np.max(DF.RMSF[DF.WTC_identifier == key]) for key in WTC_identifier]) #il minore tra i massimi delle distribuzioni delle singole dinamiche
rmsf_cutoffs = np.linspace(rmsf_min, rmsf_max, numerosity)
#BS
rmsf_BS_min = np.min(DF.RMSF[DF.Is_BS == True])
rmsf_BS_max = 0.1500#np.min([np.max(DF.RMSF[DF.Is_BS == True][DF.WTC_identifier == key]) for key in WTC_identifier])
rmsf_BS_cutoffs = np.linspace(rmsf_BS_min, rmsf_BS_max, numerosity)

rmsf_rmsf_correlators = []
rmsf_rmsf_correlations = []
rmsf_rmsf_BS_correlators = []
rmsf_rmsf_BS_correlations = []

covave_rmsf_noBS_RNA_correlators = []
covave_rmsf_noBS_RNA_correlations = []
covave_rmsf_BS_RNA_correlators = []
covave_rmsf_BS_RNA_correlations = []


rmsf_rmsf_populations = []
rmsf_rmsf_BS_populations = []
covave_rmsf_noBS_RNA_populations = []
covave_rmsf_BS_RNA_populations = []

populations_step = []
n_step = 0

for rmsf_cutoff, rmsf_BS_cutoff in zip(rmsf_cutoffs, rmsf_BS_cutoffs):

    stop = False
    n_step += 1

    rmsf_rmsf = { key :DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    rmsf_rmsf_BS = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_BS_cutoff][DF.Is_BS == True].values for key in WTC_identifier}
    covave_rmsf_noBS_RNA = {key : DF_noBS_RNA.Covariance_Mean[DF_noBS_RNA.WTC_identifier == key][DF_noBS_RNA.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    covave_rmsf_BS_RNA = {key : DF_BS_RNA.Covariance_Mean[DF_BS_RNA.WTC_identifier == key][DF_BS_RNA.RMSF > rmsf_cutoff].values for key in WTC_identifier}
    #popolazioni
    rmsf_rmsf_populations.append({ key : len(rmsf_rmsf[key]) for key in WTC_identifier})
    rmsf_rmsf_BS_populations.append({ key : len(rmsf_rmsf_BS[key]) for key in WTC_identifier})
    covave_rmsf_noBS_RNA_populations.append({ key : len(covave_rmsf_noBS_RNA[key]) for key in WTC_identifier})
    covave_rmsf_BS_RNA_populations.append({ key : len(covave_rmsf_BS_RNA[key]) for key in WTC_identifier})

    count=0
    for pop, pop_name in zip([covave_rmsf_noBS_RNA_populations[n_step-1],covave_rmsf_BS_RNA_populations[n_step-1], rmsf_rmsf_populations[n_step-1], rmsf_rmsf_BS_populations[n_step-1]], ['covave_rmsf_noBS_RNA', 'covave_rmsf_BS_RNA', 'rmsf_rmsf', 'rmsf_rmsf_BS']):

        populations_step = [pop[key] for key in WTC_identifier]
        check_idx = np.argmin(populations_step)

        if populations_step[check_idx] < population_min:
            print("{} Popolazione di {} con {} ha superato popolazione minima pari a {}\nInterrompo a step {} su {}".format(pop_name, WTC_identifier[check_idx], populations_step[check_idx], population_min, n_step-1, numerosity ))
            print("Valore dei cutoff attuali per rmsf\n{} per tutto \n{} per BS\n{} per RNA\n\n".format(rmsf_cutoff, rmsf_BS_cutoff, rmsf_RNA_cutoff))
            stop = True
            break

        count += 1

    if stop: 
        for pop in [covave_rmsf_noBS_RNA_populations,covave_rmsf_BS_RNA_populations,rmsf_rmsf_populations, rmsf_rmsf_BS_populations] :
            pop.remove(pop[-1])
        break

    else:


        #rmsf-rmsf
        rmsf_rmsf_correlators.append([np.mean(rmsf_rmsf[key]) for key in WTC_identifier])
        rmsf_rmsf_correlations.append(pearsonr(Kds_sort, rmsf_rmsf_correlators[n_step-1]))

        rmsf_rmsf_BS_correlators.append([np.mean(rmsf_rmsf_BS[key]) for key in WTC_identifier])
        rmsf_rmsf_BS_correlations.append(pearsonr(Kds_sort, rmsf_rmsf_BS_correlators[n_step-1]))
        covave_rmsf_noBS_RNA_correlators.append([np.mean(covave_rmsf_noBS_RNA[key]) for key in WTC_identifier])
        covave_rmsf_noBS_RNA_correlations.append(pearsonr(Kds_sort, covave_rmsf_noBS_RNA_correlators[n_step-1]))
        
        covave_rmsf_BS_RNA_correlators.append([np.mean(covave_rmsf_BS_RNA[key]) for key in WTC_identifier])
        covave_rmsf_BS_RNA_correlations.append(pearsonr(Kds_sort, covave_rmsf_BS_RNA_correlators[n_step-1]))




#%%
#PLOTS

#plot focus
steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
f, ax = plt.subplots()
#plt.tight_layout(rect= (0,0,2,2))

ax.plot(steps, [corr[0] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf',)
ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', )
ax.plot(steps, [corr[0] for corr in covave_rmsf_noBS_RNA_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_noBS_RNA', )
ax.plot(steps, [corr[0] for corr in covave_rmsf_BS_RNA_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_BS_RNA', )
ax.set_title('Pearson correlation with Kd further increasing cutoff\nBS treshold = {} Ang'.format(treshold), fontsize = 20, pad = 18)
ax.set_xlabel('N steps', fontsize = 16)
ax.set_ylabel('Correlation', fontsize = 16)
ax.set_ylim(-1,1)
f.set_size_inches(11, 8)
ax.legend(loc = 'lower left', title = '<population>_<cutoff>')
plt.savefig(save_path+'correlations_focus.pdf', format = 'pdf', bbox_to_inches = 'tight')
plt.show()
plt.close()


#plot focus - p-value
#%%
steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
f1, ax1 = plt.subplots()
#plt.tight_layout(rect= (0,0,2,2))

ax1.plot(steps, [corr[1] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf',)
ax1.plot(steps, [corr[1] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', )
ax1.plot(steps, [corr[1] for corr in covave_rmsf_noBS_RNA_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_noBS_RNA', )
ax1.plot(steps, [corr[1] for corr in covave_rmsf_BS_RNA_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_BS_RNA', )
ax1.hlines(0.1, steps[0], steps[-1], colors = 'grey', label = 'p = 0.1', ls = 'solid')
ax1.hlines(0.05,steps[0], steps[-1], colors = 'grey', label = 'p = 0.05', ls = 'dashed')
ax1.hlines(0.01,steps[0], steps[-1], colors = 'grey', label = 'p = 0.01', ls = 'dotted')

ax1.set_title('P - value for Pearson correlation with Kd', fontsize = 20)
ax1.set_xlabel('N steps', fontsize = 16)
ax1.set_ylabel('p-value', fontsize = 16)
ax1.set_ylim(0,1)


f1.set_size_inches(11, 8)
ax1.legend()
plt.savefig(save_path+'correlations_focus_pvalue.pdf', format = 'pdf', bbox_to_inches = 'tight')

ax1.set_ylim(0, 0.3)
plt.savefig(save_path+'correlations_focus_pvalue_zoom.pdf', format = 'pdf', bbox_to_inches = 'tight')

plt.show()
plt.close()






#%%
#GRAFICO POPOLAZIONI

rmsf_rmsf_populations = pd.DataFrame(rmsf_rmsf_populations)
rmsf_rmsf_BS_populations = pd.DataFrame(rmsf_rmsf_BS_populations)

covave_rmsf_noBS_RNA_populations = pd.DataFrame(covave_rmsf_noBS_RNA_populations)
covave_rmsf_BS_RNA_populations = pd.DataFrame(covave_rmsf_BS_RNA_populations)

# cutoff rmsf

pops = [covave_rmsf_noBS_RNA_populations, covave_rmsf_BS_RNA_populations, rmsf_rmsf_populations ]
pops_name = ['covave_rmsf_noBS_RNA', 'covave_rmsf_BS_RNA', 'rmsf_rmsf']
till = len(steps)


for pop, pop_name in zip(pops, pops_name):

    f1,ax1 = plt.subplots()
    f1.set_size_inches(11, 8)
    [ax1.plot(steps, pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]
    ax1.legend([key for key in WTC_identifier])

    ax2 = ax1.twiny()
    [ax2.plot(rmsf_cutoffs[:till], pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]


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

#%%
#cutoff rmsf_BS
pops = [rmsf_rmsf_BS_populations]
pops_name = ['rmsf_rmsf_BS']

for pop, pop_name in zip(pops, pops_name):

    f1,ax1 = plt.subplots()
    f1.set_size_inches(11, 8)
    [ax1.plot(steps, pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]
    ax1.legend([key for key in WTC_identifier])

    ax2 = ax1.twiny()
    [ax2.plot(rmsf_BS_cutoffs[:till], pop[key], 'o', ls = 'dotted', color = colors[key]) for key in WTC_identifier]


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




#cutoff covstd




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
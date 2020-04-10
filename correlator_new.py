#%%

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

now_path = '../GROMACS/'
save_path = now_path+'CORRELATIONS/'
treshold = 10
for treshold in range(8,13):
    
    print('\n\n\nTreshold = {}\n\n\n'.format(treshold))

    DF = pd.read_json(now_path+'WTC_data_frame_{}ang.json'.format(str(treshold)))
    DF_noBS_RNA = pd.concat([DF[DF.Is_BS == False], DF[DF.Is_Prot == False]], ignore_index= True)
    DF_BS_RNA = pd.concat([DF[DF.Is_BS == True], DF[DF.Is_Prot == False]], ignore_index = True)

    #DF['z_RMSF'] = (DF.RMSF.values - np.mean(DF.RMSF.values))/np.std(DF.RMSF.values)

    WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6')#'wtc5',
    Kds = [DF.Kd[DF.WTC_identifier == key].values[0] for key in WTC_identifier]
    colori = ['royalblue', 'cornflowerblue', 'forestgreen', 'goldenrod', 'orange', 'darkorchid', 'firebrick' ]
    colors = {wtc : color for (wtc, color) in zip(WTC_identifier, colori)}


    #%#%
    #1) definisco cutoffs
    numerosity = 25

    ######  RMSF
    #totale
    rmsf_min = np.min(DF.RMSF)
    rmsf_max = 0.25#np.min([np.max(DF.RMSF[DF.WTC_identifier == key]) for key in WTC_identifier]) #il minore tra i massimi delle distribuzioni delle singole dinamiche
    rmsf_cutoffs = np.linspace(rmsf_min, rmsf_max, numerosity)
    #BS
    rmsf_BS_min = np.min(DF.RMSF[DF.Is_BS == True])
    rmsf_BS_max = 0.21#np.min([np.max(DF.RMSF[DF.Is_BS == True][DF.WTC_identifier == key]) for key in WTC_identifier])
    rmsf_BS_cutoffs = np.linspace(rmsf_BS_min, rmsf_BS_max, numerosity)
    #RNA
    rmsf_RNA_min = np.min(DF.RMSF[DF.Is_Prot == False])
    rmsf_RNA_max = 0.2
    rmsf_RNA_cutoffs = np.linspace(rmsf_RNA_min, rmsf_RNA_max, numerosity)



    # %#%

    min_pop = 3

    rmsf_rmsf_correlators = []
    rmsf_rmsf_correlations = []
    rmsf_rmsf_populations = []


    rmsf_rmsf_BS_correlators = []
    rmsf_rmsf_BS_correlations = []
    rmsf_rmsf_BS_populations = []


    covave_rmsf_noBS_correlators = []
    covave_rmsf_noBS_correlations = []
    covave_rmsf_noBS_populations = []
    covave_rmsf_BS_correlators = []
    covave_rmsf_BS_correlations = []
    covave_rmsf_BS_populations = []

    populations_step = []
    n_step = 0
    stop = False

    for rmsf_cutoff, rmsf_BS_cutoff in zip(rmsf_cutoffs, rmsf_BS_cutoffs):

        n_step += 1

        #variabili da RMSF-rmsf
        rmsf_rmsf = { key :DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}
        rmsf_rmsf_populations.append({ key : len(rmsf_rmsf[key]) for key in WTC_identifier})

        rmsf_rmsf_BS = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_BS_cutoff][DF.Is_BS == True].values for key in WTC_identifier}
        rmsf_rmsf_BS_populations.append({ key : len(rmsf_rmsf_BS[key]) for key in WTC_identifier})

        covave_rmsf_noBS = {key : DF.Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == False] for key in WTC_identifier }
        covave_rmsf_noBS_populations.append({ key : len(covave_rmsf_noBS[key]) for key in WTC_identifier})

        covave_rmsf_BS = {key : DF.Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True] for key in WTC_identifier }
        covave_rmsf_BS_populations.append({ key : len(covave_rmsf_BS[key]) for key in WTC_identifier})

        pop_step = pd.DataFrame(index=WTC_identifier)
        pop_step['rmsf_rmsf'] = [rmsf_rmsf_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['rmsf_rmsf_BS'] = [rmsf_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['covave_rmsf_noBS'] = [covave_rmsf_noBS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['covave_rmsf_BS'] = [covave_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]

        for col in pop_step.columns:
            if np.min(pop_step[col].values) < min_pop:
                print('Popolazione {} di {} ha sforato il limite di {}'.format(pop_step.index[np.argmin(pop_step[col].values)], col, min_pop))
                print('Mi arresto allo step {} su {}'.format(n_step, numerosity))
                print("Valore dei cutoff attuali per rmsf\n{} per tutto \n{} per BS\n\n".format(rmsf_cutoff, rmsf_BS_cutoff))
                stop = True
                break
        
        if stop:
            rmsf_rmsf_populations.remove(rmsf_rmsf_populations[-1])
            rmsf_rmsf_BS_populations.remove(rmsf_rmsf_BS_populations[-1])
            covave_rmsf_noBS_populations.remove(covave_rmsf_noBS_populations[-1])
            covave_rmsf_BS_populations.remove(covave_rmsf_BS_populations[-1])
            break
        else:

            #rmsf-rmsf
            rmsf_rmsf_correlators.append([np.mean(rmsf_rmsf[key]) for key in WTC_identifier])
            rmsf_rmsf_correlations.append(pearsonr(Kds, rmsf_rmsf_correlators[n_step-1]))

            rmsf_rmsf_BS_correlators.append([np.mean(rmsf_rmsf_BS[key]) for key in WTC_identifier])
            rmsf_rmsf_BS_correlations.append(pearsonr(Kds, rmsf_rmsf_BS_correlators[n_step-1]))

            covave_rmsf_noBS_correlators.append([np.mean(covave_rmsf_noBS[key]) for key in WTC_identifier])
            covave_rmsf_noBS_correlations.append(pearsonr(Kds, covave_rmsf_noBS_correlators[n_step-1]))
            
            covave_rmsf_BS_correlators.append([np.mean(covave_rmsf_BS[key]) for key in WTC_identifier])
            covave_rmsf_BS_correlations.append(pearsonr(Kds, covave_rmsf_BS_correlators[n_step-1]))




    #%#%
    #PLOTS

    #plot completo
    steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
    f, ax = plt.subplots()
    #plt.tight_layout(rect= (0,0,2,2))

    ax.plot(steps, [corr[0] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf',)
    ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', )
    ax.plot(steps, [corr[0] for corr in covave_rmsf_noBS_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_noBS', )
    ax.plot(steps, [corr[0] for corr in covave_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_BS', )

    #ax.plot(steps, [corr[0] for corr in rmsf_covstd_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covstd',)
    #ax.plot(steps, [corr[0] for corr in rmsf_covave_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_covave',)
    #ax.plot(steps, [corr[0] for corr in covave_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf',)
    #ax.plot(steps, [corr[0] for corr in covave_covstd_correlations], marker = 'o', ls = 'dashed', label = 'covave_covstd',)
    #ax.plot(steps, [corr[0] for corr in covave_covstd_BS_correlations], marker = 'o', ls = 'dashed', label = 'covave_covstd_BS',)
    #ax.plot(steps, [corr[0] for corr in covstd_covave_correlations], marker = 'o', ls = 'dashed', label = 'covstd_covave',)
    #ax.plot(steps, [corr[0] for corr in covstd_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'covstd_rmsf',)
    #ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS',)

    ax.set_title('Pearson correlation with Kd increasing cutoff\nBS treshold = {} Ang'.format(treshold), fontsize = 20, pad = 18)
    ax.set_xlabel('N steps', fontsize = 16)
    ax.set_ylabel('Correlation', fontsize = 16)
    ax.set_ylim(-1,1)

    f.set_size_inches(11, 8)
    ax.legend(loc = 'upper left', title = '<population>_<cutoff>')
    plt.savefig(save_path+'complete_correlations_{}ang.pdf'.format(treshold), format = 'pdf', bbox_to_inches = 'tight')
    plt.show()
    plt.close()


# %%

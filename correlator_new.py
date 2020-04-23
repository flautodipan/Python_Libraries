#%%

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

now_path = '../GROMACS/'
save_path = now_path+'CORRELATIONS/'
for treshold in [9]:#range(6,13):
    
    print('\n\n\nTreshold = {}\n\n\n'.format(treshold))


    DF = pd.read_json(now_path+'WTC_data_frame_{}ang.json'.format(str(treshold)))
    DF['z_Covariance_Mean'] = (DF.Covariance_Mean - np.mean(DF.Covariance_Mean))/(np.std(DF.Covariance_Mean))
    DF['z_Covariance_Mean_Prot_BS'] = (DF.Covariance_Mean_Prot_BS - np.mean(DF.Covariance_Mean_Prot_BS))/(np.std(DF.Covariance_Mean_Prot_BS))
    DF['z_Covariance_Mean_RNA_BS'] = (DF.Covariance_Mean_RNA_BS - np.mean(DF.Covariance_Mean_RNA_BS))/(np.std(DF.Covariance_Mean_RNA_BS))

    #DF['z_RMSF'] = (DF.RMSF.values - np.mean(DF.RMSF.values))/np.std(DF.RMSF.values)

    WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6')#'wtc5',
    Kds = [DF.Kd[DF.WTC_identifier == key].values[0] for key in WTC_identifier]

    Kds = [4, 4, 1360, 650, 750, 1320, 1800]
    Kds_errs = [0.9, 0.9, 600, 135, 150, 350, 400]

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
    rmsf_BS_max = 0.19#np.min([np.max(DF.RMSF[DF.Is_BS == True][DF.WTC_identifier == key]) for key in WTC_identifier])
    rmsf_BS_cutoffs = np.linspace(rmsf_BS_min, rmsf_BS_max, numerosity)
    #RNA
    rmsf_RNA_min = np.min(DF.RMSF[DF.Is_Prot == False])
    rmsf_RNA_max = 0.2
    rmsf_RNA_cutoffs = np.linspace(rmsf_RNA_min, rmsf_RNA_max, numerosity)
    #covave_BS
    covave_RNA_min = np.min(DF.Covariance_Mean_RNA_BS)
    covave_RNA_max = 0.00015
    covave_RNA_cutoffs = np.linspace(covave_RNA_min, covave_RNA_max, numerosity)


    # %#%

    min_pop = 3

    rmsf_rmsf_correlators = []
    rmsf_rmsf_correlations = []
    rmsf_rmsf_populations = []


    rmsf_rmsf_BS_correlators = []
    rmsf_rmsf_BS_correlations = []
    rmsf_rmsf_BS_populations = []

    rmsf_rmsf_RNA_correlators = []
    rmsf_rmsf_RNA_correlations = []
    rmsf_rmsf_RNA_populations = []


    #covave

    pearson_rmsf_BS_BS_correlators = []
    pearson_rmsf_BS_BS_correlations = []
    pearson_rmsf_BS_BS_populations = []

    pearson_rmsf_RNA_RNA_correlators = []
    pearson_rmsf_RNA_RNA_correlations = []
    pearson_rmsf_RNA_RNA_populations = []
    
    #z-score covave
    z_covave_rmsf_noBS_correlators = []
    z_covave_rmsf_noBS_correlations = []
    z_covave_rmsf_noBS_populations = []

    z_covave_rmsf_BS_correlators = []
    z_covave_rmsf_BS_correlations = []
    z_covave_rmsf_BS_populations = []


    z_covave_rmsf_BS_BS_correlators = []
    z_covave_rmsf_BS_BS_correlations = []
    z_covave_rmsf_BS_BS_populations = []

    z_covave_rmsf_RNA_RNA_correlators = []
    z_covave_rmsf_RNA_RNA_correlations = []
    z_covave_rmsf_RNA_RNA_populations = []
    

    populations_step = []
    n_step = 0
    stop = False
    zero_pop = False

    for covave_RNA_cutoff,rmsf_cutoff, rmsf_BS_cutoff in zip(covave_RNA_cutoffs, rmsf_cutoffs, rmsf_BS_cutoffs):

        n_step += 1

        #variabili da RMSF-rmsf
        rmsf_rmsf = { key :DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff].values for key in WTC_identifier}
        rmsf_rmsf_populations.append({ key : len(rmsf_rmsf[key]) for key in WTC_identifier})

        rmsf_rmsf_BS = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_BS_cutoff][DF.Is_BS == True].values for key in WTC_identifier}
        rmsf_rmsf_BS_populations.append({ key : len(rmsf_rmsf_BS[key]) for key in WTC_identifier})

        rmsf_rmsf_RNA = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_Prot == False].values for key in WTC_identifier}
        rmsf_rmsf_RNA_populations.append({ key : len(rmsf_rmsf_RNA[key]) for key in WTC_identifier})

        z_covave_rmsf_noBS = {key : DF.z_Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == False] for key in WTC_identifier }
        z_covave_rmsf_noBS_populations.append({ key : len(z_covave_rmsf_noBS[key]) for key in WTC_identifier})

        z_covave_rmsf_BS = {key : DF.z_Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True][DF.Is_Prot == True] for key in WTC_identifier }
        z_covave_rmsf_BS_populations.append({ key : len(z_covave_rmsf_BS[key]) for key in WTC_identifier})

        z_covave_rmsf_BS_BS = {key : DF.z_Covariance_Mean_Prot_BS[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True][DF.Is_Prot == True] for key in WTC_identifier }
        z_covave_rmsf_BS_BS_populations.append({ key : len(z_covave_rmsf_BS_BS[key]) for key in WTC_identifier})
              
        z_covave_rmsf_RNA_RNA = {key : DF.z_Covariance_Mean_RNA_BS[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True][DF.Is_Prot == False] for key in WTC_identifier }
        z_covave_rmsf_RNA_RNA_populations.append({ key : len(z_covave_rmsf_RNA_RNA[key]) for key in WTC_identifier})

        pearson_rmsf_BS_BS = {key : DF.Pearson_Mean_Prot_BS[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True][DF.Is_Prot == True] for key in WTC_identifier }
        pearson_rmsf_BS_BS_populations.append({ key : len(pearson_rmsf_BS_BS[key]) for key in WTC_identifier})
       
        pearson_rmsf_RNA_RNA = {key : DF.Pearson_Mean_RNA_BS[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True][DF.Is_Prot == False] for key in WTC_identifier }
        pearson_rmsf_RNA_RNA_populations.append({ key : len(pearson_rmsf_RNA_RNA[key]) for key in WTC_identifier})
        
        pop_step = pd.DataFrame(index=WTC_identifier)

        pop_step['rmsf_rmsf'] = [rmsf_rmsf_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['rmsf_rmsf_BS'] = [rmsf_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['rmsf_rmsf_RNA'] = [rmsf_rmsf_RNA_populations[n_step-1][key] for key in WTC_identifier]

        pop_step['z_covave_rmsf_noBS'] = [z_covave_rmsf_noBS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['z_covave_rmsf_BS'] = [z_covave_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]

        pop_step['z_covave_rmsf_BS_BS'] = [z_covave_rmsf_BS_BS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['z_covave_rmsf_RNA_RNA'] = [z_covave_rmsf_RNA_RNA_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['pearson_rmsf_BS_BS'] = [pearson_rmsf_BS_BS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['pearson_rmsf_RNA_RNA'] = [pearson_rmsf_RNA_RNA_populations[n_step-1][key] for key in WTC_identifier]



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
            rmsf_rmsf_RNA_populations.remove(rmsf_rmsf_RNA_populations[-1])
            pearson_rmsf_BS_BS_populations.remove(pearson_rmsf_BS_BS_populations[-1])
            pearson_rmsf_RNA_RNA_populations.remove(pearson_rmsf_RNA_RNA_populations[-1])
            z_covave_rmsf_noBS_populations.remove(z_covave_rmsf_noBS_populations[-1])
            z_covave_rmsf_BS_populations.remove(z_covave_rmsf_BS_populations[-1])
            z_covave_rmsf_BS_BS_populations.remove(z_covave_rmsf_BS_BS_populations[-1])
            z_covave_rmsf_RNA_RNA_populations.remove(z_covave_rmsf_RNA_RNA_populations[-1])
            if n_step == 1: zero_pop = True
            break
        else:

            #rmsf-rmsf
            rmsf_rmsf_correlators.append([np.mean(rmsf_rmsf[key]) for key in WTC_identifier])
            rmsf_rmsf_correlations.append(pearsonr(Kds[1:], rmsf_rmsf_correlators[n_step-1][1:]))

            rmsf_rmsf_BS_correlators.append([np.mean(rmsf_rmsf_BS[key]) for key in WTC_identifier])
            rmsf_rmsf_BS_correlations.append(pearsonr(Kds[1:], rmsf_rmsf_BS_correlators[n_step-1][1:]))

            rmsf_rmsf_RNA_correlators.append([np.mean(rmsf_rmsf_RNA[key]) for key in WTC_identifier])
            rmsf_rmsf_RNA_correlations.append(pearsonr(Kds[1:], rmsf_rmsf_RNA_correlators[n_step-1][1:]))

            pearson_rmsf_BS_BS_correlators.append([np.mean(pearson_rmsf_BS_BS[key]) for key in WTC_identifier])
            pearson_rmsf_BS_BS_correlations.append(pearsonr(Kds[1:], pearson_rmsf_BS_BS_correlators[n_step-1][1:]))

            pearson_rmsf_RNA_RNA_correlators.append([np.mean(pearson_rmsf_RNA_RNA[key]) for key in WTC_identifier])
            pearson_rmsf_RNA_RNA_correlations.append(pearsonr(Kds[1:], pearson_rmsf_RNA_RNA_correlators[n_step-1][1:]))

            z_covave_rmsf_noBS_correlators.append([np.mean(z_covave_rmsf_noBS[key]) for key in WTC_identifier])
            z_covave_rmsf_noBS_correlations.append(pearsonr(Kds[1:], z_covave_rmsf_noBS_correlators[n_step-1][1:]))
            
            z_covave_rmsf_BS_correlators.append([np.mean(z_covave_rmsf_BS[key]) for key in WTC_identifier])
            z_covave_rmsf_BS_correlations.append(pearsonr(Kds[1:], z_covave_rmsf_BS_correlators[n_step-1][1:]))


            z_covave_rmsf_BS_BS_correlators.append([np.mean(z_covave_rmsf_BS_BS[key]) for key in WTC_identifier])
            z_covave_rmsf_BS_BS_correlations.append(pearsonr(Kds[1:], z_covave_rmsf_BS_BS_correlators[n_step-1][1:]))

            z_covave_rmsf_RNA_RNA_correlators.append([np.mean(z_covave_rmsf_RNA_RNA[key]) for key in WTC_identifier])
            z_covave_rmsf_RNA_RNA_correlations.append(pearsonr(Kds[1:], z_covave_rmsf_RNA_RNA_correlators[n_step-1][1:]))



    #%#%
    #PLOTS

    #plot completo
    steps = np.arange(n_step-1, dtype = int) if stop == True else np.arange(n_step, dtype = int)
    f, ax = plt.subplots()
    #plt.tight_layout(rect= (0,0,2,2))

    ax.plot(steps, [corr[0] for corr in rmsf_rmsf_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf',)
    ax.plot(steps, [corr[0] for corr in rmsf_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_BS', )

    #ax.plot(steps, [corr[0] for corr in covave_rmsf_noBS_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_noBS', )
    #ax.plot(steps, [corr[0] for corr in covave_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_BS', )
    ax.plot(steps, [corr[0] for corr in z_covave_rmsf_BS_BS_correlations], marker = 'o', ls = 'dashed', label = 'z_covave_rmsf_BS_BS', c = 'k' )
    ax.plot(steps, [corr[0] for corr in z_covave_rmsf_RNA_RNA_correlations], marker = 'o', ls = 'dashed', label = 'z_covave_rmsf_RNA_RNA', c = 'magenta' )

    ax.plot(steps, [corr[0] for corr in z_covave_rmsf_noBS_correlations], marker = 'o', ls = 'dashed', label = 'z_covave_rmsf_noBS', )
    ax.plot(steps, [corr[0] for corr in z_covave_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'z_covave_rmsf_BS', color = 'red')
    ax.plot(steps, [corr[0] for corr in pearson_rmsf_BS_BS_correlations], marker = 'o', ls = 'dashed', label = 'pearson_rmsf_BS_BS' )
    ax.plot(steps, [corr[0] for corr in pearson_rmsf_RNA_RNA_correlations], marker = 'o', ls = 'dashed', label = 'pearson_rmsf_RNA_RNA')
    
    ax.plot(steps, [corr[0] for corr in rmsf_rmsf_RNA_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_RNA', )
    #ax.plot(steps, [corr[0] for corr in covstd_rmsf_noBS_correlations], marker = 'o', ls = 'dashed', label = 'covstd_rmsf_noBS', )
    #ax.plot(steps, [corr[0] for corr in covstd_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'covstd_rmsf_BS', )
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

    #%#%
    # SCATTER PLOT

        #Scatterplot STEP 0 
    if not zero_pop:

        ii = 0

        f, ax = plt.subplots()
        ax.set_title('z_covave_rmsf_BS \n Step = {} Treshold = {}\n exclued : pearson = {:3.2f} p-value = {:3.2f}\n all: pearson = {:3.2f} p-value = {:3.2f}'.format(ii, treshold, *pearsonr(Kds[1:], z_covave_rmsf_BS_BS_correlators[0][1:]), *pearsonr(Kds[:], z_covave_rmsf_BS_BS_correlators[0][:])))
        ax.errorbar(Kds[:] , z_covave_rmsf_BS_correlators[0][:], fmt = 'o', xerr = Kds_errs[:], color = 'red', ecolor = 'olivedrab', label = 'correlation data')
        ax.errorbar(Kds[0] , z_covave_rmsf_BS_correlators[0][0], fmt = 'o', xerr = Kds_errs[0], color = 'olivedrab', ecolor = 'red', label = 'NMR excluded')
        ax.set_xlabel('Kd (nM)')
        ax.set_ylabel('Prot BS z Covariance with all complex ')
        ax.legend()
        plt.show()


        f, ax = plt.subplots()
        ax.set_title('z_covave_rmsf_BS_BS \n Step = {} Treshold = {}\n exclued : pearson = {:3.2f} p-value = {:3.2f}\n all: pearson = {:3.2f} p-value = {:3.2f}'.format(ii, treshold, *pearsonr(Kds[1:], z_covave_rmsf_BS_BS_correlators[0][1:]), *pearsonr(Kds[:], z_covave_rmsf_BS_BS_correlators[0][:])))
        ax.errorbar(Kds[:] , z_covave_rmsf_BS_BS_correlators[0][:], fmt = 'o', xerr = Kds_errs[:], color = 'k', ecolor = 'orange', label = 'correlation data')
        ax.errorbar(Kds[0] , z_covave_rmsf_BS_BS_correlators[0][0], fmt = 'o', xerr = Kds_errs[0], color = 'limegreen', ecolor = 'yellowgreen', label = 'NMR excluded')
        ax.set_xlabel('Kd (nM)')
        ax.set_ylabel('Prot BS z Covariance with Prot BS ')
        ax.legend()
        plt.show()

        f, ax = plt.subplots()
        ax.set_title('pearson_rmsf_BS_BS \n Step = {} Treshold = {}\n exclued : pearson = {:3.2f} p-value = {:3.2f}\n all: pearson = {:3.2f} p-value = {:3.2f}'.format(ii, treshold, *pearsonr(Kds[1:], pearson_rmsf_BS_BS_correlators[0][1:]), *pearsonr(Kds[:], pearson_rmsf_BS_BS_correlators[0][:])))
        ax.errorbar(Kds[:] , pearson_rmsf_BS_BS_correlators[0][:], fmt = 'o', xerr = Kds_errs[:], color = 'maroon', ecolor = 'orange', label = 'correlation data')
        ax.errorbar(Kds[0] , pearson_rmsf_BS_BS_correlators[0][0], fmt = 'o', xerr = Kds_errs[0], color = 'limegreen', ecolor = 'yellowgreen', label = 'NMR excluded')
        ax.set_xlabel('Kd (nM)')
        ax.set_ylabel('Prot BS z Pearson with Prot BS ')
        ax.legend()
        plt.show()

        f, ax = plt.subplots()
        ax.set_title('z_covave_rmsf_RNA_RNA \n Step = {} Treshold = {}\n exclued : pearson = {:3.2f} p-value = {:3.2f}\n all: pearson = {:3.2f} p-value = {:3.2f}'.format(ii, treshold, *pearsonr(Kds[1:], z_covave_rmsf_RNA_RNA_correlators[0][1:]), *pearsonr(Kds[:], z_covave_rmsf_RNA_RNA_correlators[0][:])))
        ax.errorbar(Kds[:] , z_covave_rmsf_RNA_RNA_correlators[0][:], fmt = 'o', xerr = Kds_errs[:], color = 'crimson', ecolor = 'darkolivegreen', label = 'correlation data')
        ax.errorbar(Kds[0] , z_covave_rmsf_RNA_RNA_correlators[0][0], fmt = 'o', xerr = Kds_errs[0], color = 'darkgoldenrod', ecolor = 'yellowgreen', label = 'NMR excluded')
        ax.set_xlabel('Kd (nM)')
        ax.set_ylabel('RNA z Covariance with RNA ')
        ax.legend()
        plt.show()

        f, ax = plt.subplots()
        ax.set_title('pearson_rmsf_RNA_RNA \n Step = {} Treshold = {}\n exclued : pearson = {:3.2f} p-value = {:3.2f}\n all: pearson = {:3.2f} p-value = {:3.2f}'.format(ii, treshold, *pearsonr(Kds[1:], pearson_rmsf_RNA_RNA_correlators[0][1:]), *pearsonr(Kds[:], pearson_rmsf_RNA_RNA_correlators[0][:])))
        ax.errorbar(Kds[:] , pearson_rmsf_RNA_RNA_correlators[0][:], fmt = 'o', xerr = Kds_errs[:], color = 'mediumpurple', ecolor = 'darkolivegreen', label = 'correlation data')
        ax.errorbar(Kds[0] , pearson_rmsf_RNA_RNA_correlators[0][0], fmt = 'o', xerr = Kds_errs[0], color = 'darkgoldenrod', ecolor = 'yellowgreen', label = 'NMR excluded')
        ax.set_xlabel('Kd (nM)')
        ax.set_ylabel('RNA z Pearson with RNA ')
        ax.legend()
        plt.show()

        f, ax = plt.subplots()
        ax.set_title('rmsf_BS \n Step = {} Treshold = {}\n exclued : pearson = {:3.2f} p-value = {:3.2f}\n all: pearson = {:3.2f} p-value = {:3.2f}'.format(ii, treshold, *pearsonr(Kds[1:], rmsf_rmsf_BS_correlators[0][1:]), *pearsonr(Kds[:], rmsf_rmsf_correlators[0][:])))
        ax.errorbar(Kds[:] , rmsf_rmsf_BS_correlators[0][:], fmt = 'o', xerr = Kds_errs[:], color = 'orange', ecolor = 'purple', label = 'correlation data')
        ax.errorbar(Kds[0] , rmsf_rmsf_BS_correlators[0][0], fmt = 'o', xerr = Kds_errs[0], color = 'purple', ecolor = 'orange', label = 'NMR excluded')
        ax.set_xlabel('Kd (nM)')
        ax.set_ylabel('RMSF BS ')
        ax.legend()
        plt.show()
# %%
#further study
# FIT

##
#per fit con  errori x,y
import scipy.odr as odr


def f(A, x):
    return A[0]*x+A[1]

linear = odr.Model(f)
# VERSIONE CON PUNTO SPERIMENTALE di wtc1 e con wtc3 (cioè con tutti)
# ERRORI presi sui 2 punti che dovrebbero essere uguali,




#ciclo sugli step
for ii in range(1):

    print('\n\n\n STEP {} \n\n\n'.format(ii))

    for to_exclude in [[3,4]]:#, [2,5]):

        max_err = np.max([z_covave_rmsf_BS_BS_correlators[ii][jj] for jj in to_exclude])
        min_err = np.min([z_covave_rmsf_BS_BS_correlators[ii][jj] for jj in to_exclude])
        err_max = (max_err - min_err)
        err = err_max/3 #procedura standard
        print(err/2)


        correlators = z_covave_rmsf_BS_BS_correlators[ii][1:]
        correlations = z_covave_rmsf_BS_BS_correlations[ii][1:]


        fig, ax = plt.subplots()
        ax.set_title('Protein BS Covariance with Protein BS  vs Kd\ncorrelation = {:3.2f} p-value = {:3.2f}'.format(*pearsonr(Kds[1:], correlators)))
        ax.errorbar(Kds[1:] , correlators, fmt = 'o', xerr = Kds_errs[1:], color = 'k', ecolor = 'orange', label = 'simulated data', mew = 0.1)
        ax.errorbar(Kds[0] , z_covave_rmsf_BS_BS_correlators[0][0], fmt = 'o', xerr = Kds_errs[0], color = 'orange', ecolor = 'k', label = 'NMR data', mew = 0.1)
        for x,y, key in zip(Kds , z_covave_rmsf_BS_BS_correlators[0], WTC_identifier):
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

    #%#%
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

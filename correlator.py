#%%

from BioAlessandria import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

now_path = '../GROMACS/'
save_path = now_path+'CORRELATIONS/'
for treshold in range(8,13):
    
    print('\n\n\nTreshold = {}\n\n\n'.format(treshold))


    DF = pd.read_json(now_path+'WTC_data_frame_{}ang.json'.format(str(treshold)))
    DF['z_Covariance_Mean'] = (DF.Covariance_Mean - np.mean(DF.Covariance_Mean))/(np.std(DF.Covariance_Mean))

    DF_noBS_RNA = pd.concat([DF[DF.Is_BS == False], DF[DF.Is_Prot == False]], ignore_index= True)
    DF_BS_RNA = pd.concat([DF[DF.Is_BS == True], DF[DF.Is_Prot == False]], ignore_index = True)

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
    rmsf_BS_max = 0.19#np.min([np.max(DF.RMSF[DF.Is_BS == True][DF.WTC_identifier == key]) for key in WTC_identifier])
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

    rmsf_rmsf_RNA_correlators = []
    rmsf_rmsf_RNA_correlations = []
    rmsf_rmsf_RNA_populations = []


    #covave
    covave_rmsf_noBS_correlators = []
    covave_rmsf_noBS_correlations = []
    covave_rmsf_noBS_populations = []

    covave_rmsf_BS_correlators = []
    covave_rmsf_BS_correlations = []
    covave_rmsf_BS_populations = []

    covave_rmsf_BS_RNA_correlators = []
    covave_rmsf_BS_RNA_correlations = []
    covave_rmsf_BS_RNA_populations = []

    #z-score covave
    z_covave_rmsf_noBS_correlators = []
    z_covave_rmsf_noBS_correlations = []
    z_covave_rmsf_noBS_populations = []

    z_covave_rmsf_BS_correlators = []
    z_covave_rmsf_BS_correlations = []
    z_covave_rmsf_BS_populations = []


    covstd_rmsf_noBS_correlators = []
    covstd_rmsf_noBS_correlations = []
    covstd_rmsf_noBS_populations = []

    covstd_rmsf_BS_correlators = []
    covstd_rmsf_BS_correlations = []
    covstd_rmsf_BS_populations = []

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

        rmsf_rmsf_RNA = { key : DF.RMSF[DF.WTC_identifier == key][DF.RMSF > rmsf_BS_cutoff][DF.Is_Prot == False].values for key in WTC_identifier}
        rmsf_rmsf_RNA_populations.append({ key : len(rmsf_rmsf_RNA[key]) for key in WTC_identifier})


        covave_rmsf_noBS = {key : DF.Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == False] for key in WTC_identifier }
        covave_rmsf_noBS_populations.append({ key : len(covave_rmsf_noBS[key]) for key in WTC_identifier})

        covave_rmsf_BS = {key : DF.Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True] for key in WTC_identifier }
        covave_rmsf_BS_populations.append({ key : len(covave_rmsf_BS[key]) for key in WTC_identifier})

        #covave_rmsf_noBS_RNA = {key : DF.Covariance_Mean_Prot_noBS[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == False][DF.Is_Prot == False] for key in WTC_identifier }
        #covave_rmsf_noBS_RNA_populations.append({ key : len(covave_rmsf_noBS_RNA[key]) for key in WTC_identifier})

        covave_rmsf_BS_RNA = {key : DF.Covariance_Mean_RNA_BS[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True][DF.Is_Prot == True] for key in WTC_identifier }
        covave_rmsf_BS_RNA_populations.append({ key : len(covave_rmsf_BS_RNA[key]) for key in WTC_identifier})


        z_covave_rmsf_noBS = {key : DF.z_Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == False] for key in WTC_identifier }
        z_covave_rmsf_noBS_populations.append({ key : len(z_covave_rmsf_noBS[key]) for key in WTC_identifier})

        z_covave_rmsf_BS = {key : DF.z_Covariance_Mean[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True] for key in WTC_identifier }
        z_covave_rmsf_BS_populations.append({ key : len(z_covave_rmsf_BS[key]) for key in WTC_identifier})


        covstd_rmsf_noBS = {key : DF.Covariance_Std[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == False] for key in WTC_identifier }
        covstd_rmsf_noBS_populations.append({ key : len(covstd_rmsf_noBS[key]) for key in WTC_identifier})
        covstd_rmsf_BS = {key : DF.Covariance_Std[DF.WTC_identifier == key][DF.RMSF > rmsf_cutoff][DF.Is_BS == True] for key in WTC_identifier }
        covstd_rmsf_BS_populations.append({ key : len(covstd_rmsf_BS[key]) for key in WTC_identifier})

        pop_step = pd.DataFrame(index=WTC_identifier)
        pop_step['rmsf_rmsf'] = [rmsf_rmsf_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['rmsf_rmsf_BS'] = [rmsf_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['rmsf_rmsf_RNA'] = [rmsf_rmsf_RNA_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['covave_rmsf_noBS'] = [covave_rmsf_noBS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['covave_rmsf_BS'] = [covave_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['covave_rmsf_BS_RNA'] = [covave_rmsf_BS_RNA_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['z_covave_rmsf_noBS'] = [z_covave_rmsf_noBS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['z_covave_rmsf_BS'] = [z_covave_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['covstd_rmsf_noBS'] = [covstd_rmsf_noBS_populations[n_step-1][key] for key in WTC_identifier]
        pop_step['covstd_rmsf_BS'] = [covstd_rmsf_BS_populations[n_step-1][key] for key in WTC_identifier]


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
            covave_rmsf_noBS_populations.remove(covave_rmsf_noBS_populations[-1])
            covave_rmsf_BS_populations.remove(covave_rmsf_BS_populations[-1])
            covave_rmsf_BS_RNA_populations.remove(covave_rmsf_BS_RNA_populations[-1])
            z_covave_rmsf_noBS_populations.remove(z_covave_rmsf_noBS_populations[-1])
            z_covave_rmsf_BS_populations.remove(z_covave_rmsf_BS_populations[-1])
            covstd_rmsf_noBS_populations.remove(covstd_rmsf_noBS_populations[-1])
            covstd_rmsf_BS_populations.remove(covstd_rmsf_BS_populations[-1])
            break
        else:

            #rmsf-rmsf
            rmsf_rmsf_correlators.append([np.mean(rmsf_rmsf[key]) for key in WTC_identifier])
            rmsf_rmsf_correlations.append(pearsonr(Kds, rmsf_rmsf_correlators[n_step-1]))

            rmsf_rmsf_BS_correlators.append([np.mean(rmsf_rmsf_BS[key]) for key in WTC_identifier])
            rmsf_rmsf_BS_correlations.append(pearsonr(Kds, rmsf_rmsf_BS_correlators[n_step-1]))

            rmsf_rmsf_RNA_correlators.append([np.mean(rmsf_rmsf_RNA[key]) for key in WTC_identifier])
            rmsf_rmsf_RNA_correlations.append(pearsonr(Kds, rmsf_rmsf_RNA_correlators[n_step-1]))

            covave_rmsf_noBS_correlators.append([np.mean(covave_rmsf_noBS[key]) for key in WTC_identifier])
            covave_rmsf_noBS_correlations.append(pearsonr(Kds, covave_rmsf_noBS_correlators[n_step-1]))
            
            covave_rmsf_BS_correlators.append([np.mean(covave_rmsf_BS[key]) for key in WTC_identifier])
            covave_rmsf_BS_correlations.append(pearsonr(Kds, covave_rmsf_BS_correlators[n_step-1]))

            covave_rmsf_BS_RNA_correlators.append([np.mean(covave_rmsf_BS_RNA[key]) for key in WTC_identifier])
            covave_rmsf_BS_RNA_correlations.append(pearsonr(Kds, covave_rmsf_BS_RNA_correlators[n_step-1]))


            z_covave_rmsf_noBS_correlators.append([np.mean(z_covave_rmsf_noBS[key]) for key in WTC_identifier])
            z_covave_rmsf_noBS_correlations.append(pearsonr(Kds, z_covave_rmsf_noBS_correlators[n_step-1]))
            
            z_covave_rmsf_BS_correlators.append([np.mean(z_covave_rmsf_BS[key]) for key in WTC_identifier])
            z_covave_rmsf_BS_correlations.append(pearsonr(Kds, z_covave_rmsf_BS_correlators[n_step-1]))

            covstd_rmsf_noBS_correlators.append([np.mean(covstd_rmsf_noBS[key]) for key in WTC_identifier])
            covstd_rmsf_noBS_correlations.append(pearsonr(Kds, covstd_rmsf_noBS_correlators[n_step-1]))
            
            covstd_rmsf_BS_correlators.append([np.mean(covstd_rmsf_BS[key]) for key in WTC_identifier])
            covstd_rmsf_BS_correlations.append(pearsonr(Kds, covstd_rmsf_BS_correlators[n_step-1]))




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
    ax.plot(steps, [corr[0] for corr in covave_rmsf_BS_RNA_correlations], marker = 'o', ls = 'dashed', label = 'covave_rmsf_BS_RNA', )
    ax.plot(steps, [corr[0] for corr in z_covave_rmsf_noBS_correlations], marker = 'o', ls = 'dashed', label = 'z_covave_rmsf_noBS', )
    ax.plot(steps, [corr[0] for corr in z_covave_rmsf_BS_correlations], marker = 'o', ls = 'dashed', label = 'z_covave_rmsf_BS', )

    #ax.plot(steps, [corr[0] for corr in rmsf_rmsf_RNA_correlations], marker = 'o', ls = 'dashed', label = 'rmsf_rmsf_RNA', )
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


# %%
#further study

##
#per fit con  errori x,y
import scipy.odr as odr
linear = odr.Model(f)
def f(A, x):
    return A[0]*x+A[1]

def f_lin(x, m, q):
    return m*x+q
##

# VERSIONE CON PUNTO SPERIMENTALE di wtc1
# ERRORI presi sui 3 punti che dovrebbero essere uguali,
# o meglio sui punti che rimangono escludend



to_exclude = [2,3, 5]
# 1) 
#ciclo sugli step
for ii in range(1):

    print('\n\n\n STEP {} \n\n\n'.format(ii))
    old_correlators = z_covave_rmsf_BS_correlators[ii]

    for kk in to_exclude:
        
        #variabili che aggiorno
        new_correlators = old_correlators.copy()
        new_correlators.remove(new_correlators[kk])
        

        new_Kds = Kds.copy()
        new_Kds.remove(new_Kds[kk])


        new_correlations = pearsonr(new_Kds, new_correlators)
        not_excluded = to_exclude.copy()
        not_excluded.remove(kk)

        max_err = np.max([old_correlators[jj] for jj in not_excluded])
        min_err = np.min([old_correlators[jj] for jj in not_excluded])
        err_max = (max_err - min_err)
        err = err_max/3 #procedura standard


        fig, ax = plt.subplots()
        ax.set_ylabel('z - Correlation')
        ax.set_xlabel('Kd (nM)')
        ax.set_title('Pearson correlation of covave_rmsf_BS at step {}\ncorrelation = {:3.2f}   p value = {:3.2f}'.format(ii, new_correlations[0], new_correlations[1]))
        ax.scatter(Kds[:], np.array(old_correlators[:]), c = 'red', label = 'data')
        ax.plot(Kds[kk], np.array(old_correlators[kk]), 'o',  color = 'orange', fillstyle = None, lw = 3, label = 'excluded from fit')
        ax.vlines(1340, min_err, max_err, linestyles = 'dashed', color = 'yellowgreen', label = 'max error bar', linewidths = 2.)
        ax.legend()
        plt.show()

    
        mydata = odr.RealData(new_Kds, new_correlators, sy = err/2,)
        myodr = odr.ODR(mydata, linear, beta0=[1.,2.])
        myoutput = myodr.run()

        fig, ax = plt.subplots()
        x = np.linspace(np.min(new_Kds), np.max(new_Kds), 1000)
        ax.set_title(r'$\bf{y = mx + q}$'+' fit for correlation of {:3.2f} \n treshold = {} Ang step = {}\n {} excluded'.format(new_correlations[0], treshold, ii, WTC_identifier[kk]), fontsize = 16)
        ax.set_xlabel('Kds (nM)')
        ax.set_ylabel('z Covariance')
        ax.errorbar(new_Kds, new_correlators, yerr = err/2, fmt = 'o', color = 'red', ecolor = 'orange', barsabove = True, lolims = True, uplims = True, label = 'data')
        ax.plot(x,f(myoutput.beta, x), color = 'yellowgreen', label = 'linear fit')
        ax.legend(title = 'm = {:3.2e}$\pm${:3.2e}\nq = {:3.2e}$\pm${:3.2e}'.format(myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1]))
        plt.show()

# %%

## VERSIONE SENZA DATO NMR

to_exclude = [2,3,5]

#ciclo step
for ii in range(1):

    print('\n\n\n STEP {} \n\n\n'.format(ii))
    old_correlators = z_covave_rmsf_BS_correlators[ii]

    for kk in to_exclude:

        #aggiornamento variabili

        new_correlators = old_correlators.copy()
        new_correlators.remove(new_correlators[kk])
        new_correlators.remove(new_correlators[0])
        new_Kds = Kds.copy()
        new_Kds.remove(new_Kds[kk])
        new_Kds.remove(new_Kds[0])

        new_correlations = pearsonr(new_Kds, new_correlators)
        not_excluded = to_exclude.copy()
        not_excluded.remove(kk)

        max_err = np.max([old_correlators[jj] for jj in not_excluded])
        min_err = np.min([old_correlators[jj] for jj in not_excluded])
        err_max = (max_err - min_err)
        err = err_max/3 #procedura standard



        fig, ax = plt.subplots()
        ax.set_ylabel('z - Correlation')
        ax.set_xlabel('Kd (nM)')
        ax.set_title('Pearson correlation of z - covave_rmsf_BS at step {}\ncorrelation = {:3.2f}   p value = {:3.2f}'.format(ii, new_correlations[0], new_correlations[1]))
        ax.scatter(Kds[1:], np.array(old_correlators[1:]), c = 'red', label = 'data')
        ax.plot(Kds[kk], np.array(old_correlators[kk]), 'o',  color = 'orange', fillstyle = None, lw = 3, label = 'excluded from fit')
        ax.plot(Kds[0], np.array(old_correlators[0]), 'o',  color = 'orange', fillstyle = None, lw = 3,)
        ax.vlines(1340, min_err, max_err, linestyles = 'dashed', color = 'yellowgreen', label = 'max error bar', linewidths = 2.)
        ax.legend()
        plt.show()

    
        mydata = odr.RealData(new_Kds, new_correlators, sy = err/2, )
        myodr = odr.ODR(mydata, linear, beta0=[1.,2.])
        myoutput = myodr.run()

        fig, ax = plt.subplots()
        x = np.linspace(np.min(new_Kds), np.max(new_Kds), 1000)
        ax.set_title(r'$\bf{y = mx + q}$'+' fit for correlation = {:3.2f}\n treshold = {} Ang step = {}\n {} excluded'.format(new_correlations[0], treshold, ii, WTC_identifier[kk]), fontsize = 16)
        ax.set_xlabel('Kds (nM)')
        ax.set_ylabel('z - Covariance')
        ax.errorbar(new_Kds, new_correlators, yerr = err/2, fmt = 'o', color = 'red', ecolor = 'orange', barsabove = True, lolims = True, uplims = True, label = 'data')
        ax.plot(x,f(myoutput.beta, x), color = 'yellowgreen', label = 'linear fit')
        ax.legend(title = 'm = {:3.2e}$\pm${:3.2e}\nq = {:3.2e}$\pm${:3.2e}'.format(myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1]))
        plt.show()

# %%

# VERSIONE TUTTI FRUTTI - errore coi 3

to_exclude = [2,5]
max_err = np.max([old_correlators[jj] for jj in to_exclude])
min_err = np.min([old_correlators[jj] for jj in to_exclude])
err_max = (max_err - min_err)
err = err_max/3 #procedura standard

fig, ax = plt.subplots()
ax.set_ylabel('z - Correlation')
ax.set_title('Pearson correlation of z_covave_rmsf_BS at step {}\ncorrelation = {:3.2f}   p value = {:3.2f}'.format(ii, z_covave_rmsf_BS_correlations[0][0], z_covave_rmsf_BS_correlations[0][1]))
ax.scatter(Kds[:], z_covave_rmsf_BS_correlators[0][:], c = 'red', label = 'data')
ax.vlines(1340, min_err, max_err, linestyles = 'dashed', color = 'yellowgreen', label = 'max error bar', linewidths = 2.)
ax.legend()
plt.show()

mydata = odr.RealData(Kds, z_covave_rmsf_BS_correlators[0], sy = err/2)
myodr = odr.ODR(mydata, linear, beta0=[1.,2.])
myoutput = myodr.run()

fig, ax = plt.subplots()
x = np.linspace(np.min(new_Kds), np.max(new_Kds), 1000)
ax.set_title(r'$\bf{y = mx + q}$'+' fit for correlation = {:3.2f} \n treshold = {} Ang step = {}\nAll included'.format(z_covave_rmsf_BS_correlations[ii][0], treshold, ii), fontsize = 16)
ax.set_xlabel('Kds (nM)')
ax.set_ylabel('z Covariance')
ax.errorbar(Kds, old_correlators, yerr = err/2, fmt = 'o', color = 'red', ecolor = 'orange', barsabove = True, lolims = True, uplims = True, label = 'data')
ax.plot(x,f(myoutput.beta, x), color = 'yellowgreen', label = 'linear fit')
ax.legend(title = 'm = {:3.2e}$\pm${:3.2e}\nq = {:3.2e}$\pm${:3.2e}'.format(myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1]))
plt.show()




# %%

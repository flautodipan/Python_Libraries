#%%
#0) preambolo

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import os
import pandas as pd
import ossaudiodev
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
from scipy.spatial.distance import euclidean
from scipy.stats import pearsonr

#exp data 
Kds_sort  = [4., 4., 700., 1320., 1350., 1360., 1800.]
Kds = [4., 4., 1360., 1350, 700, 1320, 1800]
#Kds = {'wtc{}'.format(n): kd for (n, kd) in zip(['1', '1_h', '2', '3', '4', '5', '6'], Kds_)}

Kds_new =  [4, 4, 1360, 650, 750, 1320, 1800]
Kds_errs = [0.9, 0.9, 600, 135, 150, 350, 400]

WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6')


#inizializzo dataframe

descriptor = {}
now_temp = '300 K'
scale='ns'
df = pd.DataFrame()


#%%
for treshold in [9, 11]:#range(8,13):

    for ii in ['1', '1_h', '2', '3', '4', '5', '6']:

        print('\n\n\nDinamica  wtc{} treshold = {}\n\n\n'.format(ii,treshold))

        if ii == '1_h':

            #WTC1_h
            now_path    =   '../GROMACS/WTC1_h/'
            now_name    =    'wtc1_h'
            n_frames = 20001
            time_range  =   [0, 2000000]
            time_range_eq = [1000000, 2000000]

        elif ii == '1':
            if treshold  == 9: continue
            else:
                #WTC1
                now_path    =   '../GROMACS/WTC1/'
                now_name    =    'wtc1'
                n_frames = 9700
                time_range  =   [30060, 999960]
                time_range_eq = [400060, 999960]
        elif ii == '2':
            #WTC2
            now_path    =   '../GROMACS/WTC2/'
            now_name    =    'wtc2'
            n_frames = 10001
            time_range = [0, 1000000]
            time_range_eq = [400000, 1000000]
        elif ii == '3':
            #WTC3
            now_path    =   '../GROMACS/WTC3/'
            now_name    =    'wtc3'
            n_frames = 20001
            time_range = [0, 2000000]
            time_range_eq = [600000, 2000000]
        elif ii == '4':
            #WTC4
            now_path    =   '../GROMACS/WTC4/'
            now_name    =    'wtc4'
            n_frames = 10001
            time_range = [0, 1000000]
            time_range_eq = [400000, 1000000]
        elif ii == '5':
            #WTC5
            now_path    =   '../GROMACS/WTC5/'
            now_name    =    'wtc5'
            n_frames = 30001
            time_range = [0, 3000000]
            time_range_eq = [2100000, 3000000]
        elif ii == '6':
            #WTC6
            now_path    =   '../GROMACS/WTC6/'
            now_name    =    'wtc6'
            n_frames = 10001
            time_range = [0, 1000000]
            time_range_eq = [400000, 1000000]

        now_path+='eq_frame_'+now_name+'/'
        min_CAP_dist = []

        for frame in os.listdir(now_path):


            TDP43   = BA.Protein(now_path+frame, model = False)
            TDP43.Get_CA_Coord(atom_name='CA')
            RNA     =  BA.RNA(now_path+frame, chain_id='B', model = False, initial=TDP43.initial+TDP43.CA.shape[0])
            RNA.Get_P_Coord(atom_name="P")
            RNA_start_idx = len(TDP43.CA) - 1

            print('Ok, acquisito correttamente pdb {} per {}'.format(frame, now_name))

            Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
            Dist    = BA.Dist_Matrix(Coord)
            Cont    = BA.Contacts_Matrix(Dist, treshold)
            Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))
            BS      = BA.Print_Protein_BS_old(Bonds, TDP43.lenght,prot_initial=TDP43.initial, RNA_initial= RNA.initial)
            
            distances = []
            for x_CA,y_CA,z_CA, r_number in TDP43.CA[['x_coord', 'y_coord', 'z_coord', 'residue_number']].values:
                if r_number in BS['Prot']:
                    distances.append(np.min([euclidean((x_CA, y_CA, z_CA), (x_P, y_P, z_P)) for x_P, y_P, z_P in RNA.P[['x_coord', 'y_coord', 'z_coord']].values]))

            min_CAP_dist.append(np.mean(distances))

        
        f, ax = plt.subplots()
        ax.set_title("Minum mean CA-P distace in time\n {} dynamics treshold = {}".format(now_name, treshold))
        ax.set_xlabel('Time (Number of eq frame (100 -sampled))')
        ax.set_ylabel('CA - P mean distance (Ang)')
        ax.plot(min_CAP_dist)
        plt.show()

        descriptor[now_name] = (np.mean(min_CAP_dist))

    df[str(treshold)+' ang'] = [descriptor[key] for key in WTC_identifier]



df.to_json('../GROMACS/df_CAPdist_new.json')


#%%
import scipy.odr as odr
def f_lin(A, x):
    return A[0]*x+A[1]

linear = odr.Model(f_lin)

df = pd.read_json('../GROMACS/df_CAPdist_new.json')
df = df.rename(index = { old : new for old, new in zip(range(7), WTC_identifier)})

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
    ax.errorbar(Kds_new[1:], [df[treshold][key] for key in WTC_identifier][1:], xerr = Kds_errs[1:], fmt = 'o',  color = 'green', ecolor = 'magenta',mew = 0.1, label = 'simulated data')
    ax.errorbar(Kds_new[0], [df[treshold][key] for key in WTC_identifier][0], xerr = Kds_errs[0], fmt = 'o',  color = 'magenta', ecolor = 'green', label = 'NMR', mew = 0.1)
    ax.vlines(np.mean([Kds_new[kk] for kk in to_exclude]), min_err, max_err, linestyles = 'dashed', color = 'yellowgreen', label = 'max error bar', linewidths = 2.)
    for x, key in zip(Kds_new, WTC_identifier):
        y = df[treshold][key]
        if key in ['wtc1_h']:
            key = 'wtc1'
            plt.annotate(key, (x,y), xytext = (5, -12), textcoords = 'offset points', )
        elif key == 'wtc1': continue
        elif key in ['wtc3', 'wtc6']:plt.annotate(key, (x,y), xytext = (-25, -15), textcoords = 'offset points', )

        else: plt.annotate(key, (x,y), xytext = (5, 10), textcoords = 'offset points')
    ax.legend()
    ax.set_xlabel('Kd (nM)')
    ax.set_ylabel('CA-P Mean Dist (Ang)')

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

for treshold in df.columns:

    f, ax = plt.subplots()
    ax.set_title('Old Correlation at BS treshold = {} Ang\n pearson = {:3.2f} p-value = {:3.2f}'.format(treshold, *pearsonr(Kds, [df[treshold][key] for key in WTC_identifier])))
    ax.scatter(Kds, [df[treshold][key] for key in WTC_identifier], color = 'green')
    ax.set_xlabel('Kd (nM)')
    ax.set_ylabel('CA-P Mean Dist (Ang)')


# %%

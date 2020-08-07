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
Kds_sort  = [4., 4., 70, 700., 1320., 1350., 1360., 1800.]
Kds = [4., 4., 1360., 1350, 700, 1320, 1800, 70]
#Kds = {'wtc{}'.format(n): kd for (n, kd) in zip(['1', '1_h', '2', '3', '4', '5', '6'], Kds_)}

Kds_new =  [4, 4, 1360, 650, 750, 1320, 1800, 70]
Kds_errs = [0.9, 0.9, 600, 135, 150, 350, 400, 10]

WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6', 'wtc7')

#seven reduced
red = False
wtc7_7 = True
wtc7_8 = False

#inizializzo dataframe

descriptor = {}
now_temp = '300 K'
scale='ns'
df = pd.DataFrame()


#%%
for treshold in [9]:#range(8,13):

    for ii in ['7']:#['1', '1_h', '2', '3', '4', '5', '6', '7']:

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

        elif ii == '7':
            
            if red == True:
                #WTC7
                now_path    =   '../GROMACS/WTC7/'
                now_name    =    'wtc7'
                n_frames = 10001
                #ps
                time_range = [0, 1000000]
                time_range_eq = [200000, 950000]

            elif wtc7_8 == True:

                now_path    =   '../GROMACS/WTC7_8/'
                now_name    =    'wtc7_8'
                n_frames = 3001
                time_range = [0, 300000]
                time_range_eq = time_range

            
            elif wtc7_7 == True:

                now_path    =   '../GROMACS/WTC7_7/'
                now_name    =    'wtc7_7'
                n_frames = 2301
                time_range = [0, 230000]
                time_range_eq = time_range

            else:
            
                #WTC7
                now_path    =   '../GROMACS/WTC7_red/'
                now_name    =    'wtc7'
                n_frames = 1251
                time_range = [175000, 300000]
                time_range_eq = time_range

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


if red:

    df.to_json('../GROMACS/df_CAPdist_new_red.json')

elif wtc7_7:

    df.to_json('../GROMACS/df_CAPdist_new_7.json')


elif wtc7_8:

    df.to_json('../GROMACS/df_CAPdist_new_8.json')

else:

    df.to_json('../GROMACS/df_CAPdist_new_red.json')












# %%

for treshold in df.columns:

    f, ax = plt.subplots()
    ax.set_title('Old Correlation at BS treshold = {} Ang\n pearson = {:3.2f} p-value = {:3.2f}'.format(treshold, *pearsonr(Kds, [df[treshold][key] for key in WTC_identifier])))
    ax.scatter(Kds, [df[treshold][key] for key in WTC_identifier], color = 'green')
    ax.set_xlabel('Kd (nM)')
    ax.set_ylabel('CA-P Mean Dist (Ang)')


# %%

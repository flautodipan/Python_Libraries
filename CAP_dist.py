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
red     = False
wtc7_7  = False
wtc7_8  = False
wtc7_16 = False
wtc7_24 = True

#inizializzo dataframe

descriptor = {}
now_temp = '300 K'
scale='ns'
df = pd.DataFrame()


#%%
for treshold in [9]:#range(8,13):

    for ii in ['7']:# ['1', '1_h', '2', '3', '4', '5', '6', '7']:

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

            elif wtc7_16 == True:

                now_path    =   '../GROMACS/WTC7_16/'
                now_name    =    'wtc7_16'
                n_frames = 2001
                time_range = [0, 200000]
                time_range_eq = time_range
            
            elif wtc7_24 == True:

                now_path    =   '../GROMACS/WTC7_24/'
                now_name    =    'wtc7_24'
                n_frames = 2001
                time_range = [0, 200000]
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
    #df[str(treshold)+' ang'] = [descriptor[key] for key in WTC_identifier]
    df['9 ang'] = [descriptor['wtc7_24']]

if red:

    df.to_json('../GROMACS/df_CAPdist_new_red.json')

elif wtc7_7:

    df_ = pd.read_json('../GROMACS/df_CAPdist_new.json')
    df_ = df_.drop(7)
    df = df_.append(df, ignore_index = True)
    df.to_json('../GROMACS/df_CAPdist_new_7.json')


elif wtc7_8:


    df_ = pd.read_json('../GROMACS/df_CAPdist_new.json')
    df_ = df_.drop(7)
    df = df_.append(df, ignore_index = True)
    df.to_json('../GROMACS/df_CAPdist_new_8.json')

elif wtc7_16:


    df_ = pd.read_json('../GROMACS/df_CAPdist_new.json')
    df_ = df_.drop(7)
    df = df_.append(df, ignore_index = True)
    df.to_json('../GROMACS/df_CAPdist_new_16_all.json')

elif wtc7_24:

    df_ = pd.read_json('../GROMACS/df_CAPdist_new.json')
    df_ = df_.drop(7)
    df = df_.append(df, ignore_index = True)
    df.to_json('../GROMACS/df_CAPdist_new_24_all.json')

else:

    df.to_json('../GROMACS/df_CAPdist_new_red.json')   



# %%
# PARTE per P-CA ossia studio della comparison WTC1 WTC1h

#CONFIGURAZIONE INIZIALE

treshold = 9
df_initial = pd.DataFrame()

for ii in ['1', '1_h']:

    print('\n\n\nDinamica  wtc{} treshold = {}\n\n\n'.format(ii,treshold))

    now_name    =    'wtc'+ii
    now_path = '../GROMACS/WTC'+ii+'/'

    frame = now_name+'.pdb'

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
    
    d = {}
    
    for (x_P, y_P, z_P, res_num, res_name) in RNA.P[['x_coord', 'y_coord', 'z_coord', 'residue_number', 'residue_name']].values:

        d[res_name+' '+str(res_num)] =  np.min([euclidean((x_CA, y_CA, z_CA),(x_P, y_P, z_P)) for x_CA,y_CA,z_CA, r_num in TDP43.CA[['x_coord', 'y_coord', 'z_coord', 'residue_number']].values if r_num in BS['Prot'] ])


    df_initial[now_name] = d.values()
    df_initial = df_initial.rename(index = { old : new for old, new in zip(range(11), d.keys())})

#%%
#PLOT

labels = d.keys()
x = np.arange(len(labels))
x_BS = []
for jj in x:
    cond_div1 = (jj in np.where(df_initial['wtc1'].values < 9)[0] ) & ( jj not in np.where(df_initial['wtc1_h'].values < 9)[0])
    cond_div2 = (jj in np.where(df_initial['wtc1'].values > 9)[0] ) & ( jj not in np.where(df_initial['wtc1_h'].values > 9)[0])
    if cond_div1 | cond_div2: 
        x_BS.append(jj)

labels_BS = [l for l,ii in zip(labels, range(len(labels))) if ii in x_BS]
x_BS = np.array(x_BS)
width = 0.25
print(x_BS, labels_BS)
fig, ax = plt.subplots()
ax.set_title('Initial config')
ax.hlines(9, -0.5, 10.5, color = 'yellowgreen', linestyles = 'dashed', label = 'BS 9 $\AA$')
ax.bar(x-width/2, df_initial['wtc1'], width, color = 'goldenrod')
ax.bar(x+width/2, df_initial['wtc1_h'],width , color = 'firebrick')
ax.bar(x_BS-width/2, [df_initial['wtc1'][l] for l in labels_BS], width, label = 'NMR', color = 'grey')
ax.bar(x_BS+width/2, [df_initial['wtc1_h'][l] for l in labels_BS],width, label = 'Test' , color = 'k', alpha = 0.6)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


# %%
# PARTE per P-CA ossia studio della comparison WTC1 WTC1h
#CONFIGURAZIONE all'equilibrio

treshold = 9
dfs = []

for ii in ['1', '1_h']:
    df = pd.DataFrame()
    print('\n\n\nDinamica  wtc{} treshold = {}\n\n\n'.format(ii,treshold))

    now_name    =    'wtc'+ii
    now_path = '../GROMACS/WTC'+ii+'/eq_frame_'+now_name+'/'

    lista = os.listdir(now_path)
    lista.sort()
    #sample = [lista[jj] for jj in np.arange(0, len(lista), 10)]

    for frame in lista:


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
        
        d = {}
        
        for (x_P, y_P, z_P, res_num, res_name) in RNA.P[['x_coord', 'y_coord', 'z_coord', 'residue_number', 'residue_name']].values:

            d[res_name+' '+str(res_num)] =  np.min([euclidean((x_CA, y_CA, z_CA),(x_P, y_P, z_P)) for x_CA,y_CA,z_CA, r_num in TDP43.CA[['x_coord', 'y_coord', 'z_coord', 'residue_number']].values if r_num in BS['Prot'] ])

        df[frame[10:12]] = d.values()
    
    df = df.rename(index = { old : new for old, new in zip(range(11), d.keys())})
    dfs.append(df)

DF = pd.DataFrame()
DF['wtc1'] = dfs[0].mean(axis = 1)
DF['wtc1_h'] = dfs[1].mean(axis = 1)

DF.to_json('../GROMACS/df_CPA_dist_9_ang.json', orient='index')

DF_std= pd.DataFrame()

DF_std['wtc1'] = dfs[0].std(axis = 1)
DF_std['wtc1_h'] = dfs[1].std(axis = 1)
DF_std.to_json('../GROMACS/df_CPA_dist_9_ang_std.json', orient='index')

#%%
#PLOT
DF = pd.read_json('../GROMACS/df_CPA_dist_9_ang.json').T
x_BS = []
for jj in x:
    condeq1 =  (jj in np.where(DF['wtc1'].values < 9)[0] ) & ( jj in np.where(DF['wtc1_h'].values < 9)[0])
    condeq2 =  (jj in np.where(DF['wtc1'].values > 9)[0] ) & ( jj in np.where(DF['wtc1_h'].values >  9)[0])

    if condeq1 | condeq2: x_BS.append(jj)

x_BS = np.array(x_BS)
labels = DF.index
x = np.arange(len(labels))
labels_BS = [l for l,ii in zip(labels, range(len(labels))) if ii in x_BS]
width = 0.25
fig, ax = plt.subplots()
ax.hlines(9, -0.5, 10.5, color = 'yellowgreen', linestyles = 'dashed', label = 'BS 9 $\AA$')
ax.bar(x-width/2, DF['wtc1'], width, color = 'gray')
ax.bar(x+width/2, DF['wtc1_h'],width , color = 'k', alpha = 0.6)
ax.bar(x_BS-width/2, [DF['wtc1'][l] for l in labels_BS], width, label = 'NMR', color = 'goldenrod')
ax.bar(x_BS+width/2, [DF['wtc1_h'][l] for l in labels_BS],width, label = 'Test' , color = 'firebrick')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


# %%
#%%
#PLOT con config iniziale NMR - iniziale - finale WTC
labels = DF.index
x = np.arange(len(labels))
width = 0.25
fig, ax = plt.subplots()
ax.hlines(9, -0.5, 10.5, color = 'k', linestyles = 'dashed')
ax.text(-0.75, 9.5, 'BS 9 $\AA$', )
ax.bar(x-width, [df_initial['wtc1'][l] for l in labels], width, label = 'NMR initial', color = 'goldenrod')
ax.bar(x, [df_initial['wtc1_h'][l] for l in labels],width, label = 'WTC1 initial' , color = 'olive', alpha = 0.8)
ax.bar(x+width, [DF['wtc1_h'][l] for l in labels],width, label = 'WTC1 equilibrated' , color = 'firebrick', alpha = 0.8, capsize = 2.5, yerr =  [DF_std['wtc1_h'][l] for l in labels], ecolor = 'k')

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
ax.set_xlabel('RNA$_{12}$ nucleotides')
ax.set_ylabel('Distance ($\AA$)')
fig.savefig('../Scrittek/figures/residue_NMR_comparison.pdf', format = 'pdf')

# %%
from scipy.stats import pearsonr

stats = pearsonr([DF['wtc1_h'][l] for l in labels], [df_initial['wtc1'][l] for l in labels] )
fig, ax = plt.subplots()
ax.set_xlabel('NMR')
ax.set_ylabel('WTC1_test')
ax.set_title('RNA distance from Protein NMR vs WTC1\nPearson = {:3.2f} p-value = {:3.2f}'.format(*stats))
ax.scatter([DF['wtc1_h'][l] for l in labels], [df_initial['wtc1'][l] for l in labels] )

# %%

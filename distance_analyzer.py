#%%
import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.spatial.distance import euclidean
from scipy.stats import pearsonr

path = '../GROMACS/'
wtc_keys = ['wtc1_h', 'wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
cov_keys = ['cov1', 'cov2', 'cov3', 'cov4']
mtc_keys = ['mtc1', 'mtc2']
exp_df = pd.read_excel(path+'MD_experimental_data.xlsx')
eq = 'eq'
now_keys = cov_keys[1:]

#%%

for now_name in ['wtc5_new', ]:

    now_exp_df = exp_df[exp_df.identifier == now_name]
    now_path = path + now_name.upper()  + '/'
    eq_frame_path = now_path +'/'+eq+'_frame_'+now_name+'/'

    n_frames        = int(now_exp_df.n_frames.values[0])
    time_range      = eval(now_exp_df.time_range.values[0])
    time_range_eq   = eval(now_exp_df.time_range_eq.values[0])
    treshold = now_exp_df.bs_treshold.values[0]

    min_CAP_dist = []
    min_CAO_dist = []


    for frame in os.listdir(eq_frame_path):

        if frame.endswith('.pdb'):
            #RMSF
            #acquisisco indici residui proteina da file .xvg dato da GROMACS
            res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+'_'+eq+'.xvg', now_path, skip_lines= 17)
            res = np.array(res, dtype=int)
            # faccio diventare indici dei residui RNA successivi a quelli proteina
            idx_RNA_start = np.where( res == 1.)[0][0]
            res[idx_RNA_start:] += res[idx_RNA_start-1]


            protein   = BA.Protein(eq_frame_path+frame, model = False)
            RNA     =  BA.RNA(eq_frame_path+frame, chain_id='B', model = False)
            protein.Get_Atom_Coord(atom_name = 'CA')
            protein.Get_lenght()
            protein.Get_Protein_Sequence()
            RNA.Get_Atom_Coord(atom_name="O5'")
            RNA.Get_Atom_Coord(atom_name='P')
            RNA.Get_RNA_Sequence()
            RNA.Get_lenght()
            RNA_start_idx = len(protein.CA) - 1

            print('Ok, acquisito correttamente pdb {} per {}'.format(frame, now_name))
            # 2) P
            Coord_P   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
            Dist_P    = BA.Dist_Matrix(Coord_P)
            Cont_P   = BA.Contacts_Matrix(Dist_P, treshold)
            Bonds_P   = BA.Analyze_Bond_Residues(Cont_P, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
            BS_P      = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
            BS_RNA_P = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']
        
            distances = []
            for x_CA,y_CA,z_CA, r_number in protein.CA[['x_coord', 'y_coord', 'z_coord', 'residue_number']].values:
                if r_number in BS_P:
                    distances.append(np.min([euclidean((x_CA, y_CA, z_CA), (x_P, y_P, z_P)) for x_P, y_P, z_P in RNA.P[['x_coord', 'y_coord', 'z_coord']].values]))

            min_CAP_dist.append(np.mean(distances))
    
    np.save(now_path+now_name+'_min_CAPdist_'+eq+'.npy', np.array(min_CAP_dist) )

#%%


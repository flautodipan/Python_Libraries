#%%
#0) preambolo

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

#exp data 
Kds_sort  = [4., 4., 700., 1320., 1350., 1360., 1800.]
Kds_ = [4., 4., 1360., 1350, 700, 1320, 1800]
Kds = {'wtc{}'.format(n): kd for (n, kd) in zip(['1', '1_h', '2', '3', '4', '5', '6'], Kds_)}

#inizializzo dataframe

dfs = {}

now_temp = '300 K'
scale='ns'

for treshold in range(8,13):
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

        #WTC analyzer senza figure

        WTC_traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, now_temp))
        WTC_traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )
        WTC_traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', path = now_path)
        WTC_traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq)
        filename='BS_{}_make_ndx.txt'.format(now_name)
        pdb_filename = 'average_pdb_'+now_name+'.pdb'
        #treshold = 10

        TDP43   = BA.Protein(now_path+pdb_filename, model = False)
        TDP43.Get_CA_Coord(atom_name='CA')
        RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False, initial=TDP43.initial+TDP43.CA.shape[0])
        RNA.Get_P_Coord(atom_name="P")
        print('Ok, acquisito correttamente pdb')

        Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
        Dist    = BA.Dist_Matrix(Coord)
        Cont    = BA.Contacts_Matrix(Dist, treshold)
        Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))
        BS      = BA.Print_Protein_BS_old(Bonds, TDP43.lenght,prot_initial=TDP43.initial, RNA_initial= RNA.initial)
        print(BS)


        WTC_traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path )
        WTC_traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = True, mode = 'BS', path = now_path)


        WTC_traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+'.xvg', 'RNA', path = now_path, skip_lines=17 )
        WTC_traj.Acquire_Atoms_List('rmsf_BS_'+now_name+'.xvg', 'BS', path = now_path, skip_lines = 17)
        WTC_traj.Get_RMSF(xvg_filename='rmsf_'+now_name+'.xvg', path = now_path,)


        res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+'.xvg', now_path)
        res = np.array(res, dtype=int)
        idx_RNA_start = np.where( res == 1.)[0][0]
        res[idx_RNA_start:] += res[idx_RNA_start-1]
        """
        idx_BS = np.zeros(len(BS['Prot']), dtype = int)
        for (ii,bs) in zip(range(len(BS)), BS):
            idx_BS[ii] = np.where(res == bs)[0]
        """
        cov_matrix_CAP = np.genfromtxt(now_path+now_name+'CAP_cov_matrix.txt') 
        N = cov_matrix_CAP.shape[0]
        pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP, N)

        WTC_traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA.xvg', now_path, skip_lines= 27 )


        # 2) data framing

        pearson_covariance = 'both'

        df = pd.DataFrame(res, columns=['Residue_number'])

        df['WTC_identifier'] = [now_name,]*len(res)
        df['Is_Prot'] = [True if ii < idx_RNA_start else False for ii in range(len(res))]
        df['Is_BS'] = [True if (r in BS['Prot']) | (r in BS['RNA']) else False for r in res]
        df['RMSF'] = RMSF_res

        if pearson_covariance == True:
            df['Covariance_Mean'] = [np.mean(pearson_matrix_CAP[:, ii]) for ii in range(len(res))]
            df['Covariance_Std'] = [np.std(pearson_matrix_CAP[:, ii]) for ii in range(len(res))]
        
        
        elif pearson_covariance == False:

            df['Covariance_Mean'] = [np.mean(cov_matrix_CAP[:, ii]) for ii in range(len(res))]
            df['Covariance_Std'] = [np.std(cov_matrix_CAP[:, ii]) for ii in range(len(res))]
        
        else:
            df['Pearson_Mean'] = [np.mean(pearson_matrix_CAP[:, ii]) for ii in range(len(res))]
            df['Pearson_Std'] = [np.std(pearson_matrix_CAP[:, ii]) for ii in range(len(res))]
            df['Covariance_Mean'] = [np.mean(cov_matrix_CAP[:, ii]) for ii in range(len(res))]
            df['Covariance_Std'] = [np.std(cov_matrix_CAP[:, ii]) for ii in range(len(res))]
        
        


        df['RMSD_Mean'] = [np.mean(WTC_traj.RMSD_eq), ]* len(res)
        df['RMSD_Std'] = [np.std(WTC_traj.RMSD_eq), ]*len(res)
        df['Gyradium_Mean']  = [np.mean(WTC_traj.Gyradium[WTC_traj.idx_eq:]),]*len(res)
        df['Gyradium_Std']  = [np.std(WTC_traj.Gyradium[WTC_traj.idx_eq:]),]*len(res)
        df['Kd'] = [Kds[now_name],]*len(res)

        #3) salvo

        df.to_json(now_path+now_name+'_df.json')
        dfs[now_name] = df


    
    DF = pd.concat([dfs[key] for key in dfs.keys()], ignore_index=True)
    DF.to_json('../GROMACS/WTC_data_frame_{}ang.json'.format(str(treshold)))
    DF.to_csv('../GROMACS/WTC_data_frame_{}ang.csv'.format(str(treshold)))

    # %%

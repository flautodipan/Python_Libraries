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
Kds_ = [4., 4., 1360., 1350, 700, 1320, 1800, 70]
Kds = {'wtc{}'.format(n): kd for (n, kd) in zip(['1', '1_h', '2', '3', '4', '5', '6', '7'], Kds_)}

#inizializzo dataframe

dfs = {}

now_temp = '300 K'
scale='ns'


#riferimento per lo z-score dato da studio wtc7_16


#seven reduced
red     = False
wtc7_7  = False
wtc7_8  = False
wtc7_16 = True
wtc7_24 = False
wtc7_24_1 = False

for treshold in [9]:

    for ii in ['1_h', '1','2', '3', '4', '5', '6', '7']:

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
                n_frames = 10001
                time_range = [0, 1000000]
                time_range_eq = [250000, 1000000]
            
            elif wtc7_24 == True:

                now_path    =   '../GROMACS/WTC7_24/'
                now_name    =    'wtc7_24'
                n_frames = 10001
                time_range = [0, 1000000]
                time_range_eq = [50000, 650000]

                        
            elif wtc7_24_1 == True:

                now_path    =   '../GROMACS/WTC7_24_1/'
                now_name    =    'wtc7_24'
                n_frames = 10001
                time_range = [0, 1000000]
                time_range_eq = [750000, 1000000]

            else:
            
                #WTC7
                now_path    =   '../GROMACS/WTC7_red/'
                now_name    =    'wtc7'
                n_frames = 1251
                time_range = [175000, 300000]
                time_range_eq = time_range

        #WTC analyzer senza figure

        WTC_traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, now_temp))
        WTC_traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )
        WTC_traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', path = now_path)
        WTC_traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq)


        res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+'.xvg', now_path, skip_lines= 17)
        res = np.array(res, dtype=int)
        idx_RNA_start = np.where( res == 1.)[0][0]
        res[idx_RNA_start:] += res[idx_RNA_start-1]

        filename='BS_{}_make_ndx.txt'.format(now_name)
        pdb_filename = 'average_pdb_'+now_name+'.pdb'

        TDP43   = BA.Protein(now_path+pdb_filename, model = False)
        TDP43.Get_CA_Coord(atom_name='CA')
        RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False, initial=TDP43.initial+TDP43.CA.shape[0])
        RNA.Get_P_Coord(atom_name="P")
        print('Ok, acquisito correttamente pdb')

        Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
        Dist    = BA.Dist_Matrix(Coord)
        Cont    = BA.Contacts_Matrix(Dist, treshold)
        Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))
        BS      = BA.Print_Protein_BS_old(res, Bonds, TDP43.lenght,prot_initial=TDP43.initial, RNA_initial= RNA.initial)
        print(BS)


        WTC_traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path )
        WTC_traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = True, mode = 'BS', path = now_path)


        WTC_traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+'.xvg', 'RNA', path = now_path, skip_lines=17 )
        WTC_traj.Acquire_Atoms_List('rmsf_BS_'+now_name+'.xvg', 'BS', path = now_path, skip_lines = 17)
        WTC_traj.Get_RMSF(xvg_filename='rmsf_'+now_name+'.xvg', path = now_path,)


        """
        idx_BS = np.zeros(len(BS['Prot']), dtype = int)
        for (ii,bs) in zip(range(len(BS)), BS):
            idx_BS[ii] = np.where(res == bs)[0]
        """
        print('\n\n{}\n\n'.format(now_path+now_name+'CAP_cov_matrix.txt'))
        cov_matrix_CAP = np.genfromtxt(now_path+now_name+'CAP_cov_matrix.txt')
        print('Covarianza media di {} è {}\n\n\n'.format(now_name, np.mean(cov_matrix_CAP))) 
        N = cov_matrix_CAP.shape[0]
        pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP, N)

        WTC_traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA.xvg', now_path, skip_lines= 27 )


        # 2) data framing


        df = pd.DataFrame(res, columns=['Residue_number'])

        df['WTC_identifier'] = ['wtc'+ii,]*len(res)
        df['Is_Prot'] = [True if ii < idx_RNA_start else False for ii in range(len(res))]
        df['Is_BS'] = [True if (r in BS['Prot']) | (r in BS['RNA']) else False for r in res]
        df['size(Prot_BS)'] = [len(BS['Prot'])]*len(res)
        df['RMSF'] = RMSF_res


        df['Pearson_Mean'] = [np.mean(pearson_matrix_CAP[:, ii]) for ii in range(len(res))]
        df['Pearson_Std'] = [np.std(pearson_matrix_CAP[:, ii]) for ii in range(len(res))]
        df['Covariance_Mean'] = [np.mean(cov_matrix_CAP[:, ii]) for ii in range(len(res))]
        df['Covariance_Std'] = [np.std(cov_matrix_CAP[:, ii]) for ii in range(len(res))]
 

        
        Covariance_Mean_Prot_BS = []
        Covariance_Mean_Prot_noBS = []
        Covariance_Mean_RNA_BS = []
        Covariance_Mean_RNA_noBS = []

        Pearson_Mean_Prot_BS = []
        Pearson_Mean_Prot_noBS = []
        Pearson_Mean_RNA_BS = []
        Pearson_Mean_RNA_noBS = []
    
    
        for kk in range(len(res)):
            Covariance_Mean_Prot_BS.append(np.mean([cov_matrix_CAP[jj,kk] for jj in list(np.array(BS['Prot']) - df.Residue_number[0])]))
            Covariance_Mean_Prot_noBS.append(np.mean([cov_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS['Prot'], BS['RNA']))]]))
            Covariance_Mean_RNA_BS.append(np.mean([cov_matrix_CAP[jj,kk] for jj in list(np.array(BS['RNA']) - df.Residue_number[0])]))
            Covariance_Mean_RNA_noBS.append(np.mean([cov_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS['RNA']]]))

            Pearson_Mean_Prot_BS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in list(np.array(BS['Prot']) - df.Residue_number[0])]))
            Pearson_Mean_Prot_noBS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS['Prot'], BS['RNA']))]]))
            Pearson_Mean_RNA_BS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in list(np.array(BS['RNA']) - df.Residue_number[0])]))
            Pearson_Mean_RNA_noBS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS['RNA']]]))


        df['Covariance_Mean_Prot_BS'] = Covariance_Mean_Prot_BS
        df['Covariance_Mean_Prot_noBS'] = Covariance_Mean_Prot_noBS

        df['Covariance_Mean_RNA_BS'] = Covariance_Mean_RNA_BS
        df['Covariance_Mean_RNA_noBS'] = Covariance_Mean_RNA_noBS
        
        #zscore
        df['z_Covariance_Mean'] = (df['Covariance_Mean'] - np.mean(df['Covariance_Mean']))/(np.std(df['Covariance_Mean']))
        df['z_Covariance_Mean_Prot_BS'] = (df.Covariance_Mean_Prot_BS - np.mean(df.Covariance_Mean_Prot_BS))/(np.std(df.Covariance_Mean_Prot_BS))

        #zscore con riferimento RNA12

        df['z_Covariance_Mean_12'] = (df['Covariance_Mean'] - mean_12)/std_12
        df['z_Covariance_Mean_Prot_BS_12'] = (df.Covariance_Mean_Prot_BS - mean_12)/std_12


        df['Pearson_Mean_Prot_BS'] = Pearson_Mean_Prot_BS
        df['Pearson_Mean_Prot_noBS'] = Pearson_Mean_Prot_noBS
        
        df['Pearson_Mean_RNA_BS'] = Pearson_Mean_RNA_BS
        df['Pearson_Mean_RNA_noBS'] = Pearson_Mean_RNA_noBS
        

        df['RMSD_Mean'] = [np.mean(WTC_traj.RMSD_eq), ]* len(res)
        df['RMSD_Std'] = [np.std(WTC_traj.RMSD_eq), ]*len(res)
        df['Gyradium_Mean']  = [np.mean(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right]),]*len(res)
        df['Gyradium_Std']  = [np.std(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right]),]*len(res)
        df['Kd'] = [Kds['wtc'+ii],]*len(res)

        #3) salvo

        #df.to_json(now_path+now_name+'_df.json')
        dfs[now_name] = df


    
    DF = pd.concat([dfs[key] for key in dfs.keys()], ignore_index=True)

    if wtc7_7 == True:
        DF.to_json('../GROMACS/WTC_data_frame_{}ang_7.json'.format(str(treshold)))
        DF.to_csv('../GROMACS/WTC_data_frame_{}ang_7.csv'.format(str(treshold)))

    elif wtc7_8 == True:
        DF.to_json('../GROMACS/WTC_data_frame_{}ang_8.json'.format(str(treshold)))
        DF.to_csv('../GROMACS/WTC_data_frame_{}ang_8.csv'.format(str(treshold)))
    elif wtc7_16:
        DF.to_json('../GROMACS/WTC_data_frame_{}ang_16_all.json'.format(str(treshold)))
        DF.to_csv('../GROMACS/WTC_data_frame_{}ang_16_all.csv'.format(str(treshold)))
    elif wtc7_24:
        DF.to_json('../GROMACS/WTC_data_frame_{}ang_24_all.json'.format(str(treshold)))
        DF.to_csv('../GROMACS/WTC_data_frame_{}ang_24_all.csv'.format(str(treshold)))

    elif wtc7_24_1:

        DF.to_json('../GROMACS/WTC_data_frame_{}ang_24_1_all.json'.format(str(treshold)))
        DF.to_csv('../GROMACS/WTC_data_frame_{}ang_24_1_all.csv'.format(str(treshold)))

    else:
        DF.to_json('../GROMACS/WTC_data_frame_{}ang.json'.format(str(treshold)))
        DF.to_csv('../GROMACS/WTC_data_frame_{}ang.csv'.format(str(treshold)))

    
    # %%
WTC_identifier = ('wtc1', 'wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6','wtc7')

#CHECK COVARIANZA
print('stampo le medie delle covarianze per ogni dinamica')
[np.mean(a) for a in [DF.z_Covariance_Mean_Prot_BS_12[DF.WTC_identifier == key] for key in WTC_identifier ]]
# %%


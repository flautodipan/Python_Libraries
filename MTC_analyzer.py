#%%
#0) preambolo

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


#exp data for MTC

Kds_= [529.1,]
Kds = {'MTC{}'.format(n): kd for (n, kd) in zip(['1'], Kds)}

#inizializzo dataframe

dfs = {}

# studio da wtc7_16
mean_12 = 8.552536674490916e-06
std_12 = 7.437340840337174e-05

now_temp = '300 K'
scale='ns'


for ii in ['1', ]:

    print('\n\nDinamica MTC{}\n\n'.format(ii))

    if ii == '1':

        now_path = '../GROMACS/MTC1/'
        now_name = 'MTC1'
        n_frames = 7043
        time_range = [0, 704200]
        time_range_eq = [200000, 704200]

    #WTC analyzer senza figure

        WTC_traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, now_temp))
        WTC_traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )
        WTC_traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', path = now_path)
        WTC_traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq)


        res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+'.xvg', now_path, skip_lines= 17)
        res = np.array(res, dtype=int)
        idx_RNA_start = np.where( res == 1.)[0][0]
        res[idx_RNA_start:] += res[idx_RNA_start-1]

        #prendo dati BS generati da WTC_analyzer

        BS      = {}
        BS['Prot'] = np.load(now_path+'BS.npy')
        BS['RNA'] = np.load(now_path+'BS_RNA.npy')
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

        df['MTC_identifier'] = ['MTC'+ii,]*len(res)
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
    DF.to_json('../GROMACS/MTC_data_frame.json')
    DF.to_csv('../GROMACS/MTC_data_frame.csv')


#%%
#0) preambolo

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")#inizializzo dataframe

Kds_ = [5.4, 1500, 115, 6.9 ]
Kds = {'cov{}'.format(n): kd for (n, kd) in zip(['1', '2', '3', '4'], Kds_)}

dfs = {}

now_temp = '300 K'
scale='ns'

dodici = 12

for ii in ['1', '2', '3', '4']:
    
    if ii == '1':

        now_path    =   '../GROMACS/COV1/'
        now_name    =    'cov1'
        n_frames = 2001
        time_range = [0, 200000]
        time_range_eq = time_range

        if dodici == 12:

            treshold = 12
        else:
            treshold = 9

    elif ii == '2':

        now_path    =   '../GROMACS/COV2/'
        now_name    =    'cov2'
        n_frames = 2001
        time_range = [0, 200000]
        treshold = 12

    elif ii == '3':

        now_path    =   '../GROMACS/COV3/'
        now_name    =    'cov3'
        n_frames = 2001
        time_range = [0, 200000]
        time_range_eq = time_range
        treshold = 12
    
    elif ii == '4':

        #COV4

        now_path    =   '../GROMACS/COV4/'
        now_name    =    'cov4'
        n_frames = 5001
        time_range = [0, 500000]
        time_range_eq = [180000, 500000]

        if dodici == 12 : treshold = 12
        else : treshold = 9


    
    print('\n\n\nDinamica  cov{} treshold = {}\n\n\n'.format(ii,treshold))


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
    RNA.Get_P_Coord(atom_name="05'")
    print('Ok, acquisito correttamente pdb')

    Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
    Dist    = BA.Dist_Matrix(Coord)
    Cont    = BA.Contacts_Matrix(Dist, treshold)
    Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))
    BS      = BA.Print_Protein_BS_old(res, Bonds, TDP43.lenght,prot_initial=TDP43.initial, RNA_initial= RNA.initial,)
    print(BS)               

    WTC_traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+'.xvg', 'RNA', path = now_path, skip_lines=17 )
    WTC_traj.Acquire_Atoms_List('rmsf_BS_'+now_name+'.xvg', 'BS', path = now_path, skip_lines = 17)
    WTC_traj.Get_RMSF(xvg_filename='rmsf_'+now_name+'.xvg', path = now_path,)


    print('\n\n{}\n\n'.format(now_path+now_name+'CAP_cov_matrix.txt'))
    cov_matrix_CAP = np.genfromtxt(now_path+now_name+'CAP_cov_matrix.txt')
    print('Covarianza media di {} è {}\n\n\n'.format(now_name, np.mean(cov_matrix_CAP))) 
    N = cov_matrix_CAP.shape[0]
    pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP, N)

    WTC_traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA.xvg', now_path, skip_lines= 27 )


    # 2) data framing


    df = pd.DataFrame(res, columns=['Residue_number'])

    df['cov_identifier'] = ['cov'+ii,]*(len(res))
    df['Is_Prot'] = [True if ii < idx_RNA_start else False for ii in range(len(res) )]
    df['Is_BS'] = [True if (r in BS['Prot']) | (r in BS['RNA']) else False for r in res]
    df['RMSF'] = RMSF_res


    df['Pearson_Mean'] = [np.mean(pearson_matrix_CAP[:, ii]) for ii in range(len(res))]
    df['Pearson_Std'] = [np.std(pearson_matrix_CAP[:, ii]) for ii in range(len(res) )]
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
        #Covariance_Mean_Prot_noBS.append(np.mean([cov_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS['Prot'], BS['RNA']))]]))
        Covariance_Mean_RNA_BS.append(np.mean([cov_matrix_CAP[jj,kk] for jj in list(np.array(BS['RNA']) - df.Residue_number[0])]))
        #Covariance_Mean_RNA_noBS.append(np.mean([cov_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS['RNA']]]))

        Pearson_Mean_Prot_BS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in list(np.array(BS['Prot']) - df.Residue_number[0])]))
        #Pearson_Mean_Prot_noBS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS['Prot'], BS['RNA']))]]))
        Pearson_Mean_RNA_BS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in list(np.array(BS['RNA']) - df.Residue_number[0])]))
        #Pearson_Mean_RNA_noBS.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS['RNA']]]))


    df['Covariance_Mean_Prot_BS'] = Covariance_Mean_Prot_BS
    #df['Covariance_Mean_Prot_noBS'] = Covariance_Mean_Prot_noBS


    mean_12 = 1.5320026220697674e-05
    std_12 = 8.303982395321321e-05

    df['Covariance_Mean_RNA_BS'] = Covariance_Mean_RNA_BS
    #df['Covariance_Mean_RNA_noBS'] = Covariance_Mean_RNA_noBS
    
    #zscore
    df['z_Covariance_Mean'] = (df['Covariance_Mean'] - np.mean(df['Covariance_Mean']))/(np.std(df['Covariance_Mean']))
    df['z_Covariance_Mean_Prot_BS'] = (df.Covariance_Mean_Prot_BS - np.mean(df.Covariance_Mean_Prot_BS))/(np.std(df.Covariance_Mean_Prot_BS))

    #zscore con riferimento RNA12

    df['z_Covariance_Mean_12'] = (df['Covariance_Mean'] - mean_12)/std_12
    df['z_Covariance_Mean_Prot_BS_12'] = (df.Covariance_Mean_Prot_BS - mean_12)/std_12


    df['Pearson_Mean_Prot_BS'] = Pearson_Mean_Prot_BS
    #df['Pearson_Mean_Prot_noBS'] = Pearson_Mean_Prot_noBS
    
    df['Pearson_Mean_RNA_BS'] = Pearson_Mean_RNA_BS
    #df['Pearson_Mean_RNA_noBS'] = Pearson_Mean_RNA_noBS
    

    df['RMSD_Mean'] = [np.mean(WTC_traj.RMSD_eq), ]* len(res)
    df['RMSD_Std'] = [np.std(WTC_traj.RMSD_eq), ]*len(res)
    df['Gyradium_Mean']  = [np.mean(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right]),]*len(res)
    df['Gyradium_Std']  = [np.std(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right]),]*len(res)
    df['Kd'] = [Kds['cov'+ii],]*len(res)

    #3) salvo

    df.to_json(now_path+now_name+'_df.json')
    dfs[now_name] = df



DF = pd.concat([dfs[key] for key in dfs.keys()], ignore_index=True)
DF.to_json('../GROMACS/cov_data_frame_{}.json'.format(str(dodici)))
DF.to_csv('../GROMACS/cov_data_frame_{}.csv'.format(str(dodici)))


# %%

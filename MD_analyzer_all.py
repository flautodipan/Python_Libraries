#%%
#0) preambolo

import  BioAlessandria as BA
from    Alessandria import Find_Nearest
import  numpy as np
import  warnings
import  pandas as pd
import  matplotlib.pyplot as plt
import  os
warnings.filterwarnings("ignore")


path = '../GROMACS/'
wtc_keys = ['wtc1_h', 'wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
cov_keys = ['cov1', 'cov2', 'cov3', 'cov4']
mtc_keys = ['mtc1']

exp_df = pd.read_excel(path+'MD_experimental_data.xlsx')

for now_name in wtc_keys+cov_keys[1:]+mtc_keys:

    # inputs
    now_exp_df = exp_df[exp_df.identifier == now_name]
    now_path        = path + now_name.upper() +'/'
    n_frames        = int(now_exp_df.n_frames.values[0])
    time_range      = eval(now_exp_df.time_range.values[0])
    time_range_eq   = eval(now_exp_df.time_range_eq.values[0])

    # genero traiettoria da RMSD
    traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, '300 K'))
    traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )
    traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', path = now_path)
    traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq)

    # genero residui da RMSF
    res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+'.xvg', now_path, skip_lines= 17)
    res = np.array(res, dtype=int)
    idx_RNA_start = np.where( res == 1.)[0][0]
    res[idx_RNA_start:] += res[idx_RNA_start-1]

    pdb_filename = 'average_pdb_'+now_name+'.pdb'
    protein   = BA.Protein(now_path+pdb_filename, model = False)
    RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False)
    prot_sequence = protein.pdb.amino3to1()[protein.pdb.amino3to1().chain_id == 'A'].residue_name.values
    RNA_sequence = RNA.atoms.residue_name[RNA.atoms.atom_name == "O5'"].values
    sequence = np.concatenate([prot_sequence, RNA_sequence])


    if len(sequence) != len(res):
        raise ValueError('Lunghezza residui presa da rmsf.svg {} non corrisponde a lunghezza sequenza presa da pdb {}'.format(len(res[:idx_RNA_start], len(sequence))))
    else:
        print('Ok, acquisito da pdb lista dei residui la cui lunghezza ({}) coincide con il numero di residui estratti dal pdb per {}'.format(len(res[:idx_RNA_start]), now_name))
    #prendo dati BS generati da WTC_analyzer
    if (os.path.exists(now_path+'BS.npy')) & (os.path.exists(now_path+'BS_RNA.npy')):
        BS      = {}
        BS = np.load(now_path+'BS.npy')
        BS_RNA = np.load(now_path+'BS_RNA.npy')
        print(BS)
    else : 
        print('Non hai salvato Binding Site in cartella {}'.format(now_path))
        protein.Get_CA_Coord()
        RNA.Get_P_Coord(atom_name="{}".format(RNA_main_res))
        Coord   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
        Dist    = BA.Dist_Matrix(Coord)
        Cont    = BA.Contacts_Matrix(Dist, now_exp_df.bs_treshold.values[0])
        Bonds   = BA.Analyze_Bond_Residues(Cont, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
        BS      = BA.Print_Protein_BS_old(res, Bonds, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
        BS_RNA  = BA.Print_Protein_BS_old(res, Bonds, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']
        np.save(now_path+'BS.npy', BS)
        np.save(now_path+'BS_RNA.npy', BS_RNA)


    traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path )
    traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = True, mode = 'BS', path = now_path)
    traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+'.xvg', 'RNA', path = now_path, skip_lines=17 )
    traj.Acquire_Atoms_List('rmsf_BS_'+now_name+'.xvg', 'BS', path = now_path, skip_lines = 17)
    traj.Get_RMSF(xvg_filename='rmsf_'+now_name+'.xvg', path = now_path,)

    # matrice di covarianza dei residui
    if RNA_main_res == 'O':
            
        cov_matrix_CAO = np.genfromtxt(now_path+now_name+'CAO_cov_matrix.txt')
        print('Covarianza media di {} è {}\n\n\n'.format(now_name, np.mean(cov_matrix_CAO))) 
        N = cov_matrix_CAO.shape[0]
        pearson_matrix_CAO = BA.Pearson_Matrix_from_Cov(cov_matrix_CAO, N)

    elif RNA_main_res == 'P':

        N = protein.atoms['atom_number'].size + RNA.atoms['atom_number'].size
        cov_matrix = BA.Get_Covariance_Matrix(N, 'cov_eq_'+now_name, now_path)
        # prendo gli indici atomici dei CA e P
        CA_idx = protein.CA['atom_number'].values - 1
        P_idx = RNA.P['atom_number'].values - 1 
        idx = np.concatenate((CA_idx, P_idx))
        cov_matrix_CAP = BA.CAP_Cov_Matrix(cov_matrix, idx, now_name+'CAP_cov_matrix.txt', now_path)
        pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP, N-1)


    traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA.xvg', now_path, skip_lines= 27 )

    # Salvo dataframe in ogni cartella
    # 2) data framing

    df = pd.DataFrame(res, columns=['Residue_number'])

    df['identifier'] = [now_name]*len(res)
    df['Is_Prot'] = [True if ii < idx_RNA_start else False for ii in range(len(res))]
    df['Is_BS'] = [True if (r in BS) | (r in BS_RNA) else False for r in res]
    df['size_BS'] = [len(BS)]*len(res)
    df['RMSF'] = RMSF_res
    df['residue_name'] = sequence 


    df['Pearson_Mean'] = [np.mean(pearson_matrix_CAO[:, ii]) for ii in range(len(res))]
    df['Pearson_Std'] = [np.std(pearson_matrix_CAO[:, ii]) for ii in range(len(res))]
    df['Covariance_Mean'] = [np.mean(cov_matrix_CAO[:, ii]) for ii in range(len(res))]
    df['Covariance_Std'] = [np.std(cov_matrix_CAO[:, ii]) for ii in range(len(res))]

    Covariance_Mean_Prot_BS = []
    Covariance_Mean_Prot_noBS = []
    Covariance_Mean_RNA_BS = []
    Covariance_Mean_RNA_noBS = []
    Pearson_Mean_Prot_BS = []
    Pearson_Mean_Prot_noBS = []
    Pearson_Mean_RNA_BS = []
    Pearson_Mean_RNA_noBS = []


    for kk in range(len(res)):

        Covariance_Mean_Prot_BS.append(np.mean([cov_matrix_CAO[jj,kk] for jj in list(np.array(BS) - df.Residue_number[0])]))
        Covariance_Mean_Prot_noBS.append(np.mean([cov_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS, BS_RNA))]]))
        Covariance_Mean_RNA_BS.append(np.mean([cov_matrix_CAO[jj,kk] for jj in list(np.array(BS_RNA) - df.Residue_number[0])]))
        Covariance_Mean_RNA_noBS.append(np.mean([cov_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS_RNA]]))

        Pearson_Mean_Prot_BS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in list(np.array(BS) - df.Residue_number[0])]))
        Pearson_Mean_Prot_noBS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS, BS_RNA))]]))
        Pearson_Mean_RNA_BS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in list(np.array(BS_RNA) - df.Residue_number[0])]))
        Pearson_Mean_RNA_noBS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS_RNA]]))


    df['Covariance_Mean_Prot_BS'] = Covariance_Mean_Prot_BS
    df['Covariance_Mean_Prot_noBS'] = Covariance_Mean_Prot_noBS
    df['Covariance_Mean_RNA_BS'] = Covariance_Mean_RNA_BS
    df['Covariance_Mean_RNA_noBS'] = Covariance_Mean_RNA_noBS

    """
    #zscore
    df['z_Covariance_Mean'] = (df['Covariance_Mean'] - np.mean(df['Covariance_Mean']))/(np.std(df['Covariance_Mean']))
    df['z_Covariance_Mean_Prot_BS'] = (df.Covariance_Mean_Prot_BS - np.mean(df.Covariance_Mean_Prot_BS))/(np.std(df.Covariance_Mean_Prot_BS))

    #zscore con riferimento RNA12

    df['z_Covariance_Mean_12'] = (df['Covariance_Mean'] - mean_12)/std_12
    df['z_Covariance_Mean_Prot_BS_12'] = (df.Covariance_Mean_Prot_BS - mean_12)/std_12

    """
    df['Pearson_Mean_Prot_BS'] = Pearson_Mean_Prot_BS
    df['Pearson_Mean_Prot_noBS'] = Pearson_Mean_Prot_noBS
    df['Pearson_Mean_RNA_BS'] = Pearson_Mean_RNA_BS
    df['Pearson_Mean_RNA_noBS'] = Pearson_Mean_RNA_noBS
    df['RMSD_Mean'] = [np.mean(traj.RMSD_eq), ]* len(res)
    df['RMSD_Std'] = [np.std(traj.RMSD_eq), ]*len(res)
    df['Gyradium_Mean']  = [np.mean(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]),]*len(res)
    df['Gyradium_Std']  = [np.std(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]),]*len(res)
    df['Kd'] = [now_exp_df.Kd[now_exp_df.identifier == now_name].values[0]]*len(res)
    df['Kd_err'] = [now_exp_df.Kd_err[now_exp_df.identifier == now_name].values[0]]*len(res)


    df.to_json(now_path+now_name+'_df.json')
    df.to_csv(now_path+now_name+'_df.csv')
#

# %%
# 8 )
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

skip_cov = False
path = '../GROMACS/'
wtc_keys = ['wtc1_h', 'wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
cov_keys = ['cov1', 'cov2', 'cov3', 'cov4']
mtc_keys = ['mtc1', 'mtc2', 'mtc3']

all_keys = wtc_keys[1:]+cov_keys[1:]+mtc_keys

exp_df = pd.read_excel(path+'MD_experimental_data.xlsx')
now_eqs1 = ['wtc1_h_new', 'mtc2', 'mtc3']
now_eqs2 = []

now_keys = mtc_keys[1:]

for now_name in now_keys:

    # inputs
    now_exp_df = exp_df[exp_df.identifier == now_name]
    now_path        = path + now_name.upper() +'/'
    n_frames        = int(now_exp_df.n_frames.values[0])
    time_range      = eval(now_exp_df.time_range.values[0])
    time_range_eq   = eval(now_exp_df.time_range_eq1.values[0]) if now_name in now_eqs1 else eval(now_exp_df.time_range_eq.values[0]) 
    eq = '_eq1' if now_name in now_eqs1 else '_eq'


    # 1 ) Genero traiettoria, acquisisco info RMSD e RMSF

    traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, '300 K'))
    traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )

    #RMSD 
    traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', path = now_path,)
    traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq, path = now_path,)

    #RMSF
    #acquisisco indici residui proteina da file .xvg dato da GROMACS
    res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+eq+'.xvg', now_path, skip_lines= 17)
    res = np.array(res, dtype=int)
    # faccio diventare indici dei residui RNA successivi a quelli proteina
    idx_RNA_start = np.where( res == 1.)[0][0]
    res[idx_RNA_start:] += res[idx_RNA_start-1]

    # 2) recupero info BS già pronte 
    # 2 b - ACQUISISCO PDB GENERATO da gmx trjconv e trovo BS
    #   
    filename='BS_{}_make_ndx'.format(now_name)
    pdb_filename = 'average_pdb_'+now_name+eq+'.pdb'
    treshold = exp_df.bs_treshold.values[0]

    protein   = BA.Protein(now_path+pdb_filename, model = False)
    RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False)

    print('Ok, acquisito correttamente pdb di proteina e RNA per {} eq = {}'.format(now_name, eq))

    protein.Get_Atom_Coord(atom_name = 'CA')
    protein.Get_lenght()
    protein.Get_Protein_Sequence()
    RNA.Get_Atom_Coord(atom_name="O5'")
    RNA.Get_Atom_Coord(atom_name='P')
    RNA.Get_RNA_Sequence()
    RNA.Get_lenght()
    sequence = np.concatenate([protein.sequence, RNA.sequence])

    BS_O5 = np.load(now_path+'BS_O5'+eq+'.npy')
    BS_RNA_O5 = np.load(now_path+'BS_RNA_O5'+eq+'.npy')

    Coord_P   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
    Dist_P    = BA.Dist_Matrix(Coord_P)
    Cont_P   = BA.Contacts_Matrix(Dist_P, treshold)
    Bonds_P   = BA.Analyze_Bond_Residues(Cont_P, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
    BS_P      = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
    BS_RNA_P = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']

    np.save(now_path+'BS_P'+eq+'.npy', BS_P)
    np.save(now_path+'BS_RNA_P'+eq+'.npy', BS_RNA_P)

    #2) Acquisisco gli altri RMSD generati 

    traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = False, mode = 'RNA', path = now_path,)
    traj.Get_RMSD('rmsd_'+now_name+'_BS'+eq+'.xvg',  equilibrium = False, mode = 'BS', path = now_path, )


    traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path, )
    traj.Get_RMSD('rmsd_'+now_name+'_BS'+eq+'.xvg',  equilibrium = True, mode = 'BS', path = now_path,)


    # 4 ) Gli altri RMSF


    traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+eq+'.xvg', 'RNA', path = now_path, skip_lines=17 )
    traj.Acquire_Atoms_List('rmsf_BS_'+now_name+eq+'.xvg', 'BS', path = now_path, skip_lines = 17)
    traj.Get_RMSF(xvg_filename='rmsf_'+now_name+eq+'.xvg', path = now_path, )

    # GYRADIUM
    traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA'+eq+'.xvg', now_path,)

    # 6 ) COVARIANCE ANALYSIS
    N = protein.atoms['atom_number'].size + RNA.atoms['atom_number'].size

    if os.path.exists(now_path+'cov'+eq+'_'+now_name+'.npy'):
        print('Sto prendendo una matrice covarianza di {} già salvata in file .npy\n Eq = {}'.format(now_name, eq))
        cov_matrix = np.load(now_path+'cov'+eq+'_'+now_name+'.npy')

    else:
        raise ValueError('Non esiste matrice cov{} per {}'.format(str(eq)[1:], now_name))

    #matrice covarianza residui (RNA con O5)
    CA_idx = protein.CA['atom_number'].values -1
    O5_idx = RNA.O5['atom_number'].values-1 
    idx_CAO = np.concatenate((CA_idx, O5_idx))
    cov_matrix_CAO = BA.CAP_Cov_Matrix(cov_matrix, idx_CAO, now_name+'_CAO_cov'+eq+'matrix.txt', now_path)

    #matrice covarianza residui (RNA con P)
    CA_idx = protein.CA['atom_number'].values-1
    P_idx = RNA.P['atom_number'].values-1 
    idx_CAP = np.concatenate((CA_idx, [O5_idx[0]], P_idx))
    cov_matrix_CAP = BA.CAP_Cov_Matrix(cov_matrix, idx_CAP, now_name+'_CAP_cov'+eq+'matrix.txt', now_path)

    #pearson
    print('Tutto ok ho fatto le matrici CAP e CAO cov, ora pearson per {}'.format(now_name))
    pearson_matrix_CAO = BA.Pearson_Matrix_from_Cov(cov_matrix_CAO)
    pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP)


    #
    #
    # DATA FRAMING
    #
    print('Data framing di {}'.format(now_name))
    df = pd.DataFrame(res, columns=['Residue_number'])
    df['residue_name'] = sequence 
    df['identifier'] = [now_name]*len(res)
    df['Is_Prot'] = [True if ii < idx_RNA_start else False for ii in range(len(res))]
    # binding site suddivisione P e O5
    df['Is_BS_O5'] = [True if (r in BS_O5) | (r in BS_RNA_O5) else False for r in res]
    df['size_BS_O5'] = [len(BS_O5)]*len(res)
    df['Is_BS_P'] = [True if (r in BS_P) | (r in BS_RNA_P) else False for r in res]
    df['size_BS_P'] = [len(BS_P)]*len(res)
    
    df['RMSF'] = RMSF_res

    # covarianza suddivido per O5 e P
    # O5
    size_O5 = len(res)
    df['Pearson_Mean_O5'] = [np.mean(pearson_matrix_CAO[:, ii]) for ii in range(size_O5)]
    df['Pearson_Std_O5'] = [np.std(pearson_matrix_CAO[:, ii]) for ii in range(size_O5)]
    df['Covariance_Mean_O5'] = [np.mean(cov_matrix_CAO[:, ii]) for ii in range(size_O5)]
    df['Covariance_Std_O5'] = [np.std(cov_matrix_CAO[:, ii]) for ii in range(size_O5)]
    # P
    size_P = len(res)-1
    df['Covariance_Mean_P'] = [np.mean(cov_matrix_CAP[:, ii]) for ii in range(size_O5)]
    df['Covariance_Std_P'] = [np.std(cov_matrix_CAP[:, ii]) for ii in range(size_O5)]
    df['Pearson_Mean_P'] = [np.mean(pearson_matrix_CAP[:, ii]) for ii in range(size_O5)]
    df['Pearson_Std_P'] = [np.std(pearson_matrix_CAP[:, ii]) for ii in range(size_O5)]
    
    # costruisco il futuro correlatore per O5

    Covariance_Mean_Prot_BS_O5 = []
    Covariance_Mean_Prot_noBS_O5 = []

    # costruisco il futuro correlatore per P

    Covariance_Mean_Prot_BS_P = []
    Covariance_Mean_Prot_noBS_P= []
    Pearson_Mean_Prot_BS_P = []
    Pearson_Mean_Prot_BS_O5 = []

    """
    ROBA CHE PER ORA NON SERVE PIU    
    tutto quello che sotto è commentato non serve più

    Covariance_Mean_RNA_BS = []
    Covariance_Mean_RNA_noBS = []
    Pearson_Mean_Prot_BS = []
    Pearson_Mean_Prot_noBS = []
    Pearson_Mean_RNA_BS = []
    Pearson_Mean_RNA_noBS = []
    """
    for kk in range(size_O5):

        Covariance_Mean_Prot_BS_O5.append(np.mean([cov_matrix_CAO[jj,kk] for jj in list(np.array(BS_O5) - df.Residue_number[0])]))
        Covariance_Mean_Prot_noBS_O5.append(np.mean([cov_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS_O5, BS_RNA_O5))]]))
        
        Covariance_Mean_Prot_BS_P.append(np.mean([cov_matrix_CAP[jj,kk] for jj in list(np.array(BS_P) - df.Residue_number[0])]))
        Covariance_Mean_Prot_noBS_P.append(np.mean([cov_matrix_CAP[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS_P, BS_RNA_P))]]))

        Pearson_Mean_Prot_BS_O5.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in list(np.array(BS_O5) - df.Residue_number[0])]))
        Pearson_Mean_Prot_BS_P.append(np.mean([pearson_matrix_CAP[jj,kk] for jj in list(np.array(BS_P) - df.Residue_number[0])]))

        
        """
        Covariance_Mean_RNA_BS.append(np.mean([cov_matrix_CAO[jj,kk] for jj in list(np.array(BS_RNA) - df.Residue_number[0])]))
        Covariance_Mean_RNA_noBS.append(np.mean([cov_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS_RNA]]))

        Pearson_Mean_Prot_BS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in list(np.array(BS) - df.Residue_number[0])]))
        Pearson_Mean_Prot_noBS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res if no_BS not in np.concatenate((BS, BS_RNA))]]))
        Pearson_Mean_RNA_BS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in list(np.array(BS_RNA) - df.Residue_number[0])]))
        Pearson_Mean_RNA_noBS.append(np.mean([pearson_matrix_CAO[jj,kk] for jj in [no_BS -df.Residue_number[0] for no_BS in res[idx_RNA_start:] if no_BS not in  BS_RNA]]))
        """

    df['Covariance_Mean_Prot_BS_O5'] = Covariance_Mean_Prot_BS_O5
    df['Covariance_Mean_Prot_noBS_O5'] = Covariance_Mean_Prot_noBS_O5


    df['Covariance_Mean_Prot_BS_P'] = Covariance_Mean_Prot_BS_P
    df['Covariance_Mean_Prot_noBS_P'] = Covariance_Mean_Prot_noBS_P

    df['Pearson_Mean_Prot_BS_O5'] = Pearson_Mean_Prot_BS_O5
    df['Pearson_Mean_Prot_BS_P'] = Pearson_Mean_Prot_BS_P

    df['RMSD_Mean'] = [np.mean(traj.RMSD_eq), ]* len(res)
    df['RMSD_Std'] = [np.std(traj.RMSD_eq), ]*len(res)
    df['Gyradium_Mean']  = [np.mean(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]),]*len(res)
    df['Gyradium_Std']  = [np.std(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]),]*len(res)
    df['Kd'] = [now_exp_df.Kd[now_exp_df.identifier == now_name].values[0]]*len(res)
    df['Kd_err'] = [now_exp_df.Kd_err[now_exp_df.identifier == now_name].values[0]]*len(res)


    df.to_json(now_path+now_name+'_df'+eq+'.json')
    df.to_csv(now_path+now_name+'_df'+eq+'.csv')

# %%

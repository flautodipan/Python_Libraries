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
mtc_keys = ['mtc1']

all_keys = wtc_keys+cov_keys[1:]+mtc_keys

exp_df = pd.read_excel(path+'MD_experimental_data.xlsx')
now_eqs1 = ['wtc1_h_new', 'mtc2', 'mtc3']
now_eqs2 = []

for now_name in all_keys:

    # inputs
    now_exp_df = exp_df[exp_df.identifier == now_name]
    now_path        = path + now_name.upper() +'/'
    n_frames        = int(now_exp_df.n_frames.values[0])
    time_range      = eval(now_exp_df.time_range.values[0])
    time_range_eq   = eval(now_exp_df.time_range_eq1.values[0]) if now_name in now_eqs1 else eval(now_exp_df.time_range_eq.values[0]) 
    eq = '_eq1' if now_name in now_eqs1 else '_eq'

    #colori
    color           = now_exp_df.color.values[0]
    darkcolor           = now_exp_df.darkcolor.values[0]
    brightcolor           = now_exp_df.brightcolor.values[0]
    contrastcolor           = now_exp_df.contrastcolor.values[0]
    darkcontrastcolor           = now_exp_df.darkcontrastcolor.values[0]

    # 1 ) Genero traiettoria, acquisisco info RMSD e RMSF

    traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, '300 K'))
    traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )

    #RMSD 
    traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', fig = now_name+'_RMSD', histo = now_name+'_RMSD_Histogram', bins = 50, path = now_path, color = color, ylim = (0,2))
    traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq, path = now_path, fig =  now_name+'_RMSD'+eq, alpha = 0.1, color = color, darkcolor = darkcolor, ylim = (0,2) )

    #RMSF
    #acquisisco indici residui proteina da file .xvg dato da GROMACS
    res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+eq+'.xvg', now_path, skip_lines= 17)
    res = np.array(res, dtype=int)
    # faccio diventare indici dei residui RNA successivi a quelli proteina
    idx_RNA_start = np.where( res == 1.)[0][0]
    res[idx_RNA_start:] += res[idx_RNA_start-1]

    # 2 ) TROVO CONFIGURAZIONE BINDING SITE
    # 
    # 2 a - TROVO CONFIGURAZIONE DI MEDIO RMSD

    mean_RMSD , mean_frame_idx = Find_Nearest(traj.RMSD_eq, np.mean(traj.RMSD_eq)) 
    mean_frame_time = ((mean_frame_idx+traj.idx_eq_left)*traj.timestep) - traj.initial_time
    print('Ho trovato il "centroide" della distribuzione RMSD a equilibrio\nnel frame {}\ncorrispondente al tempo {} ps\ncon RMSD = {} nm'.format(mean_frame_idx+traj.idx_eq_left, mean_frame_time, mean_RMSD))

    # 2 b - ACQUISISCO PDB GENERATO da gmx trjconv e trovo BS
    #   
    filename='BS_{}_make_ndx'.format(now_name)
    pdb_filename = 'average_pdb_'+now_name+eq+'.pdb'
    treshold = now_exp_df.bs_treshold.values[0]

    protein   = BA.Protein(now_path+pdb_filename, model = False)
    RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False)

    print('Ok, acquisito correttamente pdb')

    protein.Get_Atom_Coord(atom_name = 'CA')
    protein.Get_lenght()
    protein.Get_Protein_Sequence()
    RNA.Get_Atom_Coord(atom_name="O5'")
    RNA.Get_Atom_Coord(atom_name='P')
    RNA.Get_RNA_Sequence()
    RNA.Get_lenght()
    sequence = np.concatenate([protein.sequence, RNA.sequence])

    # Ora suddivido le BS in base all'atomo con cui rappresento 
    # 1) O5
    Coord_O5   = np.concatenate((protein.CA_Coord, RNA.O5_Coord))
    Dist_O5    = BA.Dist_Matrix(Coord_O5)
    Cont_O5   = BA.Contacts_Matrix(Dist_O5, treshold)
    Bonds_O5   = BA.Analyze_Bond_Residues(Cont_O5, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
    BS_O5      = BA.Print_Protein_BS_old(res, Bonds_O5, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
    BS_RNA_O5 = BA.Print_Protein_BS_old(res, Bonds_O5, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']
    np.save(now_path+'BS_O5.npy', BS_O5)
    np.save(now_path+'BS_RNA_O5.npy', BS_RNA_O5)
    with open  (now_path+filename+'_O5.txt', 'w') as f:
            f.write("# frame \t Protein Binding Site (BS) Residues at {} Å treshold\n".format(treshold))
            for bs in BS_O5:
                f.write('r {} | '.format(bs))
    # 2) P
    Coord_P   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
    Dist_P    = BA.Dist_Matrix(Coord_P)
    Cont_P   = BA.Contacts_Matrix(Dist_P, treshold)
    Bonds_P   = BA.Analyze_Bond_Residues(Cont_P, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
    BS_P      = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
    BS_RNA_P = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']
    np.save(now_path+'BS_P.npy', BS_P)
    np.save(now_path+'BS_RNA_P.npy', BS_RNA_P)
    with open  (now_path+filename+'_P.txt', 'w') as f:
            f.write("# frame \t Protein Binding Site (BS) Residues at {} Å treshold\n".format(treshold))
            for bs in BS_P:
                f.write('r {} | '.format(bs))
    # 3) Completo Acquisizioni e stampe degli altri RMSD e li plotto insieme

    #2) Acquisisco gli altri RMSD generati 

    traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = False, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA_tot', color = color, scale = 'ns', ylim = (0,2))
    traj.Get_RMSD('rmsd_'+now_name+'_BS'+eq+'.xvg',  equilibrium = False, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS_tot'+eq, color = color, scale = 'ns', ylim = (0,2))


    traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA', color = color, scale = 'ns')
    traj.Get_RMSD('rmsd_'+now_name+'_BS'+eq+'.xvg',  equilibrium = True, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS'+eq, color = color, scale = 'ns')


    reduction = 50
    sampling = np.arange(0, len(traj.RMSD_eq), reduction, dtype=int)
    time_steps = np.arange(traj.time_range_eq[0], traj.time_range_eq[1]+1, traj.timestep*reduction)/1000
    plt.figure()
    plt.title('Equilibrium RMSD comparison between {} groups'.format(now_name))
    plt.plot(time_steps, traj.RMSD_eq[sampling], '--', color = brightcolor, alpha = 0.8, label = 'Prot+RNA')
    plt.plot(time_steps, traj.RMSD_eq_RNA[sampling], color = color,  label = 'RNA')
    plt.plot(time_steps, traj.RMSD_eq_BS[sampling], '-.',  color = contrastcolor, alpha = 0.8, label = 'Prot BS')
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (nm)')
    plt.tight_layout()
    plt.savefig(now_path+'RMSD_comparison'+eq+'.pdf', format = 'pdf')

    # 4 ) Gli altri RMSF


    traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+eq+'.xvg', 'RNA', path = now_path, skip_lines=17 )
    traj.Acquire_Atoms_List('rmsf_BS_'+now_name+eq+'.xvg', 'BS', path = now_path, skip_lines = 17)
    traj.Get_RMSF(xvg_filename='rmsf_'+now_name+eq+'.xvg', path = now_path, fig = now_name+'_rmsf'+eq, color = color, darkcolor = contrastcolor, thirdcolor = brightcolor)

    #RMSF for residues
    text_size = 7
    for BS, atom in zip([BS_O5, BS_P], ['O5', 'P']):
        idx_BS = np.zeros(len(BS), dtype = int)
        for (ii,bs) in zip(range(len(BS)), BS):
            idx_BS[ii] = np.where(res == bs)[0]

        BS_split = []
        split = [BS[0],]

        for ii in range(len(BS)): 

            if ii < len(BS)-1 :
                if BS[ii+1] == BS[ii]+1:
                    split.append(BS[ii+1])

                else:
                    BS_split.append(split)
                    split = [BS[ii+1],]
            else:
                if BS[ii]  - 1 == BS[ii-1]:
                    split.append(BS[ii])
                    BS_split.append(split)

                else:
                    BS_split.append(split)
                    split = [BS[ii],]
                    BS_split.append(split)

        text = []   
        size = text_size
        start = 0

        for ii in range(int(BS.size/size)):
            text.append([str(v) for v in BS[start:start+size]])
            start+=size
        if BS.size%size != 0:
            final = [str(v) for v in BS[start:]]
            if len(final) != size:
                while True:
                    final.append(' ')
                    if len(final) == size:
                        break
                text.append(final)

        f, ax = plt.subplots(1,1)

        ax.stem(res[:idx_RNA_start], RMSF_res[:idx_RNA_start],  markerfmt = 'white', basefmt = 'silver' ,linefmt='silver')
        ax.stem(res[idx_RNA_start:], RMSF_res[idx_RNA_start:], markerfmt=darkcolor, basefmt=darkcolor, linefmt=color, label = 'RNA')
        for bs in BS_split:
            bs = np.array(bs, dtype = int)
            if bs[0] == BS_split[0][0]:
                ax.stem(bs, RMSF_res[np.where([res == b for b in bs])[1]], markerfmt=darkcontrastcolor, basefmt=darkcontrastcolor, linefmt=contrastcolor, label = 'BS')
            else:
                ax.stem(bs, RMSF_res[np.where([res == b for b in bs])[1]], markerfmt=darkcontrastcolor, basefmt=darkcontrastcolor, linefmt=contrastcolor,)


        ax.legend(title = 'RNA starts at res {}'.format(int(res[idx_RNA_start])))
        ax.set_title('RMSF {} for {} residues\nBS with {} atoms'.format(str(eq)[1:], now_name, atom), pad = 5)
        ax.set_xlabel('Residue number')
        ax.set_ylabel('RMSF (nm)')

        plt.tight_layout()
        f.savefig(now_path+'RMSF'+eq+'_res_'+now_name+'_{}.pdf'.format(atom), format = 'pdf', bbox_inches = 'tight')
    
    # GYRADIUM
    traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA'+eq+'.xvg', now_path, fig = now_name+'_gyradium'+eq, ylim = (0,2), alpha = 0.2, color = color, darkcolor = darkcolor, skip_lines= 27 )

    # 6 ) COVARIANCE ANALYSIS

    if not (os.path.exists(now_path+now_name+'_CAO_cov_matrix.txt') & os.path.exists(now_path+now_name+'_CAP_cov_matrix.txt')):

        #acquisisco matrice atomica
        N = protein.atoms['atom_number'].size + RNA.atoms['atom_number'].size

        if os.path.exists(now_path+'cov_eq_'+now_name+'.npy'):
            print('Sto prendendo una matrice covarianza già salvata in file .npy\n')
            cov_matrix = np.load(now_path+'cov'+eq+'_'+now_name+'.npy')
        else:
            cov_matrix = BA.Get_Covariance_Matrix(N, 'cov'+eq+'_'+now_name, now_path)
        
        BA.Print_Cov_Matrix_BS(cov_matrix, now_name,'Atoms', BS_O5, res, path = now_path, clim = (-0.005, 0.005))


        #matrice covarianza residui (RNA con O5)
        CA_idx = protein.CA['atom_number'].values -1
        O5_idx = RNA.O5['atom_number'].values-1 
        idx_CAO = np.concatenate((CA_idx, O5_idx))
        cov_matrix_CAO = BA.CAP_Cov_Matrix(cov_matrix, idx_CAO, now_name+'_CAO_cov'+eq+'matrix.txt', now_path)
        BA.Print_Cov_Matrix_BS(cov_matrix_CAO, now_name+eq, 'CAO', BS_O5, res, text = True, path = now_path, clim = (-0.005, 0.005))

        #matrice covarianza residui (RNA con P)
        P_idx = RNA.P['atom_number'].values-1 
        # barbatrucco: ci piazzo un indice in più e recupero la dimensionalità che perdo con P
        # a livello fisico/chimico il primo residuo di RNA è considerato con atomo O5, invece di P
        # però al livello della mia analisi (che riguarda prot BS), non dovrebbe cambiare troppo
        # SIINE consapevole + CAPISCI se forse è il caso di fare il barbatrucco anche da prima
        # così da non andare a definire una BS a partire da un nucleotide RNA in meno
        idx_CAP = np.concatenate((CA_idx, [O5_idx[0]], P_idx))
        cov_matrix_CAP = BA.CAP_Cov_Matrix(cov_matrix, idx_CAP, now_name+'_CAP_cov'+eq+'matrix.txt', now_path)
        BA.Print_Cov_Matrix_BS(cov_matrix_CAP, now_name+eq, 'CAP', BS_P, res, text = True, path = now_path, clim = (-0.005, 0.005))

    else:
        cov_matrix_CAO = np.genfromtxt(now_path+now_name+'_CAO_cov_'+eq+'matrix.txt')  
        cov_matrix_CAP = np.genfromtxt(now_path+now_name+'_CAP_cov_'+eq+'matrix.txt')  

    # FACCIO ANCHE LA MATRICE di PEARSON

    pearson_matrix_CAO = BA.Pearson_Matrix_from_Cov(cov_matrix_CAO)
    BA.Print_Cov_Matrix_BS(pearson_matrix_CAO, now_name, 'CAO', BS_O5, res, pearson = True, path = now_path, clim = (-1, 1))
    pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP)
    BA.Print_Cov_Matrix_BS(pearson_matrix_CAO, now_name, 'CAP', BS_P, res, pearson = True, path = now_path, clim = (-1, 1))


    #
    #
    # DATA FRAMING
    #

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


    """
    df['Covariance_Mean_RNA_BS'] = Covariance_Mean_RNA_BS
    df['Covariance_Mean_RNA_noBS'] = Covariance_Mean_RNA_noBS
    """

    """
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

    """ 

    df['RMSD_Mean'] = [np.mean(traj.RMSD_eq), ]* len(res)
    df['RMSD_Std'] = [np.std(traj.RMSD_eq), ]*len(res)
    df['Gyradium_Mean']  = [np.mean(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]),]*len(res)
    df['Gyradium_Std']  = [np.std(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]),]*len(res)
    df['Kd'] = [now_exp_df.Kd[now_exp_df.identifier == now_name].values[0]]*len(res)
    df['Kd_err'] = [now_exp_df.Kd_err[now_exp_df.identifier == now_name].values[0]]*len(res)


    df.to_json(now_path+now_name+'_df'+eq+'.json')
    df.to_csv(now_path+now_name+'_df'+eq+'.csv')

# %%

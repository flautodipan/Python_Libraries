#%%
#0) preambolo

import BioAlessandria as BA
from Alessandria import Find_Nearest, Check_Execution_Mode
import numpy as np
import os
from    os.path import join
import sys
import pandas as pd
import matplotlib.pyplot as plt

wtc_keys = ['wtc1_h', 'wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc5_new',  'wtc6', 'wtc7_16']
wtc_keys_red = [wtc_keys[1]]+wtc_keys[3:]
cov_keys = ['cov2', 'cov3', 'cov4']
mtc_keys = ['MTC1', 'mtc2', 'mtc3']
all_keys = wtc_keys+cov_keys+mtc_keys
gian_keys = [wtc_keys[1]] + wtc_keys[3:]

# IMPOSTAZIONI INIZIALI

path = '../GROMACS/'

treshold = 9
exp_df = pd.read_excel(join(path,'MD_experimental_data.xlsx'))
now_keys = ['wtc1']
now_eqs1 = []
now_eqs2 = []


#%%
for key in now_keys:

    now_path = join(path, key.upper(),)
    eq = '_eq2' if key in now_eqs2 else '_eq1' if key in now_eqs1 else '_eq' 
    eq_path = join(path, key.upper(), '{}_frame_'.format(eq[1:])+key, )
    now_exp_df = exp_df[exp_df.identifier == key]

    print('Key = {}'.format(key))
    print('Eq  = {}'.format(eq[1:]))

    #Residue_index
    res, RMSF_res = BA.Parse_xvg_skip(join(path, key.upper(), 'rmsf_res_'+key+eq+'.xvg'), path, skip_lines= 17)
    res = np.array(res, dtype=int)
    # faccio diventare indici dei residui RNA successivi a quelli proteina
    idx_RNA_start = np.where( res == 1.)[0][0]
    res = res[:idx_RNA_start]

    # PROCEDURA EQFRAME (4)

    #1) acquisisco info generali Prot
    protein   = BA.Protein(join(now_path,  'average_pdb_'+key+eq+'.pdb'), model = False)
    protein.Get_Atom_Coord(atom_name = 'CA')
    protein.Get_Protein_Sequence()
    if len(protein.sequence) != len(res): raise ValueError ('Lughezza proteina non combaciante # residui')
    
    df = pd.DataFrame()
    df['res_number'] = res
    df['res_name'] = protein.sequence
    
    # 2) analizzo contatti
    contacts = []
    print('Tot number of equilibrium frame for {} is {}'.format(key, len(os.listdir(eq_path))))
    print("Comincio procedura analisi contatti all'equilibrio per {} ".format(key))

    for frame in os.listdir(eq_path):

        if frame.endswith('.pdb'):
            protein   = BA.Protein(join(eq_path,frame), model = False)
            RNA     =  BA.RNA(join(eq_path,frame), chain_id='B', model = False)
            protein.Get_Atom_Coord(atom_name = 'CA')
            protein.Get_lenght()
            protein.Get_Protein_Sequence()
            RNA.Get_Atom_Coord(atom_name='P')
            RNA.Get_Atom_Coord(atom_name="O5'")
            RNA.Get_RNA_Sequence()
            RNA.Get_lenght()
            RNA_start_idx = len(protein.CA) - 1
            print('Ok, acquisito correttamente pdb {} per {}'.format(frame, key))

            Coord_P   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
            Dist_P    = BA.Dist_Matrix(Coord_P)
            Cont_P   = BA.Contacts_Matrix(Dist_P, treshold)
            Bonds_P   = BA.Analyze_Bond_Residues(Cont_P, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
            BS_P      = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
            BS_RNA_P = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']
        
            contacts.append(BS_P)

    # 3) rendo contatti istogramma

        counts = np.zeros(res.size)
        for r, ii in zip(res, range(len(res))):
            how_many = 0
            for c in contacts:
                if r in c: how_many += 1
                else: pass
            counts[ii] = how_many
            
        df['contacts'] = counts


        # 4) salvo in rispettiva cartella

        df.to_json(join(now_path, 'df_{}_contacts_{}.json'.format(eq[1:],key)))
        df.to_csv(join(now_path, 'df_{}_contacts_{}.csv'.format(eq[1:],key)))
        print('Ok ho salvato i risultati in {}'.format(now_path))


        # 5) visualizzo 

        f, ax = plt.subplots()

        ax.plot(res, df.contacts.values, c = now_exp_df.color.values[0], label = '# contacts')

        ax.set_xticks(label)
        ax.set_title('Number of contacts vs residue number')

        # 6) identifico e salvo i residui più importanti


#%%
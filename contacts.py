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

treshold = 9
path = '../GROMACS/'
now_keys = ['wtc1']


#%%

res = np.arange(96, 270, 1)
how_many_contacts = {}
numbers = {}
df = pd.DataFrame(index=res)


for key in WTC_identifier:

    contacts = []
    now_path = join(path, key.upper(), 'eq_frame_'+key, )
    #now_path = path+key.upper()+'/' if key != 'wtc1_h' else '../GROMACS/WTC1_h/eq_frame_wtc1_h/'
    numbers[key] = len(os.listdir(now_path))
    print('Tot number of equilibrium frame for {} is {}'.format(key, len(os.listdir(now_path))))

    for frame in os.listdir(now_path):

        protein   = BA.Protein(join(now_path,frame), model = False)
        RNA     =  BA.RNA(join(now_path,frame), chain_id='B', model = False)
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

    counts = np.zeros(res.size)
    for r, ii in zip(res, range(len(res))):
        how_many = 0
        for c in contacts:
            if r in c: how_many += 1
            else: pass
        counts[ii] = how_many
        
    how_many_contacts[key] = counts
    df[key] = how_many_contacts[key]/numbers[key]

df.to_json(join(path,'df_contacts_9ang_16.json'))
df.to_csv(join(path, 'df_contacts_9ang_16.csv'))
# %%
f, ax = plt.subplots(1,1)
for key in WTC_identifier:
    ax.plot(res, df[key], color = 'firebrick' if '7' in key else 'gray', label = key if '7' in key else None)
ax.set_xlabel('Chain element')
ax.set_ylabel('Contact frequency')
ax.set_title('Comparison of contact frequency')
plt.legend()
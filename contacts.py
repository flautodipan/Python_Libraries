#%%
#0) preambolo

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

terminal = False if sys.argv[1] == '-f' else True

if terminal:
    print(sys.argv[1])
    print('\n\n Eseguo programma in modalità da terminale\n\n')
else:
    print(sys.argv[1])
    print('\n\n Eseguo programma in modalità interattiva jupyter\n\n')
treshold = 9
WTC_identifier = ('wtc1_h', 'wtc2', 'wtc3', 'wtc4',  'wtc5', 'wtc6','wtc7_16')
if terminal:
    path = sys.argv[1]
else:
    path = '../GROMACS/'


#%%

res = np.arange(96, 270, 1)
how_many_contacts = {}
numbers = {}
df = pd.DataFrame(index=res)


for key in WTC_identifier:

    contacts = []
    
    now_path = path+key.upper()+'/eq_frame_'+key+'/' if key != 'wtc1_h' else '../GROMACS/WTC1_h/eq_frame_wtc1_h/'
    numbers[key] = len(os.listdir(now_path))
    print('Tot number of equilibrium frame for {} is {}'.format(key, len(os.listdir(now_path))))

    for frame in os.listdir(now_path):

        TDP43   = BA.Protein(now_path+frame, model = False)
        TDP43.Get_CA_Coord(atom_name='CA')
        RNA     =  BA.RNA(now_path+frame, chain_id='B', model = False, initial=TDP43.initial+TDP43.CA.shape[0])
        RNA.Get_P_Coord(atom_name="P")
        RNA_start_idx = len(TDP43.CA) - 1

        print('Ok, acquisito correttamente pdb {} per {}'.format(frame, key))

        Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
        Dist    = BA.Dist_Matrix(Coord)
        Cont    = BA.Contacts_Matrix(Dist, treshold)
        Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))
        BS      = BA.Print_Protein_BS_old(Bonds, TDP43.lenght,prot_initial=TDP43.initial, RNA_initial= RNA.initial)
        
        contacts.append(BS['Prot'])
    counts = np.zeros(res.size)
    for r, ii in zip(res, range(len(res))):
        how_many = 0
        for c in contacts:
            if r in c: how_many += 1
            else: pass
        counts[ii] = how_many
        
    how_many_contacts[key] = counts
    df[key] = how_many_contacts[key]/numbers[key]

df.to_json(path+'df_contacts_9ang_16.json')
df.to_csv(path+'df_contacts_9ang_16.csv')
# %%
f, ax = plt.subplots(1,1)
for key in WTC_identifier:
    ax.plot(res, df[key], color = 'firebrick' if '7' in key else 'gray', label = key if '7' in key else None)
ax.set_xlabel('Chain element')
ax.set_ylabel('Contact frequency')
ax.set_title('Comparison of contact frequency')
plt.legend()
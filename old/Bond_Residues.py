#%%
# Script per capire regione di interazione tra domini RRMs della TDP43 e RNA1 

import os
import sys
import BioAlessandria as BA
import numpy as np
import warnings
warnings.filterwarnings("ignore")

#%%
pdb_filename    =   sys.argv[1]
dir_filename    =   'WTC1'
txt_filename    =   'BS_'+pdb_filename

os.system('mkdir zerk')
path            =   '../MD/'+dir_filename+'/zerk/'
save_path       =   path+'BS/'



#%%
#1) importo il pdb
#path        =   '../MD/WTC1/'
#pdb_filename = 'zerk_initial1.pdb'
TDP43 = BA.Protein(pdb_filename)
RNA1  =  BA.RNA(pdb_filename, chain_id='B')
print('Ok, acquisito correttamente pdb')


# %%

# costruisco matrice distanze

TDP43.Get_CA_Coord()
RNA1.Get_P_Coord()

Coord = np.concatenate((TDP43.CA_Coord, RNA1.P_Coord))
Dist = BA.Dist_Matrix(Coord)

# %%
#vado coi contatti in threshold a 1 angstorm
# e identifico 

for treshold in range(6,12,1):
    Cont = BA.Contacts_Matrix(Dist, treshold)#,  fig='Protein-RNA1')
    Bonds = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA1.lenght), ("TDP43", "RNA1"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))#, fig = True, verbose = True)

    print("La proteina e l'RNA a distanza %3.2f fanno %d legami" %(treshold, len(Bonds)))

# %%

BA.Print_Bonds_HDOCK(Bonds, 'contact.txt', 11., TDP43.initial)


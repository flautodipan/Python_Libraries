#%%
# Script per capire regione di interazione tra domini RRMs della TDP43 e RNA1 

import BioAlessandria as BA
import numpy as np
import warnings
warnings.filterwarnings("ignore")


#%%
#1) importo il pdb

TDP43 = BA.Protein('4bs2.pdb')
RNA1  =  BA.RNA('4bs2.pdb', chain_id='B')

# %%
# costruisco matrice distanze

TDP43.Get_CA_Coord()
RNA1.Get_P_Coord()

Coord = np.concatenate((TDP43.CA_Coord, RNA1.P_Coord))
Dist = BA.Dist_Matrix(Coord)

# %%
#vado coi contatti in threshold a 1 angstorm
# e identifico 


Cont = BA.Contacts_Matrix(Dist, 11.,  fig='Protein-RNA1')
Bonds = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA1.lenght), ("TDP43", "RNA1"), TDP43.initial)


# %%


# %%

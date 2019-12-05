#%%
# Script per capire regione di interazione tra domini RRMs della TDP43 e RNA1 

import BioAlessandria as BA
import numpy as np

#%%
#1) importo il pdb

TDP43 = BA.Protein('4bs2.pdb')

RNA1  =  BA.RNA('4bs2.pdb', chain_id='B')


# %%

BA.BioStructure('"4bs2.pdb"', 2, ('TDP43', 'RNA1'), ('Protein', 'RNA'), ('"A"', '"B"'), (1,1))

# %%
# costruisco matrice distanze

TDP43.Get_CA_Coord()
RNA1.Get_P_Coord()

Coord = np.concatenate((TDP43.CA_Coord, RNA1.P_Coord))
Dist = BA.Dist_Matrix(Coord)
# %%
#vado coi contatti in threshold a 9 angstorm

Cont = BA.Contacts_Matrix(Dist, 9.)

# %%

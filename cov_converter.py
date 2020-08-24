import BioAlessandria as BA
import numpy as np




filename='BS_{}_make_ndx.txt'.format(now_name)
pdb_filename = 'average_pdb_'+now_name+'.pdb'
treshold = 9

TDP43   = BA.Protein(now_path+pdb_filename, model = False)
RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False)

print('Ok, acquisito correttamente pdb')

TDP43.Get_CA_Coord()
RNA.Get_P_Coord()

Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
Dist    = BA.Dist_Matrix(Coord)

Cont    = BA.Contacts_Matrix(Dist, treshold)
Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))

BS      = BA.Print_Protein_BS_old(Bonds, TDP43.lenght, prot_initial=TDP43.initial, RNA_initial=RNA.initial)['Prot']
BS_RNA  = BA.Print_Protein_BS_old(Bonds, TDP43.lenght, prot_initial=TDP43.initial, RNA_initial=RNA.initial)['RNA']

print('Ok analizzato correttamente BA')

N = TDP43.atoms['atom_number'].size + RNA.atoms['atom_number'].size
cov_matrix = BA.Get_Covariance_Matrix(N, 'cov_eq_'+now_name, now_path)
BA.Print_Cov_Matrix_BS(cov_matrix, now_name,'Atoms', BS, res, path = now_path, clim = (-0.005, 0.005))

CA_idx = TDP43.CA['atom_number'].values -1
P_idx = RNA.P['atom_number'].values -1 
idx = np.concatenate((CA_idx, P_idx))

cov_matrix_CAP = BA.CAP_Cov_Matrix(cov_matrix, idx, now_name+'CAP_cov_matrix.txt', now_path)

BA.Print_Cov_Matrix_BS(cov_matrix_CAP, now_name, 'CA', BS, res, text = True, path = now_path, clim = (-0.1, 0.1))


print('Ok analizzato correttamente BA')

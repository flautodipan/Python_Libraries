#%%
import BioAlessandria as BA
import numpy as np
import warnings
warnings.filterwarnings("ignore")

now_path    =   '../MD/WTC1/zerk/'
filename    =   'Protein_BS_WTC1.txt'
treshold    =   10 #angstrom
#%%

#f   =   open(now_path+filename, 'w')
with open  (now_path+filename, 'w') as f:

    f.write("# frame \t Protein Binding Site (BS) Residues\n")

    for ii in range (0,600,1):
        
        pdb_filename =   'frame_'+str(ii)+'.pdb'
        TDP43   = BA.Protein(now_path+pdb_filename)
        RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B')
        
        print('Ok, acquisito correttamente pdb n %d'%(ii))
        

        TDP43.Get_CA_Coord()
        RNA.Get_P_Coord()

        Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
        Dist    = BA.Dist_Matrix(Coord)

        Cont    = BA.Contacts_Matrix(Dist, treshold)
        Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))

        BS      = BA.Print_Protein_BS_old(Bonds, TDP43.lenght, initial=TDP43.initial)

        line    = "{}\t{}\n".format(ii,np.array2string(BS, max_line_width= 100))
        f.write(line)

#f.close()





# %%

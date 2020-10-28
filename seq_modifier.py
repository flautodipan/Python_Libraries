#%%
# SCRIPT utile per le mutazioni di TDP43
# 0) crea in /home/Desktop/GROMACS/ una cartella MDCi (MoDified Complex dove i varia)

# 1) prende il pdb di riferimento
# 2) ne estrae la proteina TDP43 
# 3) modifica la sequenza in base a volontà 
# 4) fa girare il programma Scwrl4 per ottenere pdb della proteina mutata

import BioAlessandria as BA
import os

pdb_path = '../GROMACS/'
pdb_filename = 'TDP43.pdb'

now_name = 'MTC4'


#check: già esiste 
if not os.path.exists(pdb_path+now_name):

    os.system('cd {} && mkdir {}'.format(pdb_path, now_name))
else:
    raise ValueError('Attenzione il nome {} per cartella già esiste, cambia!'.format(now_name))

now_path = pdb_path + now_name+'/'


print('PDB inserito da prendere come input per la mutazione: {}'.format(pdb_filename[:-4]))

#%%


protein = BA.Protein(pdb_path+pdb_filename, model = False)

sequence = protein.pdb.amino3to1()
#lista di char
sequence = sequence[sequence.chain_id == 'A']
# stringa unica con tutte maiuscole
up_sequence = ''.join(sequence.residue_name)
# stringa unica con tutte minuscole
low_sequence = up_sequence.lower()


# %%
# effettuo le mutazioni in sequenza 

residues_to_mutate= [105, 254]
old_residues = ['d', 's']
mutations = ['a', 'a']

new_sequence = low_sequence

for idx, r, r_old in zip(residues_to_mutate, mutations, old_residues):

    new_sequence = new_sequence[:idx-protein.initial]+r+new_sequence[idx-protein.initial+1:]

    old = low_sequence[idx-protein.initial]


    # check: lettera presa da pdb è quella che mi aspetto
    if old == r_old:
        print('Ok il residuo {} acquisito dal pdb è effettivamente {}'.format(idx, old))
    new = new_sequence[idx-protein.initial]

    #check: ho sostituito bene

    if new_sequence[idx-protein.initial] != r:
        raise ValueError('Non sono riuscito a sostituire in modo appropriato')
    else:
        print('Ho sostituito il residuo {} di tipo {} con {}\n{}'.format(idx, old, new, old+str(idx)+new))

    
    # (specifico per software SCWRL4) rendo maiuscole le lettere vicine alle mutazioni e delle mutazione, 
    # significa che iol software farà fluttuare le posizioni più facilmente

    new_sequence = new_sequence[:idx-protein.initial-1] + new_sequence[idx-protein.initial-1].upper() + new_sequence[idx-protein.initial-1+1:]
    new_sequence = new_sequence[:idx-protein.initial] + new_sequence[idx-protein.initial].upper() + new_sequence[idx-protein.initial+1:]
    new_sequence = new_sequence[:idx-protein.initial+1] + new_sequence[idx-protein.initial+1].upper() + new_sequence[idx-protein.initial+1+1:]

with open(now_path+'seq_modified_{}.txt'.format(now_name), 'w') as out:
    out.write('{}'.format(new_sequence))

print('Ok, stampata la stringa contentente nuova sequenza nel file "seq_modified.txt in {}'.format(now_path))

#%%
#check visivo

for ii in range(len(low_sequence)):
    print('{} {} {}'.format(ii+protein.initial, low_sequence[ii], new_sequence[ii]))


# %%
# lancio il programma Scwrl4 


os.system('cd {} && Scwrl4 -i ../TDP43.pdb -o TDP43_{}.pdb -s seq_modified_{}.txt'.format(now_path, now_name, now_name))

#%%
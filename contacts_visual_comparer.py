#%%
#0) preambolo

import  BioAlessandria      as BA
from    Alessandria         import Find_Nearest, Check_Execution_Mode
import  numpy               as np
import  os
from    os.path             import join
import  sys
import  pandas              as pd
import  matplotlib.pyplot   as plt
import  seaborn             as sns

wtc_keys = ['wtc1','wtc1_h_new', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
cov_keys = ['cov2', 'cov3', 'cov4']
mtc_keys = ['MTC1', 'mtc2', 'mtc3', 'mtc4']
all_keys = wtc_keys+cov_keys+mtc_keys
gian_keys = [wtc_keys[1]] + wtc_keys[3:]

# IMPOSTAZIONI INIZIALI
BS_treshold = 12
path = join('..', 'GROMACS')
now_keys = wtc_keys +mtc_keys

now_eqs1 = ['wtc1_h_new', 'MTC1',  'mtc3']
now_eqs2 = [ 'mtc2']
#%%
# 1) Acquisisco DataFrame con dati sui contatti e binding site
exp_df = pd.read_excel(join(path,'MD_experimental_data.xlsx'))
    #gacquisizione

dfs = {}
BSs = {}
eq_len = {}

for key in now_keys:  

    now_path    = join(path, key.upper())
    now_exp_df  = exp_df[exp_df.identifier == key]


    eq = '_eq2' if key in now_eqs2 else '_eq1' if key in now_eqs1 else '_eq' 
    if ((key in now_eqs1) & (key in now_eqs2)):
        raise ValueError('{} in both eq1 and eq2'.format(key))
    # /1000 perché frame ogni 100 ps, eq_sampling ogni 10
    eq_len[key]  = int((eval(getattr(now_exp_df, 'time_range'+eq).values[0])[1] - eval(getattr(now_exp_df, 'time_range'+eq).values[0])[0])/1000)
    if key == 'wtc1' : exp_eq_len = eq_len[key]

    BSs[key] = np.load(join(now_path, 'BS_P'+eq+'.npy'))
    dfs[key] = pd.read_json(join(now_path, 'df'+eq+'_contacts_'+key+'.json'))
    print('Acquisisco {} con eq = {}'.format(key, eq, ))


# %%

#2) RESIDUI PIU IMPORTANTI DI WTC1 e confronto con pdb iniziali 

treshold = 300 #valore basato su pdf a partire da distribuzione valori




truth = {}
#if 'wtc1' in now_keys : now_keys.remove('wtc1')

for key in now_keys:

    if key == 'wtc1':

        exp_contacts        = dfs[key].res_number[dfs[key].contacts > treshold]
        exp_contacts_name   = dfs[key].res_name[dfs[key].contacts > treshold]
        print('Ho trovato {} residui per la dinamica NMR che hanno contatti in più di {} frame su {}\n\n'.format(len(exp_contacts), treshold, int(exp_eq_len)))

    else:
        truth[key] = [True if bs in exp_contacts.values else False for bs in BSs[key]]
        print('La configurazione iniziale di {} ritrova\n{} contatti sperimentali su {}\n'.format(key, sum(truth[key]), len(exp_contacts)))


#%%
#3) stampo figure stemplot

for key in now_keys:  

    now_path    = join(path, key.upper())
    now_exp_df  = exp_df[exp_df.identifier == key]

    eq = '_eq2' if key in now_eqs2 else '_eq1' if key in now_eqs1 else '_eq' 
    if ((key in now_eqs1) & (key in now_eqs2)):
        raise ValueError('{} in both eq1 and eq2'.format(key))

    #acquisizione
    f, ax = plt.subplots()
    
    ax.stem(dfs[key].res_number.values, dfs[key].contacts.values, markerfmt = now_exp_df.color.values[0], basefmt = now_exp_df.color.values[0], linefmt = now_exp_df.color.values[0])
    y_lims = ax.get_ylim()
    [ax.vlines(bs, y_lims[0], y_lims[1], colors = now_exp_df.contrastcolor.values[0], alpha = 0.4, label = 'Experimental BS' if bs == exp_contacts.values[0] else None) for bs in exp_contacts]

    ax.set_title('{} Number of contacts vs residue number\n Treshold = {} ang   # tot frame = {} '.format(key, BS_treshold, eq_len[key]))
    ax.set_xlabel('Protein Residue Number')
    ax.set_ylabel('# of contacts')

    plt.legend()
    f.savefig(join(now_path,key+'_contacts'+eq+'_{}.pdf'.format(BS_treshold)), format = 'pdf')

# %%

# 4) Matrice dei contatti
exp_key = 'wtc1'
now_keys_discarded_ = ['wtc1_h_new', 'wtc7_16', 'MTC1', 'mtc2', 'mtc3', 'mtc4']
now_keys_discarded = [ key + '_discarded' for key in now_keys_discarded_]

prot_res =  dfs[exp_key].res_number.values
df = pd.DataFrame(index =prot_res )
df_dist = pd.DataFrame(index =prot_res )

NOW_KEYS = now_keys+now_keys_discarded
NOW_KEYS.sort()

for key in NOW_KEYS:

    # differenzio key real key, perché i discarded prendono file da real key
    real_key = key[:-10] if '_discarded' in key else key
    now_path = join(path, real_key.upper())
    now_pdb = key+'.pdb'
    eq = '_eq2' if key in now_eqs2 else '_eq1' if real_key in now_eqs1 else '_eq' 
    if ((real_key in now_eqs1) & (real_key in now_eqs2)):
        raise ValueError('{} in both eq1 and eq2'.format(real_key))
    
    print('\n\nNuovo ciclo path = {} \nchiave = {}, chiave reale = {}\n\n'.format(now_path, key, real_key))


    # acquisisco RMSF per res e frame iniziale 
    res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+real_key+eq+'.xvg', now_path, skip_lines= 17)
    res = np.array(res, dtype=int)
    # faccio diventare indici dei residui RNA successivi a quelli proteina
    idx_RNA_start = np.where( res == 1.)[0][0]
    res[idx_RNA_start:] += res[idx_RNA_start-1]

    protein   = BA.Protein(os.path.join(now_path, now_pdb), model = False)
    RNA     =  BA.RNA(os.path.join(now_path, now_pdb), chain_id='B' if real_key == 'wtc1' else '', model = False)

    print('Ok, acquisito frame iniziale {} da {}\n'.format('discarded' if 'disc' in key else '', join(now_path, now_pdb)))


    protein.Get_Atom_Coord(atom_name = 'CA')
    protein.Get_lenght()
    #protein.Get_Protein_Sequence()
    RNA.Get_Atom_Coord(atom_name='P')
    #RNA.Get_RNA_Sequence()
    RNA.Get_lenght()
    #sequence = np.concatenate([protein.sequence, RNA.sequence])
    Coord_P   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
    Dist_P    = BA.Dist_Matrix(Coord_P)
    Cont_P   = BA.Contacts_Matrix(Dist_P, BS_treshold)
    Bonds_P   = BA.Analyze_Bond_Residues(Cont_P, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
    BS_P      = BA.Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']

    if key == exp_key:
        df[key] = [1 if r in exp_contacts.values else 0 for r in prot_res]
        df[key+'_ini'] = [1 if r in BS_P else 0 for r in prot_res]
        df_dist[key] = [np.min([ d for d in Dist_P[ii][idx_RNA_start:]]) for ii in range(len(prot_res))]

        protein   = BA.Protein(os.path.join(now_path, 'average_pdb_wtc1_eq.pdb'), model = False)
        RNA     =  BA.RNA(os.path.join(now_path, now_pdb), chain_id='B' if real_key == 'wtc1' else '', model = False)
        protein.Get_Atom_Coord(atom_name = 'CA')
        protein.Get_lenght()
        #protein.Get_Protein_Sequence()
        RNA.Get_Atom_Coord(atom_name='P')
        #RNA.Get_RNA_Sequence()
        RNA.Get_lenght()
        #sequence = np.concatenate([protein.sequence, RNA.sequence])
        Coord_P   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
        Dist_P    = BA.Dist_Matrix(Coord_P)
        df_dist[key+'_ini'] = [np.min([ d for d in Dist_P[ii][idx_RNA_start:]]) for ii in range(len(prot_res))]

    else:
        df[key] = [1 if r in BS_P else 0 for r in prot_res]
    
        df_dist[key] = [np.min([ d for d in Dist_P[ii][idx_RNA_start:]]) for ii in range(len(prot_res))]

df = df[['wtc1', 'wtc1_ini', 'wtc1_h_new', 'wtc1_h_new_discarded', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16', 'wtc7_16_discarded', 'MTC1', 'MTC1_discarded', 'mtc2', 'mtc2_discarded', 'mtc3', 'mtc3_discarded',  'mtc4', 'mtc4_discarded']]
df_dist = df_dist[['wtc1', 'wtc1_ini', 'wtc1_h_new', 'wtc1_h_new_discarded', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16', 'wtc7_16_discarded', 'MTC1', 'MTC1_discarded', 'mtc2', 'mtc2_discarded', 'mtc3', 'mtc3_discarded',  'mtc4', 'mtc4_discarded']]

df.to_csv(join(path, 'df_contacts_comparison_{}ang.csv'.format(BS_treshold)))
df_dist.to_csv(join(path, 'df_distance_comparison.csv'))

# %% 4) FIGURA dei Contatti

f, ax = plt.subplots(figsize=(19, 8))
ax.pcolor(df.values, cmap = 'hot', edgecolors = 'white', linewidth = .1)

x_ticks = np.arange(0, df.values.shape[1], step = 1)
y_labels = [str(int(a)) for a in ax.get_yticks() + prot_res[0] ]
x_labels = ['NMR', 'NMR_ini','wtc1', 'wtc1_d', 'wtc2','wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7', 'wtc7_d', 'mtc1', 'mtc1_d', 'mtc2', 'mtc2_d', 'mtc3', 'mtc3_d', 'mtc4', 'mtc4_d' ]
ax.set_xticks(x_ticks + 0.5, minor = False)
#ax.set_yticks(y_ticks)
ax.set_xticklabels(x_labels)
ax.set_yticklabels(y_labels)
ax.set_ylabel('Residue number')
ax.set_xlabel('Structure identifier')
ax.set_title('Contacts comparison of structure with NMR\nBS treshold = {} ang'.format(BS_treshold))
plt.xticks(fontsize='large')
# %% 5) FIGURA delle minime distanze

f, ax = plt.subplots(figsize=(19, 8))
mappa = ax.pcolor(df_dist.values, cmap = 'hot', edgecolors = 'white', linewidth = .1)
plt.colorbar(mappa)

x_ticks = np.arange(0, df_dist.values.shape[1], step = 1)
y_labels = [str(int(a)) for a in ax.get_yticks() + prot_res[0] ]
x_labels = ['NMR', 'NMR_ini','wtc1', 'wtc1_d', 'wtc2','wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7', 'wtc7_d', 'mtc1', 'mtc1_d', 'mtc2', 'mtc2_d', 'mtc3', 'mtc3_d', 'mtc4', 'mtc4_d' ]
ax.set_xticks(x_ticks + 0.5, minor = False)
#ax.set_yticks(y_ticks)
ax.set_xticklabels(x_labels)
ax.set_yticklabels(y_labels)
ax.set_ylabel('Residue number')
ax.set_xlabel('Structure identifier')
ax.set_title('Minimum C_a - P distance (angstrom) comparison of structures with NMR = {} ang'.format(BS_treshold))
plt.xticks(fontsize='large')

# %%

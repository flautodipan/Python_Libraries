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
treshold = 9
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

    ax.set_title('{} Number of contacts vs residue number\n Treshold = {} ang   # tot frame = {} '.format(key, treshold, eq_len[key]))
    ax.set_xlabel('Protein Residue Number')
    ax.set_ylabel('# of contacts')

    plt.legend()
    f.savefig(join(now_path,key+'_contacts'+eq+'.pdf'), format = 'pdf')

# %%

# 4) Matrice dei contatti

now_keys_discarded = ['wtc1_h_new', 'wtc7_16', 'MTC1', 'mtc2', 'mtc3', 'mtc4']
now_keys_discarded = [ key + '_discarded' for key in now_keys_discarded]

res =  dfs['wtc1'].res_number.values
df = pd.DataFrame(index =res )

NOW_KEYS = now_keys+now_keys_discarded
NOW_KEYS.sort()

for key in NOW_KEYS:

    # acquisisco frame iniziale 
    now_path = join(path, key.upper() if 'discarded' not in key else key[:-10].upper())
    now_pdb = key if 'discarded' not in key else key[:-10]
    now_pdb += '.pdb'

    protein   = BA.Protein(os.path.join(now_path, now_pdb), model = False)
    RNA     =  BA.RNA(os.path.join(now_path, now_pdb), chain_id=' ', model = False)

    print('Ok, acquisito frame iniziale {} da {}\n'.format('discarded' if 'disc' in key else ' ', now_path))

    #BS_P, _ = 
    
    
    
    temp = np.zeros(len(res))
    df[key] = [1 if r in exp_contacts else 0 for r in res]

df
# %%
# %%
    f, ax = plt.subplots()
    ax.pcolor(df.values, cmap = 'PiYG', edgecolors = 'white', linewidth = .1)

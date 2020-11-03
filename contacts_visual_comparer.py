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

wtc_keys = ['wtc1_h_new', 'wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5', 'wtc6', 'wtc7_16']
cov_keys = ['cov2', 'cov3', 'cov4']
mtc_keys = ['MTC1', 'mtc2', 'mtc3']
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

for now_name in now_keys:  

    now_path    = join(path, now_name.upper())

    eq = '_eq2' if now_name in now_eqs2 else '_eq1' if now_name in now_eqs1 else '_eq' 
    if ((now_name in now_eqs1) & (now_name in now_eqs2)):
        raise ValueError('{} in both eq1 and eq2'.format(now_name))

    BSs[now_name] = np.load(join(now_path, 'BS_P'+eq+'.npy'))
    dfs[now_name] = pd.read_json(join(now_path, 'df'+eq+'_contacts_'+now_name+'.json'))
    print('Acquisisco {} con eq = {}'.format(now_name, eq, ))

#%%
#2) stampo figure

for now_name in now_keys:  

    now_path    = join(path, now_name.upper())
    now_exp_df  = exp_df[exp_df.identifier == now_name]
    eq = '_eq2' if now_name in now_eqs2 else '_eq1' if now_name in now_eqs1 else '_eq' 
    if ((now_name in now_eqs1) & (now_name in now_eqs2)):
        raise ValueError('{} in both eq1 and eq2'.format(now_name))
    # /1000 perché frame ogni 100 ps, eq_sampling ogni 10
    eq_len      = (eval(getattr(now_exp_df, 'time_range'+eq).values[0])[1] - eval(getattr(now_exp_df, 'time_range'+eq).values[0])[0])/1000
    if now_name == 'wtc1' : exp_eq_len = eq_len
    #acquisizione
    f, ax = plt.subplots()

    ax.stem(dfs[now_name].contacts.values, markerfmt = now_exp_df.color.values[0], basefmt = now_exp_df.color.values[0], linefmt = now_exp_df.color.values[0])

    ax.set_title('{} Number of contacts vs residue number\n Treshold = {} ang   # tot frame = {} '.format(now_name, treshold, int(eq_len)))
    ax.set_xlabel('Protein Residue Number')
    ax.set_ylabel('# of contacts')

    f.savefig(join(now_path,now_name+'_contacts'+eq+'.pdf'), format = 'pdf')
# %%

#3) RESIDUI PIU IMPORTANTI DI WTC1 e confronto con pdb iniziali 

treshold = 300 #valore basato su pdf a partire da distribuzione valori

exp_contacts        = dfs['wtc1'].res_number[dfs['wtc1'].contacts > treshold]
exp_contacts_name   = dfs['wtc1'].res_name[dfs['wtc1'].contacts > treshold]

print('Ho trovato {} residui per la dinamica NMR che hanno contatti in più di {} frame su {}\n\n'.format(len(exp_contacts), treshold, int(exp_eq_len)))

truth = {}
if 'wtc1' in now_keys : now_keys.remove('wtc1')

for key in now_keys:

    truth[key] = [True if bs in exp_contacts.values else False for bs in BSs[key]]

    print('La configurazione iniziale di {} ritrova\n{} contatti sperimentali su {}\n'.format(key, sum(truth[key]), len(exp_contacts)))

# %%

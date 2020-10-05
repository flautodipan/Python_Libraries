#%%

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt

now_path = '../GROMACS/'
skip_cov = False


# %%
# VERSIONE single ---> figure e tutto


# 0 ) INPUT
#
# prendo info sperimentali e di dinamica MD da file excel in 
# ../GROMACS/MD_experimental_data.xlsx

now_name = 'cov1'
exp_df = pd.read_excel(now_path+'MD_experimental_data.xlsx')
exp_df = exp_df[exp_df.identifier == now_name]

now_path        += now_name.upper() +'/'
n_frames        = int(exp_df.n_frames.values[0])
time_range      = eval(exp_df.time_range.values[0])
time_range_eq   = eval(exp_df.time_range_eq.values[0])

#colori
color           = exp_df.color.values[0]
darkcolor           = exp_df.darkcolor.values[0]
brightcolor           = exp_df.brightcolor.values[0]
contrastcolor           = exp_df.contrastcolor.values[0]
darkcontrastcolor           = exp_df.darkcontrastcolor.values[0]
# %%
# 1 ) Genero traiettoria, acquisisco info RMSD e RMSF

traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, '300 K'))
traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )

#RMSD 
traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', fig = now_name+'_RMSD', histo = now_name+'_RMSD_Histogram', bins = 50, path = now_path, color = color, ylim = (0,2))
traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq, path = now_path, fig =  now_name+'_RMSD_eq', alpha = 0.1, color = color, darkcolor = darkcolor, ylim = (0,2) )

#RMSF
#acquisisco indici residui proteina da file .xvg dato da GROMACS
res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+'.xvg', now_path, skip_lines= 17)
res = np.array(res, dtype=int)
# faccio diventare indici dei residui RNA successivi a quelli proteina
idx_RNA_start = np.where( res == 1.)[0][0]
res[idx_RNA_start:] += res[idx_RNA_start-1]


# %%
# 2 ) TROVO CONFIGURAZIONE BINDING SITE
# 
# 2 a - TROVO CONFIGURAZIONE DI MEDIO RMSD

mean_RMSD , mean_frame_idx = Find_Nearest(traj.RMSD_eq, np.mean(traj.RMSD_eq)) 
mean_frame_time = ((mean_frame_idx+traj.idx_eq_left)*traj.timestep) - traj.initial_time
print('Ho trovato il "centroide" della distribuzione RMSD a equilibrio\nnel frame {}\ncorrispondente al tempo {} ps\ncon RMSD = {} nm'.format(mean_frame_idx+traj.idx_eq_left, mean_frame_time, mean_RMSD))

#%%
# 2 b - ACQUISISCO PDB GENERATO da gmx trjconv e trovo BS
#   
filename='BS_{}_make_ndx.txt'.format(now_name)
pdb_filename = 'average_pdb_'+now_name+'.pdb'
treshold = exp_df.bs_treshold.values[0]

protein   = BA.Protein(now_path+pdb_filename, model = False)
RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False)

print('Ok, acquisito correttamente pdb')

protein.Get_CA_Coord()
RNA.Get_P_Coord(atom_name="O5'")

Coord   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
Dist    = BA.Dist_Matrix(Coord)

Cont    = BA.Contacts_Matrix(Dist, treshold)
Bonds   = BA.Analyze_Bond_Residues(Cont, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))

BS      = BA.Print_Protein_BS_old(res, Bonds, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
BS_RNA  = BA.Print_Protein_BS_old(res, Bonds, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']

with open  (now_path+filename, 'w') as f:
        f.write("# frame \t Protein Binding Site (BS) Residues at {} Å treshold\n".format(treshold))
        for bs in BS:
            f.write('r {} | '.format(bs))

# me lo salvo
np.save(now_path+'BS.npy', BS)
np.save(now_path+'BS_RNA.npy', BS_RNA)

#%%
# 3) Completo Acquisizioni e stampe degli altri RMSD e li plotto insieme

#2) Acquisisco gli altri RMSD generati 

traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = False, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA_tot', color = color, scale = 'ns', ylim = (0,2))
traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = False, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS_tot', color = color, scale = 'ns', ylim = (0,2))


traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA', color = color, scale = 'ns')
traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = True, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS', color = color, scale = 'ns')


reduction = 50
sampling = np.arange(0, len(traj.RMSD_eq), reduction, dtype=int)
time_steps = np.arange(traj.time_range_eq[0], traj.time_range_eq[1]+1, traj.timestep*reduction)/1000
plt.figure()
plt.title('Equilibrium RMSD comparison between {} groups'.format(now_name))
plt.plot(time_steps, traj.RMSD_eq[sampling], '--', color = brightcolor, alpha = 0.8, label = 'Prot+RNA')
plt.plot(time_steps, traj.RMSD_eq_RNA[sampling], color = color,  label = 'RNA')
plt.plot(time_steps, traj.RMSD_eq_BS[sampling], '-.',  color = contrastcolor, alpha = 0.8, label = 'Prot BS')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.tight_layout()
plt.savefig(now_path+'RMSD_comparison.pdf', format = 'pdf')

#%%
# 4 ) Gli altri RMSF

traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+'.xvg', 'RNA', path = now_path, skip_lines=17 )
traj.Acquire_Atoms_List('rmsf_BS_'+now_name+'.xvg', 'BS', path = now_path, skip_lines = 17)
traj.Get_RMSF(xvg_filename='rmsf_'+now_name+'.xvg', path = now_path, fig = now_name+'_rmsf', color = color, darkcolor = contrastcolor, thirdcolor = brightcolor)

#RMSF for residues
text_size = 7

idx_BS = np.zeros(len(BS), dtype = int)
for (ii,bs) in zip(range(len(BS)), BS):
    idx_BS[ii] = np.where(res == bs)[0]

BS_split = []
split = [BS[0],]

for ii in range(len(BS)): 

    if ii < len(BS)-1 :
        if BS[ii+1] == BS[ii]+1:
            split.append(BS[ii+1])

        else:
            BS_split.append(split)
            split = [BS[ii+1],]
    else:
        if BS[ii]  - 1 == BS[ii-1]:
            split.append(BS[ii])
            BS_split.append(split)

        else:
            BS_split.append(split)
            split = [BS[ii],]
            BS_split.append(split)


text = []   
size = text_size
start = 0

for ii in range(int(BS.size/size)):
    text.append([str(v) for v in BS[start:start+size]])
    start+=size
if BS.size%size != 0:
    final = [str(v) for v in BS[start:]]
    if len(final) != size:
        while True:
            final.append(' ')
            if len(final) == size:
                break
        text.append(final)

f, ax = plt.subplots(1,1)

ax.stem(res[:idx_RNA_start], RMSF_res[:idx_RNA_start],  markerfmt = 'white', basefmt = 'silver' ,linefmt='silver')
ax.stem(res[idx_RNA_start:], RMSF_res[idx_RNA_start:], markerfmt=darkcolor, basefmt=darkcolor, linefmt=color, label = 'RNA')
for bs in BS_split:
    bs = np.array(bs, dtype = int)
    if bs[0] == BS_split[0][0]:
        ax.stem(bs, RMSF_res[np.where([res == b for b in bs])[1]], markerfmt=darkcontrastcolor, basefmt=darkcontrastcolor, linefmt=contrastcolor, label = 'BS')
    else:
        ax.stem(bs, RMSF_res[np.where([res == b for b in bs])[1]], markerfmt=darkcontrastcolor, basefmt=darkcontrastcolor, linefmt=contrastcolor,)


ax.legend(title = 'RNA starts at res {}'.format(int(res[idx_RNA_start])))
ax.set_title('RMSF for {} residues'.format(now_name), pad = 5)
ax.set_xlabel('Residue number')
ax.set_ylabel('RMSF (nm)')

plt.tight_layout()
f.savefig(now_path+'RMSF_res_'+now_name+'.pdf', format = 'pdf', bbox_inches = 'tight')



# %%
# 5 ) Gyration Radium

traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA.xvg', now_path, fig = now_name+'_gyradium', ylim = (0,2), alpha = 0.2, color = color, darkcolor = darkcolor, skip_lines= 27 )


# %%

# 6 ) COVARIANCE ANALYSIS


# MATRICE COVARIANZA ATOMICA
if not skip_cov:
    N = protein.atoms['atom_number'].size + RNA.atoms['atom_number'].size
    cov_matrix = BA.Get_Covariance_Matrix(N, 'cov_eq_'+now_name, now_path)
    BA.Print_Cov_Matrix_BS(cov_matrix, now_name,'Atoms', BS, res, path = now_path, clim = (-0.005, 0.005))

# prendo gli indici atomici dei CA e P
CA_idx = protein.CA['atom_number'].values - 1
P_idx = RNA.P['atom_number'].values - 1 
idx = np.concatenate((CA_idx, P_idx))

if not skip_cov:
    cov_matrix_CAP = BA.CAP_Cov_Matrix(cov_matrix, idx, now_name+'CAP_cov_matrix.txt', now_path)
else:
    cov_matrix_CAP = np.genfromtxt(now_path+now_name+'CAP_cov_matrix.txt')  

BA.Print_Cov_Matrix_BS(cov_matrix_CAP, now_name, 'CA', BS, res, text = True, path = now_path, clim = (-0.005, 0.005))

# FACCIO ANCHE LA MATRICE di PEARSON
N = idx.size
pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP, N)
BA.Print_Cov_Matrix_BS(pearson_matrix_CAP, now_name, 'CA', BS, res, pearson = True, path = now_path, clim = (-1, 1))

#%%
# 7 ) STAMPO CARATTERISTICHE MEDIE
print("caratteristiche per {}".format(now_name))
print("Tempo totale = {} ns".format(np.array(time_range)/1000))
print("Range equilibrio = {} ns".format(np.array(time_range_eq)/1000))
print("RMSD totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(traj.RMSD_eq), np.std(traj.RMSD_eq)))
print("RMSD RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(traj.RMSD_eq_RNA), np.std(traj.RMSD_eq_RNA)))
print("RMSD BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(traj.RMSD_eq_BS), np.std(traj.RMSD_eq_BS)))
print("RMSF totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(RMSF_res), np.std(RMSF_res)))
print("RMSF RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(RMSF_res[idx_RNA_start:]), np.std(RMSF_res[idx_RNA_start:])))
print("RMSF BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(RMSF_res[idx_BS]), np.std(RMSF_res[idx_BS])))
print("Gyration radium medio a eq per BS-RNA = {:3.2f}\nstdev = {:3.2f}".format(np.mean(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]), np.std(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right])))
print('La BS della porteina consta di {} residui'.format(len(BS)))

with open(now_path+now_name+'_finals.txt', 'w') as f:

    f.write(("caratteristiche per {}\n\n".format(now_name)))
    f.write("Tempo totale = {} ns\n0.80(0.02)".format(np.array(time_range)/1000))
    f.write("Range equilibrio = {} ns\n".format(np.array(time_range_eq)/1000))
    f.write("RMSD totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(traj.RMSD_eq), np.std(traj.RMSD_eq)))
    f.write("RMSD RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(traj.RMSD_eq_RNA), np.std(traj.RMSD_eq_RNA)))
    f.write("RMSD BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(traj.RMSD_eq_BS), np.std(traj.RMSD_eq_BS)))
    f.write("RMSF totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(RMSF_res), np.std(RMSF_res)))
    f.write("RMSF RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(RMSF_res[idx_RNA_start:]), np.std(RMSF_res[idx_RNA_start:])))
    f.write("RMSF BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(RMSF_res[idx_BS]), np.std(RMSF_res[idx_BS])))
    f.write("Gyration radium medio a eq per BS-RNA = {:3.2f}\nstdev = {:3.2f}\n".format(np.mean(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right]), np.std(traj.Gyradium[traj.idx_eq_left:traj.idx_eq_right])))
    f.write('La BS della porteina consta di {} residui'.format(len(BS)))

# %%

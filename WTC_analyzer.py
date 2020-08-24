"""
PROGRAMMA ANALISI delle traiettorie dei complessi
--> RMSD, RMSF, Normal mode, etc...
"""
#%%
#0) preambolo

import BioAlessandria as BA
from Alessandria import Find_Nearest
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")



now_temp = '300 K'
scale='ns'
#WTC7_24

now_path    =   '../GROMACS/WTC7_24/'
now_name    =    'wtc7_24'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [50000, 650000]
color = 'black'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastcolor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)

"""
#WTC1_h
now_path    =   '../GROMACS/WTC1_h/'
now_name    =    'wtc1_h'
n_frames = 20001
time_range  =   [0, 2000000]
time_range_eq = [1000000, 2000000]
color = 'cornflowerblue'
darkcolor = 'midnightblue'
brightcolor='magenta'
contrastcolor = 'orange'
darkcontrastcolor = 'orangered'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)

#WTC1
now_path    =   '../GROMACS/WTC1/'
now_name    =    'wtc1'
n_frames = 9700
time_range  =   [30060, 999960]
time_range_eq = [400060, 999960]
color = 'royalblue'
darkcolor = 'navy'
brightcolor = 'magenta'
contrastcolor='orange'
darkcontrastcolor='orangered'
ylim = (0,1)
gyrad_ylim = (1.4, 2.2)


#WTC2
now_path    =   '../GROMACS/WTC2/'
now_name    =    'wtc2'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [400000, 1000000]
color       =   'forestgreen'
darkcolor   =   'darkgreen'
contrastcolor = 'goldenrod'
contrastcolor = 'crimson'
darkcontrastcolor = 'darkred'
ylim = (0,1)
gyrad_ylim = (1.3, 1.9)


#WTC3
now_path    =   '../GROMACS/WTC3/'
now_name    =    'wtc3'
n_frames = 20001
time_range = [0, 2000000]
time_range_eq = [600000, 2000000]
color       =   'goldenrod'
darkcolor   =   'darkgoldenrod'
brightcolor = 'mediumslateblue'
contrastcolor = 'indianred'
darkcontrastcolor = 'darkred'
ylim = (0,1.4)
gyrad_ylim = (1.65, 2.3)


#WTC4
now_path    =   '../GROMACS/WTC4/'
now_name    =    'wtc4'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [400000, 1000000]
color = 'orange'
darkcolor = 'darkorange'
brightcolor = 'yellowgreen'
contrastcolor = 'darkorchid'
darkcontrastcolor = 'darkslateblue'
ylim = (0,1.2)
gyrad_ylim = (0.8, 1.5))

#WTC5
now_path    =   '../GROMACS/WTC5/'
now_name    =    'wtc5'
n_frames = 30001
time_range = [0, 3000000]
time_range_eq = [2100000, 3000000]
color = 'darkorchid'
darkcolor = 'indigo'
brightcolor = 'crimson'
contrastcolor = 'chartreuse'
darkcontrastcolor = 'yellowgreen'
ylim = (0,2)
gyrad_ylim = (1.1, 1.7)


#WTC6
now_path    =   '../GROMACS/WTC6/'
now_name    =    'wtc6'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [400000, 1000000]
color = 'firebrick'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastclor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)


#WTC7
now_path    =   '../GROMACS/WTC7/'
now_name    =    'wtc7'
n_frames = 10001
#ps
time_range = [0, 1000000]
time_range_eq = [200000, 950000]
color = 'firebrick'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastcolor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)

#WTC7_red

now_path    =   '../GROMACS/WTC7_red/'
now_name    =    'wtc7_red'
n_frames = 1251
time_range = [175000, 300000]
time_range_eq = time_range
color = 'black'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastcolor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)

#WTC7_7

now_path    =   '../GROMACS/WTC7_7/'
now_name    =    'wtc7_7'
n_frames = 2301
time_range = [0, 230000]
time_range_eq = time_range
color = 'black'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastcolor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)

#WTC7_8

now_path    =   '../GROMACS/WTC7_8/'
now_name    =    'wtc7_8'
n_frames = 3001
time_range = [0, 300000]
time_range_eq = time_range
color = 'black'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastcolor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)


#WTC7_16

now_path    =   '../GROMACS/WTC7_16/'
now_name    =    'wtc7_16'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [250000, 1000000]
color = 'black'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastcolor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)

#WTC7_24

now_path    =   '../GROMACS/WTC7_24/'
now_name    =    'wtc7_24'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [50000, 650000]
color = 'black'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'
darkcontrastcolor = 'darkgoldenrod'
ylim = (0,1)
gyrad_ylim = (1.1, 1.8)

"""





#%%
#1) acquisizione dati spazio essenziale e RMSD



WTC_traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, now_temp))
WTC_traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )
WTC_traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', fig = now_name+'_RMSD', histo = now_name+'_RMSD_Histogram', bins = 50, path = now_path, color = color, scale = 'ns', ylim = ylim)
WTC_traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq, path = now_path, fig =  now_name+'_RMSD_eq', alpha = 0.1, color = color, darkcolor = darkcolor, scale = 'ns', ylim = ylim)


#%%

#1) Stampa per gmx make_ndx i residui del binding site
#   data una treshold e un frame


#1a) indentifico frame medio dell'RMSD equilibrio

mean_RMSD , mean_frame_idx = Find_Nearest(WTC_traj.RMSD_eq, np.mean(WTC_traj.RMSD_eq)) 
mean_frame_time = ((mean_frame_idx+WTC_traj.idx_eq_left)*WTC_traj.timestep) - WTC_traj.initial_time
print('Ho trovato il "centroide" della distribuzione RMSD a equilibrio\nnel frame {}\ncorrispondente al tempo {} ps\ncon RMSD = {} nm'.format(mean_frame_idx+WTC_traj.idx_eq_left, mean_frame_time, mean_RMSD))


# %%
#1b) prendo la BS di quel frame e te la stampo
filename='BS_{}_make_ndx.txt'.format(now_name)
pdb_filename = 'average_pdb_'+now_name+'.pdb'
treshold = 9

TDP43   = BA.Protein(now_path+pdb_filename, model = False)
RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B', model = False)

print('Ok, acquisito correttamente pdb')

#%%
#stampa

TDP43.Get_CA_Coord()
RNA.Get_P_Coord()

Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
Dist    = BA.Dist_Matrix(Coord)

Cont    = BA.Contacts_Matrix(Dist, treshold)
Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))

BS      = BA.Print_Protein_BS_old(Bonds, TDP43.lenght, prot_initial=TDP43.initial, RNA_initial=RNA.initial)['Prot']
BS_RNA  = BA.Print_Protein_BS_old(Bonds, TDP43.lenght, prot_initial=TDP43.initial, RNA_initial=RNA.initial)['RNA']

with open  (now_path+filename, 'w') as f:
        f.write("# frame \t Protein Binding Site (BS) Residues\n")
        for bs in BS:
            f.write('r {} | '.format(bs))
      
#%%
#2) Acquisisco gli altri RMSD generati 

WTC_traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = False, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA_tot', color = color, scale = scale, ylim = ylim)
WTC_traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = False, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS_tot', color = color, scale = scale, ylim = ylim)

print('ok')
#%%

WTC_traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA', color = color, scale = scale)
WTC_traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = True, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS', color = color, scale = scale)



# %%
#3) Plotto insieme, versione campionata
reduction = 50
sampling = np.arange(0, len(WTC_traj.RMSD_eq), reduction, dtype=int)
time_steps = np.arange(WTC_traj.time_range_eq[0], WTC_traj.time_range_eq[1]+1, WTC_traj.timestep*reduction)/1000
plt.figure()
plt.title('Equilibrium RMSD comparison between {} groups'.format(now_name))
plt.plot(time_steps, WTC_traj.RMSD_eq[sampling], '--', color = brightcolor, alpha = 0.8, label = 'Prot+RNA')
plt.plot(time_steps, WTC_traj.RMSD_eq_RNA[sampling], color = color,  label = 'RNA')
plt.plot(time_steps, WTC_traj.RMSD_eq_BS[sampling], '-.',  color = contrastcolor, alpha = 0.8, label = 'Prot BS')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.tight_layout()
plt.savefig(now_path+'RMSD_comparison.pdf', format = 'pdf')

# %%
# RMSF

WTC_traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+'.xvg', 'RNA', path = now_path, skip_lines=17 )
WTC_traj.Acquire_Atoms_List('rmsf_BS_'+now_name+'.xvg', 'BS', path = now_path, skip_lines = 17)
WTC_traj.Get_RMSF(xvg_filename='rmsf_'+now_name+'.xvg', path = now_path, fig = now_name+'_rmsf', color = color, darkcolor = contrastcolor, thirdcolor = brightcolor)

#%%
#RMSF for residues
text_size = 7

res, RMSF_res = BA.Parse_xvg_skip('rmsf_res_'+now_name+'.xvg', now_path)
res = np.array(res, dtype=int)
#scambio RNA e PROT per coerenza con altre figure 
idx_RNA_start = np.where( res == 1.)[0][0]
res[idx_RNA_start:] += res[idx_RNA_start-1]
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
"""
ax.table(text)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
"""
ax.set_xlabel('Residue number')
ax.set_ylabel('RMSF (nm)')
plt.tight_layout()
f.savefig(now_path+'RMSF_res_'+now_name+'.pdf', format = 'pdf', bbox_inches = 'tight')

#%%
#GYRADIUM
WTC_traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA.xvg', now_path, fig = now_name+'_gyradium', ylim = gyrad_ylim, alpha = 0.2, color = color, darkcolor = darkcolor, skip_lines= 27 )


# %%
# COVARIANCE ANALYSIS

#versione senza figure
"""
WTC_traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, now_temp))
WTC_traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )
WTC_traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', path = now_path)
WTC_traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq)
"""



#%%
# MATRICE COVARIANZA ATOMICA
skip_cov = True
if not skip_cov:
    N = TDP43.atoms['atom_number'].size + RNA.atoms['atom_number'].size
    cov_matrix = BA.Get_Covariance_Matrix(N, 'cov_eq_'+now_name, now_path)
    BA.Print_Cov_Matrix_BS(cov_matrix, now_name,'Atoms', BS, res, path = now_path, clim = (-0.005, 0.005))
# MATRICE COVARIANZA CAP
#%%
# prendo gli indici atomici dei CA e P
CA_idx = TDP43.CA['atom_number'].values -1
P_idx = RNA.P['atom_number'].values -1 
idx = np.concatenate((CA_idx, P_idx))
if not skip_cov:

    cov_matrix_CAP = BA.CAP_Cov_Matrix(cov_matrix, idx, now_name+'CAP_cov_matrix.txt', now_path)
else:
    cov_matrix_CAP = np.genfromtxt(now_path+now_name+'CAP_cov_matrix.txt')  


BA.Print_Cov_Matrix_BS(cov_matrix_CAP, now_name, 'CA', BS, res, text = True, path = now_path, clim = (-0.1, 0.1))

#%%
#trasformo in matrice di pearson


N = idx.size
pearson_matrix_CAP = BA.Pearson_Matrix_from_Cov(cov_matrix_CAP, N)

BA.Print_Cov_Matrix_BS(pearson_matrix_CAP, now_name, 'CA', BS, res, pearson = True, path = now_path, clim = (-1, 1))


# %%
print("caratteristiche per {}".format(now_name))
print("Tempo totale = {} ns".format(np.array(time_range)/1000))
print("Range equilibrio = {} ns".format(np.array(time_range_eq)/1000))
print("RMSD totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(WTC_traj.RMSD_eq), np.std(WTC_traj.RMSD_eq)))
print("RMSD RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(WTC_traj.RMSD_eq_RNA), np.std(WTC_traj.RMSD_eq_RNA)))
print("RMSD BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(WTC_traj.RMSD_eq_BS), np.std(WTC_traj.RMSD_eq_BS)))
print("RMSF totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(RMSF_res), np.std(RMSF_res)))
print("RMSF RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(RMSF_res[idx_RNA_start:]), np.std(RMSF_res[idx_RNA_start:])))
print("RMSF BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}".format(np.mean(RMSF_res[idx_BS]), np.std(RMSF_res[idx_BS])))
print("Gyration radium medio a eq per BS-RNA = {:3.2f}\nstdev = {:3.2f}".format(np.mean(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right]), np.std(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right])))

with open(now_path+now_name+'_finals.txt', 'w') as f:

    f.write(("caratteristiche per {}\n\n".format(now_name)))
    f.write("Tempo totale = {} ns\n0.80(0.02)".format(np.array(time_range)/1000))
    f.write("Range equilibrio = {} ns\n".format(np.array(time_range_eq)/1000))
    f.write("RMSD totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(WTC_traj.RMSD_eq), np.std(WTC_traj.RMSD_eq)))
    f.write("RMSD RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(WTC_traj.RMSD_eq_RNA), np.std(WTC_traj.RMSD_eq_RNA)))
    f.write("RMSD BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(WTC_traj.RMSD_eq_BS), np.std(WTC_traj.RMSD_eq_BS)))
    f.write("RMSF totale medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(RMSF_res), np.std(RMSF_res)))
    f.write("RMSF RNA medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(RMSF_res[idx_RNA_start:]), np.std(RMSF_res[idx_RNA_start:])))
    f.write("RMSF BS medio all'equilibrio = {:3.2f}\n Con stdev = {:3.2f}\n".format(np.mean(RMSF_res[idx_BS]), np.std(RMSF_res[idx_BS])))
    f.write("Gyration radium medio a eq per BS-RNA = {:3.2f}\nstdev = {:3.2f}\n".format(np.mean(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right]), np.std(WTC_traj.Gyradium[WTC_traj.idx_eq_left:WTC_traj.idx_eq_right])))
200    # %%


# %%

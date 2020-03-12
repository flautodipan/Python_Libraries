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
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


#WTC2
now_path    =   '../GROMACS/WTC2/'
now_name    =    'wtc2'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [400000, 1000000]
color       =   'darkolivegreen'
darkcolor   =   'crimson'
thirdcolor = 'mediumslateblue'
contrastcolor = 'darkred'
ylim = (0,1)
gyrad_ylim = (1.3, 1.9)


now_temp = '350 K'
scale='ns'


"""
#WTC1_h
now_path    =   '../GROMACS/WTC1_h/'
now_name    =    'wtc1_h'
n_frames = 10001
time_range  =   [0, 1000000]
time_range_eq = [400000, 1000000]
color = 'seagreen'
darkcolor = 'forestgreen'
thirdcolor = 'orange'
contrastcolor = 'darkred'
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
contrastcolor='orange'
ylim = (0,1)
gyrad_ylim = (1.4, 2.2)

#WTC2
now_path    =   '../GROMACS/WTC2/'
now_name    =    'wtc2'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [400000, 1000000]
color       =   'darkolivegreen'
darkcolor   =   'crimson'
thirdcolor = 'mediumslateblue'
contrastcolor = 'darkred'
ylim = (0,1)
gyrad_ylim = (1.3, 1.9)


#WTC3
now_path    =   '../GROMACS/WTC3/'
now_name    =    'wtc3'
n_frames = 20001
time_range = [0, 2000000]
time_range_eq = [600000, 2000000]
color       =   'darkgoldenrod'
darkcolor   =   'chocolate'
thirdcolor = 'mediumslateblue'
contrastcolor = 'darkred'
ylim = (0,1.4)
gyrad_ylim = (1.65, 2.3)


#WTC4
now_path    =   '../GROMACS/WTC4/'
now_name    =    'wtc4'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [400000, 1000000]
color = 'darkorange'
darkcolor = 'slategray'
thirdcolor = 'yellowgreen'
contrastcolor = 'darkred'
ylim = (0,1.2)
gyrad_ylim = (0.8, 1.5)

#WTC5
now_path    =   '../GROMACS/WTC5/'
now_name    =    'wtc5'
n_frames = 20001
time_range = [0, 2000000]
time_range_eq = [1400000, 2000000]
color = 'darkmagenta'
darkcolor = 'chartreuse'
thirdcolor = 'palevioletred'
contrastcolor = 'lawngreen'
ylim = (0,2)
gyrad_ylim = (1.1, 1.7)

#WTC6

now_path    =   '../GROMACS/WTC6/'
now_name    =    'wtc6'
n_frames = 10001
time_range = [0, 1000000]
time_range_eq = [400000, 1000000]
color = 'firebrick'
darkcolor = 'gold'
ylim = (0,1)


"""

#%%
#1) acquisizione dati spazio essenziale e RMSD

WTC_traj = BA.Trajectory(bio_name='{} in water (T = {})'.format(now_name, now_temp))
WTC_traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100 )
WTC_traj.Get_RMSD(xvg_filename = 'rmsd_'+now_name+'.xvg', fig = now_name+'_RMSD', histo = now_name+'_RMSD_Histogram', bins = 50, path = now_path, color = color, scale = 'ns', ylim = ylim)
#%%
WTC_traj.Define_Equilibrium_by_RMSD(time_range_eq = time_range_eq, path = now_path, fig =  now_name+'_RMSD_eq', alpha = 0.1, color = color, darkcolor = darkcolor, scale = 'ns', ylim = ylim)


#%%

#1) Stampa per gmx make_ndx i residui del binding site
#   data una treshold e un frame


#1a) indentifico frame medio dell'RMSD equilibrio

mean_RMSD , mean_frame_idx = Find_Nearest(WTC_traj.RMSD_eq, np.mean(WTC_traj.RMSD_eq)) 
mean_frame_time = ((mean_frame_idx+WTC_traj.idx_eq)*WTC_traj.timestep) - WTC_traj.initial_time
print('Ho trovato il "centroide" della distribuzione RMSD a equilibrio\nnel frame {}\ncorrispondente al tempo {} ps\ncon RMSD = {} nm'.format(mean_frame_idx+WTC_traj.idx_eq, mean_frame_time, mean_RMSD))


# %%
#1b) prendo la BS di quel frame e te la stampo
filename='BS_{}_make_ndx.txt'.format(now_name)
pdb_filename = 'average_pdb.pdb'
treshold = 10

TDP43   = BA.Protein(now_path+pdb_filename)
RNA     =  BA.RNA(now_path+pdb_filename, chain_id='B')

print('Ok, acquisito correttamente pdb')

#%%
#stampa

TDP43.Get_CA_Coord()
RNA.Get_P_Coord()

Coord   = np.concatenate((TDP43.CA_Coord, RNA.P_Coord))
Dist    = BA.Dist_Matrix(Coord)

Cont    = BA.Contacts_Matrix(Dist, treshold)
Bonds   = BA.Analyze_Bond_Residues(Cont, (TDP43.lenght, RNA.lenght), ("TDP43", "RNA"), first=  ('RNA', 1), second = ('Proteina', TDP43.initial))

BS      = BA.Print_Protein_BS_old(Bonds, TDP43.lenght, initial=TDP43.initial)

with open  (now_path+filename, 'w') as f:
        f.write("# frame \t Protein Binding Site (BS) Residues\n")
        for bs in BS:
            f.write('r {} | '.format(bs))
      
#%%
#2) Acquisisco gli altri RMSD generati 

WTC_traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = False, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA_tot', color = color, scale = scale, ylim = ylim)
WTC_traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = False, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS_tot', color = color, scale = scale, ylim = ylim)



WTC_traj.Get_RMSD('rmsd_'+now_name+'_RNA.xvg', equilibrium = True, mode = 'RNA', path = now_path, fig = now_name+'_RMSD_RNA', color = color, scale = scale)
WTC_traj.Get_RMSD('rmsd_'+now_name+'_BS.xvg',  equilibrium = True, mode = 'BS', path = now_path, fig = now_name+'_RMSD_BS', color = color, scale = scale)



# %%
#3) Plotto insieme, versione campionata
reduction = 50
sampling = np.arange(0, len(WTC_traj.RMSD_eq), reduction, dtype=int)
time_steps = np.arange(WTC_traj.time_range_eq[0], WTC_traj.time_range_eq[1]+1, WTC_traj.timestep*reduction)
plt.figure()
plt.title('Equilibrium RMSD comparison between {} groups'.format(now_name))
plt.plot(time_steps, WTC_traj.RMSD_eq[sampling], '--', color = color, alpha = 0.8, label = 'Prot+RNA')
plt.plot(time_steps, WTC_traj.RMSD_eq_RNA[sampling], color = darkcolor, alpha = 0.8, label = 'RNA')
plt.plot(time_steps, WTC_traj.RMSD_eq_BS[sampling], '-.',  color = thirdcolor, alpha = 0.8, label = 'Prot BS')
plt.legend()
plt.savefig(now_path+'RMSD_comparison.pdf', format = 'pdf')

# %%
# Acquisizione Raggio di Girazione
WTC_traj.Acquire_Atoms_List('rmsf_RNA_'+now_name+'.xvg', 'RNA', path = now_path, skip_lines=17 )
WTC_traj.Acquire_Atoms_List('rmsf_BS_'+now_name+'.xvg', 'BS', path = now_path, skip_lines = 17)
WTC_traj.Get_RMSF(xvg_filename='rmsf_'+now_name+'.xvg', path = now_path, fig = now_name+'_rmsf', color = color, darkcolor = contrastcolor, thirdcolor = thirdcolor)
#%%
WTC_traj.Get_Gyradium('gyration_'+now_name+'_BS_RNA.xvg', now_path, fig = now_name+'_gyradium', ylim = gyrad_ylim, alpha = 0.2, color = color, darkcolor = darkcolor, skip_lines= 27 )

# %%

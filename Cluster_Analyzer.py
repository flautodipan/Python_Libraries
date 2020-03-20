#%%
import BioAlessandria as BA
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot



#RNA6
now_name    =    'RNA6'
now_temp    =   '500K'
color = 'firebrick'
darkcolor = 'darkred'
brightcolor = 'limegreen'
contrastcolor='gold'

n_frames = 10001
time_range = [0, 1000000]
now_path    =   '../GROMACS/'+now_name+'/analysis/'

"""

#RNA1
now_name    =    'RNA1'
now_temp    =   '400 K'
color = 'royalblue'
darkcolor = 'navy'
brightcolor = 'magenta'
contrastcolor='orange'

#RNA2
now_name    =    'RNA2'
now_temp    =   '350 K'
color       =   'lightseagreen'
darkcolor   =   'darkslategray'
brightcolor = 'mediumslateblue'
contrastcolor = 'crimson'

#RNA3
now_name    =    'RNA3'
now_temp    =   '350 K'
color       =   'goldenrod'
darkcolor   =   'darkgoldenrod'
brightcolor = 'limegreen'
contrastcolor = 'indianred'
darkcontrastcolor = 'darkred'

#RNA4
now_name    =    'RNA4'
now_temp    =   '400K'
color = 'orange'
darkcolor = 'darkorange'
brightcolor = 'yellowgreen'
contrastcolor = 'darkorchid'
darkcontrastcolor = 'darkslateblue'

#RNA5
now_name    =    'RNA5'
now_temp    =   '350K'
color = 'darkorchid'
darkcolor = 'indigo'
brightcolor = 'crimson'
contrastcolor = 'chartreuse'



"""
#%%


#1) acquisisco e analizzo autovalori 

RNA_traj = BA.Trajectory(bio_name = now_name+' in water (T = '+now_temp+')')
RNA_traj.Set_Time_Info(n_frames = n_frames, time_range = time_range, timestep = 100)
RNA_traj.Get_EigenValues(path=now_path)
RNA_traj.Analyze_Eig_Variance_Ratio(n_eig = 2)
RNA_traj.Print_EigenValues(path=now_path,  fig = 'Eigen_Values', color = darkcolor, contrastcolor = contrastcolor)

#%%
#2) acquisisco proiezione 2D sui due autovettori principali e RMSD della traiettoria

RNA_traj.Get_RMSD(xvg_filename = 'rmsd_md_'+now_name+'.xvg', fig = now_name+'_RMSD', scale = 'ns', histo = now_name+'_RMSD_Histogram', bins = 50, path = now_path, color = color)
RNA_traj.Get_2D_Traj(path=now_path, fig = 'Essential')
RNA_traj.Get_Terminals_Dist(xvg_filename = 'ter_dist.xvg', skip_lines = 17, fig = now_name+'_Ter_dist', histo = now_name+'_ter_dist_Histogram', bins = 50, path = now_path, color = color)


#%%
# 3A) faccio il fit della distribuzione con N  gaussiane 
# ordine params (mu, sigma, A)

p0 = np.genfromtxt(now_path+now_name+'_Nmodalfit.txt', skip_header=1)[:, 1]
N_gauss = 2 if len(p0)==6 else 3

img_kwargs  =   {'path': now_path, 'fig_fit' : 'Gauss_fit_ter_dist', 'color_fold' : darkcolor, 'fit_color' : brightcolor, 'color_unfold': contrastcolor}
RNA_traj.Analyze_Traj_by_Descriptor(descriptor = 'terminals distance', N_gauss=N_gauss,  p0=p0, bins=50, **img_kwargs )

with open (now_path+now_name+'_Nmodalfit.txt', 'w') as f:
    f.write('Fit_Params = %s'%(str(RNA_traj.df)))


#%%
#4A) impongo il limite mu - sigma*factor all'rmsd e in base a quello divido in due lo spazio essenziale

img_kwargs  =   {'path': now_path, 'fig' : 'ter_dist_division', 'color': darkcolor, 'darkcolor': contrastcolor}
RNA_traj.Divide_Traj_by_Descriptor(descriptor = 'ter_dist',  sigma_factor = 0., **img_kwargs)

#%%
#5A) faccio analisi dell'UNFOLD
RNA_Unfold     =   BA.Cluster_2DAnalysis('KMeans', now_name+' Unfold States')
xy              =   np.array([RNA_traj.x_unfold, RNA_traj.y_unfold]).T
RNA_Unfold.Get_by_Passing(xy,  timestep = RNA_traj.timestep)
RNA_Unfold.Silhouette_KMeans(kmax= 15, path = now_path, fig = 'KMeans_Silhouette_Unfold' )
RNA_Unfold.Clusterize(verbose = True, descriptor=RNA_traj.ter_dist, fig = 'RNA_unfold_clusters', path = now_path)
RNA_traj.Which_Part_of_Traj(RNA_Unfold.clusterchiefs)

#%%
#6)
#Stampe finali
print('Caratteristiche di {}'.format(now_name))
print("Tempo tot = {}".format(RNA_traj.time_range))
print("Temperature = {}".format(now_temp))
print("Mean RMSD = {:3.2f}\n StdDev = {:3.2f}".format(np.mean(RNA_traj.RMSD), np.std(RNA_traj.RMSD)))
print("Mean terminals distance = {:3.2f}\n StdDev = {:3.2f}".format(np.mean(RNA_traj.ter_dist), np.std(RNA_traj.ter_dist)))
print('Unfold population = {}'.format(RNA_traj.unfold_percentage))
print('Explained variance = {}'.format(RNA_traj.perceigen))


#%%


# %%

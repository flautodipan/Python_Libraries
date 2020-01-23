#%%
import BioAlessandria as BA
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot



now_path    =   'RNA4/400K/analysis/'
now_name    =    'RNA4'
now_temp    =   '400K'
#%%


#1) acquisisco e analizzo autovalori 

RNA_traj = BA.Trajectory(bio_name = now_name+' in water (T = '+now_temp+')')
RNA_traj.Get_EigenValues(path=now_path,  fig = 'Eigen_Values')
RNA_traj.Analyze_Eig_Variance_Ratio(n_eig = 2)

#%%
#2) acquisisco proiezione 2D sui due autovettori principali e RMSD della traiettoria

RNA_traj.Get_2D_Traj(time_range = [0, 1000000], path=now_path, fig = 'Essential', verbose = True)
RNA_traj.Get_RMSD(xvg_filename = 'rmsd.xvg', fig = 'RMSD', skip_lines = 18, histo = 'RMSD_Histogram',  bins = 100, path = now_path)
RNA_traj.Get_Terminals_Dist(xvg_filename = 'ter_dist.xvg', skip_lines = 17, fig = 'Ter_dist', histo = 'ter_dist_Histogram',  bins = 100, path = now_path)


#%%
# 3A) faccio il fit della distribuzione con N  gaussiane 

p0 = (1.37, .4, 400, 3.7, .7, 70)
#p0 = (1,1,1,1,1,1, 1, 1, 1)
img_kwargs  =   {'path': now_path, 'fig_fit' : 'Gauss_fit_ter_dist'}
RNA_traj.Analyze_Traj_by_Descriptor(descriptor = 'ter_dist', N_gauss=2,  p0=p0, bins=100, **img_kwargs )

#%%
#4A) impongo il limite mu - sigma*factor all'rmsd e in base a quello divido in due lo spazio essenziale

img_kwargs  =   {'path': now_path, 'fig' : 'ter_dist_division'}
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
#3B) clusterizzo su tutta la traiettoria

RNA_Clust = BA.Cluster_2DAnalysis('KMeans', now_name+' in water')
RNA_Clust.Get_by_Traj(RNA_traj)
RNA_Clust.Silhouette_KMeans(kmax=15, path=now_path,  fig = 'KMeans_Silhouette_All' )
#RNA_Clust.Clusterize(by_silhouette = False, n_clust = 3, verbose=True, fig = 'MD_clusters', path = now_path)
RNA_Clust.Clusterize(verbose=True, fig = 'MD_clusters', path = now_path)


#%%
#4B) cluster RMSD division : acquisisco solo cluster il cui centroide ha RMSD più alto 

RNA_Clust.Cluster_RMSD_Division(RNA_traj.RMSD, bins = 100, path=now_path)#, verbose = True)


#%%
#5B) faccio analisi del cluster UNFOLD

RNA_Unfold     =   BA.Cluster_2DAnalysis('KMeans', 'RNA2 Unfold States')
RNA_Unfold.Get_by_Passing(RNA_Clust.xy[RNA_Clust.Cluster.labels_ == RNA_Clust.unfoldest_idx],  timestep = RNA_Clust.timestep)
RNA_Unfold.Silhouette_KMeans(kmax= 15, path = now_path, fig = 'KMeans_Silhouette_Unfold' )
RNA_Unfold.Clusterize(verbose = True, fig = 'RNA2Unfold_clusters', path = now_path)



# %%

a =  (1499, 246175, 305569)
b = (105689, 157984, 244375)


# %%
#provo metodo A con altro descrittore
RNA_traj = BA.Trajectory(bio_name = now_name+' in water (T = 350K)')
RNA_traj.Get_EigenValues(path=now_path,  fig = 'Eigen_Values')
RNA_traj.Analyze_Eig_Variance_Ratio(n_eig = 2)

#%%
#2) acquisisco proiezione 2D sui due autovettori principali e RMSD della traiettoria

RNA_traj.Get_2D_Traj(time_range = [0, 1000000], path=now_path, fig = 'Essential', verbose = True)
RNA_traj.Get_RMSD(xvg_filename = 'rmsd.xvg', fig = 'RMSD', skip_lines = 18, histo = 'RMSD_Histogram',  bins = 100, path = now_path)
RNA_traj.Get_Terminals_Dist(xvg_filename = 'ter_dist.xvg', skip_lines = 17, fig = 'Ter_dist', histo = 'ter_dist_Histogram',  bins = 100, path = now_path)

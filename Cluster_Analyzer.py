#%%
import BioAlessandria as BA
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot

now_path    =   'RNA2/350K/analysis/'

#%%
#1) acquisisco e analizzo autovalori 
RNA2_traj = BA.Trajectory(bio_name = 'RNA2 in water (T = 350K)')
RNA2_traj.Get_EigenValues(path=now_path,  fig = 'Eigen_Values')
RNA2_traj.Analyze_Eig_Variance_Ratio(n_eig = 2)

#%%
#2) acquisisco proiezione 2D sui due autovettori principali e RMSD

RNA2_traj.Get_2D_Traj(time_range = [0, 1000000], path=now_path, fig = 'Essential', verbose = True)
RNA2_traj.Get_RMSD(xvg_filename = 'rmsd.xvg',path = now_path, fig = 'RMSD',  skip = 10, skip_lines = 18)

#%%
#3) clusterizzo
RNA2_Clust = BA.Cluster_2DAnalysis(RNA2_traj,'KMeans', 'RNA2 in water')
RNA2_Clust.Silhouette_KMeans(kmax=15, path=now_path,  fig = 'KMeans_Silhouette' )
RNA2_Clust.Clusterize(verbose=True, fig = 'MD_clusters', path = now_path)






#%%
import BioAlessandria as BA
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot

#%%
#RNA 2
EigVal_RNA2 = BA.Eigen_Trj(path=path, fig = 'Eigen_Values')
#EigVal_RNA2.Analyze_Variance_Explained_Ratio(0.65)


# %%
#creare un old path
x,y =  BA.Parse_xvg('2dproj.xvg', path= '1micros/')

#per un motivo stupiddo ci sta un frame di troppo
x   =   x[np.arange(x.size-1)]
y   =   y[np.arange(y.size-1)]

print('Taglia array pre reduction: %d' %(x.size))
#faccio sampling ogni 100 frame
indices = np.arange(0, x.size, 100)
x   =   x[indices]
y   =   y[indices]

print('Taglia array post reduction: %d' %(x.size))


fig = plt.figure()
c = range(x.size)
plt.scatter(x,y, c=c, cmap= 'viridis', s = 10)
plt.colorbar()
plt.xlabel(' PC 1 ')
plt.ylabel(' PC 2 ')
plt.title('Projection of MD trajectory in Essential Space')
plt.show()
fig.savefig('Essential.png')

#%%
xy = np.array([x,y]).T
RNA2_Clust = BA.Cluster_2DAnalysis(xy, 1000000)
RNA2_Clust.Silhouette_KMeans(15)

# %%

RNA2_Clust.Clusterize(verbose=True)
RNA2_Clust.Figure_Clusters()

#%%
#parte2) prendo il cluster con RMSD più alto e lo ri-clusterizzo
#Je do de sotto cluster
#1) divisione by RMSD e identificazione cluster con RMSD più alto

RNA2_Clust.Cluster_RMSD_Division('rmsd_red.xvg', path = '1micros/', skip = 100, fig_RMSD='RMSD')


# %%


RNA2_Unfold  =   BA.Cluster_2DAnalysis(RNA2_Clust.xy[RNA2_Clust.Cluster.labels_ == RNA2_Clust.unfoldest_idx],  final_time = 1000000)
RNA2_Unfold.Silhouette_KMeans(15)

# %%

RNA2_Unfold.Clusterize(verbose=True)
RNA2_Unfold.Figure_Clusters()

# %%

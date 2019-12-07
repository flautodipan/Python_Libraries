#%%
import time
start = time.process_time()
from operator import methodcaller
import BioAlessandria as BA
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot
from sklearn.cluster import KMeans


#%%
#RNA 2
EigVal_RNA2 = BA.Eigen_Trj(path = '3 micros/', fig = True)
EigVal_RNA2.Analyze_Variance_Explained_Ratio(0.65)

# %%

x,y =  BA.Parse_xvg('1micros/2dproj.xvg')
#per un motivo stupiddo ci sta un frame di troppo
x   =   x[np.arange(x.size-1)]
y   =   y[np.arange(y.size-1)]
#faccio sampling ogni 100 frame
indices = np.arange(0, x.size, 100)
x   =   x[indices]
y   =   y[indices]




fig = plt.figure()
c = range(x.size)
plt.scatter(x,y, c=c, cmap= 'viridis', s = 10)
plt.colorbar()
plt.xlabel(' PC 1 ')
plt.ylabel(' PC 2 ')
plt.title('Projection of MD trajectory in Essential Space')
plt.show()
fig.savefig('3 micros/Essential.png')

#%%
xy = np.array([x,y]).T
RNA2_Clust = BA.Cluster_2DAnalysis(xy)
RNA2_Clust.Silhouette_KMeans(15)

# %%

RNA2_Clust.Clusterize(verbose=True)
RNA2_Clust.Figure_Clusters()

#%%

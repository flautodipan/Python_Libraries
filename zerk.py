#%%
import BioAlessandria as BA
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

#%%


#fig =   plt.figure()
dist        =   ()
means       =   np.zeros(5)

plt.figure()
plt.title('Distribution for distance between Zernike momenta')
plt.xlim(0.1, 0.5)
dynamics = ['wtc1', 'wtc2', 'wtc3', 'wtc4', 'wtc5']
colors = ['blue', 'red', 'yellow', 'green', 'black']
i_order =   [1,4,5,3,2]

for (n, i, color) in zip(i_order, range(5), colors):
    
    now_path    =   '../MD/GROMACS/ZERNiKE/'
    file_RNA    =   'df_zernike_shape_RNA_zerk'+str(n)+'.csv'
    file_PROT   =   'df_zernike_shape_protein_zerk'+str(n)+'.csv'

    #1) leggo i file e li salvo in DataFrame

    PROT_df = pd.read_csv(now_path+file_PROT, delimiter = ' ')#, dtype=np.float64, usecols=[1:])
    RNA_df  = pd.read_csv(now_path+file_RNA, delimiter = ' ')# dtype=np.float64, usecols=[1:])

    print("Ho letto correttamente i file in input per WTC"+str(n))

    #2) calcolo vettore distanza tra i vettori di zernike di PROT e RNA
    #   ad ogni frame

    dist = dist + (np.sum(np.abs(PROT_df - RNA_df), axis = 1),)
    means[i]=np.mean(dist[i].values)
    #3) plotto
    
    _, _, _ =  plt.hist(dist[i].values, bins = 150, alpha = 0.6, label = 'WTC'+str(n), histtype='barstacked')
    i+=1

plt.legend()
plt.savefig(now_path+'zerk_dist_distribution.png')


# %%
#  plot correlazione

K_d     =   [4, 700, 1320, 1350, 1360]
i_order =   [1,4,5,3,2]

plt.figure()
plt.title('Correlation mean_zernike vs K_d exp')
plt.xlabel('K_d exp')
plt.ylabel('K_d zerk')
plt.scatter(K_d, means)

for (i,dyn) in zip(range(5), dynamics):
    plt.annotate(str(i_order[i]), (K_d[i], means[i]))
plt.savefig(now_path+'zerk_correlation.png')


# %%
np.array(dist)

# %%

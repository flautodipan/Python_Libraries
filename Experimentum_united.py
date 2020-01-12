#%%
import      os
import      numpy               as      np
import      matplotlib.pyplot   as      plt
import      lib_Experimentum    as      Exp
from        Alessandria         import  *
import      time
# %%
#0) Acquisisco dati e inizializzo oggetti Spectrum per ognuno su una matrice (n_rows, n_cols)
#   e compio alcune operazioni di sistema utili per salvataggio dati

spectra_path        =   '../Claudia/DaticellBoniPuntiDoppi/'
spectra_filename    =   '20191218_K27M'

VIPA_path           =   '../Claudia/DaticellBoniPuntiDoppi/picchi_elastici_con_filtro_100msexp/Pos0/'
VIPA_filename       =   'img_000000000_Default_000'

os.system('cd .. & mkdir '+ spectra_filename+'_analysis')
now_path            =   '../'+spectra_filename+'_analysis/'

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, spectra_path, var_name = 'y')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
matrix = Exp.Initialize_Matrix(n_rows,n_cols)


# %%
#1) Acquisisco VIPA e Spettri
start = time.process_time()
for ii in range(n_rows):
    for jj in range(n_cols):
        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        matrix[ii][jj].Get_Spectrum(how_to_get = 'by_passing', y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183.)

print('tempo impiegato per acquisizione spettri: %f s'%(time.process_time()-start))

plt.figure()

plt.plot(matrix[9][8].y)
plt.plot(matrix[74][43].y)
plt.plot(matrix[31][81].y)
plt.savefig(now_path+'esempio_spettri.png')
plt.close()


# %%

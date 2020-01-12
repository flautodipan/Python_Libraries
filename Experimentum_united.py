#%%
import      os
spectra_path        =   '../Claudia/DaticellBoniPuntiDoppi/'
spectra_filename    =   '20191218_K27M'

VIPA_path           =   '../Claudia/DaticellBoniPuntiDoppi/picchi_elastici_con_filtro_100msexp/Pos0/'
VIPA_filename       =   'img_000000000_Default_000.tif'

os.system('cd .. & mkdir '+ spectra_filename+'_analysis')
now_path            =   '../'+spectra_filename+'_analysis/'


#%%
import      numpy               as      np
import      matplotlib.pyplot   as      plt
import      lib_Experimentum    as      Exp
from        Alessandria         import  *
import      time
# %%
#0) Acquisisco dati e inizializzo oggetti Spectrum per ognuno su una matrice (n_rows, n_cols)
#   e compio alcune operazioni di sistema utili per salvataggio dati

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, spectra_path, var_name = 'y')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
matrix = Exp.Initialize_Matrix(n_rows,n_cols)

#definisco quantità di interesse

saturated   =   ()
excluded    =   () 
fitted      =   ()  
non_fitted  =   ()
failed      =   ()



# %%
#1) Acquisisco VIPA e Spettri
start = time.process_time()
matrix[0][0].Get_VIPA_tif(VIPA_filename, VIPA_path, offset = 183., )

for ii in range(n_rows):
    for jj in range(n_cols):
        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183.)
        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA
        
        if matrix[ii][jj].Is_Saturated(treshold = 40000 ):
            saturated   =   saturated   +   ((ii,jj), )


print('tempo impiegato per acquisizione spettri: %f s'%(time.process_time()-start))
print('Totale spettri saturati : %d'%(len(saturated)/2), saturated)


# %%
#2) Faccio operazioni di modifica spettro

matrix[0][0].How_Many_Peaks_To_VIPA(treshold = 30)
matrix[0][0].Fit_Pixel2GHz(fig = True)
matrix[0][0].VIPA_Pix2GHz(fig=True)
matrix[0][0].Spectrum_Pix2GHz(fig=True)
matrix[0][0].Fit_VIPA_Gaussian()

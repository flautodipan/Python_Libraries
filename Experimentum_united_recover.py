#%%
######### VERSIONE EXP_UNITED CHE RECUPERA INFORMAZIONI FIT  #############

# idea è la parte di acquisizione e modifica è uguale, sono le fit info ad essere recuperate invece che calcolate
# così si può strutturare una versione da terminale e poi un recupero delle info da qua


#%%

#%%
import      os
spectra_path        =   '../BRILLOUIN/Claudia/DaticellBoniPuntiDoppi/'
spectra_filename    =   '20191218_K27M'

VIPA_path           =   '../BRILLOUIN/Claudia/DaticellBoniPuntiDoppi/picchi_elastici_con_filtro_100msexp/Pos0/'
VIPA_filename       =   'img_000000000_Default_000.tif'

os.system('cd .. & mkdir '+ spectra_filename+'_analysis')
now_path            =   '../'+spectra_filename+'_analysis/'


cols      = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')

#%%
import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Experimentum    import  *
from        Alessandria         import  *
import      time

super_start         =   time.process_time()
tempo               =   ()

#%%
#0) Acquisisco dati e inizializzo oggetti Spectrum per ognuno su una matrice (n_rows, n_cols)
#   e compio alcune operazioni di sistema utili per salvataggio dati

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, spectra_path, var_name = 'y')
n_rows  =   5#len(dati)
n_cols  =   5#len(dati[0])
dim     =   n_cols*n_rows
matrix = Initialize_Matrix(n_rows,n_cols)

#definisco quantità di interesse

invisible           =   () 
saturated           =   () 
brillouin_higher    =   ()
boni                =   ()


syg_kwargs   =   {'height': 20, 'distance': 20, 'width': 5.}

#
"""

1) implementare semplicissima recover della tupla fit con eval() dell'intero testo recuperato dalla lettura di fit.txt
2) fare Unpack_Fit() di questa tupla 
3) implementare semplice ciclo da qua utilizzando per ogni classe il metodo Recover_Fit_Params()

e a questo punto hai recuperato tutte le info importanti

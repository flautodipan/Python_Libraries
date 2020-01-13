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

excluded            =   () 
fitted              =   ()  
non_fitted          =   ()
failed              =   ()
saturated           =   () 
brillouin_higher    =   ()


syg_kwargs   =   {'height': 20, 'distance': 25, 'width': 5.}
# %%
#1) Acquisisco VIPA e Spettri
start = time.process_time()
matrix[0][0].Get_VIPA_tif(VIPA_filename, VIPA_path, offset = 183.)


for ii in range(n_rows):
    for jj in range(n_cols):
        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183., cut = True, cut_range = (200, 600))
        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs)

        check   =   matrix[ii][jj].Check_Spectrum(saturation_width = 15.)
        if (check == 1):

            saturated   =   saturated   +   ((ii,jj), )

        elif (check == 2):

            brillouin_higher    =   brillouin_higher    +   ((ii,jj),)




print('tempo impiegato per acquisizione spettri: %f s'%(time.process_time()-start))
print('Totale spettri saturati : %d\n'%(len(saturated)/2), saturated)
print('Totale spettri con Brillouin strani : %d\n'%(len(brillouin_higher)/2), brillouin_higher)

# %%
#2) Faccio operazioni di modifica spettro
start = time.process_time()

matrix[0][0].How_Many_Peaks_To_VIPA(treshold = 30)
matrix[0][0].Fit_Pixel2GHz(fig = True)
matrix[0][0].VIPA_Pix2GHz(fig=True)

matrix[0][0].Spectrum_Pix2GHz(fig=True)
matrix[0][0].How_Many_Peaks_To()
matrix[0][0].Cut_n_Estimate_Spectrum(estimate = True)
matrix[0][0].Fit_VIPA_Gaussian()


for ii in range(n_rows):
    for jj in range(n_cols):

        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))

        if ((ii,jj) != (0,0) ) & ((ii,jj) not in saturated):
                    
            matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA_freq
            matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA
            
            matrix[ii][jj].Poly2GHz =   matrix[0][0].Poly2GHz
            matrix[ii][jj].Spectrum_Pix2GHz()

            if (ii,jj) in i_know:
                matrix[ii][jj].How_Many_Peaks_To(i_know_it_is = True)
            else:
                matrix[ii][jj].How_Many_Peaks_To()
            matrix[ii][jj].Cut_n_Estimate_Spectrum()


print('tempo impiegato per modifica spettri: %f s'%(time.process_time()-start))
#%%
for (ii,jj) in brillouin_higher:
    plt.figure()
    plt.plot(matrix[ii][jj].y)
    plt.title(str((ii,jj)))

# %%
plt.plot(matrix[6][45].y)

# %%

for ii in range(2):
    for jj in range(10):
        pass
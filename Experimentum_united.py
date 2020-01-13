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
dim     =   n_cols*n_rows
matrix = Exp.Initialize_Matrix(n_rows,n_cols)

#definisco quantità di interesse

invisible           =   () 
saturated           =   () 
brillouin_higher    =   ()
boni                =   ()


syg_kwargs   =   {'height': 20, 'distance': 20, 'width': 5.}
# %%
#1) Acquisisco VIPA e Spettri
start = time.process_time()
matrix[0][0].Get_VIPA_tif(VIPA_filename, VIPA_path, offset = 183.)


for ii in range(n_rows):
    for jj in range(n_cols):
        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183., cut = True, cut_range = (200, 600))
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs)
        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA

        check   =   matrix[ii][jj].Check_Spectrum(saturation_width = 13.5)
        
        if (check == 1):

            saturated   =   saturated   +   ((ii,jj), )

        elif (check == 2):

            brillouin_higher    =   brillouin_higher    +   ((ii,jj),)
            boni                =   boni                +   ((ii,jj),)

        elif (check == 3):

            invisible           =   invisible           +   ((ii,jj),)
        
        else:

            boni                =   boni                +   ((ii,jj),)


print('tempo impiegato per acquisizione spettri: %f s'%(time.process_time()-start))
print('Totale spettri saturati : %d\n'%(len(saturated)), saturated)
print('Totale spettri con Brillouin più alti : %d\n'%(len(brillouin_higher)), brillouin_higher)
print('Totale spettri con Brillouin invisibili: %d\n'%(len(invisible)), invisible)
print('Totale spettri boni :   %d'%(len(boni)))
print('Totale   spettri : %d\ndi cui %d inutilizzabili'%(dim, len(invisible)+len(saturated)))
print('ossia il %3.2f percento'%(float(len(invisible)+len(saturated))*100/dim))
# %%
#2) Faccio operazioni di modifica spettro
excluded    =   saturated + invisible

start = time.process_time()

matrix[0][0].How_Many_Peaks_To_VIPA(treshold = 30)
matrix[0][0].Fit_Pixel2GHz(fig = True)
matrix[0][0].VIPA_Pix2GHz(fig=True)

matrix[0][0].Spectrum_Pix2GHz(fig=True)
matrix[0][0].Get_Spectrum_4_Peaks_by_Height()
matrix[0][0].Cut_n_Estimate_Spectrum(estimate = True, distanza = 0.25)
matrix[0][0].Fit_VIPA_Gaussian()


for ii in range(n_rows):
    for jj in range(n_cols):

        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))

        if ((ii,jj) != (0,0)) & ((ii,jj) not in excluded):
                    
            matrix[ii][jj].x_VIPA_freq   =   matrix[0][0].x_VIPA_freq
            matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA
            
            matrix[ii][jj].Poly2GHz =   matrix[0][0].Poly2GHz
            matrix[ii][jj].Spectrum_Pix2GHz()

            if (ii,jj) in brillouin_higher:
            
                matrix[ii][jj].Get_Spectrum_Peaks(height = 5., distance = 50, width = 5)
                matrix[ii][jj].Get_Spectrum_4_Peaks_by_Order()
            
            else :
                
                matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                

            matrix[ii][jj].Cut_n_Estimate_Spectrum(distanza = 0.25)


print('tempo impiegato per modifica spettri: %f s'%(time.process_time()-start))


#%%
#3) faccio il primo fit markoviano

fit                 =   ()
fitted              =   ()  
non_fitted          =   ()
failed              =   ()
isolated = Exp.Get_Isolated_Elements(excluded)

#prima riga, stimo da sx, eccetto il primo
start = time.process_time()

ii = 0
for jj in range(n_cols):
    print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
    if (ii,jj) in boni:
        if (ii,jj) in isolated:
            fit = fit + (matrix[ii][jj].Estimate_Initial_Parameters(matrix[0][0].p0.values[0], 1000), (ii,jj),)
        else:
            fit = fit + (matrix[ii][jj].Estimate_Initial_Parameters(matrix[ii][jj-1].p0.values[0], 1000), (ii,jj),)

#tutto il resto
#%%
for ii in range(1, n_rows,1 ):
    for jj in range(n_cols):
        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        if (ii,jj) in boni:
            if (ii,jj) in isolated:
                fit =   fit + (matrix[ii][jj].Estimate_Initial_Parameters(matrix[0][0].p0.values[0], 1000), (ii,jj),)
            else: 
                fit =   fit + (matrix[ii][jj].Estimate_Initial_Parameters(matrix[ii-1][jj].p0.values[0], 1000), (ii,jj),)

print('tempo impiegato per fit markoviani: %f s'%(time.process_time()-start))






#%%
#VERIFICHE DA QUI IN POI






















#%%

#stima dei recuperati, ma anche scialla: problema è come li distingui
# idee  : sicuramente quelli con un picco a sx dell'elastico
#         e sicuramente quelli con picchi Brillouin troppo vicini a elastici, o troppo vicini tra loro

for (ii,jj) in invisible:
    print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
    matrix[ii][jj].How_Many_Peaks_To(width = 1)
    plt.figure()
    plt.plot(matrix[ii][jj].x_pix, matrix[ii][jj].y)
    plt.title(str((ii,jj)))
    idx     = find_peaks(matrix[ii][jj].y, height = matrix[ii][jj].spectrum_cut_height, distance = matrix[ii][jj].spectrum_peaks_dist, width = 1.)[0]
    plt.plot(matrix[ii][jj].x_pix[idx],matrix[ii][jj].y[idx], '*' )

# %%
for (ii,jj) in brillouin_higher[0:20]:
    matrix[ii][jj].Get_Spectrum_Peaks(height = 5, distance = 50, width = 5.)
    plt.figure()
    plt.plot(matrix[ii][jj].x_pix, matrix[ii][jj].y)
    plt.plot(matrix[ii][jj].x_pix[matrix[ii][jj].peaks[0]], matrix[ii][jj].y[matrix[ii][jj].peaks[0]], '*')
    plt.title(str((ii,jj))+str(matrix[ii][jj].peaks[1]['widths'][matrix[ii][jj].peaks[1]['peak_heights'].argmax()]))
# %%
a = ()
for (ii,jj) in brillouin_higher:

    matrix[ii][jj].Get_Spectrum_Peaks(height = 5., distance = 50, width = 5)
    a = a + (matrix[ii][jj].n_peaks, )
    
    if matrix[ii][jj].n_peaks == 7:
        print(str((ii,jj)))



# %%

for (ii,jj) in boni[640:660]:
    plt.figure()
    plt.plot(matrix[ii][jj].x_freq, matrix[ii][jj].y)
    plt.title(str((ii,jj)))
    #plt.plot(matrix[ii][jj].x_freq[matrix[ii][jj].peaks['peaks_idx']], matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx']], '*')

# %%

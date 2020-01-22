#%%
import      os
spectra_path        =   '../Claudia/DaticellBoniPuntiDoppi/'
spectra_filename    =   '20191218_K27M'

VIPA_path           =   '../Claudia/DaticellBoniPuntiDoppi/picchi_elastici_con_filtro_100msexp/Pos0/'
VIPA_filename       =   'img_000000000_Default_000.tif'

os.system('cd .. & mkdir '+ spectra_filename+'_analysis')
now_path            =   '../'+spectra_filename+'_analysis/'


cols        = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')

#%%
import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Experimentum    import  *
from        Alessandria         import  *
import      time

super_start         =   time.process_time()

# %%
#0) Acquisisco dati e inizializzo oggetti Spectrum per ognuno su una matrice (n_rows, n_cols)
#   e compio alcune operazioni di sistema utili per salvataggio dati

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, spectra_path, var_name = 'y')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
dim     =   n_cols*n_rows
matrix = Initialize_Matrix(n_rows,n_cols)

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
                #boni
                matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                

            matrix[ii][jj].Cut_n_Estimate_Spectrum(distanza = 0.25)


print('tempo impiegato per modifica spettri: %f s'%(time.process_time()-start))


#%%
#3) faccio il primo fit markoviano
fit                 =   ()
start = time.process_time()
recover = False

if recover:

    params  =   Parse_Parameter_Save(path=now_path+'save_15_01/')
    accomplished  =   ((0, 0),(0, 15),(0, 17),(0, 19),(0, 27),(0, 29),(0, 31),(0, 35),(0, 37),(0, 39),(0, 45),(0, 49),(0, 51),(0, 52),(0, 53),(0, 59),(0, 64),(0, 71),(0, 73),(0, 77),(0, 79),(0, 83),(1, 10),(1, 36),(1, 37),(1, 51),(1, 55),(1, 56),(2, 35),(2, 36),(2, 37),(2, 38),(2, 39),(2, 52),(2, 59),(3, 35),(3, 37),(3, 39),(3, 48),(3, 51),(3, 60),(4, 15),(4, 16),(4, 35),(4, 36),(4, 41),(4, 42),(4, 51),(5, 35),(5, 36),(5, 38),(5, 48),(5, 52),(5, 55),(6, 34),(6, 35),(6, 36),(6, 42),(6, 51),(6, 52),(6, 55),(6, 56),(7, 10),(7, 30),(7, 34),(7, 35),(7, 36),(7, 37),(7, 38),(7, 44),(7, 51),(7, 52),(7, 56),(8, 31),(8, 36),(8, 37),(8, 41),(8, 45),(8, 47),(8, 52),(9, 15),(9, 21),(9, 30),(9, 34),(9, 35),(9, 36),(9, 37),(9, 38),(9, 40),(9, 41),(9, 42),(9, 45),(9, 47),(9, 48),(9, 52),(9, 56),(10, 35),(10, 36),(10, 37),(10, 42),(10, 45),(10, 48),(10, 52),(11, 34),(11, 40),(11, 41),(11, 43),(11, 44),(11, 47),(11, 48),(11, 51),(11, 52),(12, 36),(12, 37),(12, 39),(12, 40),(12, 41),(12, 42),(12, 44),(12, 51),(12, 52),(13, 34),(13, 35),(13, 37),(13, 38),(13, 40),(13, 41),(13, 42),(13, 52),(14, 31),(14, 33),(14, 34),(14, 35),(14, 36),(14, 38),(14, 40),(14, 41),(14, 42),(14, 47),(14, 48),(14, 51),(14, 52),(14, 53),(14, 60),(15, 31),(15, 33),(15, 35),(15, 37),(15, 39),(15, 42),(15, 45),(15, 47),(15, 60),(16, 28),(16, 36),(16, 37),(16, 40),(16, 41),(16, 42),(16, 45),(16, 48),(16, 52),(16, 53),(17, 28),(17, 34),(17, 35),(17, 37),(17, 38),(17, 39),(17, 40),(17, 41),(17, 42),(17, 45),(17, 48),(17, 51),(17, 52),(17, 53),(18, 28),(18, 34),(18, 36),(18, 40),(18, 41),(18, 42),(18, 45),(18, 48),(18, 52),(19, 11),(19, 34),(19, 36),(19, 38),(19, 39),(19, 40),(19, 51),(19, 52),(19, 53),(19, 54),(20, 28),(20, 32),(20, 35),(20, 36),(20, 38),(20, 39),(20, 40),(20, 45),(20, 47),(20, 48),(20, 51),(20, 52),(20, 53),(21, 28),(21, 31),(21, 32),(21, 34),(21, 35),(21, 37),(21, 38),(21, 41),(21, 42),(21, 45),(21, 47),(21, 48),(21, 51),(21, 52),(22, 29),(22, 31),(22, 35),(22, 36),(22, 38),(22, 40),(22, 41),(22, 43),(22, 45),(22, 47),(22, 48),(23, 27),(23, 35),(23, 36),(23, 37),(23, 41),(23, 45),(23, 47),(23, 48),(23, 50),(23, 52),(24, 29),(24, 30),(24, 31),(24, 34),(24, 35),(24, 38),(24, 39),(24, 40),(24, 41),(24, 42),(24, 43),(24, 45),(24, 51),(24, 52),(25, 10),(25, 24),(25, 28),(25, 29),(25, 30),(25, 32),(25, 34),(25, 35),(25, 36),(25, 37),(25, 39),(25, 41),(25, 44),(25, 45),(25, 50),(25, 52),(25, 55),(25, 56),(25, 63),(26, 28),(26, 31),(26, 32),(26, 33),(26, 34),(26, 35),(26, 39),(26, 40),(26, 42),(26, 44),(26, 48),(26, 50),(26, 51),(26, 54),(26, 55),(26, 56),(27, 9),(27, 33),(27, 34),(27, 35),(27, 37),(27, 38),(27, 39),(27, 45),(27, 48),(27, 51),(27, 53),(27, 55),(27, 72),(28, 4),(28, 14),(28, 29),(28, 30),(28, 32),(28, 34),(28, 36),(28, 40),(28, 42),(28, 43),(28, 50),(28, 51),(28, 52),(28, 54),(28, 59),(29, 15),(29, 24),(29, 25),(29, 29),(29, 30),(29, 31),(29, 33),(29, 36),(29, 38),(29, 40),(29, 41),(29, 43),(29, 45),(29, 47),(29, 50),(29, 59),(29, 60),(30, 10),(30, 11),(30, 21),(30, 25),(30, 30),(30, 32),(30, 41),(30, 42),(30, 47),(30, 48),(30, 50),(30, 53),(30, 56),(30, 59),(31, 8),(31, 9),(31, 10),(31, 11),(31, 15),(31, 16),(31, 25),(31, 29),(31, 30),(31, 32),(31, 35),(31, 36),(31, 38),(31, 41),(31, 43),(31, 44),(31, 47),(31, 48),(31, 50),(31, 51),(31, 52),(31, 59),(32, 8),(32, 9),(32, 10),(32, 15),(32, 16),(32, 21),(32, 22),(32, 27),(32, 31),(32, 32),(32, 33),(32, 34),(32, 37),(32, 38),(32, 41),(32, 45),(32, 47),(32, 50),(32, 51),(32, 59),(32, 65),(33, 5),(33, 7),(33, 8),(33, 10),(33, 11),(33, 13),(33, 14),(33, 15),(33, 19),(33, 22),(33, 24),(33, 25),(33, 27),(33, 28),(33, 29),(33, 31),(33, 34),(33, 36),(33, 37),(33, 38),(33, 39),(33, 41),(33, 42),(33, 43),(33, 44),(33, 47),(33, 50),(33, 51),(33, 52),(33, 53),(33, 56),(33, 60),(33, 62),(34, 4),(34, 5),(34, 7),(34, 11),(34, 13),(34, 14),(34, 15),(34, 16),(34, 17),(34, 18),(34, 20),(34, 22),(34, 23),(34, 25),(34, 26),(34, 28),(34, 29),(34, 30),(34, 31),(34, 35),(34, 36),(34, 37),(34, 40),(34, 41),(34, 43),(34, 44),(34, 45),(34, 47),(34, 50),(34, 51),(34, 58),(34, 59),(34, 69),(35, 5),(35, 9),(35, 11),(35, 12),(35, 13),(35, 14),(35, 17),(35, 20),(35, 21),(35, 22),(35, 23),(35, 24),(35, 25),(35, 26),(35, 28),(35, 29),(35, 30),(35, 37),(35, 38),(35, 40),(35, 41),(35, 42),(35, 43),(35, 45),(35, 48),(35, 51),(35, 52),(35, 56),(35, 58),(35, 59),(35, 65),(36, 4),(36, 5),(36, 9),(36, 10),(36, 11),(36, 12),(36, 13),(36, 14),(36, 15),(36, 18),(36, 19),(36, 21),(36, 23),(36, 24),(36, 27),(36, 28),(36, 29),(36, 30),(36, 32),(36, 33),(36, 34),(36, 36),(36, 37),(36, 39),(36, 47),(36, 51),(36, 53),(36, 56),(36, 58),(36, 59),(36, 60),(36, 61),(36, 69),(36, 70),(37, 8),(37, 10),(37, 13),(37, 14),(37, 15),(37, 17),(37, 18),(37, 19),(37, 20),(37, 21),(37, 24),(37, 26),(37, 28),(37, 29),(37, 32),(37, 33),(37, 34),(37, 35),(37, 36),(37, 37),(37, 41),(37, 42),(37, 43),(37, 44),(37, 45),(37, 46),(37, 47),(37, 49),(37, 50),(37, 52),(37, 53),(37, 56),(37, 59),(37, 60),(37, 62),(37, 70),(38, 8),(38, 9),(38, 14),(38, 17),(38, 18),(38, 20),(38, 21),(38, 22),(38, 23),(38, 25),(38, 26),(38, 27),(38, 28),(38, 31),(38, 32),(38, 34),(38, 35),(38, 36),(38, 39),(38, 43),(38, 44),(38, 46),(38, 49),(38, 52),(38, 53),(38, 56),(38, 59),(38, 60),(38, 62),(38, 71),(38, 72),(39, 10),(39, 14),(39, 15),(39, 17),(39, 18),(39, 21),(39, 24),(39, 25),(39, 27),(39, 29),(39, 31),(39, 32),(39, 34),(39, 38),(39, 39),(39, 42),(39, 43),(39, 45),(39, 46),(39, 47),(39, 48),(39, 52),(39, 56),(39, 62),(39, 71),(39, 72),(39, 77),(40, 9),(40, 14),(40, 17),(40, 21),(40, 22),(40, 24),(40, 25),(40, 27),(40, 28),(40, 29),(40, 30),(40, 32),(40, 34),(40, 37),(40, 38),(40, 39),(40, 40),(40, 41),(40, 42),(40, 44),(40, 46),(40, 47),(40, 48),(40, 50),(40, 53),(40, 56),(40, 60),(41, 13),(41, 14),(41, 15),(41, 17),(41, 18),(41, 21),(41, 25),(41, 26),(41, 28),(41, 30),(41, 33),(41, 37),(41, 38),(41, 39),(41, 40),(41, 41),(41, 44),(41, 45),(41, 46),(41, 48),(41, 51),(41, 53),(41, 55),(41, 56),(41, 72),(42, 10),(42, 11),(42, 13),(42, 14),(42, 15),(42, 17),(42, 18),(42, 25),(42, 26),(42, 28),(42, 29),(42, 30),(42, 33),(42, 34),(42, 35),(42, 37),(42, 43),(42, 44),(42, 47),(42, 48),(42, 49),(42, 50),(42, 52),(42, 53),(42, 54),(42, 56),(43, 10),(43, 12),(43, 13),(43, 14),(43, 15),(43, 16),(43, 17),(43, 20),(43, 24),(43, 25),(43, 29),(43, 30),(43, 31),(43, 36),(43, 39),(43, 43),(43, 44),(43, 45),(43, 50),(43, 55),(43, 56),(43, 58),(44, 9),(44, 14),(44, 16),(44, 17),(44, 19),(44, 20),(44, 21),(44, 23),(44, 24),(44, 27),(44, 28),(44, 29),(44, 30),(44, 34),(44, 35),(44, 38),(44, 39),(44, 41),(44, 42),(44, 43),(44, 44),(44, 45),(44, 51),(44, 53),(44, 56),(45, 10),(45, 11),(45, 12),(45, 13),(45, 14),(45, 15),(45, 21),(45, 23),(45, 24),(45, 29),(45, 31),(45, 35),(45, 36),(45, 39),(45, 40),(45, 41),(45, 42),(45, 43),(45, 44),(45, 45),(45, 46),(45, 54),(45, 60),(45, 78),(46, 4),(46, 10),(46, 11),(46, 16),(46, 18),(46, 19),(46, 22),(46, 23),(46, 24),(46, 25),(46, 27),(46, 28),(46, 29),(46, 31),(46, 33),(46, 35),(46, 37),(46, 39),(46, 40),(46, 43),(46, 44),(46, 45),(46, 46),(46, 47),(46, 49),(46, 50),(46, 53),(46, 60),(47, 4),(47, 15),(47, 17),(47, 18),(47, 21),(47, 25),(47, 26),(47, 27),(47, 28),(47, 30),(47, 32),(47, 34),(47, 35),(47, 37),(47, 39),(47, 40),(47, 41),(47, 42),(47, 43),(47, 45),(47, 49),(47, 51),(47, 52),(47, 53),(48, 9),(48, 14),(48, 15),(48, 16),(48, 17),(48, 18),(48, 26),(48, 28),(48, 29),(48, 30),(48, 31),(48, 32),(48, 33),(48, 35),(48, 41),(48, 42),(48, 43),(48, 44),(48, 45),(48, 47),(48, 52),(48, 53),(48, 54),(48, 56),(49, 14),(49, 15),(49, 17),(49, 18),(49, 20),(49, 21),(49, 22),(49, 24),(49, 25),(49, 26),(49, 30),(49, 31),(49, 32),(49, 35),(49, 39),(49, 40),(49, 41),(49, 43),(49, 49),(49, 52),(49, 54),(49, 55),(49, 56),(49, 60),(49, 72),(50, 9),(50, 10),(50, 14),(50, 15),(50, 18),(50, 20),(50, 22),(50, 24),(50, 26),(50, 27),(50, 28),(50, 29),(50, 30),(50, 31),(50, 33),(50, 34),(50, 37),(50, 39),(50, 41),(50, 42),(50, 44),(50, 45),(50, 46),(50, 48),(50, 49),(50, 50),(50, 51),(50, 55),(50, 60),(51, 15),(51, 17),(51, 19),(51, 21),(51, 25),(51, 26),(51, 27),(51, 29),(51, 31),(51, 33),(51, 34),(51, 37),(51, 41),(51, 42),(51, 43),(51, 45),(51, 47),(51, 50),(51, 55),(51, 59),(51, 71),(52, 15),(52, 17),(52, 18),(52, 21),(52, 24),(52, 25),(52, 27),(52, 28),(52, 30),(52, 33),(52, 35),(52, 36),(52, 37),(52, 38),(52, 39),(52, 40),(52, 44),(52, 47),(52, 49),(52, 51),(52, 52),(52, 53),(52, 59),(53, 15),(53, 17),(53, 18),(53, 24),(53, 25),(53, 26),(53, 27),(53, 33),(53, 34),(53, 41),(53, 43),(53, 44),(53, 46),(53, 51),(53, 52),(53, 64),(54, 13),(54, 14),(54, 15),(54, 18),(54, 21),(54, 23),(54, 25),(54, 26),(54, 27),(54, 30),(54, 31),(54, 32),(54, 33),(54, 35),(54, 36),(54, 37),(54, 40),(54, 41),(54, 43),(54, 45),(54, 47),(54, 48),(54, 51),(54, 52),(54, 53),(54, 59),(55, 15),(55, 17),(55, 18),(55, 28),(55, 31),(55, 32),(55, 33),(55, 34),(55, 35),(55, 36),(55, 37),(55, 42),(55, 43),(55, 44),(55, 47),(55, 48),(55, 56),(55, 59),(56, 13),(56, 14),(56, 17),(56, 18),(56, 21),(56, 23),(56, 24),(56, 25),(56, 26),(56, 28),(56, 30),(56, 32),(56, 35),(56, 36),(56, 37),(56, 39),(56, 40),(56, 42),(56, 43),(56, 45),(56, 47),(56, 51),(56, 52),(56, 54),(56, 56),(56, 59),(57, 15),(57, 18),(57, 19),(57, 21),(57, 23),(57, 25),(57, 28),(57, 31),(57, 32),(57, 33),(57, 34),(57, 35),(57, 39),(57, 41),(57, 42),(57, 43),(57, 44),(57, 45),(57, 47),(57, 48),(57, 49),(57, 50),(57, 51),(57, 53),(57, 59),(58, 4),(58, 18),(58, 20),(58, 21),(58, 23),(58, 24),(58, 27),(58, 28),(58, 32),(58, 35),(58, 37),(58, 39),(58, 40),(58, 41),(58, 43),(58, 44),(58, 47),(58, 49),(58, 50),(58, 51),(58, 52),(58, 53),(58, 62),(59, 14),(59, 17),(59, 24),(59, 26),(59, 27),(59, 28),(59, 29),(59, 30),(59, 31),(59, 33),(59, 36),(59, 37),(59, 39),(59, 41),(59, 42),(59, 43),(59, 44),(59, 45),(59, 47),(59, 48),(59, 49),(59, 50),(59, 52),(59, 54),(59, 72),(60, 14),(60, 17),(60, 18),(60, 21),(60, 26),(60, 28),(60, 29),(60, 30),(60, 31),(60, 33),(60, 34),(60, 36),(60, 37),(60, 43),(60, 44),(60, 46),(60, 51),(60, 52),(60, 56),(60, 62),(60, 72),(61, 17),(61, 18),(61, 20),(61, 21),(61, 22),(61, 27),(61, 28),(61, 30),(61, 31),(61, 33),(61, 34),(61, 36),(61, 37),(61, 38),(61, 41),(61, 42),(61, 45),(61, 49),(61, 50),(61, 52),(61, 53),(61, 72),(62, 21),(62, 22),(62, 24),(62, 28),(62, 30),(62, 32),(62, 33),(62, 34),(62, 36),(62, 37),(62, 38),(62, 40),(62, 42),(62, 44),(62, 45),(62, 46),(62, 47),(62, 48),(62, 49),(62, 50),(62, 52),(62, 53),(63, 18),(63, 19),(63, 20),(63, 22),(63, 23),(63, 24),(63, 25),(63, 26),(63, 28),(63, 30),(63, 31),(63, 32),(63, 33),(63, 34),(63, 35),(63, 37),(63, 41),(63, 44),(63, 45),(63, 47),(63, 48),(63, 49),(63, 52),(64, 19),(64, 22),(64, 23),(64, 24),(64, 25),(64, 26),(64, 30),(64, 31),(64, 33),(64, 34),(64, 35),(64, 36),(64, 37),(64, 38),(64, 42),(64, 47),(64, 48),(64, 52),(64, 53),(64, 71),(65, 19),(65, 21),(65, 23),(65, 24),(65, 25),(65, 28),(65, 31),(65, 33),(65, 34),(65, 35),(65, 37),(65, 40),(65, 42),(65, 44),(65, 47),(65, 48),(65, 49),(65, 52),(66, 21),(66, 23),(66, 24),(66, 25),(66, 28),(66, 30),(66, 32),(66, 33),(66, 34),(66, 35),(66, 37),(66, 39),(66, 40),(66, 42),(66, 44),(66, 45),(66, 48),(67, 18),(67, 22),(67, 23),(67, 25),(67, 28),(67, 32),(67, 34),(67, 35),(67, 37),(67, 39),(67, 40),(67, 41),(67, 42),(67, 44),(67, 53),(67, 71),(68, 21),(68, 23),(68, 25),(68, 26),(68, 28),(68, 30),(68, 32),(68, 33),(68, 34),(68, 35),(68, 37),(68, 38),(68, 41),(68, 42),(68, 44),(68, 45),(68, 47),(69, 14),(69, 23),(69, 24),(69, 25),(69, 26),(69, 27),(69, 28),(69, 33),(69, 34),(69, 38),(69, 39),(69, 41),(69, 42),(69, 45),(70, 26),(70, 27),(70, 31),(70, 34),(70, 35),(70, 37),(70, 38),(70, 42),(70, 43),(70, 44),(70, 45),(70, 46),(70, 47),(70, 48),(70, 73),(71, 19),(71, 23),(71, 25),(71, 26),(71, 31),(71, 33),(71, 34),(71, 37),(71, 38),(71, 39),(71, 46),(71, 48),(72, 21),(72, 24),(72, 25),(72, 26),(72, 28),(72, 30),(72, 32),(72, 33),(72, 35),(72, 39),(72, 42),(73, 19),(73, 23),(73, 26),(73, 27),(73, 33),(73, 34),(73, 36),(73, 38),(73, 39),(73, 71),(74, 15),(74, 19),(74, 20),(74, 21),(74, 25),(74, 27),(74, 29),(74, 31),(74, 32),(74, 33),(74, 35),(74, 36),(74, 37),(74, 38),(74, 39),(74, 73),(75, 19),(75, 21),(75, 23),(75, 24),(75, 25),(75, 26),(75, 27),(75, 29),(75, 30),(75, 33),(75, 34),(75, 35),(75, 37),(75, 38),(75, 39),(75, 42),(75, 48),(76, 20),(76, 22),(76, 23),(76, 25),(76, 27),(76, 28),(76, 30),(76, 34),(76, 38),(76, 41),(76, 42),(77, 22),(77, 24),(77, 26),(77, 28),(77, 30),(77, 33),(77, 34),(77, 36),(77, 39),(77, 41),(77, 42),(77, 82),(78, 14),(78, 20),(78, 22),(78, 23),(78, 25),(78, 26),(78, 27),(78, 29),(78, 31),(78, 32),(78, 33),(78, 35),(78, 38),(78, 39),(78, 42),(78, 49),(79, 17),(79, 19),(79, 23),(79, 24),(79, 25),(79, 26),(79, 27),(79, 28),(79, 29),(79, 32),(79, 33),(79, 37),(80, 19),(80, 20),(80, 21),(80, 22),(80, 25),(80, 26),(80, 28),(80, 29),(80, 31),(80, 32),(80, 33),(80, 34),(80, 35),(80, 38),(80, 42),(81, 20),(81, 21),(81, 25),(81, 26),(81, 27),(81, 32),(81, 34),(81, 35),(81, 39),(82, 16),(82, 17),(82, 19),(82, 21),(82, 24),(82, 25),(82, 30),(82, 34),(82, 44),(83, 16),(83, 17),(83, 18),(83, 21),(83, 28),(83, 30),(83, 31),(83, 34),(83, 35),(83, 37),(84, 21),(84, 29),(84, 30),(84, 31),(84, 33),(84, 35),(85, 30),(85, 32),(85, 33),(85, 34),)
    
    if len(params) != len(accomplished):
        raise ValueError('Problema dim parametri recuperati %d != dim accomplished recuperati %d'%(len(params),len(accomplished)))

    for ((ii,jj), p) in zip(accomplished, params):

        matrix[ii][jj].Recover_Initial_Parameters(p)

else:

    isolated = Get_Isolated_Elements(excluded)
    #prima riga, stimo da sx, eccetto il prim

    ii = 0
    for jj in range(n_cols):
        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        if (ii,jj) in boni:
            if (ii,jj) in isolated:
                fit = fit + ((matrix[ii][jj].Estimate_Initial_Parameters(matrix[0][0].p0.values[0], 1000), (ii,jj)),)
            else:
                fit = fit + ((matrix[ii][jj].Estimate_Initial_Parameters(matrix[ii][jj-1].p0.values[0], 1000), (ii,jj)),)

    #tutto il resto

    for ii in range(1,n_rows,1):
        for jj in range(n_cols):
            print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
            if (ii,jj) in boni:
                
                if (ii,jj) in isolated:
                    fit =   fit + ((matrix[ii][jj].Estimate_Initial_Parameters(matrix[0][0].p0.values[0], 1000), (ii,jj)),)
                else: 
                    fit =   fit + ((matrix[ii][jj].Estimate_Initial_Parameters(matrix[ii-1][jj].p0.values[0], 1000), (ii,jj)),)

tot_time    = time.process_time()-start
print('tempo impiegato per fit markoviani: %f s'%(tot_time))
print('tempo impiegato ore = %3.2f'%(tot_time/3600))

#%%

# 4) after - fit markoviano

non_fitted, accomplished, exceded = Unpack_Fit(fit)

Fit_Map     =   Get_Fit_Map(n_rows, n_cols, non_fitted, exceded, excluded, fig = 'Markov_Fit_Map', path = now_path)
Omega_Map   =   Get_Parameter_Map('Omega', matrix, n_rows, n_cols, accomplished+exceded, excluded ,fig = 'Markov_Omega_Map', path = now_path)
Gamma_Map   =   Get_Parameter_Map('Gamma', matrix, n_rows, n_cols, accomplished+exceded, excluded ,fig = 'Markov_Gamma_Map', path = now_path)

#%%
# 5)fit completo

fit_tot                    =   ()
percents                    =   (0.2, 0.1, 0.15, 'positive', 'positive', 0.15, 0.15, 0.05, 0.05, 0.05, np.inf, np.inf)
start = time.process_time()

for ii in range(n_rows):
    for jj in range(n_cols):
        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        if (ii,jj) in accomplished:
            matrix[ii][jj].Get_Fit_Bounds(percents, columns = cols)
            fit_tot =   fit_tot + (((matrix[ii][jj].Non_Linear_Least_Squares(matrix[ii][jj].p0.values[0], bound = (matrix[ii][jj].bounds['down'].values,matrix[ii][jj].bounds['up'].values), ftol = None)), (ii,jj)),)

tot_time    = time.process_time()-start
print('tempo impiegato per fit totali: %f s'%(tot_time))
print('tempo impiegato ore = %3.2f'%(tot_time/3600))

#%%
#6) after fit

non_fitted_tot, accomplished_tot, exceded_tot = Unpack_Fit(fit_tot)

Fit_Map     =   Get_Fit_Map(n_rows, n_cols, non_fitted_tot, exceded_tot, excluded, fig = 'Total_Fit_Map', path = now_path)
Omega_Map   =   Get_Parameter_Map('Omega', cols_red, matrix, n_rows, n_cols, accomplished_tot+exceded_tot, excluded ,fig = 'Omega_Map', path = now_path)
Gamma_Map   =   Get_Parameter_Map('Gamma', cols_red,colmatrix, n_rows, n_cols, accomplished_tot+exceded_tot, excluded ,fig = 'Gamma_Map', path = now_path)

Save_Params(now_path, n_rows, n_cols, matrix, accomplished_tot+exceded_tot)
Save_Fit_Info(now_path, n_rows, n_cols, fit_tot)

super_time    = time.process_time()-super_start
print('tempo impiegato per esecuzione dello script ore = %3.2f'%(tot_time/3600))
#%%
#pre accomplished cambia accomplished in boni
#VERIFICHE DA QUI IN POI
""" 
#ESCLUSIONE A MANO
to_add  = Whose_Gamma_Too_High(5., matrix, accomplished, exceded)
excluded = Escludi_a_Mano(to_add, excluded)
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
"""
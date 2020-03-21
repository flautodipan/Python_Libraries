#%%

# 1) SETTINGS and INPUTS

#libraries
import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Experimentum    import  *
from        Alessandria         import  *
import      time
import      os


#I/O 
now_path        =   '../BRILLOUIN/TDP43/NO_ARS_12_02/'
spectra_filename    =   'NO_ARS_12_02'
VIPA_filename       =   'NO_ARS_12_02_VIPA_quasisat.tif'

#operatives
syg_kwargs_test          =   {'height': 20, 'distance': 31, 'width': 2.}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}

#elementi da eliminare a mano
too_add              =   []


# %%
#2) Acquisisco dati e inizializzo oggetti Spectrum per ognuno su una matrice (n_rows, n_cols)
#   e compio alcune operazioni di sistema utili per salvataggio dati

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, now_path, var_name = 'y3')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
matrix, rows, cols = Initialize_Matrix(0,0, n_rows, n_cols)
dim     =   len(rows)*len(cols)
#%%

four = []
more_than_four = []
less_than_four = []
bizarre = []

#3) Riempio oggetti
matrix[0][0].Get_VIPA_tif(VIPA_filename, now_path, offset = 183.)

for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183., cut = False, cut_range = (200, 600))
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_test)
        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA

        if matrix[ii][jj].n_peaks == 4:
            four.append((ii,jj))
        elif matrix[ii][jj].n_peaks > 4:
            more_than_four.append((ii,jj))
        elif matrix[ii][jj].n_peaks < 4 :
            less_than_four.append((ii,jj))


# tolgo i saturati
not_saturated, saturated = Get_Saturated_Elements(matrix, len(rows), len(cols))

#varie aggiunte a mano
excluded = saturated.copy()
excluded = Escludi_a_Mano(too_add, excluded)

print("Lunghezza di spettri con 4 picchi: {}".format(len(four)))
print('Lunghezza di spettri con più di 4 picchi: {}'.format(len(more_than_four)))
print('Lunghezza di spettri con meno di 4 picchi: {}'.format(len(less_than_four)))

#%%
#VERIFICHIAMO IL VIPA
matrix[0][0].How_Many_Peaks_To_VIPA(treshold = 0, **syg_kwargs_VIPA, fig = True, verbose = True)


#%%
# Prendo per tutti 4 picchi

lowest_height = ()
print('Lunghezza iniziale four = {}'.format(len(four)))
print('Lunghezza iniziale less_than_four = {}'.format(len(less_than_four)))
print('Lunghezza iniziale more_than_four = {}'.format(len(more_than_four)))

for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

        if (ii,jj) in excluded:
            if (ii,jj) in more_than_four :  more_than_four.remove((ii,jj))
            elif (ii,jj) in less_than_four: less_than_four.remove((ii,jj))
            elif (ii,jj) in four: four.remove((ii,jj))

        
        else:
            if ((ii,jj) in more_than_four):
                matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                more_than_four.remove((ii,jj))
                four.append((ii,jj))

            elif ((ii,jj) in less_than_four):
                matrix[ii][jj].How_Many_Peaks_To(width = syg_kwargs_test['width'], distance = syg_kwargs_test['distance'])
                if hasattr(matrix[ii][jj], 'spectrum_cut_height'):
                    matrix[ii][jj].Get_Spectrum_Peaks(height =  matrix[ii][jj].spectrum_cut_height, width = syg_kwargs_test['width'], distance = syg_kwargs_test['distance'])

                if matrix[ii][jj].n_peaks == 4:
                    lowest_height += ((matrix[ii][jj].spectrum_cut_height, (ii,jj)),)
                    less_than_four.remove((ii,jj))
                    four.append((ii,jj))
                else: 
                    bizarre.append((ii,jj))
                    print('\nSpettro {} ha evidentemente meno di quattro picchi, check it\n'.format(str((ii,jj))))

print('\n\nSpettri tot = {}\nSpettri esclusi = {} \n--> Spettri disponibili = {}'.format(dim, len(excluded), (dim -len(excluded))))
print('Ho trovato {} spettri su {} con 4 picchi\n'.format(len(four),  (dim -len(excluded))))
print('Ho trovato {} spettri su {} con meno di 4 picchi\n'.format(len(less_than_four),  (dim -len(excluded))))
print('Ho trovato {} spettri su {} con più  4 picchi\n'.format(len(more_than_four),  (dim -len(excluded))))

if lowest_height:
    arg_min = np.argmin([heights[0] for heights in lowest_height])
    min_elastic_height = lowest_height[arg_min][0]
    min_elastic_name   = lowest_height[arg_min][1]
    print("Il minimo dell'altezza per trovare 4 picchi è stata {} nello spettro {}".format(min_elastic_height, min_elastic_name))

else:
    print('Tutto ok, non vi erano spettri less_than_four non esclusi')

if len(four) != (dim -len(excluded)):
    print('\nHo trovato {} spettri di cui non sono riuscito a trovare 4 picchi\n'.format(len(bizarre)))
    print(bizarre)
    

#%%
# STATSTICHE:     costruisco variabili che mi permettono di studiare le posizioni relative
#                 e le altezze dei vari picchi -> individuo gli anomali e i parametri operativi
#                 syg_kwargs
#                 il senso è che (ipotizzando posizioni elastic, brill_stockes, brill_anti, elastic)
#                 syg_kwargs['height'] = min(min(height_second), min(height_third)) 
#                 concentro presa dati su picchi brillouin
#
#                 syg_kwargs['distance'] = min(min(dist_01), min(dist_23))
#                 in quanto i picchi vicini sono elastico-stokes e antistokes-elastic
#           
#                 DA ELEMENTI STRANI IN QUESTE TUPLE INDIVIDUO CHI È CHE NON E' STATO VALUTATO 
#                 CORRETTAMENTE (mi aspetto distribuzioni normali intorno alla media)
#                 --> la statistica va centrata su quelli che già sono quattro

four.sort()


height_first   = ()
height_second  = ()
height_third   = ()
height_fourth  = ()

width_first    = ()
width_second   = ()
width_third    = ()
width_fourth   = ()
dist_01        = ()
dist_12        = ()
dist_23        = ()


for (ii,jj) in four:

        height_first  += (matrix[ii][jj].peaks['heights'][0],)
        height_second += (matrix[ii][jj].peaks['heights'][1],)
        height_third  += (matrix[ii][jj].peaks['heights'][2],)
        height_fourth += (matrix[ii][jj].peaks['heights'][3],)

        width_first  += (matrix[ii][jj].peaks['widths'][0],)
        width_second += (matrix[ii][jj].peaks['widths'][1],)
        width_third  += (matrix[ii][jj].peaks['widths'][2],)
        width_fourth += (matrix[ii][jj].peaks['widths'][3],)

        dist_01       += (np.abs(matrix[ii][jj].x[matrix[ii][jj].peaks['idx'][0]] - matrix[ii][jj].x[matrix[ii][jj].peaks['idx'][1]]), )
        dist_12       += (np.abs(matrix[ii][jj].x[matrix[ii][jj].peaks['idx'][1]] - matrix[ii][jj].x[matrix[ii][jj].peaks['idx'][2]]), )
        dist_23       += (np.abs(matrix[ii][jj].x[matrix[ii][jj].peaks['idx'][2]] - matrix[ii][jj].x[matrix[ii][jj].peaks['idx'][3]]), )




stats = (height_first, height_second, height_third, height_fourth, width_first, width_second, width_third, width_fourth, dist_01, dist_12, dist_23)


for (s, what) in zip(stats, ('height picco 1', 'Picco 2 height', 'Picco 3 height', 'Picco 4 height', 'Picco 1 width', 'Picco 2 width', 'Picco 3 width', 'Picco 4 width', 'dist_01', 'dist12', 'dist23')):
    print('\n{}\nMedia = {:.3}\nMin = {:3}\nMax = {:3}\n'.format(what, np.mean(s), np.min(s), np.max(s)))

#%%
# faccio pure i plot

for (s, what) in zip(stats, ('Picco 1 height', 'Picco 2 height', 'Picco 3 height', 'Picco 4 height', 'Picco 1 width', 'Picco 2 width', 'Picco 3 width', 'Picco 4 width','dist_01', 'dist12', 'dist23')):
    if (what == 'Picco 1 height') | (what == 'Picco 4 height'):
        log = True
    else : log = False

    plt.figure()
    plt.title('Scatter di {} vs indice tupla'.format(what))
    plt.xlabel('Idx of {} tuple'.format(what))
    plot(s, '.', color = 'darkblue', label = what)
    plt.legend()
    plt.show
    plt.figure()
    plt.title('Histogram for {}'.format(what))
    plt.hist(s, bins = 100, rwidth=0.85, alpha = 0.7, color = 'darkblue', label = what, log = log)
    plt.xlabel('{} value'.format(what))
    plt.legend()
    plt.show()


#%%
#Suggestions

how_many = 15
for tup, what in zip([height_first, height_second, height_third, height_fourth, dist_01, dist_12, dist_23], ['picco el sx', 'picco brillouin sx', 'picco brillouin dx', 'picco el dx', 'dist 01', 'dist 12','dist2']):
        
    print('\n\nPrimi {} elementi più bassi di {} :\n'.format(how_many, what), np.sort(tup)[:how_many],'\n')
    #print('{:3.2f}\t'.format([float(t) for t in tup[:how_many]]))



syg_kwargs_height = np.min([np.min(height_third), np.min(height_second)])
print('Ti suggerisco di usare come height di syg_kwargs quella minore del brillouin minore: {}\n'.format(syg_kwargs_height))

syg_kwargs_dist = np.min([np.min(dist_23), np.min(dist_01)])
print('Ti suggerisco di usare come dist di syg_kwargs quella minore tra brillouin ed elastico vicino: {}\n'.format(syg_kwargs_dist))

syg_kwargs_brill_height = np.min([np.min(height_first), np.min(height_fourth)])
print('Ti suggerisco di usare come height di syg_kwargs_brill quella minore tra i due picchi elastici: {}\n'.format(syg_kwargs_brill_height))

syg_kwargs_width = np.min([np.min(width_first), np.min(width_fourth)])
print('Ti suggerisco di usare come width di syg_kwargs, quella minore tra i due picchi elastici: {}\n'.format(syg_kwargs_width))

print('Ti suggerisco di usare come medie distanze da elastici per taglio sulle x_freq:\n {:3.2} per 01\n{:3.2} per 23\n'.format(np.mean(dist_01), np.mean(dist_23)))




#%%

"""
Per trovare il (ii,jj) di un elemento di queste tuple che (per es) dai plot risulta asim
metrico, uas np.where() TRASFORMANDO TUPLA IN ARRAY PERÒ
--> esce degli indici
es. 4288
poi si fa un ciclo sui more_than_four per avere quale (ii,jj) 

count = 0
for (ii,jj) in four:
    if (count == 3989):
        print(str((ii,jj)))
        
    count+=1

--> poi passi a exp_single
"""
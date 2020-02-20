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

now_path        =   '../BRILLOUIN/TDP43/ARS_13_02/'
spectra_filename    =   'ARS_13_02'
VIPA_filename       =   'NO_ARS_13_02_VIPA_not_sat.tif'

#variables
invisible           =   [] 
brillouin_higher    =   []
brillouin_highest   =   []
boni                =   []
excluded            =   []

cols        = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position',  'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_width', 'delta_position', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')
# %%
#2) Acquisisco dati e inizializzo oggetti Spectrum per ognuno su una matrice (n_rows, n_cols)
#   e compio alcune operazioni di sistema utili per salvataggio dati

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, now_path, var_name = 'y3')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
#matrix, rows, cols = Initialize_Matrix(1,8, 2+1, 10)
matrix, rows, cols = Initialize_Matrix(0,0, n_rows, n_cols)
dim     =   len(rows)*len(cols)
#%%

#3) Riempio oggetti
#prendo i 4 più alti
matrix[0][0].Get_VIPA_tif(VIPA_filename, now_path, offset = 183.)
syg_kwargs_test          =   {'height': 20, 'distance': 31, 'width': 3.}

for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183., cut = False, cut_range = (200, 600))
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_test)
        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA

not_saturated, saturated = Get_Saturated_Elements(matrix, len(rows), len(cols))

#%%
#VERIFICHIAMO IL VIPA
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
matrix[0][0].How_Many_Peaks_To_VIPA(treshold = 0, **syg_kwargs_VIPA, fig = True, verbose = True)


#%%

#varie aggiunte a mano
too_add = [(66,3),]
excluded = saturated.copy()
excluded = Escludi_a_Mano(too_add, excluded)

#%%
#4) statistiche : costruisco variabili che mi permettono di studiare le posizioni relative
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

syg_kwargs_test          =   {'height': 20, 'distance': 31, 'width': 3.}

while True:

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


    four           = ()
    three          = ()
    two            = ()
    less_than_four = ()
    more_than_four = ()


    for ii in range(len(rows)):
        for jj in range(len(cols)):
            if (ii,jj) not in excluded:

                if matrix[ii][jj].n_peaks >=4:
                    if matrix[ii][jj].n_peaks == 4:
                        four += ((ii,jj), )
                    more_than_four += ((ii,jj),)
                    matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                elif matrix[ii][jj].n_peaks  == 3:
                    three += ((ii,jj), )
                elif matrix[ii][jj].n_peaks  == 2:
                    two += ((ii,jj), )
                else:
                    raise ValueError('BOH?')

    if len(more_than_four) == (dim-len(excluded)):
        break

    elif syg_kwargs_test['height'] <= 0:
        print('Superato lo zero, qualcosa non va in \n {} spettri con due picchi\n'.format(len(two)), two)
        print('Superato lo zero, qualcosa non va in \n {} spettri con tre picchi\n'.format(len(three)), three)
        
        break
    else:
        syg_kwargs_test['height'] -= 1
        for (ii,jj) in two: matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_test)
        for (ii,jj) in three: matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_test)

print('Spettri tot = {}\nSpettri esclusi = {} \n--> Spettri disponibili = {}\nHo trovato {} spettri su {} con 4 picchi\n'.format(dim, len(excluded), (dim -len(excluded)),len(four),  (dim -len(excluded))))
print(syg_kwargs_test['height'])
    
#%%

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


# %%

#I/O 

now_path        =   '../BRILLOUIN/TDP43/ARS_13_02/'
spectra_filename    =   'ARS_13_02'
VIPA_filename       =   'NO_ARS_13_02_VIPA_not_sat.tif'

#operatives

syg_kwargs          =   {'height': 80, 'distance': 31, 'width': 3.}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
syg_kwargs_brill    =  {'height': 18, 'distance': 31, 'width': 3.}
VIPA_treshold       =   6
sat_height          =   50000
sat_width           =   13.5

"""
Per trovare il (ii,jj) di un elemento di queste tuple che (per es) dai plot risulta asim
metrico, uas np.where() TRASFORMANDO TUPLA IN ARRAY PERÒ
--> esce degli indici
es. 4288
poi si fa un ciclo sui more_than_four per avere quale (ii,jj) 

count = 0
for (ii,jj) in more_than_four:
    if (count == 3989):
        print(str((ii,jj)))
        
    count+=1

--> poi passi a exp_single
"""
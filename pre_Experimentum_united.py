#%%

#####
###### 1) SETTINGS and INPUTS
#####

### EXTERNAL PACKAGES

import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Experimentum    import  *
from        Alessandria         import  *
import      time
import      os
from        os.path             import  join
import      configparser

# nomi dei parametri modelli viscoelastici
cols_mark           = ('Co', 'Omega', 'Gamma', 'delta_position', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark_nodelta   = ('Co', 'Omega', 'Gamma', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real           = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_real_nodelta   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'shift', 'offset')

### MANUAL INPUTS

#I/O 
spectra_filename    =   'ARS_10_02'
now_path            =   join('..', 'BRILLOUIN', 'TDP43', spectra_filename)
VIPA_filename       =   'ARS_10_02_VIPA_notsat'

#to modify
syg_kwargs_test     =   {'height': 10, 'distance': 32, 'width': 2.1}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
syg_kwargs_brill    =   {'height': 18, 'distance': 31, 'width': 2.1}

# OPERATIVES FOR ACQUISITION
transpose           =   False
data_offset         =   183.
almost_treshold     =   150000
saturation_height   =   50000
saturation_width    =   15.
pre_cut             =   True
pre_cut_range       =   (30,210) if pre_cut else None
initial             =   'right'
to_add              =   []

#OPERATIVES FOR FIT
fit_algorithm       =   'trf'
cut                 =   True
exclude_delta       =   True
recover_markov      =   False
rules_markov_bounds =   {'Co' : 'positive', 'Omega' : 0.2, 'Gamma' :'positive', 'delta_position' : [-2,2], 'delta_width': 'positive', 'delta_amplitude' : 'positive', 'A' : 'positive', 'mu' : 0.01, 'sigma' : 0.001,  'shift' : 'inf', 'offset' :'inf'}
skip_tot            =   True
rules_tot_bounds    =   {'Co' : 0.2, 'Omega' : 0.01, 'Gamma' : 0.01, 'Delta': 'positive', 'tau' : 'positive', 'delta_position': [-2,2], 'delta_width' : 0.01, 'delta_amplitude' : 0.01, 'shift':'inf', 'offset' : 'inf'}


# %%

#####
######  2)  DATA ACQUISITION AND OBJECTS INITIALIZATION
#####

dati                =   Import_from_Matlab(spectra_filename, now_path, transpose = transpose)
n_rows              =   len(dati)
n_cols              =   len(dati[0])
matrix, rows, cols  =   Initialize_Matrix(0,0, n_rows, n_cols)
dim                 =   len(rows)*len(cols)
#%%

#####
######  3A)  DATA CLASSIFICATION: FASE I
#####

four            = []
more_than_four  = []
less_than_four  = []
bizarre         = []

# Riempio oggetti e classifico in chi ha 4 picchi, chi ne ha di più, chi di meno

for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('row = {}/{} col = {}/{}'.format(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = data_offset, cut = pre_cut, cut_range = pre_cut_range)
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_test)


        if matrix[ii][jj].n_peaks == 4:
            four.append((ii,jj))
        elif matrix[ii][jj].n_peaks > 4:
            more_than_four.append((ii,jj))
        elif matrix[ii][jj].n_peaks < 4 :
            less_than_four.append((ii,jj))


# classifico i saturati, in base all'input saturation_height
not_saturated, saturated = Get_Saturated_Elements(matrix, len(rows), len(cols), saturation_height = saturation_height, saturation_width=saturation_width)

# formo la categoria degli esclusi
excluded = (saturated + to_add)

print("\n\nNumero di spettri con 4 picchi:\n{}".format(len(four)))
print('Numero di spettri con più di 4 picchi:\n{}'.format(len(more_than_four)))
print('Numero di spettri con meno di 4 picchi:\n{}'.format(len(less_than_four)))
print('Numero di spettri non saturati:\n{}/{}'.format(len(not_saturated), dim))
print('Numero di spettri saturati:\n{}/{}'.format(len(saturated), dim))
print('Numero di spettri che saranno esclusi da analisi:\n{}/{}'.format(len(excluded), dim))


#%%
#####
######  3B)  DATA CLASSIFICATION: FASE II
#####        tutti con 4 picchi

lowest_height = ()

for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('row = {}/{} col = {}/{}'.format(ii,len(rows)-1, jj, len(cols)-1))

        if (ii,jj) in excluded:
            if (ii,jj) in more_than_four :  more_than_four.remove((ii,jj))
            elif (ii,jj) in less_than_four: less_than_four.remove((ii,jj))
            elif (ii,jj) in four: four.remove((ii,jj))
        
        else:
            if ((ii,jj) in more_than_four):
                if matrix[ii][jj].n_peaks == 5:
                    matrix[ii][jj].Exclude_Glass_Peak()
                else:
                    matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                more_than_four.remove((ii,jj))
                four.append((ii,jj))

            elif ((ii,jj) in less_than_four):
                print('Per lo spettro {} sto usando How_Many_Peaks_To\n'.format(str((ii,jj))))
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

Plot_Elements_Spectrum(matrix, bizarre, peaks = True, pix = True)
to_add += bizarre
#%%
#####
###### 4) STATS:  costruisco variabili che mi permettono di studiare le posizioni relative
#####             e le altezze dei vari picchi -> individuo gli anomali e i parametri operativi
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

#salvo medie distanze elastico-brillouin
mean_dist_01 = np.mean(dist_01)
mean_dist_23 = np.mean(dist_23)

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
# stampo i syg_kwargs consigliati risultanti delle analisi statistiche

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

syg_kwargs_width = np.min([np.min(width_first), np.min(width_fourth), np.min(width_second), np.min(width_third)])

print('Ti suggerisco di usare come width di syg_kwargs, quella minore tra i due picchi elastici: {}\n'.format(syg_kwargs_width))

print('Ti suggerisco di usare come medie distanze da elastici per taglio sulle x_freq:\n {:3.2} per 01\n{:3.2} per 23\n'.format(np.mean(dist_01), np.mean(dist_23)))


#%%
#####
######  5)  FIND PROBLEMS (solo se necessario)
#####       

where_is = Find_Problems(four, which = dist_23, what_is = 'greater', value = 36)

Plot_Elements_Spectrum(matrix, where_is[:10], pix = True, peaks = True, x_range= (), y_range=())
#%%
#####
######  6)  FIRST NORMAL and FIRST ALMOST and VIPA check
#####       fit per initial conditions 

for ii,jj in serpentine_range(len(rows), len(cols), initial):

    if matrix[ii][jj].y.max() > 15000:
        first_almost = (ii,jj)
        break

first_normal = serpentine_range(len(rows), len(cols), initial)[0]

print(" Primo normale è {}, Primo almost è {}\nOra faccio i fit per stimare le migliori condizioni iniziali\n".format(first_normal, first_almost))

p0s = []
for (ii,jj) in [first_normal, first_almost]:
    print('Faccio il fit per {}\n'.format(str((ii,jj))))

    matrix[ii][jj].Get_VIPA_tif(VIPA_filename, now_path, offset = data_offset)
    matrix[ii][jj].How_Many_Peaks_To_VIPA(treshold = 6, **syg_kwargs_VIPA, fig = True, verbose = False)
    matrix[ii][jj].Fit_Pixel2GHz(fig = True if (ii,jj) == first_normal else False , data_color = 'maroon', fit_color = 'yellowgreen', )
    matrix[ii][jj].VIPA_Pix2GHz(fig = True if (ii,jj) == first_normal else False )
    matrix[ii][jj].Poly2GHz      =   matrix[ii][jj].Poly2GHz
    matrix[ii][jj].Align_Spectrum(alignment = matrix[ii][jj].alignment)
    matrix[ii][jj].Spectrum_Pix2GHz()
    matrix[ii][jj].Cut_n_Estimate_Spectrum(estimate = True, cut = cut, mean_dist01 = mean_dist_01, mean_dist23 = mean_dist_23)
    matrix[ii][jj].Fit_VIPA_Gaussian(fig = True if (ii,jj) == first_normal else False )
    columns = cols_mark
    rules_bounds = rules_markov_bounds

    matrix[ii][jj].Get_VIPA_for_fit('interpolate', interpolation_density = 500)
    matrix[ii][jj].Get_Best_p0([matrix[ii][jj].p0[list(columns)].values[0]], columns)
    matrix[ii][jj].Get_Fit_Bounds(rules_bounds, columns)
    if fit_algorithm == 'trf':
        matrix[ii][jj].Non_Linear_Least_Squares_Markov(columns, bounds = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values),  max_nfev = 100, ) 
    else: matrix[ii][jj].Non_Linear_Least_Squares_Markov(columns, method = 'lm',  max_nfev = 100, )
    p0s.append(list(matrix[ii][jj].Markov_Fit_Params.T.Values.values))
    print('Ok ho stimato i parametri iniziali per {}\n'.format(str((ii,jj)),))

    for col, val in zip(columns, matrix[ii][jj].Markov_Fit_Params.T.Values.values):
        print('{} = {:3.2f}'.format(col, val))

p0_normal = p0s[0]
p0_almost = p0s[1]
#%%
#####
######  7)  GENERO FILE .ini
#####       

config = configparser.ConfigParser()

config['I/O'] = {'spectra_filename' : spectra_filename, 'VIPA_filename' : VIPA_filename}
config['syg_kwargs'] = { 'height' : Get_Around(syg_kwargs_height, 0.01)[0], 'width' : Get_Around(syg_kwargs_width, 0.01)[0], 'distance' : Get_Around(syg_kwargs_dist, 0.01)[0]}
config['syg_kwargs_brill'] = {'height' : Get_Around(syg_kwargs_brill_height, 0.01)[0], 'width' : Get_Around(syg_kwargs_width, 0.01)[0], 'distance' : Get_Around(syg_kwargs_dist, 0.01)[0]}
config['syg_kwargs_VIPA'] = {'width' : Get_Around(syg_kwargs_width, 0.01)[0], 'distance' : Get_Around(syg_kwargs_dist, 0.01)[0]}
config['Operatives'] = {'transpose' : transpose, 'data_offset' : data_offset, 'almost_treshold':almost_treshold, 'sat_height': saturation_height, 'sat_width':saturation_width, 'pre_cut' : pre_cut, 'pre_cut_range' : pre_cut_range, 'initial': initial, 'to_add' : to_add,  'fit_algorithm' : fit_algorithm, 'cut': cut,  'exclude_delta' : exclude_delta,'mean_dist_01' : mean_dist_01, 'mean_dist_23' : mean_dist_23, }
config['Markov'] = {'recover_markov': recover_markov, 'first_normal' : first_normal, 'p0_normal' : p0_normal, 
'first_almost': first_almost, 'p0_almost' :p0_almost, 'rules_markov_bounds':  rules_markov_bounds }
config['Tot'] = {'skip_tot' : skip_tot, 'rules_tot_bounds' : rules_tot_bounds}

with open(join(now_path,'config.ini'), 'w') as fin:
    config.write(fin)

#%%

"""
Per trovare il (ii,jj) di un elemento di queste tuple che (per es) dai plot risulta asim
metrico, uas np.where() TRASFORMANDO TUPLA IN ARRAY PERÒ
--> esce degli indici
es. 4288
poi si fa un ciclo sui more_than_four per avere quale (ii,jj) 

count = 0
for (ii,jj) in four:
    if (count == 8557):
        print(str((ii,jj)))
        
    count+=1

--> poi passi a exp_single
"""

# %%


# %%

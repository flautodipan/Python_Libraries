
#########################                                                   #######################
#########################           EXPERIMENTUM UNITED                     #######################

#%%
######################################################################################################################
#######   ||||    PREAMBOLO: - acquisiamo i file di INPUT
#######    ||                - eseguiamo operazioni su cartelle
#######    ||   
#######    ||    
#######   ||||


#libraries
import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Experimentum    import  *
from        Alessandria         import  *
import      time
import      os

#I/O 

now_path            =   '../BRILLOUIN/TDP43/ARS_13_02/'
spectra_filename    =   'ARS_13_02'
VIPA_filename       =   'NO_ARS_13_02_VIPA_quasisat.tif'
log_file            =   'log_'+spectra_filename
analysis_dir       =   'analysis_best/'

#operatives

#esclusi a mano
to_add              =   [(66, 3),]

syg_kwargs          =   {'height': 80, 'distance': 31, 'width': 3.}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
syg_kwargs_brill    =  {'height': 18, 'distance': 31, 'width': 3.}
VIPA_treshold       =   6
sat_height          =   50000
sat_width           =   13.5
#quanto mi allontano dal VIPA
cut_distance        =   0.25
#markov_fit
recover_markov = False
percents_markov     =   ('positive', 0.2, 'positive', np.inf, 'positive', 'positive', 0.2, 0.2, 0.2,  np.inf, np.inf)
#tot fit
percents_tot        = (0.2, 0.1, 0.5, 'positive', 'positive', 0.15,  0.15, 0.15, np.inf, np.inf)

#variables

invisible           =   []
brillouin_higher    =   []
brillouin_highest   =   []
boni                =   []
excluded            =   []

cols        = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position',  'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_position','delta_width',  'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')
# %%
#2) ######################################################################################################################

#######   ||||  ||||    DATA ACQUISITION AND TREATMENT : - acquisiamo i file di INPUT
#######    ||    ||                                  - eseguiamo operazioni su cartelle
#######    ||    ||
#######    ||    ||
#######   ||||  ||||

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, now_path, var_name = 'y3')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
#matrix, rows, cols = Initialize_Matrix(0,0, 1+1, 1+1)
matrix, rows, cols = Initialize_Matrix(0,0, n_rows, n_cols)
dim     =   len(rows)*len(cols)


os.system('cd ' + now_path +' & mkdir ' + now_path+analysis_dir)
analysis_path            =   now_path +analysis_dir

with open(analysis_path+log_file, 'w') as f_log:
    f_log.write('#This is a log file: you will find info on script run for {}\n'.format(spectra_filename))
    f_log.write('\n\nHo inizializzato una matrice {}x{}, per un totale di {} spettri'.format(len(matrix), len(matrix[0]), len(matrix)*len(matrix[0])))
    
# %%
### riempio oggetti spectrum

tempo = ()
super_start = time.process_time()
start = super_start

matrix[0][0].Get_VIPA_tif(VIPA_filename, now_path, offset = 183.)

for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183., cut = False, cut_range = (200, 600))
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs)

        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA

### catalogo la natura degli spettri

not_saturated, saturated = Get_Saturated_Elements(matrix, len(rows), len(cols), saturation_height = sat_height, saturation_width = sat_width)
excluded        = saturated.copy()
excluded        = Escludi_a_Mano(to_add, excluded)

for ii in range(len(rows)):
    for jj in range(len(cols)):
            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
            if (ii,jj) not in excluded:
                #Spettri normali , quattro picchi di cui i picchi più alti sono elastici
                
                if matrix[ii][jj].n_peaks >= 4:
                    if matrix[ii][jj].Check_Brillouin_Distances(average = 70, stdev = 70/10):
                        invisible += [(ii,jj), ]
                    else: boni += [(ii,jj),]
                elif (matrix[ii][jj].n_peaks == 2):
                    matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_brill)
                    brillouin_highest += [(ii,jj), ]
                    boni += [(ii,jj),]
                elif  (matrix[ii][jj].n_peaks == 3):
                    matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_brill)
                    brillouin_higher += [(ii,jj), ]
                    boni += [(ii,jj),]
                else:
                    raise ValueError("Numero di picchi non previsti dal codice per spettro {}".format(str((ii,jj))))

acq_time    =   time.process_time()-start
tempo       =   tempo + (('acquisizione', acq_time),)
print('\nTempo impiegato per acquisizione spettri: {} s\n'.format(acq_time))

with open(analysis_path+log_file, 'a') as f_log:

    f_log.write('\n\n###############################INFO ON DATA ACQUISITION#######################################\n\n')
    f_log.write('\nTotale spettri saturati : {}\n'.format(len(saturated)))
    f_log.write(str(saturated))
    f_log.write('\nTotale spettri con Brillouin più alti : {}\n'.format(len(brillouin_higher)))
    f_log.write(str(brillouin_higher))
    f_log.write('\nTotale spettri con Brillouin invisibili: {}\n'.format(len(invisible)))
    f_log.write(str(invisible))
    f_log.write('\nTotale spettri boni :  {}'.format(len(boni)))
    f_log.write('\n\nTotale spettri : {}\ndi cui {} inutilizzabili, '.format(dim, len(invisible)+len(saturated)))
    f_log.write('ossia il %3.2f percento\n'%(float(len(invisible)+len(saturated))*100/dim))
    f_log.write('Di cui il {} invisibili e il {} saturati\n\n'.format(len(invisible)*100/dim, len(saturated)*100/dim))

# %%
#2) Faccio operazioni di modifica spettro

start = time.process_time()

matrix[0][0].How_Many_Peaks_To_VIPA(treshold = VIPA_treshold, **syg_kwargs_VIPA)
matrix[0][0].Fit_Pixel2GHz()
matrix[0][0].VIPA_Pix2GHz()
matrix[0][0].Spectrum_Pix2GHz()
matrix[0][0].Align_Spectrum()
matrix[0][0].Cut_n_Estimate_Spectrum(estimate = True, distanza = cut_distance)
matrix[0][0].Fit_VIPA_Gaussian()

for ii in range(len(rows)):
    for jj in range(len(cols)):
            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
            if ((ii,jj) not in excluded) & ((ii,jj) != (0,0)):
                matrix[ii][jj].x_VIPA_freq   =   matrix[0][0].x_VIPA_freq
                matrix[ii][jj].y_VIPA        =   matrix[0][0].y_VIPA
                matrix[ii][jj].Poly2GHz      =   matrix[0][0].Poly2GHz
                matrix[ii][jj].Spectrum_Pix2GHz()
                matrix[ii][jj].Align_Spectrum(alignment = matrix[0][0].alignment)
                matrix[ii][jj].Cut_n_Estimate_Spectrum(distanza = cut_distance)
            elif ((ii,jj) in excluded):
                matrix[ii][jj].Poly2GHz      =   matrix[0][0].Poly2GHz
                matrix[ii][jj].Spectrum_Pix2GHz()

mod_time    =   time.process_time()-start
tempo       =   tempo + (('modifica',mod_time), )
print('tempo impiegato per modifica spettri: %f s'%(mod_time))


# salvo info spettri e VIPA
Save_XY_position(matrix, len(rows), len(cols), path = analysis_path)
Save_XY_VIPA(matrix[0][0].x_VIPA_freq, matrix[0][0].y_VIPA, path = analysis_path)
print('\n I saved xy info on xy.txt and xy_VIPA.txt in your analysis directory {}\n\n'.format(analysis_path))

#%%
#######################################################################################################################

#######   ||||  ||||  ||||   MARKOVIAN FIT: - opero fit markoviano con tutti i parametri
#######    ||    ||    ||                     tranne Delta e tau (quindi anche Gauss)
#######    ||    ||    ||
#######    ||    ||    ||
#######   ||||  ||||  ||||


if recover_markov == False:
        
    print('\n\n You chose to do the markovian fit\n\n')
    fit                 =   ()
    start = time.process_time()
    isolated = Get_Isolated_Elements(excluded)
    
    for (ii,jj) in boni:

        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

        if (ii,jj) in isolated:

            matrix[ii][jj].Get_p0(matrix[0][0].p0.values[0], cols_mark)
        
        else:

            if ii == 0:
                matrix[ii][jj].Get_p0(matrix[ii][jj-1].p0.values[0], cols_mark)
            else:
                matrix[ii][jj].Get_p0(matrix[ii-1][jj].p0.values[0], cols_mark)

        matrix[ii][jj].Get_cost_markov(matrix[ii][jj].p0.values[0])
        print('Cost before fitting = {}'.format(matrix[ii][jj].cost_markov))
        matrix[ii][jj].Get_Fit_Bounds(percents_markov, cols_mark)
        fit = fit + ((matrix[ii][jj].Non_Linear_Least_Squares_Markov(bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values),  max_nfev = 100),(ii,jj)),)
        matrix[ii][jj].Get_cost_markov(matrix[ii][jj].Markov_Fit_Params.values[0])
        print('Cost after fitting = {}\n'.format(matrix[ii][jj].cost_markov))

        del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].y_markov_convolution
        #print(fit)
else:

    print('\n\n You chose to SKIP the markovian fit and recover info \n\n')

    ############## faccio er ricovery

    with open(analysis_path+'markov_fit.txt', 'r') as fin:
        fit     =   eval(fin.read())

    non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

    with open(analysis_path+'markov_fit_params.txt', 'r') as fin:
        lines   =   fin.readlines()

    if (len(fitted) != len(lines)):
        raise ValueError("Incompatibilità tra lunghezza file parametri ({}) e informazioni fit ({})".format(len(fitted), len(lines)))

    for (line, (ii,jj)) in zip(lines, fitted) :
        matrix[ii][jj].Recover_Markov_Fit_Params(line)

markov_time     =   time.process_time()-start
tempo           =   tempo + (('fit markoviano', markov_time),)

print('tempo impiegato per fit markoviani: %f s'%(markov_time))
print('tempo impiegato ore = %3.2f'%(markov_time/3600))


# 4) after - fit markoviano

non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

#too_markov         =   Whose_Param_Too_High('Gamma', 2., matrix, fitted)
Save_Fit_Info(fit, filename = 'markov_fit.txt', path=analysis_path)
Save_Markov_Fit_Parameters(matrix, fitted, out_filename = 'markov_fit_params.txt', path = analysis_path)
Save_y_markov_fit(matrix, boni, path = analysis_path)
Save_cost_markov(matrix, boni, path = analysis_path)

#%%
##################################################################################################################################
# fit tot

fit_tot = ()

start = time.process_time()
for (ii,jj) in boni:

    print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
    p_gauss = matrix[ii][jj].Markov_Fit_Params[list(cols_gauss)].values[0]
    matrix[ii][jj].Initials_Parameters_from_Markov(matrix[ii][jj].Markov_Fit_Params.T['Values'].values)
    matrix[ii][jj].Get_Fit_Bounds(percents_tot, columns = cols_real)
    matrix[ii][jj].Get_cost_tot(matrix[ii][jj].p0.values[0], p_gauss)
    print('\nCost before fitting = {}\n'.format(matrix[ii][jj].cost_tot))
    fit_tot =   fit_tot + (((matrix[ii][jj].Non_Linear_Least_Squares(p_gauss, cols_real, bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values), max_nfev = 35)), (ii,jj)),)
    matrix[ii][jj].Get_cost_tot(matrix[ii][jj].Tot_Fit_Params.values[0], p_gauss)
    print('\nCost after fitting = {}\n'.format(matrix[ii][jj].cost_tot))
    #del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].res_lsq, matrix[ii][jj].bounds



tot_time     =   time.process_time()-start
tempo           =   tempo + (('fit tot', tot_time),)

#after fit
non_fitted_tot, accomplished_tot, exceded_tot, fitted_tot = Unpack_Fit(fit_tot)

Save_Fit_Info(fit_tot, filename = 'tot_fit.txt', path=analysis_path)
Save_Tot_Fit_Parameters(matrix, fitted_tot, out_filename = 'tot_fit_params.txt', path = analysis_path)
Save_y_fit(matrix, boni, path = analysis_path)
Save_cost_tot(matrix, boni, path = analysis_path)


super_time = super_start - time.process_time()
   
print('tempo impiegato per esecuzione dello script ore = %3.2f\n '%(super_time/3600))

for (what,t) in tempo:

    print('di cui %f secondi =  %f  ore in %s \n' %(t, t/3600, what))


# %%

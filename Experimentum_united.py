
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


############

#I/O 
now_path        =   '../BRILLOUIN/TDP43/NO_ARS_12_02/'
spectra_filename    =   'NO_ARS_12_02'
VIPA_filename       =   'NO_ARS_12_02_VIPA_quasisat.tif'
log_file            =   'log_'+spectra_filename
analysis_dir        =   'dabuttare/'

#operatives

#esclusi a mano
to_add              =   []

syg_kwargs          =   {'height': 119, 'distance': 31, 'width': 3.}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
syg_kwargs_brill    =  {'height': 23, 'distance': 31, 'width': 3.}
VIPA_treshold       =   6
sat_height          =   50000
sat_width           =   13.5
almost_treshold     =   15000

#quanto mi allontano dal VIPA
pre_cut             =   False
cut                 =   True

mean_dist_01 = 37
mean_dist_23 = 34
#markov_fit
p0_normal = np.array([ 1.07378474e-01,  7.57148558e+00,  1.49128813e-01,  1.19015861e-01,
        1.448930518e-01,  8.34614271,  4.79747192e+03, -1.00904973e+01,
        1.58007162e+01,  2.11019859e-01, -3.10388495e-01])
p0_brillouin = np.array([ 1.07378474e-01,  7.57148558e+00,  1.49128813e-01,  1.19015861e-01,
        1.48930518e-01,  2.34614271e-01,  4.79747192e+03, -1.00904973e+01,
        1.58007162e+01,  2.11019859e-01, -3.10388495e-01])
p0_almost = np.array([ 1.07186924e-01,  7.63051819e+00,  1.33280055e-01,  1.97510814e+00,
        5.09986043e-01,  1.66616101e+00,  4.33362727e+03, -1.00496864e+01,
        1.59365161e+01,  2.77695117e-01,  6.43211621e+00])

recover_markov = False
rules_markov_bounds     =   ('positive', 0.2, 'positive', [-2,2] , 'positive', 'positive', 0.2, 0.01, 0.001,  'inf', [-2,2])
#tot fit
skip_tot = False
rules_tot_bounds                   =   (0.2, 0.01, 0.01, 'positive', 'positive', [-2,2], 0.01, 0.01, 'inf', 0.5)
############

#variables

invisible           =   []
brillouin_higher    =   []
brillouin_highest   =   []
boni                =   []
excluded            =   []
almost_height       =   []
normals             =   []

cols_smart  =  ('Co', 'Omega', 'Gamma', 'delta_position',  'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
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
#matrix, rows, cols = Initialize_Matrix(0,0,3,3)
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
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183., cut = pre_cut, cut_range = (200, 600))
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
                    matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                    if matrix[ii][jj].Check_Brillouin_Distances(average = 70, stdev = 70/10):
                        invisible += [(ii,jj), ]
                    else: 

                        if matrix[ii][jj].y.max() > almost_treshold:
                            almost_height += [(ii,jj),]
                        else:
                            normals += [(ii,jj),]

                        boni += [(ii,jj),]
                elif (matrix[ii][jj].n_peaks == 2):
                    matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_brill)
                    brillouin_highest += [(ii,jj), ]
                    boni += [(ii,jj),]
                elif  (matrix[ii][jj].n_peaks == 3):
                    matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_brill)
                    if (matrix[ii][jj].y.argmax() == matrix[ii][jj].peaks['idx'][1]) | (matrix[ii][jj].y.argmax() == matrix[ii][jj].peaks['idx'][2]):
                        brillouin_highest += [(ii,jj), ]
                    else:
                        brillouin_higher += [(ii,jj), ]
                    boni += [(ii,jj),]
                else:
                    raise ValueError("Numero di picchi non previsti ({} )dal codice per spettro {}".format(matrix[ii][jj].n_peaks, str((ii,jj))))

excluded    += invisible
excluded.sort()

acq_time    =   time.process_time()-start
tempo       =   tempo + (('acquisizione', acq_time),)
print('\nTempo impiegato per acquisizione spettri: {} s\n'.format(acq_time))

with open(analysis_path+log_file, 'a') as f_log:

    f_log.write('\n\n###############################INFO ON DATA ACQUISITION#######################################\n\n')
    f_log.write('\nTotale spettri saturati : {}\n'.format(len(saturated)))
    f_log.write(str(saturated))
    f_log.write('\nTotale spettri con Brillouin più alti di un elastico : {}\n'.format(len(brillouin_higher)))
    f_log.write(str(brillouin_higher))
    f_log.write('\nTotale spettri con Brillouin più alti in assoluto : {}\n'.format(len(brillouin_highest)))
    f_log.write(str(brillouin_highest))
    f_log.write('\nTotale spettri con Brillouin invisibili (= aggiunti dopo ulteriore controllo): {}\n'.format(len(invisible)))
    f_log.write(str(invisible))
    f_log.write('\nTotale spettri boni :  {} di cui\n'.format(len(boni)))
    f_log.write('\nTotale spettri con elastici più alti di {} : {}\n'.format(almost_treshold, len(almost_height)))
    f_log.write('\nTotale spettri normali : {}\n'.format(len(normals)))
    f_log.write('\n\nTotale spettri : {}\ndi cui {} inutilizzabili, '.format(dim, len(invisible)+len(saturated)))
    f_log.write('ossia il %3.2f percento\n'%(float(len(invisible)+len(saturated))*100/dim))
    f_log.write('Di cui il {} invisibili e il {} saturati\n\n'.format(len(invisible)*100/dim, len(saturated)*100/dim))

# %%
#2) Faccio operazioni di modifica spettro

start = time.process_time()

matrix[0][0].How_Many_Peaks_To_VIPA(treshold = VIPA_treshold, **syg_kwargs_VIPA)
matrix[0][0].Fit_Pixel2GHz()
matrix[0][0].VIPA_Pix2GHz()
matrix[0][0].Align_Spectrum()
matrix[0][0].Spectrum_Pix2GHz()
matrix[0][0].Cut_n_Estimate_Spectrum(estimate = True, cut = cut, mean_dist01 = mean_dist_01, mean_dist23 = mean_dist_23)
matrix[0][0].Fit_VIPA_Gaussian()

for ii in range(len(rows)):
    for jj in range(len(cols)):
            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
            if ((ii,jj) not in excluded) & ((ii,jj) != (0,0)):
                matrix[ii][jj].x_VIPA_freq   =   matrix[0][0].x_VIPA_freq
                matrix[ii][jj].y_VIPA        =   matrix[0][0].y_VIPA
                matrix[ii][jj].Poly2GHz      =   matrix[0][0].Poly2GHz
                matrix[ii][jj].Align_Spectrum(alignment = matrix[0][0].alignment)
                matrix[ii][jj].Spectrum_Pix2GHz()
                matrix[ii][jj].Cut_n_Estimate_Spectrum(cut = cut, mean_dist01 = mean_dist_01, mean_dist23 = mean_dist_23)
            elif ((ii,jj) in excluded):
                matrix[ii][jj].Poly2GHz      =   matrix[0][0].Poly2GHz
                matrix[ii][jj].Spectrum_Pix2GHz()

mod_time    =   time.process_time()-start
tempo       =   tempo + (('modifica',mod_time), )

#%%
with open(analysis_path+log_file, 'a') as f_log:
    f_log.write('\n\n######### MODIFICA SPETTRI #####################\n\n')
    f_log.write('\nTempo impiegato per modifica spettri: {} s\nTaglio spettri è {}\n\n'.format(mod_time, cut))

# salvo info spettri e VIPA
Save_XY_position(matrix, len(rows), len(cols), path = analysis_path)
Save_XY_VIPA(matrix[0][0].x_VIPA_freq, matrix[0][0].y_VIPA, path = analysis_path)

with open(analysis_path+log_file, 'a') as f_log:
    f_log.write('\n\nI saved xy info on xy.txt and xy_VIPA.txt in your analysis directory {}\n\n'.format(analysis_path))

#%%
#######################################################################################################################

#######   ||||  ||||  ||||   MARKOVIAN FIT: - opero fit markoviano con tutti i parametri
#######    ||    ||    ||                     tranne Delta e tau (quindi anche Gauss)
#######    ||    ||    ||
#######    ||    ||    ||
#######   ||||  ||||  ||||

print("\nI'm beginning markovian fit\n")

if recover_markov == False:
        
    print('\n\n You chose to do the markovian fit\n\n')
    fit                 =   ()
    start = time.process_time()
    isolated = Get_Isolated_Elements(excluded)
    
    # je faccio fa un primo giro perchè no, così lo controllo e miglioro la mia stima di p0
    """
    matrix[0][0].Get_Fit_Bounds(percents_markov, cols_mark)
    _ =  matrix[0][0].Non_Linear_Least_Squares_Markov(cols_mark, bound = (matrix[0][0].bounds['down'].values, matrix[0][0].bounds['up'].values),  max_nfev = 100)
    Plot_Elements_Spectrum(matrix, [(0,0)], fit = 'markov')
    """
    for (ii,jj) in serpentine_range(len(rows), len(cols), 'right'):

        if (ii,jj) in boni:

            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

            matrix[ii][jj].Get_VIPA_for_fit('interpolate', interpolation_density = 500)
            p0s = Get_p0_by_Neighbours(matrix, ii, jj, len(rows), len(cols))
        
            if (ii,jj) in almost_height:
                p0s.append(p0_almost)
            elif ((ii,jj) in brillouin_higher) | ((ii,jj) in brillouin_highest):
                p0s.append(p0_brillouin)
            else:#normals
                p0s.append(p0_normal)

            matrix[ii][jj].Get_Best_p0(p0s, cols_mark)

            matrix[ii][jj].Get_cost_markov(matrix[ii][jj].p0[list(cols_mark)].values[0], cols_mark)
            print('Cost before fitting = {}'.format(matrix[ii][jj].cost_markov))
            matrix[ii][jj].Get_Fit_Bounds(rules_markov_bounds, cols_mark)
            fit = fit + ((matrix[ii][jj].Non_Linear_Least_Squares_Markov(cols_mark, bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values),  max_nfev = 100),(ii,jj)),)
            matrix[ii][jj].Get_cost_markov(matrix[ii][jj].Markov_Fit_Params.values[0], cols_mark)
            print('Cost after fitting = {}\n'.format(matrix[ii][jj].cost_markov))

            if (ii,jj) in almost_height:
                p0_almost = matrix[ii][jj].Markov_Fit_Params.values[0]
            elif ((ii,jj) in brillouin_higher) | ((ii,jj) in brillouin_highest):
                p0_brillouin =  matrix[ii][jj].Markov_Fit_Params.values[0]
            else:#normals
                p0_normal =  matrix[ii][jj].Markov_Fit_Params.values[0]


            del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].y_markov_convolution

        #afterfit

    markov_time     =   time.process_time()-start
    tempo           =   tempo + (('fit markoviano', markov_time),)
    
    with open(analysis_path+log_file, 'a') as f_log:

        f_log.write('\n\n#####################   MARKOVIAN     FIT     ##########################\n\n')
        f_log.write('\nTempo impiegato per fit markoviani: {:3.2} s'.format(markov_time))
        f_log.write('\n\Tempo impiegato ore = {:3.2}\n'.format(markov_time/3600))

    # 4) after - fit markoviano

    non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

    #too_markov         =   Whose_Param_Too_High('Gamma', 2., matrix, fitted)
    Save_Fit_Info(fit, filename = 'markov_fit.txt', path=analysis_path)
    Save_Markov_Fit_Parameters(matrix, fitted, out_filename = 'markov_fit_params.txt', path = analysis_path)
    Save_y_markov_fit(matrix, fitted, path = analysis_path)
    Save_cost_markov(matrix, fitted, path = analysis_path)
    print('\nHo salvato informazioni fit markoviano su {}\n'.format(analysis_path))


else:

    print('\n\n You chose to SKIP the markovian fit and recover info from {}\n\n'.format(now_path+analysis_dir))

    ############## faccio er ricovery

    with open(analysis_path+'markov_fit.txt', 'r') as fin:

        fit     =   eval(fin.read())

    non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

    with open(analysis_path+'markov_fit_params.txt', 'r') as fin:
        lines   =   fin.readlines()

    if (len(boni) != len(lines)):
        raise ValueError("Incompatibilità tra lunghezza file parametri ({}) e lunghezza boni fit ({})\nFILE SBAGLIATO".format(len(lines), len(boni)))

    for (line, (ii,jj)) in zip(lines, fitted) :
        matrix[ii][jj].Recover_Markov_Fit_Params(line)

    print("\n\nI've correctely recovered markovian fit info \n\n")



#%%
##################################################################################################################################
######################################################################################################################


#######   ||||   ||      ||   TOT FIT: - opero fit totale con tutti i parametri
#######    ||     ||    ||               viscoelastici ma senza più i gaussiani
#######    ||      ||  ||
#######    ||       ||||   
#######   ||||       ||

if not skip_tot:
        
    fit_tot = ()
    print("\n\nI'm beginning total fit\n\n")
    start = time.process_time()

    for (ii,jj) in serpentine_range(len(rows), len(cols), 'right'):

        if (ii,jj) in boni:

            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

            p_gauss = matrix[ii][jj].Markov_Fit_Params[list(cols_gauss)].values[0]
            kernel   = matrix[ii][jj].VIPA_w_j/(p_gauss[0]*(np.exp(-((matrix[ii][jj].w_j_VIPA-p_gauss[1])**2)/(2*(p_gauss[2]**2)))))

            matrix[ii][jj].Initials_Parameters_from_Markov(matrix[ii][jj].Markov_Fit_Params, cols_mark)
            matrix[ii][jj].Get_Fit_Bounds(rules_tot_bounds, columns = cols_real)

            matrix[ii][jj].Get_cost_tot(matrix[ii][jj].p0.values[0], p_gauss, kernel)
            print('\nCost before fitting = {}\n'.format(matrix[ii][jj].cost_tot))
            fit_tot =   fit_tot + (((matrix[ii][jj].Non_Linear_Least_Squares(p_gauss, cols_real, bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values), max_nfev = 35)), (ii,jj)),)
            matrix[ii][jj].Get_cost_tot(matrix[ii][jj].Tot_Fit_Params.values[0], p_gauss, kernel)
            print('\nCost after fitting = {}\n'.format(matrix[ii][jj].cost_tot, kernel))
            #del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].res_lsq, matrix[ii][jj].bounds



    tot_time     =   time.process_time()-start
    tempo           =   tempo + (('fit tot', tot_time),)

    #after fit
    non_fitted_tot, accomplished_tot, exceded_tot, fitted_tot = Unpack_Fit(fit_tot)

    Save_Fit_Info(fit_tot, filename = 'tot_fit.txt', path=analysis_path)
    Save_Tot_Fit_Parameters(matrix, fitted_tot, out_filename = 'tot_fit_params.txt', path = analysis_path)
    Save_y_fit(matrix, fitted_tot, path = analysis_path)
    Save_cost_tot(matrix, fitted_tot, path = analysis_path)

    with open(analysis_path+log_file, 'a') as f_log:
        f_log.write("\n\n################## TOTAL FIT ####################\n\n")

else:
    print("\n\nSkippato fit tot\n\n")
    with open(analysis_path+log_file, 'a') as f_log:

        f_log.write("\n\n################## TOTAL FIT ####################\n\n")
        f_log.write("\n\nskipped\n\n")

super_time = super_start - time.process_time()
### FINAL
#%%
with open(analysis_path+log_file, 'a') as f_log:
    f_log.write("\n\n################### FINAL PERFORMACES ########################\n\n")
    f_log.write('tempo impiegato per esecuzione dello script ore = %3.2f\n '%(super_time/3600))
    for (what,t) in tempo:
        f_log.write('di cui %f secondi =  %f  ore in %s \n' %(t, t/3600, what))


# %%


# %%

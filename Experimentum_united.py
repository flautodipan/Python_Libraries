
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
import      configparser
import      sys


#%%
if sys.argv[1] != '-f':

    print('Sto in modalità terminale\n')
    spectra_filename = sys.argv[1]
    now_path = '../BRILLOUIN/TDP43/'+spectra_filename+'/'
    print("By default setting, I'm taking data from directory {}, hope it's correct\n".format(now_path))

    while True:

        print('Insert analysis directory name.\n Considera che per la presa dati {} ci sono le cartelle:'.format(spectra_filename))
        os.system('cd '+now_path +' && ls -lh -lt -d */' )
        analysis_name = input()

        if os.path.exists(now_path +analysis_name+'/'):

            print("\nDirectory with such name already exists.\nPress 'o' to overwrite it, or 'n' to generate a new directory for analysis\n")
            ans = input()
            if (ans == 'n'):

                print("Insert name of new directory\n")
                new_name = input()
                os.system('cd '+now_path +' && mkdir '+new_name)
                analysis_path = now_path  + new_name +'/'
                break

            elif (ans == 'o'):

                os.system('cd '+now_path +' && rm -r '+ analysis_name +'_backup/')
                os.system('cd '+now_path +' && cp -r '+ analysis_name+' '+analysis_name+'_backup/')
                print('I backed up your directory, for any inconvenience...\n')
                analysis_path = now_path + analysis_name+ '/'

                break
            
            else:
                print('\nValue inserted not correct\n Try again motherfucker\n')
        else:

            os.system('cd '+now_path +' && mkdir '+analysis_name+'/')
            analysis_path = now_path + analysis_name+ '/'
            break


elif sys.argv[1] == '-f': 

    print('Sto in modalità interattiva')
    spectra_filename = 'ARS_13_02'
    now_path = '../BRILLOUIN/TDP43/'+spectra_filename+'/'
    analysis_name = 'dabuttare'
    os.system('cd '+now_path +' && mkdir '+analysis_name+'/')
    analysis_path = now_path + analysis_name +'/'

else:

    raise ValueError('argv[1] = {} not recognized by software. Check \n'.format(sys.argv[1]))


#%%

inputs = configparser.ConfigParser()
with open('../BRILLOUIN/TDP43/ARS_13_02/config.ini', 'r') as f:
    inputs.read_file(f)

############
#I/O 

VIPA_filename       =   inputs['I/O']['VIPA_filename']
log_file            =   'log_'+spectra_filename

#operatives
initial             =   inputs['Operatives']['initial']
to_add              =   eval(inputs['Operatives']['to_add'])
exclude_delta       =   inputs.getboolean('Operatives', 'exclude_delta')
syg_kwargs          =  {item[0] : float(item[1]) for item in inputs.items('syg_kwargs')}
syg_kwargs_VIPA     =  {item[0] : float(item[1]) for item in inputs.items('syg_kwargs_VIPA')}
syg_kwargs_brill    =  {item[0] : float(item[1]) for item in inputs.items('syg_kwargs_brill')}
VIPA_treshold       =  inputs.getfloat('Operatives','VIPA_treshold')
sat_height          =  inputs.getfloat('Operatives','sat_height')
sat_width           =  inputs.getfloat('Operatives','sat_width')
almost_treshold     =  inputs.getfloat('Operatives','almost_treshold')
pre_cut             =  inputs.getboolean('Operatives','pre_cut')
cut                 =  inputs.getboolean('Operatives','cut')
mean_dist_01        =  inputs.getfloat('Operatives','mean_dist_01')
mean_dist_23        =  inputs.getfloat('Operatives','mean_dist_23')
#markov_fit

recover_markov      = inputs.getboolean('Markov', 'recover_markov')
p0_normal           = np.array(eval(inputs['Markov']['p0_normal']))
p0_brillouin        = np.array(eval(inputs['Markov']['p0_brillouin']))
p0_almost           = np.array(eval(inputs['Markov']['p0_almost']))

rules_markov_bounds =   eval(inputs['Markov']['rules_markov_bounds'])
#tot fit
skip_tot            =  inputs.getboolean('Tot', 'skip_tot')
rules_tot_bounds    =   eval(inputs['Tot']['rules_tot_bounds'])

############

if sys.argv[1] != '-f':

    if recover_markov:
        print('You decided to recover markov fit from {} \nIt is correct? Enter "ok" if so, any other key to change this option'.format(analysis_path))  
    else:
        print('You will perform markov fit. Enter "ok" to continue, any other key to modify this opt\n')
    if input() == 'ok':
        pass  
    else:
        while(True):
            print('Inserire "yes" to perform Markov, "no" to not perfom it')
            inp = input()
            if inp == 'yes':
                print('You will perform markov fit')
                time.sleep(2.)
                break
            elif inp == 'no':
                
                print('You will NOT perform markov fit and recover it from {}\nPress any key too continue otherwise change you analysis path running again script'.format(analysis_path))
                if input(): break
            else:
                print('Did not understand. Retry')

else:

    recover_markov = recover_markov
    skip_tot = skip_tot

##

with open(analysis_path+'config.ini', 'w') as fin:
    inputs.write(fin)
print('Ho salvato configurazione iniziale in file config.ini in directory {}'.format(analysis_path))


#%%

#variables

invisible           =   []
brillouin_higher    =   []
brillouin_highest   =   []
boni                =   []
excluded            =   []
almost_height       =   []
normals             =   []

cols        = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position',  'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_position','delta_width',  'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark_nodelta  = ('Co', 'Omega', 'Gamma', 'A', 'mu', 'sigma', 'shift', 'offset')
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
#matrix, rows, cols = Initialize_Matrix(11,34,13,36)
matrix, rows, cols = Initialize_Matrix(0,0, n_rows, n_cols)
dim     =   len(rows)*len(cols)


with open(analysis_path+log_file, 'w') as f_log:
    f_log.write('#This is a log file: you will find info on script run for {}\n'.format(spectra_filename))
    f_log.write('\n\nHo inizializzato una matrice {}x{}, per un totale di {} spettri'.format(len(matrix), len(matrix[0]), len(matrix)*len(matrix[0])))
    
# %%
### riempio oggetti spectrum

tempo = ()
super_start = time.process_time()
start = super_start

matrix[0][0].Get_VIPA_tif(VIPA_filename, now_path, offset = 183.)

for (ii, ii_true) in zip(range(len(rows)), rows):  
    for (jj, jj_true) in zip(range(len(cols)), cols):

        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii_true][jj_true],np.max(dati[ii_true][jj_true].shape)) , offset = 183., cut = pre_cut, cut_range = (200, 600))
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs)

        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA

### catalogo la natura degli spettri

not_saturated, saturated = Get_Saturated_Elements(matrix, len(rows), len(cols), saturation_height = sat_height, saturation_width = sat_width)
excluded        = saturated.copy()
excluded        = Escludi_a_Mano(to_add, excluded)

for (ii, ii_true) in zip(range(len(rows)), rows):  
    for (jj, jj_true) in zip(range(len(cols)), cols):
            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

            if (ii_true,jj_true) not in excluded:
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
    for (ii,jj) in serpentine_range(len(rows), len(cols), initial):

        if (ii,jj) in boni:

            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
                        
            if (exclude_delta) & ((ii,jj) not in almost_height):

                columns = cols_mark_nodelta
                rules_bounds = rules_markov_bounds[0:3]+rules_markov_bounds[6:]

            else:

                columns = cols_mark
                rules_bounds = rules_markov_bounds


            matrix[ii][jj].Get_VIPA_for_fit('interpolate', interpolation_density = 500)

            p0s = Get_p0_by_Neighbours(matrix, columns,  ii, jj, len(rows), len(cols))
        
            if (ii,jj) in almost_height:
                p0s.append(p0_almost)
            elif ((ii,jj) in brillouin_higher) | ((ii,jj) in brillouin_highest):
                p0s.append(p0_brillouin)
            else:#normals
                p0s.append(p0_normal)

            matrix[ii][jj].Get_Best_p0(p0s, columns)
            matrix[ii][jj].Get_Fit_Bounds(rules_bounds, columns)
            matrix[ii][jj].Get_cost_markov(matrix[ii][jj].p0[list(columns)].values[0], columns)
            print('Cost before fitting = {}'.format(matrix[ii][jj].cost_markov))
            fit = fit + ((matrix[ii][jj].Non_Linear_Least_Squares_Markov(columns, bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values),  max_nfev = 100),(ii,jj)),)
            matrix[ii][jj].Get_cost_markov(matrix[ii][jj].Markov_Fit_Params.values[0], columns)
            print('Cost after fitting = {}\n'.format(matrix[ii][jj].cost_markov))

            if (ii,jj) in almost_height:
                p0_almost = matrix[ii][jj].Markov_Fit_Params.values[0]
            elif ((ii,jj) in brillouin_higher) | ((ii,jj) in brillouin_highest):
                p0_brillouin =  matrix[ii][jj].Markov_Fit_Params.values[0]
            else:
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

    for (ii,jj) in serpentine_range(len(rows), len(cols), initial):

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

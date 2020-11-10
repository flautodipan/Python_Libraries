#########################                                                   #######################
#########################           EXPERIMENTUM UNITED                     #######################

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
from        datetime            import date
from        subprocess          import run
import      os
from        os.path             import join
import      configparser
import      sys

# viscoelastic models
cols        = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position',  'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_position','delta_width',  'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark_nodelta  = ('Co', 'Omega', 'Gamma', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_real_nodelta =  ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')

# dho models
cols_dho            = ('dho_position', 'dho_width', 'dho_amplitude', 'delta_position', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_dho_nodelta    = ('dho_position', 'dho_width', 'dho_amplitude','A', 'mu', 'sigma',  'shift', 'offset')

#%%
#ANALYSIS PATH 

mode = Check_Execution_Mode(sys.argv)

if mode == 'terminal':

    print("I'm executing by command prompt \n")
    spectra_filename = sys.argv[1]
    now_path = '../BRILLOUIN/TDP43/'+spectra_filename+'/'
    print("By default setting, I'm taking data from directory {}, hope it's correct\n".format(now_path))
    analysis_path = Get_Analysis_Path_From_Terminal(now_path, spectra_filename)

elif mode == 'interactive': 

    print("Executing python interactive window using Jupyter")
    spectra_filename = 'ARS_10_02'
    now_path = '../BRILLOUIN/TDP43/'+spectra_filename+'/'
    analysis_name = 'prova_{}'.format(date.today().strftime("%d_%m_%Y"))

    copy = 0
    backup = analysis_name
    
    while True:

        if not os.path.exists(now_path +analysis_name+'/'):
            #run(['cd', now_path, '&&', 'mkdir', analysis_name], shell = True)
            os.system('cd '+now_path +' && mkdir '+analysis_name+'/')
            break

        else: 
            copy += 1
            analysis_name = backup
            analysis_name += ('_' + str(copy))

        

    analysis_path = now_path + analysis_name +'/'


else:

    raise ValueError('argv[1] = {} not recognized by software. Check \n'.format(sys.argv[1]))


#%%
# 0) ACQUISISCO DATI DI INPUT 

##############
####   config.ini
##############

inputs = configparser.ConfigParser()

with open(now_path+'config.ini', 'r') as f:
    inputs.read_file(f)

############
#I/O 

VIPA_filename       =   inputs['I/O']['VIPA_filename']
log_file            =   'log_'+spectra_filename
inputs.set('I/O', 'analysis_path', analysis_path)
inputs.set('I/O', 'now_path', now_path)


# OPERATIVES for ACQUISITION
syg_kwargs          =   {item[0] : float(item[1]) for item in inputs.items('syg_kwargs')}
syg_kwargs_VIPA     =   {item[0] : float(item[1]) for item in inputs.items('syg_kwargs_VIPA')}
syg_kwargs_brill    =   {item[0] : float(item[1]) for item in inputs.items('syg_kwargs_brill')}

transpose           =   inputs.getboolean('Operatives', 'transpose')
data_offset         =   inputs.getfloat('Operatives', 'data_offset')
alignment           =   inputs['Operatives']['alignment']
almost_treshold     =   inputs.getfloat('Operatives','almost_treshold')
sat_height          =   inputs.getfloat('Operatives','sat_height')
sat_width           =   inputs.getfloat('Operatives','sat_width')
pre_cut             =   inputs.getboolean('Operatives','pre_cut')
pre_cut_range       =   eval(inputs['Operatives']['pre_cut_range'])
initial             =   inputs['Operatives']['initial']
to_add              =   eval(inputs['Operatives']['to_add'])

# OPERATIVES for FITs
model               =   inputs['Operatives']['model']
fit_algorithm       =   inputs['Operatives']['fit_algorithm']
cut                 =   inputs.getboolean('Operatives','cut')
exclude_delta       =   inputs.getboolean('Operatives', 'exclude_delta')
mean_dist_01        =   inputs.getfloat('Operatives','mean_dist_01')
mean_dist_12        =   inputs.getfloat('Operatives','mean_dist_12')
mean_dist_23        =   inputs.getfloat('Operatives','mean_dist_23')

#markov_fit
recover_markov      = inputs.getboolean('Markov', 'recover_markov')
ii_0,jj_0           = eval(inputs.get('Markov', 'first_normal'))
first_almost        = inputs.get('Markov', 'first_almost')
p0_normal           = np.array(eval(inputs['Markov']['p0_normal']))
p0_almost           = np.array(eval(inputs['Markov']['p0_almost']))
rules_markov_bounds =   eval(inputs['Markov']['rules_markov_bounds'])

#tot fit
skip_tot            =  inputs.getboolean('Tot', 'skip_tot')
rules_tot_bounds    =   eval(inputs['Tot']['rules_tot_bounds'])

############

# check delle cose operazionali

if mode == 'terminal':

    partial = False
    print('p0s lenght is {} for normal, {} for almost\n'.format(len(p0_normal), len(p0_almost)))
    recover_markov, skip_tot, exclude_delta, fit_algorithm = Check_Settings_From_Terminal(recover_markov, skip_tot, exclude_delta, fit_algorithm)
    inputs.set('Operatives', 'exclude_delta', str(exclude_delta))
    inputs.set('Markov', 'recover_markov', str(recover_markov))
    inputs.set('Tot', 'skip_tot', str(skip_tot) )
   
else:
    partial = False
    log_file            =   'inter_log_'+spectra_filename
    recover_markov = recover_markov
    skip_tot = skip_tot
    exclude_delta = exclude_delta
    initial = initial
    fit_algorithm = fit_algorithm
    print('Rec Mark = {}\nSkip Tot = {}\nexlude_delta={}\ninitial={}'.format(recover_markov, skip_tot, exclude_delta, initial))

inputs.set('Markov', 'fit_algorithm', fit_algorithm)


####################
#######  classes.ini
####################

classes = configparser.ConfigParser()

with open(join(now_path, 'classes.ini'), 'r') as f_in:
    classes.read_file(f_in)

normals             = eval(classes.get('tuples', 'normals'))
brillouin_higher    = eval(classes.get('tuples', 'brillouin_higher'))
almost_height       = eval(classes.get('tuples', 'almost_height'))
excluded            = eval(classes.get('tuples', 'excluded'))

# aggiuntivi
boni      = []
invisible = []
now_excluded = []

# %%
#2) ######################################################################################################################

#######   ||||  ||||    DATA ACQUISITION AND TREATMENT : - acquisiamo i file di INPUT
#######    ||    ||                                      - eseguiamo operazioni su cartelle
#######    ||    ||                                      - converto e taglio lo spettro
#######    ||    ||
#######   ||||  ||||

############
#  1) ACQUISITION  
############

dati    =   Import_from_Matlab(spectra_filename, now_path, transpose = transpose)
n_rows  =   len(dati)
n_cols  =   len(dati[0])
if partial:
    matrix, rows, cols = Initialize_Matrix(0,0,3,3)
    ii_0, jj_0 = serpentine_range(len(rows), len(cols), initial)[0]

else:
    matrix, rows, cols = Initialize_Matrix(0,0, n_rows, n_cols)
dim     =   len(rows)*len(cols)
inputs.set('I/O','n_rows', str(len(rows)))
inputs.set('I/O','n_cols', str(len(cols)))

with open(analysis_path+'config.ini', 'w') as fin:
    inputs.write(fin)
del inputs

with open(analysis_path+log_file, 'w') as f_log:
    f_log.write('#This is a log file: you will find info on script run for {}\n'.format(spectra_filename))
    f_log.write('\n\nHo inizializzato una matrice {}x{}, per un totale di {} spettri'.format(len(matrix), len(matrix[0]), len(matrix)*len(matrix[0])))
    
### riempio oggetti spectrum

tempo = ()
super_start = time.process_time()
start = super_start

matrix[ii_0][jj_0].Get_VIPA_tif(VIPA_filename, now_path, offset = 183.)

for (ii, ii_true) in zip(range(len(rows)), rows):  
    for (jj, jj_true) in zip(range(len(cols)), cols):

        print('row = {}/{} col = {}/{}'.format(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii_true][jj_true],np.max(dati[ii_true][jj_true].shape)) , offset = 183., cut = pre_cut, cut_range = pre_cut_range)
        matrix[ii][jj].Find_Spectrum_Peaks(**syg_kwargs)
    
        if (ii_true,jj_true) not in excluded:

            if matrix[ii][jj].n_peaks == 4:
                pass
            elif matrix[ii][jj].n_peaks < 4:
                matrix[ii][jj].Find_Spectrum_Peaks(**syg_kwargs_brill)
            elif matrix[ii][jj].n_peaks == 5:
                matrix[ii][jj].Exclude_Glass_Peak()
            else:
                matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()

            matrix[ii][jj].Align_Spectrum(alignment = alignment)

            # check ulteriore
            if matrix[ii][jj].Check_Brillouin_Distances(mean_dist_01 = mean_dist_01, mean_dist_12 = mean_dist_12, mean_dist_23 = mean_dist_23, stdev = 0.15):
                print('Acc! Something went wrong with spectrum (ii,jj) = {}.\nFor safety, I will exclude it'.format(str((ii,jj))))
                invisible.append((ii,jj),)
                now_excluded.append((ii.jj),)
            else:
                boni.append((ii,jj), )
            # acquisisco VIPA
            matrix[ii][jj].x_VIPA   =   matrix[ii_0][jj_0].x_VIPA
            matrix[ii][jj].y_VIPA   =   matrix[ii_0][jj_0].y_VIPA

            # set_nature
            if (ii,jj) in normals: matrix[ii][jj].Set_Spectrum_Nature('normal')
            elif (ii,jj) in almost_height: matrix[ii][jj].Set_Spectrum_Nature('almost_height')
            elif (ii,jj) in brillouin_higher: matrix[ii][jj].Set_Spectrum_Nature('brillouin_higher')
        
        else: now_excluded.append((ii,jj))

# aggiorno esclusi
if len(invisible) > 0:
    print('Attenzione: rispetto agli esclusi previsti da pre_Experimentum.py, ho aggiunto {} elementi'.format(len(invisible)))
    excluded    += invisible
    excluded.sort()
else: print("Ok, tutto bene nell'acquisizione spettri")

# check 
if (len(now_excluded) + len(boni)) != dim:
    raise ValueError('ATTENZIONE: excluded elements ({}) and ready_for_fit(boni) ({}) elements are not complementary to total dimension {}'.format(len(excluded), len(boni), dim))

acq_time    =   time.process_time()-start
tempo       =   tempo + (('acquisizione', acq_time),)
print('\nTempo impiegato per acquisizione spettri: {} s\n'.format(acq_time))

#%%
# scrivo su log file

with open(analysis_path+log_file, 'a') as f_log:

    f_log.write('\n\n###############################INFO ON DATA ACQUISITION#######################################\n\n')
    f_log.write('\nTotale spettri boni :  {} di cui\n'.format(len(boni)))
    f_log.write('\nTotale spettri con elastici più alti di {} : {}\n'.format(almost_treshold, len(almost_height)))
    f_log.write('\nTotale spettri normali : {}\n'.format(len(normals))) 
    f_log.write('\nTotale spettri almost height : {}\n'.format(len(almost_height)))
    f_log.write('\nTotale spettri con Brillouin più alti di un elastico : {}\n'.format(len(brillouin_higher)))
    f_log.write(str(brillouin_higher))
    f_log.write('\nTotale spettri con Brillouin invisibili (= aggiunti dopo ulteriore controllo): {}\n'.format(len(invisible)))
    f_log.write(str(invisible))
    f_log.write('\nTotale spettri esclusi : {}\n'.format(len(excluded)))
    f_log.write(str(excluded))   
    f_log.write('\n\nTotale spettri : {}\ndi cui {} inutilizzabili, '.format(dim, len(invisible)+len(excluded)))
    f_log.write('ossia il %3.2f percento\n'%(float(len(invisible)+len(excluded))*100/dim))

# %%
############
#  1) TREATMENT  
############
start = time.process_time()

matrix[ii_0][jj_0].How_Many_Peaks_To_VIPA(treshold = 6, **syg_kwargs_VIPA)
matrix[ii_0][jj_0].Fit_Pixel2GHz()
matrix[ii_0][jj_0].VIPA_Pix2GHz()
matrix[ii_0][jj_0].Spectrum_Pix2GHz()
matrix[ii_0][jj_0].Estimate_Spectrum_Parameters(verbose = True)
matrix[ii_0][jj_0].Cut_Spectrum(mean_dist01 = mean_dist_01, mean_dist23 = mean_dist_23)
matrix[ii_0][jj_0].Fit_VIPA_Gaussian()

for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('row = {}/{} col = {}/{}'.format(ii,len(rows)-1, jj, len(cols)-1))

        if ((ii,jj) not in excluded) & ((ii,jj) != (ii_0,  jj_0)):

            matrix[ii][jj].x_VIPA_freq   =   matrix[ii_0][jj_0].x_VIPA_freq
            matrix[ii][jj].y_VIPA        =   matrix[ii_0][jj_0].y_VIPA
            matrix[ii][jj].Poly2GHz      =   matrix[ii_0][jj_0].Poly2GHz
            matrix[ii][jj].Spectrum_Pix2GHz()
            matrix[ii][jj].Cut_Spectrum(mean_dist01 = mean_dist_01, mean_dist23 = mean_dist_23)

        elif ((ii,jj) in excluded):
            matrix[ii][jj].Poly2GHz      =   matrix[ii_0][jj_0].Poly2GHz
            matrix[ii][jj].Spectrum_Pix2GHz()

mod_time    =   time.process_time()-start
tempo       =   tempo + (('modifica',mod_time), )

#%%
# salvataggi

with open(analysis_path+log_file, 'a') as f_log:
    f_log.write('\n\n######### MODIFICA SPETTRI #####################\n\n')
    f_log.write('\nTempo impiegato per modifica spettri: {} s\nTaglio spettri è {}\n\n'.format(mod_time, cut))

# salvo info spettri e VIPA
Save_XY_position(matrix, len(rows), len(cols), path = analysis_path)
Save_XY_VIPA(matrix[ii_0][jj_0].x_VIPA_freq, matrix[ii_0][jj_0].y_VIPA, path = analysis_path)

print('\n\nI saved xy info on xy.txt and xy_VIPA.txt in your analysis directory {}\n\n'.format(analysis_path))
with open(analysis_path+log_file, 'a') as f_log:
    f_log.write('\n\nI saved xy info on xy.txt and xy_VIPA.txt in your analysis directory {}\n\n'.format(analysis_path))

#%%
#######################################################################################################################

#######   ||||  ||||  ||||   MARKOVIAN FIT: - opero fit markoviano con tutti i parametri
#######    ||    ||    ||                     tranne Delta e tau (quindi anche Gauss)
#######    ||    ||    ||
#######    ||    ||    ||
#######   ||||  ||||  ||||

print("\nI'm beginning markovian fit\nFit algorithm = {}\n\n".format(fit_algorithm))
time.sleep(3)
if recover_markov == False:
        
    print('\n\n You chose to do the markovian fit\n\n')
    fit                 =   ()
    start = time.process_time()
    isolated = Get_Isolated_Elements(excluded)

    for (ii,jj) in serpentine_range(len(rows), len(cols), initial):

        if (ii,jj) in boni:

            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

            if (exclude_delta) & ((ii,jj) not in almost_height):
                columns = cols_mark_nodelta
                rules_bounds = {key : rules_markov_bounds[key] for key in columns}
                if len(p0_normal) != len(rules_bounds):
                    p0_normal = np.concatenate([p0_normal[0:3], p0_normal[6:]])
            else:
                columns = cols_mark
                rules_bounds = rules_markov_bounds

            matrix[ii][jj].Get_VIPA_for_fit('interpolate', interpolation_density = 500)

            p0s = Get_p0_by_Neighbours(matrix, columns, ii, jj, len(rows), len(cols))
        
            if (ii,jj) in almost_height:

                p0s.append(p0_almost)

            else:#normals
                p0s.append(p0_normal)

            matrix[ii][jj].Get_Best_p0(p0s, columns)
            matrix[ii][jj].Get_cost_markov(matrix[ii][jj].p0[list(columns)].values[0], columns)
            print('Cost before fitting = {:3.2f}'.format(matrix[ii][jj].cost_markov))
            matrix[ii][jj].Get_Fit_Bounds(rules_bounds, columns)
            if fit_algorithm == 'trf':
                fit = fit + ((matrix[ii][jj].Non_Linear_Least_Squares(columns, bounds = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values),  max_nfev = 100, ),(ii,jj)),)
            elif fit_algorithm == 'lm':
                fit = fit + ((matrix[ii][jj].Non_Linear_Least_Squares(columns,  max_nfev = 100, method = 'lm'),(ii,jj)),)
            else: raise ValueError('Possibile non si è scelto un valore per fit_algorithm = {} ?'.format(fit_algorithm))
            matrix[ii][jj].Get_cost_markov(matrix[ii][jj].Markov_Fit_Params.values[0], columns)
            print('Cost after fitting = {:3.2f}\n'.format(matrix[ii][jj].cost_markov))

            if (ii,jj) in almost_height:
                p0_almost = matrix[ii][jj].Markov_Fit_Params.values[0]

            else:#normals
                p0_normal =  matrix[ii][jj].Markov_Fit_Params.values[0]


            del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].y_markov_convolution

        else:
            fit += ((0, (ii, jj)), )


        #afterfit

    markov_time     =   time.process_time()-start
    tempo           =   tempo + (('fit markoviano', markov_time),)
    
    with open(analysis_path+log_file, 'a') as f_log:

        f_log.write('\n\n#####################   MARKOVIAN     FIT     ##########################\n\n')
        f_log.write('\nTempo impiegato per fit markoviani: {:3.2} s'.format(markov_time))
        f_log.write('\nTempo impiegato ore = {:3.2}\n'.format(markov_time/3600))

    # 4) after - fit markoviano

    non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

    #too_markov         =   Whose_Param_Too_High('Gamma', 2., matrix, fitted)
    Save_Fit_Info(fit, filename = 'markov_fit.txt', path=analysis_path)
    Save_Markov_Fit_Parameters(matrix, fitted, out_filename = 'markov_fit_params.txt', path = analysis_path)
    Save_y_markov_fit(matrix, fitted, path = analysis_path)
    Save_cost_markov(matrix, fitted, path = analysis_path)
    print('\nHo salvato informazioni fit markoviano su {}\n'.format(analysis_path))


else:

    print('\n\n You chose to SKIP the markovian fit and recover info from {}\n\n'.format(analysis_path))

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


            if 'delta_position' not in matrix[ii][jj].Markov_Fit_Params.keys():
                columns = cols_real_nodelta
                rules_bounds = rules_tot_bounds[0:5]+rules_tot_bounds[8:]
            else:
                columns = cols_real
                rules_bounds = rules_tot_bounds


            p_gauss = matrix[ii][jj].Markov_Fit_Params[list(cols_gauss)].values[0]
            kernel   = matrix[ii][jj].VIPA_w_j/(p_gauss[0]*(np.exp(-((matrix[ii][jj].w_j_VIPA-p_gauss[1])**2)/(2*(p_gauss[2]**2)))))

            matrix[ii][jj].Initials_Parameters_from_Markov()
            matrix[ii][jj].Get_Fit_Bounds(rules_bounds, columns = columns)
            matrix[ii][jj].Get_cost_tot(matrix[ii][jj].p0.values[0], p_gauss, kernel, columns)
            print('\nCost before fitting = {:3.2f}\n'.format(matrix[ii][jj].cost_tot))
            fit_tot =   fit_tot + (((matrix[ii][jj].Non_Linear_Least_Squares_NoGauss(p_gauss, cols_real, bounds = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values), max_nfev = 35)), (ii,jj)),)
            matrix[ii][jj].Get_cost_tot(matrix[ii][jj].Tot_Fit_Params.values[0], p_gauss, kernel, columns)
            print('\nCost after fitting = {:3.2f}\n'.format(matrix[ii][jj].cost_tot))
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


# %%

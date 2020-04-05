"""
VERSIONE FINALE PROGRAMMA DI FIT per dati BRILLOUIN(Aprile 2020)
TESI MAGISTRALE

Funzionamento schematico:

    Premessa: il programma è differente nei primi passi se viene eseguito da finestra interattiva o da terminale, ma solo a livello di 
              controllo diretto che si può fare o meno sulle variabili operative e sul nome della cartella di analisi

    I)  PREAMBOLO : 
        - importo librerie utilizzate nell'esecuzione
        - inizializzo variabili globali che mi serviranno per tutto il corso dell'esecuzione
        - programma gira con argv[1] = <nome della cartella dati nel percorso ../BRILLOUIN/TDP43/> := "cartella madre"
        - mi viene chiesto da terminale nome della cartella da creare in quel percorso dove si salvano tutti i risultati 
        - eseguo operazioni di I/O da file config.ini che sta in cartella madre (generato da pre_Experimentum_united.py, v.)
        - ho la possibilità di controllare e cambiare opzioni operative (esecuzione/recupero fit markoviano, esecuzione fit tot, 
          esculsione delta dai fit, algoritmo di fit)
        - salvo la configurazione iniziale di input direttamente legata all'esecuzione nella cartella di analisi generata

    II) DATA ACQUISITION AND TREATMENT:
        - importo i dati dai file nella cartella madre e li metto in una matrice di oggetti Spectrum (creati in lib_Brillouin.py, v.)
        - classifico gli spettri a classi di appartenenza 


"""

#%%
######################################################################################################################
#######   ||||    PREAMBOLO: - acquisiamo i file di INPUT
#######    ||                - eseguiamo operazioni su cartelle
#######    ||   
#######    ||    
#######   ||||

################################################################
###########################################################
##################################################

#   I.I importo librerie e definisco variabili globali

#libraries

import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Brillouin       import  *
from        Alessandria         import  *
import      time
import      os
import      configparser
import      sys

#variables immutable

invisible           =   []
boni                =   []
excluded            =   []
almost_height       =   []
normals             =   []

cols        = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position',  'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_position','delta_width',  'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark_nodelta  = ('Co', 'Omega', 'Gamma', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_real_nodelta =  ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')

#%%
#   I.II Creo cartella di analisi

#ANALYSIS PATH 
if sys.argv[1] != '-f':

    print('Sto in modalità terminale\n')
    spectra_filename = sys.argv[1]
    now_path = '../BRILLOUIN/TDP43/'+spectra_filename+'/'
    print("By default setting, I'm taking data from directory {}, hope it's correct\n".format(now_path))
    
    analysis_path = Get_Analysis_Path_From_Terminal(now_path, spectra_filename)

elif sys.argv[1] == '-f': 

    print('Sto in modalità interattiva')
    spectra_filename = 'NO_ARS_13_02'
    now_path = '../BRILLOUIN/TDP43/'+spectra_filename+'/'
    analysis_name = 'dabuttare'
    if not os.path.exists(now_path +analysis_name+'/'):
        os.system('cd '+now_path +' && mkdir '+analysis_name+'/')
    analysis_path = now_path + analysis_name +'/'

else:

    raise ValueError('argv[1] = {} not recognized by software. Check \n'.format(sys.argv[1]))

#%%
#   I.III Importo il file config.ini in cui sono contenute le variabili necessarie a esecuzione programma

inputs = configparser.ConfigParser()

with open(now_path+'config.ini', 'r') as f:
    inputs.read_file(f)

############

#I/O 

VIPA_filename       =   inputs['I/O']['VIPA_filename']
log_file            =   'log_'+spectra_filename
transpose           =   inputs.getboolean('I/O', 'transpose')

inputs.set('I/O', 'analysis_path', analysis_path)
inputs.set('I/O', 'now_path', now_path)


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
first_normal        = inputs.get('Markov', 'first_normal')
first_almost         = inputs.get('Markov', 'first_almost')
p0_normal           = {first_normal : pd.Series(np.array(eval(inputs['Markov']['p0_normal'])), list(cols_mark_nodelta) if exclude_delta else list(cols_mark))}
p0_almost           = {first_almost : pd.Series(np.array(eval(inputs['Markov']['p0_almost'])), list(cols_mark_nodelta) if exclude_delta else list(cols_mark))}

rules_markov_bounds =   eval(inputs['Markov']['rules_markov_bounds'])
#tot fit
skip_tot            =  inputs.getboolean('Tot', 'skip_tot')
rules_tot_bounds    =   eval(inputs['Tot']['rules_tot_bounds'])

############

if sys.argv[1] != '-f':

    print('p0s lenght is {} for normal, {} for almost\n'.format(len(p0_normal), len(p0_almost)))
    recover_markov, skip_tot, exclude_delta, method = Check_Settings_From_Terminal(recover_markov, skip_tot, exclude_delta)
    inputs.set('Operatives', 'exclude_delta', str(exclude_delta))
    inputs.set('Markov', 'recover_markov', str(recover_markov))
    inputs.set('Tot', 'skip_tot', str(skip_tot) )
else:

    log_file            =   'inter_log_'+spectra_filename
    recover_markov = False
    skip_tot = False
    exclude_delta = False
    initial = initial
    method = 'trf'
    print('Rec Mark = {}\nSkip Tot = {}\nexlude_delta={}\ninitial={}'.format(recover_markov, skip_tot, exclude_delta, initial))

inputs.set('Markov', 'method', method)

# %%
#2) ######################################################################################################################

#######   ||||  ||||    DATA ACQUISITION AND TREATMENT : - acquisiamo i file di INPUT
#######    ||    ||                                  - eseguiamo operazioni su cartelle
#######    ||    ||
#######    ||    ||
#######   ||||  ||||

###################################################################
##########################################################
###################################

#   II.I    Importo dati da file .mat e inizializzo la matrice di oggetti e la riempio

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, now_path, transpose = transpose, var_name = 'y3')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
#matrix, rows, cols = Initialize_Matrix(0,9,2,11)
matrix, rows, cols = Initialize_Matrix(0,0, n_rows, n_cols)
dim     =   len(rows)*len(cols)
ii_0,jj_0 = eval(first_normal)


### riempio oggetti spectrum

tempo = ()
super_start = time.process_time()
start = super_start

matrix[ii_0][jj_0].Get_VIPA_tif(VIPA_filename, now_path, offset = 183.)

for (ii, ii_true) in zip(range(len(rows)), rows):  
    for (jj, jj_true) in zip(range(len(cols)), cols):

        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii_true][jj_true],np.max(dati[ii_true][jj_true].shape)) , offset = 183., cut = pre_cut, cut_range = (200, 600))
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs)

        matrix[ii][jj].x_VIPA   =   matrix[ii_0][jj_0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[ii_0][jj_0].y_VIPA


#           A questo punto salvo file di config legato a esecuzione (mi serve info righe colonne)
#           Genero file di log

inputs.set('I/O','n_rows', str(len(rows)))
inputs.set('I/O','n_cols', str(len(cols)))
with open(analysis_path+'config.ini', 'w') as fin:
    inputs.write(fin)
del inputs

print('Ho salvato configurazione iniziale in file config.ini in directory {}'.format(analysis_path))

with open(analysis_path+log_file, 'w') as f_log:
    f_log.write('#This is a log file: you will find info on script run for {}\n'.format(spectra_filename))
    f_log.write('\n\nHo inizializzato una matrice {}x{}, per un totale di {} spettri'.format(len(matrix), len(matrix[0]), len(matrix)*len(matrix[0])))
    
#%%
#   II.II Catalogo la natura degli spettri

# Prima divisione: non saturati, esclusi = saturati + spettri che per vari motivi (v. pre_Exp.py) sono stati scartati

not_saturated, saturated = Get_Saturated_Elements(matrix, len(rows), len(cols), saturation_height = sat_height, saturation_width = sat_width)
excluded    = (saturated+to_add).sort()

# Seconda suddivisione: spettri normali (in cui sono stati inglobati anche i brillouin_higher di un tempo)
#                       e spettri almost_height che sono quelli con l'elastico molto alto, che avrà un peso nel fit, 
#                       nonostante il taglio

# OSS: (ii_true, jj_true) servono per quando scelgo di prendere un sottoinsieme dei dati
#       ed evitare che prenda sempre i primi elementi ma realmente quelli che voglio
for (ii, ii_true) in zip(range(len(rows)), rows):  
    for (jj, jj_true) in zip(range(len(cols)), cols):
            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

            if (ii_true,jj_true) not in excluded:

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

                elif (matrix[ii][jj].n_peaks < 4) & (matrix[ii][jj].n_peaks > 1):

                    matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_brill)
                    normals.append((ii,jj))
                    boni += [(ii,jj),]

                else:
                    raise ValueError("Numero di picchi non previsti ({} )dal codice per spettro {}".format(matrix[ii][jj].n_peaks, str((ii,jj))))

excluded    += invisible
excluded.sort()

acq_time    =   time.process_time()-start
tempo       =   tempo + (('acquisizione', acq_time),)
print('\nTempo impiegato per acquisizione spettri: {} s\n'.format(acq_time))

# Salvo queste info
if sys.argv[1] != '-f':
        
    with open(analysis_path+log_file, 'a') as f_log:

        f_log.write('\n\n###############################INFO ON DATA ACQUISITION#######################################\n\n')
        f_log.write('\nTotale spettri saturati : {}\n'.format(len(saturated)))
        f_log.write(str(saturated))
        f_log.write('\nTotale spettri con Brillouin invisibili (= aggiunti dopo ulteriore controllo): {}\n'.format(len(invisible)))
        f_log.write(str(invisible))
        f_log.write('\nTotale spettri boni :  {} di cui\n'.format(len(boni)))
        f_log.write('\nTotale spettri con elastici più alti di {} : {}\n'.format(almost_treshold, len(almost_height)))
        f_log.write('\nTotale spettri normali : {}\n'.format(len(normals))) 
        f_log.write('\nTotale spettri almost height : {}\n'.format(len(almost_height)))
        f_log.write('\n\nTotale spettri : {}\ndi cui {} inutilizzabili, '.format(dim, len(invisible)+len(saturated)))
        f_log.write('ossia il %3.2f percento\n'%(float(len(invisible)+len(saturated))*100/dim))
        f_log.write('Di cui il {} invisibili e il {} saturati\n\n'.format(len(invisible)*100/dim, len(saturated)*100/dim))

#%%
#   II.III  Converto spettri in GHz e li taglio
#   
start = time.process_time()

matrix[ii_0][jj_0].How_Many_Peaks_To_VIPA(treshold = VIPA_treshold, **syg_kwargs_VIPA)
matrix[ii_0][jj_0].Fit_Pixel2GHz()
matrix[ii_0][jj_0].VIPA_Pix2GHz()
matrix[ii_0][jj_0].Align_Spectrum()
matrix[ii_0][jj_0].Spectrum_Pix2GHz()
matrix[ii_0][jj_0].Cut_n_Estimate_Spectrum(estimate = True, cut = cut, mean_dist01 = mean_dist_01, mean_dist23 = mean_dist_23)
matrix[ii_0][jj_0].Fit_VIPA_Gaussian()

for ii in range(len(rows)):
    for jj in range(len(cols)):
            print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
            if ((ii,jj) not in excluded) & ((ii,jj) != (0,0)):
                matrix[ii][jj].x_VIPA_freq   =   matrix[ii_0][jj_0].x_VIPA_freq
                matrix[ii][jj].y_VIPA        =   matrix[ii_0][jj_0].y_VIPA
                matrix[ii][jj].Poly2GHz      =   matrix[ii_0][jj_0].Poly2GHz
                matrix[ii][jj].Align_Spectrum(alignment = matrix[ii_0][jj_0].alignment)
                matrix[ii][jj].Spectrum_Pix2GHz()
                matrix[ii][jj].Cut_n_Estimate_Spectrum(cut = cut, mean_dist01 = mean_dist_01, mean_dist23 = mean_dist_23)
            elif ((ii,jj) in excluded):
                matrix[ii][jj].Poly2GHz      =   matrix[ii_0][jj_0].Poly2GHz
                matrix[ii][jj].Spectrum_Pix2GHz()

mod_time    =   time.process_time()-start
tempo       =   tempo + (('modifica',mod_time), )

# %%

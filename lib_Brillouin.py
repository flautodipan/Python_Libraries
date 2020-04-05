"""
VERSIONE FINALE LIBRERIA PER PROGRAMMA DI FIT per dati BRILLOUIN(Aprile 2020)
TESI MAGISTRALE

Struttura schematica:
    I) Librerie importate e variabili globali
    II) Classe Spectrum per il trattamento del singolo spettro
    III) Funzioni esterne alla classe per gestione I/O e della matrice di oggetti spectrum
         dell'acquisizione completa

"""

# PARTE I : librerie importate e variabili globali

import  time
import  numpy               as np
import  pandas              as pd
from    matplotlib.pyplot   import plot
import  matplotlib.pyplot   as plt
from    scipy.signal        import find_peaks
from    scipy.optimize      import leastsq
from    scipy.optimize      import least_squares
from    scipy.optimize      import curve_fit
from    scipy.io            import loadmat
import  json
import  os

from    Models              import S_Dynamical_Form_Factor_2, S_2_Generate, S_Dynamical_Form_Factor_0, S_0_Generate, S_Dynamical_Form_Factor_0_nodelta, S_Dynamical_Form_Factor_2_nodelta     
from    Alessandria         import *
from    lmfit               import Model

free_spectral_range =   29.9702547 #GHz

cols      = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_position', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark_nodelta  = ('Co', 'Omega', 'Gamma', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_real_nodelta =  ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')


#   PARTE II : Classe Spectrum


class Spectrum  :

    #   Parte 0 : costruttore

    def __init__(self, name):

        self.name       =   name

    def __str__(self):

        return '{}'.format(self.name)

    #   Parte 1 : esecuzione del programma

    def Get_Spectrum(self, y, offset = 183., cut = False, cut_range = None, fig = False):

        """

        Acquisisco lo spettro da un array di segnale y - offset
        i pixel sono generati come range di interi da 0 a len(y)
        associo errore sperimentale come radice di y, processo poissoniano

        """

        self.x      =   np.arange(1, len(y)+1, 1)
        self.offset     =   offset     
        self.y          =   y - self.offset

        ###levo elementi <0.01 problematici per le incertezze poissoniane (ne escono alcuni negativi)
        ### (W la fi(si)ca)

        self.x          =   self.x[self.y >= 0.01]
        self.y          =   self.y[self.y >= 0.01]

        if cut:

            self.x      =   self.x[cut_range[0]:cut_range[1]]
            self.y          =   self.y[cut_range[0]:cut_range[1]]

        self.y_err      =   np.sqrt(self.y)

        if fig:

            plt.figure()
            plt.plot(self.x, self.y)
            plt.plot(self.x[self.peaks[0]], self.y[self.peaks[0]], '*')
            plt.xlabel('Pixels')
            plt.title('Experimental spectrum for '+self.__str__())
            plt.savefig(fig+'.png')
            plt.show()


    def Get_Spectrum_Peaks(self, **syg_kwargs):

        """
        Funzione base wrapper di scipy.find_peaks() 
        Trova picchi dello spettro (array y al momento della chiamata) e genera attributo
        con le caratteristiche date da syg_kwargs (altezza, distanza, ampiezza)
        """

        pk =  find_peaks(self.y, **syg_kwargs)

        self.peaks      =   {'idx': pk[0], 'heights': pk[1]['peak_heights'], 'widths' : pk[1]['widths']}       
        self.n_peaks    =   self.peaks['idx'].size
    

    def Get_Spectrum_4_Peaks_by_Height(self):

        """
        Funzione che in ogni caso mi ritorna un dict con informazioni sui 4 picchi dello spettro, disposti in ordine 
        elastico, brillouin stockes, brillouin antistockes, elastico

        """
        self.peaks      =   Find_Highest_n_peaks(self.peaks, 4)
        self.n_peaks    =   self.peaks['idx'].size

    def Get_VIPA_tif(self, tif_filename, path='', offset = 'same', **img_kwargs):

            """
            Acquisisco i dati della funzione di risposta del VIPA, necessari poi per convoluzione
            Questi sono in formato di immagine .tiff
            """

            VIPA                =   Import_TIF(path+tif_filename)
            
            if (offset != 'same'):

                self.offset =   offset
                
            self.y_VIPA         =   np.mean(VIPA,1) - self.offset
            self.x_VIPA         =   np.arange(1, len(self.y_VIPA)+1, 1)

            if 'fig' in img_kwargs:
                
                plt.figure()
                plt.plot(self.x_VIPA, self.y_VIPA)
                plt.xlabel('Pixels')
                plt.ylabel('VIPA transfer function')
                plt.title('Funzione di trasferimento VIPA in pixel da '+tif_filename)
                plt.savefig(img_kwargs['save_path']+img_kwargs['fig']+'.png')
                plt.show()

    
    def Check_Brillouin_Distances(self, average, stdev):
        """
        Function that is supposed to be called when peaks are already four
        --> check if the distance of Brillouin peaks is in the average 
            otherwise it means they are not brillouin
        """
        if (np.abs(self.x[self.peaks['idx'][1]] - self.x[self.peaks['idx'][2]]) > (average*stdev)):
            return True
        else:
            return False

    def How_Many_Peaks_To_VIPA(self, treshold, n_GHz_peaks = 5, n_gauss_peaks = 3, delta = 1., fig = False, verbose = False, **syg_kwargs):

        h_save = ()

        for n_peaks in (n_GHz_peaks, n_gauss_peaks):
                
            pk      =   find_peaks(self.y_VIPA, height = treshold, **syg_kwargs)
            height  =   np.max(pk[1]['peak_heights'])/2
    
            while True:
                pk      =   find_peaks(self.y_VIPA, height= height, **syg_kwargs)
                if (height > treshold) & (pk[0].size == n_peaks):
                    h_save  =   h_save + (height,)
                    if verbose:
                        print("Ho trovato valore dell'altezza per avere %d picchi: %f\n"%(n_peaks, height), pk)
                        _ = Analyze_Peaks(self.x_VIPA, self.y_VIPA, 'GHz', fig = fig, verbose = verbose, height = height, **syg_kwargs)
                    break
                elif (height <= treshold):
                    print(pk)
                    raise ValueError('Errore: superata altezza minima %f\nQualcosa è andato storto'%(treshold))
                else: 
                    height-=delta
    
    def Align_Spectrum(self, alignment = False):

        """
        Funzione che allinea correttamente gli spettri con gli elastici di convoluzione VIPA 
        """
        if alignment:
            pass
        else:
            alignment = getattr(self, 'alignment')
        if alignment == 'dx':
            self.x      = self.x - self.x[self.peaks['idx'][0]]

        elif alignment == 'sx':
            self.x      = self.x - self.x[self.peaks['idx'][3]]


    def Spectrum_Pix2GHz (self, fig = False):

        #modifico in GHz le ascisse dello spettro

        self.x_freq     =   ((self.x**3)*self.Poly2GHz[0])+ ((self.x**2)*self.Poly2GHz[1]) + (self.x*self.Poly2GHz[2]) + (self.Poly2GHz[3])
        
        if fig:

            plt.figure()
            plt.xlabel('GHz')
            plt.plot(self.x_freq, self.y)
            plt.title('Spettro Exp in GHz')
            plt.show()
            plt.close()

    def Cut_n_Estimate_Spectrum(self, cut = True, mean_dist01 = 37, mean_dist23 = 34, estimate = False, verbose = False):
        
        """

        Funzione che esegue 
        #
        # taglio      :     trovo i valori dei picchi elastici e delle relative ampiezze 
        #                   (grazie scipy) e banalmente mi sposto di distanza rispetto 
        #                   alle posizioni dei picchi sull'array delle frequenze
        #                   il taglio sulle y è consequenziale, il tutto grazie agli indici
        #
        # stima dei parametri iniziali :
        #                    riesco a stimare qualche parametro e genero un p0 tarato per fit markoc
        #                    che in ogni caso è il primo che faccio
        #  
        #   OSS: tutto funziona perchè mi aspetto due picchi elastici ad aprire e
        #        chiudere lo spettro
        #      

        """
        if self.n_peaks == 2:
            
            #per chi non ha elastici
            idx_min               =   self.peaks['idx'][0] - int(mean_dist01/2)
            idx_max               =   self.peaks['idx'][1] + int(mean_dist23/2)

        else:
            idx_min               =   self.peaks['idx'][1] - int(mean_dist01/2)
            idx_max               =   self.peaks['idx'][2] + int(mean_dist23/2)
        
        # STIMA PARAMETRI INIZIALI della funzione teorica
        # (quelli che posso, e devo farlo  prima di tagliare)
        
        self.p0  =   pd.DataFrame({}, columns = cols, index = ['Values'])

        # 1)    stima dei parametri dai dati
        #       Omega come la media della posizione dei massimi brillouin
        #       Gamma (=Delta)come la media dell'ampiezza restituita da find_peaks per i Brillouin /10 (mah..)
        #       offset come la media dei valori dello spettro nella zona tra i due picchi

        if  estimate:
                
            self.p0['Omega']            =   [np.absolute(self.x_freq[self.peaks['idx'][2]] - self.x_freq[self.peaks['idx'][1]])*0.5]
            self.p0['Gamma']            =   [0.1]
            self.p0['offset']           =   [0]
            self.p0['Co']               =   [0.01]
            self.p0['shift']            =   [0.]
            self.p0['delta_position']   =   [0.]
            self.p0['delta_amplitude']  =   [1]
            self.p0['delta_width']      =   [0.1]
            self.p0['Delta']            =   self.p0['Gamma']
            self.p0['tau']              =   [1.]

        
        # procedo con il taglio se è voluto

        if cut:           
            
            self.x_freq         =   self.x_freq[idx_min:idx_max]
            self.x              =   self.x[idx_min:idx_max]
            self.y              =   self.y[idx_min:idx_max]
            self.y_err          =   self.y_err[idx_min:idx_max]


    #   Parte 2 : analisi dei risultati, funzioni di recovery 



#   PARTE III: Funzioni esterne di gestione I/O e della matrice di oggetti Spectrum

def Get_Analysis_Path_From_Terminal(now_path, spectra_filename):

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
    
    return analysis_path

def Check_Settings_From_Terminal(recover_markov, skip_tot, exclude_delta ):

    delay = 3. #sec

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
                time.sleep(delay)
                break
            elif inp == 'no':
                
                print('You will NOT perform markov fit and recover it from {}\nPress any key too continue otherwise change you analysis path running again script'.format(analysis_path))
                if input(): break
            else:
                print('Did not understand. Retry')    

    if skip_tot:
        print('Skipping fit_tot is active. Press "ok" to confirm, any other to change')
        if input() == 'ok':
            pass
        else:
            skip_tot = False
            print('You will perform fit tot')
            time.sleep(delay)
    else:
        print('Fit tot will be performed. Press "ok" to confirm, any other to change')
        if input() == 'ok':
            pass
        else:
            skiptot = True
            print('You will exclude delta in all fits')
            time.sleep(delay)

    
    if exclude_delta:
        print('Exclude delta from fit is active. Press "ok" to confirm, any other to change')
        if input() == 'ok':
            pass
        else:
            exclude_delta = False
            print('You will include delta in all fits')
            time.sleep(delay)
    else:
        print('Fit will all be performed with delta. Press "ok" to confirm, any other to change')
        if input() == 'ok':
            pass
        else:
            exclude_delta = True
            print('You will exclude delta in all fits')
            time.sleep(delay)
    
    print('Insert fit algorithm: lm for Levenberg-Marqadrart, trf per trust region')

    while True:
        method = input()
        if method == 'lm':
            print('You choose lm')
            time.sleep(delay)
            break
        elif method == 'trf':
            print('You choose trf')
            time.sleep(delay)
            break
        else: print('Did not understand. Retry.\n')

    return recover_markov, skip_tot, exclude_delta, method

    
def Initialize_Matrix(ii_0, jj_0, ii_stop, jj_stop):

    matrix = ()
    
    rows = np.arange(ii_0, ii_stop, 1)
    cols = np.arange(jj_0, jj_stop, 1)
    
    for ii in rows:
        riga = ()
        for jj in cols:
            riga  = riga + (Spectrum('Element ('+str(ii)+','+str(jj)+')'),)
        matrix = matrix + (riga,)

    print('Ho inizializzato una matrice %dx%d, per un totale di %d spettri'%(len(matrix), len(matrix[0]), len(matrix)*len(matrix[0])  ))

    return (matrix, rows, cols)



def Get_Saturated_Elements(matrix, n_rows, n_cols, saturation_height = 40000, saturation_width = 15.):

    saturated = []
    not_saturated = []

    for ii in range(n_rows):
        for jj in range(n_cols):
            pk_max_idx  =   np.argmax(matrix[ii][jj].peaks['heights'])
            if (matrix[ii][jj].y.max() >= saturation_height)  |  (matrix[ii][jj].peaks['widths'][pk_max_idx] > saturation_width):
                saturated+= [(ii,jj),]
            else:
                not_saturated+= [(ii,jj),]

    print('Ho trovato {} elementi saturati sul totale di {}\n'.format(len(saturated), n_cols*n_rows))

    return not_saturated, saturated
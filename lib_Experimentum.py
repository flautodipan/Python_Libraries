import  time
import  numpy               as np
import  pandas              as pd
from    matplotlib.pyplot   import plot
import  matplotlib.pyplot   as plt
from matplotlib_scalebar import scalebar
from    scipy.signal        import find_peaks
from    scipy.optimize      import leastsq
from    scipy.optimize      import least_squares
from    scipy.optimize      import curve_fit
from    scipy.io            import loadmat
import  json
import  os
from    os.path             import join
from    subprocess          import run

from    Models              import S_Dynamical_Form_Factor_2, S_2_Generate, S_Dynamical_Form_Factor_0, S_0_Generate, S_Dynamical_Form_Factor_0_nodelta, S_Dynamical_Form_Factor_2_nodelta     
from    Alessandria         import *
from    lmfit               import Model

free_spectral_range =   29.9702547 #GHz

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

#       PARTE   2       Classi 

class Spectrum  :

    #   Oggetto che sarà alla base di ogni spettro sperimentale

    def __init__(self, name):

        """
        Definisci il nome dell'oggetto per le stampe e il numero di parametri del tuo modello

        """
        self.name       =   name

    def __str__(self):

        return '{}'.format(self.name)
        
    def Get_Spectrum(self, y, offset, cut = False, cut_range = None, fig = False):

        """

        #  Acquisisce spettro tramite passaggio dell'array delle ordinate
            'cut' = True è per gli spettri troppo lunghi e va associata a cut_range (min, max)

        """

        self.x          =   np.arange(1, len(y)+1, 1)
        self.offset     =   offset     
        self.y          =   y - self.offset

        ###levo elementi <0.01 problematici per le incertezze poissoniane (ne escono alcuni negativi)
        ### (W la fi(si)ca)

        self.x          =   self.x[self.y >= 0.01]
        self.y          =   self.y[self.y >= 0.01]

        if cut:

            self.x      =   self.x[cut_range[0]:cut_range[1]]
            self.y      =   self.y[cut_range[0]:cut_range[1]]

        self.y_err      =   np.sqrt(self.y)

        if fig:

            plt.figure()
            plt.plot(self.x, self.y)
            plt.plot(self.x[self.peaks[0]], self.y[self.peaks[0]], '*')
            plt.xlabel('Pixels')
            plt.title('Experimental spectrum for '+self.__str__())
            plt.savefig(fig+'.png')
            plt.show()

    def Recover_Spectrum(self, x, y):

        self.x_freq = x
        self.y = y
        self.y_err = np.sqrt(self.y)
    
    def Recover_VIPA(self, x_VIPA, y_VIPA):

        self.x_VIPA_freq = x_VIPA
        self.y_VIPA = y_VIPA

    def Find_Spectrum_Peaks(self, **syg_kwargs):

        pk =  find_peaks(self.y, **syg_kwargs)

        self.peaks      =   {'idx': pk[0], 'heights': pk[1]['peak_heights'], 'widths' : pk[1]['widths'], 'prominences' : pk[1]['prominences']}       
        self.n_peaks    =   self.peaks['idx'].size
    
    def Set_Spectrum_Peaks(self, peaks):

        setattr(self, 'peaks', peaks)
        setattr(self, 'n_peaks', len(self.peaks['idx']))

    def Get_Spectrum_4_Peaks_by_Height(self):

        """
        Funzione che in ogni caso mi ritorna un dict con informazioni sui 4 picchi dello spettro, disposti in ordine 
        elastico, brillouin stockes, brillouin antistockes, elastico

        """
        self.peaks      =   Find_Highest_n_peaks(self.peaks, 4)
        self.n_peaks    =   self.peaks['idx'].size
        
    def Get_Spectrum_4_Peaks_by_Order(self):

        if self.n_peaks ==  6:
            
            self.peaks  =   Find_First_n_peaks(self.peaks, 5, exclude = [3])

        elif    (self.n_peaks == 5) | (self.n_peaks == 4):

            self.peaks  =   Find_First_n_peaks(self.peaks, 4)

        elif    self.n_peaks == 7:

            self.peaks  =   Find_First_n_peaks(self.peaks, 6, exclude = [1, 3])
        
        else:

            raise ValueError("Problema: numero di picchi non previsto dal codice: %d"%(self.n_peaks))


        self.n_peaks    =   self.peaks['idx'].size

    def Exclude_Glass_Peak(self):

        self.peaks      =   Find_First_n_peaks(self.peaks, 4, exclude = [3])
        self.n_peaks    =   self.peaks['idx'].size

    def Set_Spectrum_Nature(self, nature):
        
        setattr(self, 'nature', nature)

    def Get_VIPA_mat(self, mat_filename, path='./', tunable = None, offset = 'same', fig = False):
        
        """
        Prende dati funzione di trasferimento da file MATLAB
        OSS AGGIORNARE LEVANDO QUELLA MERDA DI H5PY e usare loadmat, soltanto capendo come 

        Parametri

        tunable         indica che la struttura è una matrice e va usato un certo tipo di presa

        offset          impostato su 'same' perchè tendenzialmente prendo prima lo spettro della VIPA, 
                        ma per altri scopi basta settarlo e si può usare indipendentemente il dato VIPA

        path            su WINDOWS bisogna stare attenti, l'argomento di default './' non è accettato

        """
    
        f                   =   h5py.File(path+mat_filename, 'r')
        VIPA                =   f.get('data')
        VIPA                =   np.array(VIPA)

        if (offset != 'same'):

            self.offset =   offset


        if (tunable is not None):

            self.nu0            =   tunable
            y_VIPA         =   VIPA[:,self.nu0] - self.offset
            self.x_VIPA         =   np.arange(1, len(self.y_VIPA)+1, 1)        

        else:

            raise ValueError ("Cojone non hai scritto il codice per file matlab non tunable, non sei de che stamo a parla'\n\n")

        if fig:
            
            plt.figure()
            plt.plot(self.x_VIPA, self.y_VIPA)
            plt.xlabel('Pixels')
            plt.ylabel('VIPA transfer function')
            plt.title('Funzione di trasferimento VIPA in pixel da %s' %(mat_filename))
            plt.savefig(fig+'.png')

    def Get_VIPA_tif(self, tif_filename, path='./', offset = 'same', **img_kwargs):

        """
        Importo il VIPA direttamente dall'immagine .tiff
        """
        VIPA                =   Import_TIF(join(path,tif_filename+'.tif'))
        
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
            plt.savefig(join(img_kwargs['save_path'],img_kwargs['fig']+'.png'))
            plt.show()

    def Check_Spectrum_Saturation(self, saturation_height = 40000, saturation_width = 15.):

        pk_max_idx  =   np.argmax(self.peaks['heights'])
        condition_peaks_height  =   (self.peaks['heights'] < 1000).all()

        if (self.y.max() >= saturation_height)  |  (self.peaks['widths'][pk_max_idx] > saturation_width):
                print('spettro saturato')
                return          1


        elif (self.n_peaks >= 4) & (self.n_peaks <= 7):

            condition_peaks_pos     =   ((self.y[self.peaks['idx'][0]] < self.y[self.peaks['idx'][1]]) | (self.y[self.peaks['idx'][3]] < self.y[self.peaks['idx'][2]]))
            
            if condition_peaks_pos:
                
                if condition_peaks_height:
                    print('Spettro con Brillouin più alti')
                    return          2
                
                else:
                    print('Spettro invisibile')
                    return  3
            
    def How_Many_Peaks_To(self, n_peaks = 4, delta = 1., treshold = 5, fig = False, verbose = False, i_know_it_is = False, **syg_kwargs):

        pk      =   find_peaks(self.y, height = treshold, **syg_kwargs)
        height  =   np.max(pk[1]['peak_heights'])/2


        
        while True:

            pk      =   find_peaks(self.y, height = height, **syg_kwargs)

            if (height > treshold) & (pk[0].size == n_peaks):

                if verbose:
                    print("Ho trovato valore dell'altezza per avere %d picchi: %f\n"%(n_peaks, height), pk)
                    _ = Analyze_Peaks(self.x_freq, self.y, 'GHz', fig = fig, verbose = verbose, height= height, **syg_kwargs)
                self.spectrum_cut_height        =   height
                break
            
            elif (height <= treshold):

                print(pk)
                print('Errore: superata altezza minima {}\nQualcosa è andato storto per il picco {}'.format(treshold, self.__str__()))
                break
            
            else: 

                height-=delta
        
        """
        #check che i picchi Brillouin siano dentro agli elastici, premesso che so che possano essere più alti degli elastici

        condition_peaks =   ((self.y[pk[0][0]] < self.y[pk[0][1]]) | ((self.y[pk[0][3]] < self.y[pk[0][2]])))
        
        if condition_peaks &  (i_know_it_is == False):
            raise ValueError("Picchi Brillouin non sono interni o sono più alti degli elastici, controllare")
        """

    
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

        self.GHz_fit_height     =   h_save[0]
        self.gauss_fit_height   =   h_save[1]
        self.VIPA_peaks_dist    =   syg_kwargs['distance']


    def Get_Spectrum_Alignment_From_Normal(self,):
        
        # supponendo che la nature normal dello spettro sia corretta, stimo alignment

        if self.peaks['heights'][0] < self.peaks['heights'][3]:
            print('La presa dati è allineata a sinistra, nel senso che il picco elastico del secondo ordine è a sinistra del principale (v. immagine VIPA)')
            return 'sx'
        elif self.peaks['heights'][0] > self.peaks['heights'][3]:
            print('La presa dati è allineata a destra, nel senso che il picco elastico del secondo ordine è a destra del principale (v. immagine VIPA)')
            return 'dx'
        else : raise ValueError('Picchi a ordini differenti uguali?\npicco sx = {}\npicco dx = {}'.format(self.peaks['heights'][0],self.peaks['heights'][3]))


    def Fit_VIPA_Gaussian(self, fig = False, verbose = False):

        peaks_idx   =   find_peaks(self.y_VIPA, height = self.gauss_fit_height, distance = self.VIPA_peaks_dist, width = 1)[0]
        print(peaks_idx)
        print(self.x_VIPA_freq[peaks_idx])
        gmod = Model(gaussian)
        result = gmod.fit(self.y_VIPA[peaks_idx], x = self.x_VIPA_freq[peaks_idx], A = 1., mu = 10, sigma = 10)
        
        A = result.values['A']
        mu = result.values['mu']
        sigma = result.values['sigma']

        #mi salvo anche i valori della gaussiana come valori iniziali per il fit

        self.p0['A']        =   [A]
        self.p0['mu']       =   [mu]
        self.p0['sigma']    =   [sigma]   

        if verbose:

            print ('Ho stimato i parametri della gaussiana come A = %3.2f\tmu  = %3.2f\tsigma = %3.2f' %(A, mu, sigma))
            print ('E li ho aggiunti ai parametri iniziali per il fit. Ora conosco %d parametri su %d \n' %(self.p0.size, 12))

        if fig:

            x = np.linspace(self.x_VIPA_freq[peaks_idx].min()*3/2, self.x_VIPA_freq[peaks_idx].max()*3/2, 1000)
            y = gaussian(x, A, mu, sigma)
            plot(self.x_VIPA_freq[peaks_idx], self.y_VIPA[peaks_idx], '*', label = 'VIPA Max')
            plot(x,y,label = 'Gauss_fit for VIPA \nA = %3.2f mu = %3.2f sigma = %3.2f' %(A, mu, sigma))
            plot(self.x_VIPA_freq, self.y_VIPA, label = 'VIPA transfer function')
            plt.plot(self.x_freq, self.y, label = 'Exp Spectrum')
            plt.legend(loc='upper center', bbox_to_anchor=(1.5, 1.05), fancybox=True, shadow=True)
            plt.xlim(-50,50)
            plt.show()
   
    def Fit_Pixel2GHz(self,  fig = False, savepath = None, **img_kwargs):

        """
        Dai dati VIPA calcolo funzione di conversione in GHz
        Parametri
        (wrap di scipy.signal.find_peaks())
        altezza         altezza minima per essere ammesso nel novero dei picchi
                        default su 0
        dist            distanza tra un picco e l'altro, onde evitare assunzioni strane
                        di default il 'min del free spectral range' in  pixel è 100
        """
 

        #1)     trovo picchi elastici dello spettro e li salvo

        peaks               =   find_peaks(self.y_VIPA, height=self.GHz_fit_height, distance = self.VIPA_peaks_dist, width = 1)
        peaks_idx           =   peaks[0]
        peaks_pix           =   self.x_VIPA[peaks_idx]
        peaks_counts        =   self.y_VIPA[peaks_idx]
        peaks_ref           =   peaks_pix[peaks_counts.argmax()]

        #2)     costruisco
        #       Delta Pixel rispetto il picco più alto
        # 
        #       Delta Ghz so che è a passi di 30
        #       faccio  in modo che passi per 0,0


        Delta_Pixel         =   peaks_pix - peaks_ref
        zero_idx            =   np.argwhere(Delta_Pixel==0.)
        Delta_GHz           =   np.arange(-zero_idx*free_spectral_range, (len(Delta_Pixel)-(zero_idx+1))*(free_spectral_range)+1, free_spectral_range) 

        #3)     faccio il fit polinomiale terzo grado e mi salvo i coefficienti

        self.Poly2GHz           =   np.polyfit(Delta_Pixel, Delta_GHz, 3)
        self.Poly2Pix           =   np.polyfit(Delta_GHz, Delta_Pixel, 3)

        if fig:

            f, ax = plt.subplots()

            x_fit = np.linspace(Delta_Pixel[0], Delta_Pixel[len(Delta_Pixel)-1], 100)
            y_fit = self.Poly2GHz[0]*(x_fit**3) + self.Poly2GHz[1]*(x_fit**2) + self.Poly2GHz[2]*x_fit + self.Poly2GHz[3]
           
            ax.plot(Delta_Pixel, Delta_GHz, '*', label='Peaks Data', c = img_kwargs['data_color'] if 'data_color' in img_kwargs else None )
            ax.plot(x_fit, y_fit, label='Polynomial fit', alpha = .8, c = img_kwargs['fit_color'] if 'fit_color' in img_kwargs else None  )
            ax.legend()
            ax.set_title('Conversion GHz vs Pixel')
            ax.xaxis.set_label_text('Delta Pixels')
            ax.yaxis.set_label_text('Delta Frequencies (GHz)')
            plt.tight_layout()
            if savepath:

                f.savefig(savepath+'GHz_conversion.pdf', format = 'pdf')
            #plt.savefig('figure/'+fig+'.png')
            plt.show()

            plt.close()
        
        
    def VIPA_Pix2GHz (self, fig = False):
            
        #modifico in GHz le ascisse della funz trasf 
        # --> prima devo convertire in DeltaPixel

        self.x_VIPA     =   self.x_VIPA - self.x_VIPA[self.y_VIPA.argmax()]
        self.x_VIPA_freq =   ((self.x_VIPA**3)*self.Poly2GHz[0])+ ((self.x_VIPA**2)*self.Poly2GHz[1]) + (self.x_VIPA*self.Poly2GHz[2]) + (self.Poly2GHz[3])

        if fig:

            fig  = plt.figure()
            plt.plot(self.x_VIPA_freq, self.y_VIPA)
            plt.title('Funzione di trasf VIPA in GHz')
            plt.xlabel('GHz')
            plt.show()
            plt.close()
    
    def Check_Brillouin_Distances(self, mean_dist_01, mean_dist_12, mean_dist_23,  stdev):

        """
        Function that is supposed to be called when peaks are already four
        --> check if the distance of Brillouin peaks is in the average 
            otherwise it means they are not brillouin
        """
        if ((self.n_peaks == 4) | ((self.n_peaks == 3) & (self.alignment == 'dx'))):
            if (np.abs(self.x[self.peaks['idx'][0]] - self.x[self.peaks['idx'][1]]) > (mean_dist_01+mean_dist_01*stdev)):
                print(1) 
                return True
            elif (np.abs(self.x[self.peaks['idx'][2]] - self.x[self.peaks['idx'][3]]) > (mean_dist_23+mean_dist_23*stdev)):
                print(2) 
                return True
            elif (np.abs(self.x[self.peaks['idx'][1]] - self.x[self.peaks['idx'][2]]) > (mean_dist_12+mean_dist_12*stdev)):
                print(3) 
                return True
            else: return False

        elif ((self.n_peaks == 3) & (self.alignment == 'sx')):
            if (np.abs(self.x[self.peaks['idx'][0]] - self.x[self.peaks['idx'][1]]) > (mean_dist_12+mean_dist_01*stdev)):
                print(4) 
                return True
            elif (np.abs(self.x[self.peaks['idx'][1]] - self.x[self.peaks['idx'][2]]) > (mean_dist_23+mean_dist_12*stdev)):
                print(5) 
                return True
            else: return False


        else: raise ValueError('Numero picchi {} non riconosciuto'.format(self.n_peaks))

    def Align_Spectrum(self, alignment):

        """
        Funzione che allinea correttamente gli spettri con gli elastici di convoluzione VIPA 
        """
        setattr(self, 'alignment', alignment)

        if self.alignment == 'dx':
            if (self.n_peaks == 4) | (self.n_peaks == 3):            
                self.x      = self.x - self.x[self.peaks['idx'][0]]
            else: raise ValueError('Il picco non ha nè 4, né 3 picchi, ma {}. Non so che fare'.format(self.n_peaks))

        elif self.alignment == 'sx':
            if self.n_peaks == 4:
                self.x      = self.x - self.x[self.peaks['idx'][3]]
            elif self.n_peaks == 3:
                self.x      = self.x - self.x[self.peaks['idx'][2]]
            else: raise ValueError('Il picco non ha nè 4, né 3 picchi, ma {}. Non so che fare'.format(self.n_peaks))

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

    def Cut_Spectrum(self, mean_dist01, mean_dist23,  factor = 0.5, ):

        """
        Esegue il taglio degli spettri, eliminando gli elastici dall'analisi.
        (mi allargo a partire dai Brillouin di un fattore pari a factor*meandist(brillouin, elastico))
        Se la natura è normal o almost, si aspetta di avere tutti e quattro i picchi
        Se la natura è brillouin_higher, si aspetta di avere 3 o 2 picchi addirittura
        Nel caso di 3, 

        """
        # 4 picchi
        if ((self.nature == 'normal') | (self.nature == 'almost_height')):

            idx_min               =   self.peaks['idx'][1] - int(mean_dist01*factor)
            idx_max               =   self.peaks['idx'][2] + int(mean_dist23*factor)

        
        elif self.nature == 'brillouin_higher':

            # 3 picchi
            if self.n_peaks == 3:
                if self.alignment == 'sx':
                    idx_min               =   self.peaks['idx'][0] - int(mean_dist01*factor)
                    idx_max               =   self.peaks['idx'][1] + int(mean_dist23*factor)
                else:
                    idx_min               =   self.peaks['idx'][1] - int(mean_dist01*factor)
                    idx_max               =   self.peaks['idx'][2] + int(mean_dist23*factor)
            # 2 picchi
            elif self.n_peaks == 2:
                    idx_min               =   self.peaks['idx'][0] - int(mean_dist01*factor)
                    idx_max               =   self.peaks['idx'][1] + int(mean_dist23*factor)

            else: raise ValueError('Numero picchi non riconosciuto come valido: {}'.format(self.n_peaks))


        else: raise ValueError('Natura dello spettro non riconosciuta: {}'.format(self.nature))

        #eseguo il taglio
        self.x_freq         =   self.x_freq[idx_min:idx_max]
        self.x              =   self.x[idx_min:idx_max]
        self.y              =   self.y[idx_min:idx_max]
        self.y_err          =   self.y_err[idx_min:idx_max]

            
    def Estimate_Spectrum_Parameters(self, model, verbose = False):

        # STIMA PARAMETRI INIZIALI della funzione teorica
        # (quelli che posso, e devo farlo  prima di tagliare)
        
        self.p0  =   pd.DataFrame({}, columns = cols if model == 'viscoelastic' else cols_dho, index = ['Values'])

        # 1)    stima dei parametri dai dati
        #       Omega come la media della posizione dei massimi brillouin
        #       Gamma (=Delta)come la media dell'ampiezza restituita da find_peaks per i Brillouin /10 (mah..)
        #       offset come la media dei valori dello spettro nella zona tra i due picchi  
        #   
        if model == 'viscoelastic':
            #visco
            self.p0['Omega']            =   [np.absolute(self.x_freq[self.peaks['idx'][2]] - self.x_freq[self.peaks['idx'][1]])*0.5]
            self.p0['Gamma']            =   [0.1]
            self.p0['Co']               =   [0.01]
            self.p0['Delta']            =   self.p0['Gamma']
            self.p0['tau']              =   [1.]            
            #delta
            self.p0['delta_position']   =   [0.]
            self.p0['delta_amplitude']  =   [1]
            self.p0['delta_width']      =   [0.1]
            #gen
            self.p0['offset']           =   [0]
            self.p0['shift']            =   [0.]
        else:
            #dho
            self.p0['dho_position']     =   [np.absolute(self.x_freq[self.peaks['idx'][2]] - self.x_freq[self.peaks['idx'][1]])*0.5]
            self.p0['dho_width']        =   [0.1]
            self.p0['dho_amplitude']    =   [0.1]
            #delta
            self.p0['delta_position']   =   [0.]
            self.p0['delta_amplitude']  =   [1]
            self.p0['delta_width']      =   [0.1]
            #gen
            self.p0['offset']           =   [0]
            self.p0['shift']            =   [0.]


    def Get_Cost_DHO(self, p0, columns):
        
        if columns == cols_dho_nodelta:
            attribute = 'Residuals_DHO_nodelta'
        elif columns == cols_dho:
            attribute = 'Residuals_DHO'
        else : raise ValueError('Choose columns for cost')

        self.cost_DHO = 0.5*np.sum(getattr(self, attribute)(p0, self.y)**2)

    def Get_cost_markov(self, p0, columns):
        
        if columns == cols_mark_nodelta:
            attribute = 'Residuals_Markov_nodelta'
        elif columns == cols_mark:
            attribute = 'Residuals_Markov'
        else : raise ValueError('Choose columns for cost')

        self.cost_markov        =   0.5*np.sum(getattr(self, attribute)(p0, self.y)**2)

    def Get_cost_tot(self, p0, p_gauss, kernel, columns):

        if columns == cols_real_nodelta:
            attribute = 'Residuals_NoGauss_nodelta'
        elif columns == cols_real:
            attribute = 'Residuals_NoGauss'
        else : raise ValueError('Choose columns for fit')

        self.cost_tot           =   0.5*np.sum(getattr(self, attribute)(p0, self.y, p_gauss, kernel)**2)
    
    
    def Get_p0_by_Markov(self, p0, treshold, **kwargs):

        """
        Funzione che dato un p0 e una treshold, stima il valore della cost function usata dall'algoritmo di fit
        -   se il valore è < treshold, quello è p0
        -   se il valore è >= treshold, esegue fit markoviano, stima i parametri e poi aggiunge
            Delta = Gamma e tau = 100.
            
        """
        
        self.Get_p0(p0, cols)
        self.cost           =   0.5*np.sum(self.Residuals(p0, self.y)**2)
        print('costo = '+str(self.cost))
        
        if self.cost >= treshold:
            
            percents        =   ('positive', 0.2, 'positive', 'positive', 'positive', 0.1, 0.1, 0.1,  np.inf, np.inf)
            self.Get_p0(self.p0[list(cols_mark)].values[0], cols_mark)
            self.Get_Fit_Bounds(percents, columns = cols_mark)
            self.Non_Linear_Least_Squares(bound = (self.bounds['down'].values, self.bounds['up'].values), **kwargs)
            self.Get_p0(np.concatenate((self.Markov_Fit_Params.T['Values'].values[:3], (self.p0['Gamma']['Values'], 1.), self.Markov_Fit_Params.T['Values'].values[3:])), cols)
            self.cost           =   0.5*np.sum(self.Residuals(self.p0.values[0], self.y)**2)

            print('costo dopo fit = '+str(self.cost))
            
            if self.res_lsq['success']:
                return  1
            else:
                return  2

        else:
            return  0


    def Initials_Parameters_from_Markov(self):
        
        if 'delta_position' not in self.Markov_Fit_Params.keys():             columns = cols_mark_nodelta
        else:   columns = cols_mark

        self.p0 =  pd.DataFrame({}, columns = cols, index = ['Values'])
        self.p0.T.Values[list(columns)] = [value for value in self.Markov_Fit_Params[list(columns)].values[0]]
        self.p0['Delta']            =   self.p0['Gamma']
        self.p0['tau']              =   [1.]


    def Gauss_Convolve_Theoretical_Response (self, p, fantoccio = False, fig = False):

        #       Funzione che convolve la mia funzione teorica con la funzione di risposta
        #       del VIPA --> è la funzione su cui dovrò fare il fit
        #       
        #
        #       fantoccio deve essere del formato doppio es. [-15, 15], (3, 3), etc...
        #       ed è fantoccio perchè allarga rispetto a self.xrange
        #
        #       OSS fa anche la divisione per la gaussiana 
        #       p è il vettore dei parametri per il fit
        #
        #       nomenclatura:   p[0] = Co
        #                       p[1] = Omega
        #                       p[2] = Gamma
        #                       p[3] = Delta
        #                       p[4] = tau
        #
        #                       p[5] = delta position
        #   param delta         p[6] = delta amplitude
        #                       p[7] = delta factor
        #
        #   param gauss         p[8] = A
        #                       p[9] = mu 
        #                       p[10] = sigma 
        #
 
        #   questi ultimi       p[11] = shift
        #   per forza           p[12] = offset       

        if (p.size != 13):
            raise ValueError("Calcola stai a usa modello Real = Markov + Exp\n p deve avere 13 elementi! Ne ha {}".format(p.size))

        if fantoccio :

            conv_range = np.linspace(fantoccio[0], fantoccio[1], 200)

        else:
            
            conv_range = self.x_freq

        self.y_convolution      =       np.zeros(conv_range.size)
        w_j_VIPA                =       np.linspace(-35,35, self.x_VIPA_freq.size)
        kernel                  =       self.Interpolate_VIPA(w_j_VIPA)
        kernel                  =       kernel/(p[8]*(np.exp(-((w_j_VIPA-p[9])**2)/(2*(p[10]**2)))))
        
        for  ii in range(len(conv_range)):

            delta_w                 =   conv_range[ii] -   w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_2(delta_w-p[11], *p[0:8])
            self.y_convolution[ii]  =   np.sum(theor*kernel)

        if fig:

                plt.figure()
                plt.plot(conv_range, self.y_convolution)
                plt.xlabel('GHz')
                plt.title('convoluzione dati VIPA con funz S2 for '+self.__str__())
                plt.show
        
 

        self.y_Gauss_convolution   =   p[12] + self.y_convolution*p[8]*np.exp(-((conv_range - p[9])**2)/(2*(p[10]**2)))

        if fig:

            plt.figure()
            plt.plot(conv_range, self.y_Gauss_convolution)
            plt.title('convoluzione dati VIPA con funz S2 più moltiplicazione Gauss')

    def Get_VIPA_for_fit(self, mode, **kwargs):

        if mode == 'natural':
                
            _ , idx_min                 =       Find_Nearest(self.x_VIPA_freq, -35.)
            _ , idx_max                 =       Find_Nearest(self.x_VIPA_freq, 35.)
            self.w_j_VIPA                    =       self.x_VIPA_freq[idx_min:idx_max]
            self.VIPA_w_j                    =       self.y_VIPA[idx_min:idx_max-1]
            self.Delta_w_j_VIPA         =       Get_Delta_Between_Array_Elements(self.w_j_VIPA)
            self.w_j_VIPA               =       self.w_j_VIPA[:-1]#-1 per azzeccare dim coi bins


        elif mode == 'interpolate':

            self.w_j_VIPA               =       np.linspace(-35., 35, kwargs['interpolation_density'])
            self.VIPA_w_j               =       self.Interpolate_VIPA(self.w_j_VIPA)             



    def Convolve_Theoretical_Response_Fast (self, p, p_gauss, kernel):

        """
        Versione per non fittare parametri della gaussiana nel fit totale
        Funziona con modello Real (Markov + Rilassamento exp) --> 12 parametri
        """

        self.y_convolution      =       np.zeros(self.x_freq.size)
        
        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   self.w_j_VIPA
            theor                   =   lorentian(delta_w, *p[5:8]) + S_Dynamical_Form_Factor_2_nodelta(delta_w-p[8], *p[0:5])
            self.y_convolution[ii]  =   np.sum(theor*kernel)

        self.y_Gauss_convolution   =   p[9] + self.y_convolution*p_gauss[0]*np.exp(-((self.x_freq - p_gauss[1])**2)/(2*(p_gauss[2]**2)))

        return self.y_Gauss_convolution
    
    
    def Gauss_Convolve_Theoretical_Response_Fast (self, p):

        """
        Versione veloce della funzione Gauss_Convolve
        Funziona con modello Real (Markov + Rilassamento exp) --> 12 parametri
        """
        self.y_convolution      =       np.zeros(self.x_freq.size)
        _ , idx_min             =       Find_Nearest(self.x_VIPA_freq, -35.)
        _ , idx_max             =       Find_Nearest(self.x_VIPA_freq, 35.)
        w_j_VIPA                =       self.x_VIPA_freq[idx_min:idx_max]#-1 per azzeccare dim coi bins
        VIPA_w_j                =       self.y_VIPA[idx_min:idx_max-1]
        Delta_w_j_VIPA          =       Get_Delta_Between_Array_Elements(w_j_VIPA)
        w_j_VIPA                =       w_j_VIPA[:w_j_VIPA.size-1]
        kernel                  =       VIPA_w_j*Delta_w_j_VIPA/(p[8]*(np.exp(-((w_j_VIPA-p[9])**2)/(2*(p[10]**2)))))
        
        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_2(delta_w-p[11], *p[0:8])
            self.y_convolution[ii]  =   np.sum(theor*kernel)
 
        self.y_Gauss_convolution   =   p[12] + self.y_convolution*p[8]*np.exp(-((self.x_freq - p[9])**2)/(2*(p[10]**2)))

        return self.y_Gauss_convolution

    def Gauss_Convolve_Markovian_Response (self, p, fantoccio = False, fig = False, compare = False, zoom = False):
        
        #       nomenclatura:   p[0] = Co
        #                       p[1] = Omega
        #                       p[2] = Gamma

        #   param delta         p[3] = delta position
        #                       p[4] = delta width
        #                       p[5] = delta amplitude
        #
        #   param gauss         p[6] = A
        #                       p[7] = mu 
        #                       p[8] = sigma 
        #
        #   questi ultimi       p[9] = shift
        #   per forza           p[10] = offset       

    
        if fantoccio :

            conv_range = np.linspace(fantoccio[0], fantoccio[1], 200)

        else:
            
            conv_range = self.x_freq

        self.y_markov_convolution      =       np.zeros(conv_range.size)
        w_j_VIPA                       =       np.linspace(-35,35, 3000)
        VIPA_w_j                        =       self.Interpolate_VIPA(w_j_VIPA)
        kernel                         =       VIPA_w_j/(p[6]*(np.exp(-((w_j_VIPA-p[7])**2)/(2*(p[8]**2)))))
        
        for  ii in range(len(conv_range)):

            delta_w                         =   conv_range[ii] -   w_j_VIPA
            theor                           =   delta_function(delta_w, *p[3:6]) + S_Dynamical_Form_Factor_0_nodelta(delta_w-p[8], *p[0:3])
            self.y_markov_convolution[ii]   =   np.sum(theor*kernel)

        self.y_Gauss_markov_convolution   =   p[10] + self.y_markov_convolution*p[6]*np.exp(-((conv_range - p[7])**2)/(2*(p[8]**2)))
        
        
        if fig:

            plt.figure()
            plt.plot(conv_range, self.y_Gauss_markov_convolution, label = 'initial_guess')
            plt.title('convoluzione dati VIPA con funz S0 più moltiplicazione Gauss')
            if compare:
                plt.plot(self.x_freq, self.y, '+', label = 'data')
            plt.legend()
            plt.show
            if zoom:
                if self.alignment == 'dx':
                    plt.xlim(3, 27)
                else:
                    plt.xlim(-27, -3)
                plt.ylim(-10, 250)

        return self.y_Gauss_markov_convolution

    def Interpolate_VIPA (self, freq):

        # funzione che interpola il valore della funzione di trasf
        # richiesta nel valore freq a partire dai dati sperim
        # funziona per valore singolo e per array di valori


        freq = np.array(freq)

        if (freq.size == 1):
        
            _ , idx             =   Find_Nearest(self.x_VIPA_freq, freq)
            x_fittino           =   np.array([self.x_VIPA_freq[idx-1], self.x_VIPA_freq[idx], self.x_VIPA_freq[idx+1]])
            y_fittino           =   np.array([self.y_VIPA[idx-1], self.y_VIPA[idx], self.y_VIPA[idx+1]])
            Parameters          =   np.polyfit(x_fittino, y_fittino, 2)

            interpolate         =   ((freq**2) *(Parameters[0])) + (freq * (Parameters[1])) + (Parameters[2])

        else:

            interpolate         =   np.zeros(np.size(freq))
            _ , idx             =   Find_Nearest_Array(self.x_VIPA_freq, freq)

            for ii in range(np.size(freq)): 

                x_fittino           =   np.array([self.x_VIPA_freq[idx[ii]-1], self.x_VIPA_freq[idx[ii]], self.x_VIPA_freq[idx[ii]+1]])
                y_fittino           =   np.array([self.y_VIPA[idx[ii]-1], self.y_VIPA[idx[ii]], self.y_VIPA[idx[ii]+1]])
                Parameters          =   np.polyfit(x_fittino, y_fittino, 2)
                interpolate[ii]     =   ((freq[ii]**2) *(Parameters[0])) + (freq[ii] * (Parameters[1])) + (Parameters[2])

        return              interpolate


    def Gauss_Convolve_Markovian_Response_Fast_nodelta (self, p):

        """
        Versione veloce della funzione Gauss_Convolve_Markovian
        Funziona con modello Markov --> 8 parametri
        """
        self.y_markov_convolution      =       np.zeros(self.x_freq.size)

        if hasattr(self, 'Delta_w_j_VIPA'):

            kernel                  =       self.VIPA_w_j*self.Delta_w_j_VIPA/(p[3]*(np.exp(-((self.w_j_VIPA-p[4])**2)/(2*(p[5]**2)))))

        else:

            kernel                  =       self.VIPA_w_j/(p[3]*(np.exp(-((self.w_j_VIPA-p[4])**2)/(2*(p[5]**2)))))


        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   self.w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_0_nodelta(delta_w-p[6], *p[0:3])
            self.y_markov_convolution[ii]  =   np.sum(theor*kernel)
 
        self.y_Gauss_markov_convolution   =   p[7] + self.y_markov_convolution*p[3]*np.exp(-((self.x_freq - p[4])**2)/(2*(p[5]**2)))

        return self.y_Gauss_markov_convolution 


    def Gauss_Convolve_Markovian_Response_Fast (self, p):

        """
        Versione veloce della funzione Gauss_Convolve_Markovian
        Funziona con modello Markov --> 11 parametri
        """
        self.y_markov_convolution      =       np.zeros(self.x_freq.size)
        if hasattr(self, 'Delta_w_j_VIPA'):

            kernel                  =       self.VIPA_w_j*self.Delta_w_j_VIPA/(p[6]*(np.exp(-((self.w_j_VIPA-p[7])**2)/(2*(p[8]**2)))))

        else:

            kernel                  =       self.VIPA_w_j/(p[6]*(np.exp(-((self.w_j_VIPA-p[7])**2)/(2*(p[8]**2)))))


        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   self.w_j_VIPA
            theor                   =   lorentian(delta_w, *p[3:6])  +  S_Dynamical_Form_Factor_0_nodelta(delta_w-p[9], *p[0:3])
            self.y_markov_convolution[ii]  =   np.sum(theor*kernel)
 
        self.y_Gauss_markov_convolution   =   p[10] + self.y_markov_convolution*p[6]*np.exp(-((self.x_freq - p[7])**2)/(2*(p[8]**2)))

        return self.y_Gauss_markov_convolution 
    

    def Gauss_Convolve_DHO_Response_Fast (self, p):

        """
        Funziona con modello DHO --> 11 parametri
        """

        self.y_DHO_convolution      =       np.zeros(self.x_freq.size)

        if hasattr(self, 'Delta_w_j_VIPA'):

            kernel                  =       self.VIPA_w_j*self.Delta_w_j_VIPA/(p[6]*(np.exp(-((self.w_j_VIPA-p[7])**2)/(2*(p[8]**2)))))

        else:

            kernel                  =       self.VIPA_w_j/(p[6]*(np.exp(-((self.w_j_VIPA-p[7])**2)/(2*(p[8]**2)))))


        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   self.w_j_VIPA
            theor                   =   lorentian(delta_w, *p[3:6])  +  DHO(delta_w-p[9], *p[0:3])
            self.y_DHO_convolution[ii]  =   np.sum(theor*kernel)
 
        self.y_Gauss_DHO_convolution   =   p[10] + self.y_DHO_convolution*p[6]*np.exp(-((self.x_freq - p[7])**2)/(2*(p[8]**2)))

        return self.y_Gauss_DHO_convolution 

    def Gauss_Convolve_Markovian_Response_Fast_nodelta (self, p):

        """
        Funziona con modello DHO nodelta --> 8 parametri
        """
        self.y_DHO_convolution      =       np.zeros(self.x_freq.size)

        if hasattr(self, 'Delta_w_j_VIPA'):

            kernel                  =       self.VIPA_w_j*self.Delta_w_j_VIPA/(p[3]*(np.exp(-((self.w_j_VIPA-p[4])**2)/(2*(p[5]**2)))))

        else:

            kernel                  =       self.VIPA_w_j/(p[3]*(np.exp(-((self.w_j_VIPA-p[4])**2)/(2*(p[5]**2)))))


        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   self.w_j_VIPA
            theor                   =   DHO(delta_w-p[6], *p[0:3])
            self.y_DHO_convolution[ii]  =   np.sum(theor*kernel)
 
        self.y_Gauss_DHO_convolution   =   p[7] + self.y_DHO_convolution*p[3]*np.exp(-((self.x_freq - p[4])**2)/(2*(p[5]**2)))

        return self.y_Gauss_DHO_convolution 

    def Residuals(self, p, y):
        
        return (self.Gauss_Convolve_Theoretical_Response_Fast(p) - y)/self.y_err

    
    def Residuals_Markov_nodelta(self, p, y):
        
        return (self.Gauss_Convolve_Markovian_Response_Fast_nodelta(p) - y)/self.y_err
    
    def Residuals_DHO(self, p, y):
        
        return (self.Gauss_Convolve_DHO_Response_Fast(p) - y)/self.y_err
    
    def Residuals_DHO_nodelta(self, p, y):
        
        return (self.Gauss_Convolve_DHO_Response_Fast_nodelta(p) - y)/self.y_err

    def Residuals_NoGauss(self, p, y, p_gauss, kernel):
        
        return (self.Convolve_Theoretical_Response_Fast(p, p_gauss, kernel) - y)/self.y_err

    def Non_Linear_Least_Squares_NoGauss (self, p_gauss, columns, p0 = 'auto', fig = False, **kwargs):

        """
        OSS le kwargs sono quelle di least_squares() di scipy.optimize
        """       
        attribute = 'Residuals_NoGauss'
 

        start                   =    time.process_time()
        kernel                  =       self.VIPA_w_j/(p_gauss[0]*(np.exp(-((self.w_j_VIPA-p_gauss[1])**2)/(2*(p_gauss[2]**2)))))


        if p0 == 'auto':

            self.res_lsq     =    least_squares(getattr(self, attribute), self.p0[list(columns)].values[0], args= ([self.y, p_gauss, kernel]),  **kwargs)    
        else:
            self.res_lsq     =    least_squares(getattr(self, attribute), p0, args= ([self.y, p_gauss, kernel]),  **kwargs)    
            
        print("s impiegati a fare il fit totale ", time.process_time()-start, '\n')

        Parameters       =    self.res_lsq.x

        try:
            J                =    self.res_lsq.jac
            cov              =    np.linalg.inv(J.T.dot(J))
            Delta_Parameters =    np.sqrt(np.diagonal(cov))

        except  np.linalg.LinAlgError as err:

            if 'Singular matrix' in str(err):
                print('Ho trovato matrice singolare')
                Delta_Parameters    =   np.empty(len(self.p0[list(columns)].values[0]))
                Delta_Parameters[:]    =   np.nan
            else :
                raise


        self.y_fit  = self.y_Gauss_convolution

        if fig:

            plt.figure()
            plt.title('Fit for '+self.__str__())
            plot(self.x_freq, self.y, '+', label='Data')
            plot(self.x_freq, self.y_fit, label= 'Fit')
            plt.legend()
        
        
        self.Tot_Fit_Params = pd.DataFrame((Parameters, Delta_Parameters, self.p0[list(columns)].values[0]), index = ('Values', 'StdErrs', 'Initials'), columns = columns)
        if (self.res_lsq['success'] == False):
            return 2
        else: return 1


    def Get_Fit_Bounds(self, rules, columns): 

        """
        Funzione che definisce un attributo self.bounds, un dataframe con le colonne up and down
        ogni valore di self.p0 ha due elementi associati che corrispondono ai limiti che si vogliono imporre nel fit

        il calcolo dei bounds si basa sulla funzione Get_Around e sul parametro da passare percents, il quale DEVE essere di lunghezza uguale al 
        p0 di riferimento, con quell'ordine di paramentri

        se un elemento di percents è np.inf, la funzione ritorna l'asse reale come limiti per il parametro corrispondente
        se un elemento di percents è 'positive', la funzione ritorna [0, np.inf] come limiti #

        la funzione funziona sia per Markov che per Real

        """

        if len(rules) != len(columns):

            raise ValueError ("Lunghezza della tupla percents = %d errata, deve essere uguale alla dimensione delle colonne passate = %d"%(len(rules), len(columns)))
                
        self.bounds = pd.DataFrame( {}, index= columns, columns=('down', 'up'))

        for (col, (_, rule)) in zip(columns, rules.items()):
    
            if  rule    ==  'inf':

                self.bounds._set_value(col, 'down', -np.inf)   
                self.bounds._set_value(col, 'up', +np.inf)   
            
            elif    rule == 'positive':

                self.bounds._set_value(col, 'down', 0)   
                self.bounds._set_value(col, 'up', +np.inf) 

            elif (type(rule) == list):
                
                self.bounds._set_value(col, 'down', rule[0])   
                self.bounds._set_value(col, 'up', rule[1]) 

            else:
                    
                bound   =   Get_Around(self.p0[col]['Values'], rule)
                self.bounds._set_value(col, 'down', bound[0])   
                self.bounds._set_value(col, 'up', bound[1]) 


    def Residuals_Markov(self, p, y):
        
        return (self.Gauss_Convolve_Markovian_Response_Fast(p) - y)/self.y_err

    
    def Non_Linear_Least_Squares (self, columns, p0 = 'auto', fig = False, zoom = False, **kwargs):
        
        
        if columns == cols_mark_nodelta:
            attribute = 'Residuals_Markov_nodelta'
        elif columns == cols_mark:
            attribute = 'Residuals_Markov'
        elif columns == cols_dho:
            attribute = 'Residuals_DHO'
        elif columns == cols_dho_nodelta:
            attribute = 'Residuals_DHO_nodelta'

        else : raise ValueError('Choose columns for fit')

        start            =    time.process_time()

        if p0 == 'auto':
            
            self.res_lsq     =    least_squares(getattr(self, attribute), x0 = self.p0[list(columns)].values[0], args= ([self.y]),  **kwargs)    
        else:
            self.res_lsq     =    least_squares(getattr(self, attribute), p0, args= ([self.y]),  **kwargs)    

        print("s impiegati a fare il fit ", time.process_time()-start, '\n')
        Parameters       =    self.res_lsq.x
        
        try:
            J                =    self.res_lsq.jac
            cov              =    np.linalg.inv(J.T.dot(J))
            Delta_Parameters =    np.sqrt(np.diagonal(cov))

        except  np.linalg.LinAlgError as err:

            if 'Singular matrix' in str(err):
                print('Ho trovato matrice singolare')
                Delta_Parameters        =   np.empty(len(self.p0[list(columns)].values[0]))
                Delta_Parameters[:]     =   np.nan
            else :
                raise

        self.y_markov_fit  = self.y_Gauss_markov_convolution
      
        if fig:
            plt.figure()
            plt.title('Fit for '+self.__str__())
            plot(self.x_freq, self.y, '+', label='Data')
            plot(self.x_freq, self.y_markov_fit, label= 'Fit')
            if zoom:
        
                if self.alignment == 'dx':
                    plt.xlim(3, 27)
                else:
                    plt.xlim(-27, -3)
                plt.ylim(-10, 250)

            plt.legend()
           
        self.Markov_Fit_Params = pd.DataFrame((Parameters, Delta_Parameters, self.p0[list(columns)].values[0]), index = ('Values', 'StdErrs', 'Initials'), columns = columns)
        
        if (self.res_lsq['success'] == False):
            return 2
        else: return 1
    
    def Non_Linear_Least_Squares_DHO (self, columns, p0 = 'auto', fig = False, zoom = False, **kwargs):
        
        
        if columns == cols_dho_nodelta:
            attribute = 'Residuals_DHO_nodelta'
        elif columns == cols_dho:
            attribute = 'Residuals_DHO'

        else : raise ValueError('Choose columns for fit')

        start            =    time.process_time()

        if p0 == 'auto':
            
            self.res_lsq     =    least_squares(getattr(self, attribute), x0 = self.p0[list(columns)].values[0], args= ([self.y]),  **kwargs)    
        else:
            self.res_lsq     =    least_squares(getattr(self, attribute), p0, args= ([self.y]),  **kwargs)    

        print("s impiegati a fare il fit ", time.process_time()-start, '\n')
        Parameters       =    self.res_lsq.x
        
        try:
            J                =    self.res_lsq.jac
            cov              =    np.linalg.inv(J.T.dot(J))
            Delta_Parameters =    np.sqrt(np.diagonal(cov))

        except  np.linalg.LinAlgError as err:

            if 'Singular matrix' in str(err):
                print('Ho trovato matrice singolare')
                Delta_Parameters        =   np.empty(len(self.p0[list(columns)].values[0]))
                Delta_Parameters[:]     =   np.nan
            else :
                raise

        self.y_DHO_fit  = self.y_Gauss_DHO_convolution
      
        if fig:
            plt.figure()
            plt.title('Fit for '+self.__str__())
            plot(self.x_freq, self.y, '+', label='Data')
            plot(self.x_freq, self.y_DHO_fit, label= 'Fit')
            if zoom:
        
                if self.alignment == 'dx':
                    plt.xlim(3, 27)
                else:
                    plt.xlim(-27, -3)
                plt.ylim(-10, 250)

            plt.legend()
           
        self.DHO_Fit_Params = pd.DataFrame((Parameters, Delta_Parameters, self.p0[list(columns)].values[0]), index = ('Values', 'StdErrs', 'Initials'), columns = columns)
        
        if (self.res_lsq['success'] == False):
            return 2
        else: return 1
    
    def Recover_Markov_Fit_Params(self, dictio_string):

        """
        Funzione che tramite libreria json ricostruisce DataFrame Markov_Fit_Params per l'oggetto
        dictio_string contiene sotto forma di stringa le info in formato json (che sono un dizionario)

        """
        self.Markov_Fit_Params =   pd.DataFrame(json.loads(dictio_string))

    def Recover_y_markov_fit(self, y_markov_fit):

        self.y_markov_fit = y_markov_fit

    def Recover_y_fit(self, y_tot_fit):

        self.y_fit = y_tot_fit

    def Recover_Tot_Fit_Params(self, dictio_string):

        self.Tot_Fit_Params =   pd.DataFrame(json.loads(dictio_string))

    def Recover_Gauss_Parameter(self, dictio_string):

        df              =   pd.DataFrame(json.loads(dictio_string))

        self.Get_p_gauss(df[list(cols_gauss)].values[0])

    def Get_p0(self,p0, columns):

        if not hasattr(self, 'p0'): 
            self.p0  =   pd.DataFrame({}, columns = columns, index = ['Values'])
        self.p0.T.Values[list(columns)] = [value for value in p0]

        #self.p0      =   pd.DataFrame({idx : value for (idx, value) in zip(columns, p0)}, index = ['Values'])

    def Get_p_gauss(self, p_gauss):

        self.p_gauss =   pd.DataFrame({idx : value for (idx, value) in zip(cols_gauss, p_gauss)}, index = ['Values'])
   
    def Recover_cost_markov(self, cost):
        self.cost_markov = cost

    def Recover_cost_tot(self, cost):
        self.cost_tot = cost

    def Get_Best_p0(self, p0s, columns, model):

        costs = np.zeros(len(p0s))

        for (p0, kk) in zip(p0s, range(costs.size)):
            if model == 'viscoelastic':
                self.Get_cost_markov(p0, columns)
                costs[kk] = self.cost_markov
            elif model == 'DHO':
                self.Get_cost_DHO(p0, columns)
                costs[kk] = self.cost_DHO

        self.Get_p0(p0s[np.argmin(costs)], columns)

def Initialize_Matrix(ii_0, jj_0, ii_stop, jj_stop):

    matrix = ()
    
    rows = np.arange(ii_0, ii_stop, 1)
    cols = np.arange(jj_0, jj_stop, 1)
    
    for ii in rows:
        riga = ()
        for jj in cols:
            riga  = riga + (Spectrum('Element ('+str(ii)+','+str(jj)+')'),)
        matrix = matrix + (riga,)

    print('Initialized %dx%d matrix, for a total number of spectra %d'%(len(matrix), len(matrix[0]), len(matrix)*len(matrix[0])  ))

    return (matrix, rows, cols)

def Get_Isolated_Elements(excluded):

    """

    Essendo il fit basato su prendere condizioni iniziali da elemento a sx per la prima riga, elemento sopra per tutti gli altri
    Questa funzione ritorna una tupla con tutti gli (ii,jj) affetti da mancanza di questo vicino

    (0,0) isolato per def -> di conseguenza 

    """

    isolated = ((0,0),)

    for (ii,jj) in excluded:

        if (ii == 0):

            isolated    =   isolated    +   ((ii, jj+1),)
            isolated    =   isolated    +   ((ii+1, jj),)

        else :

            isolated    =   isolated    +   ((ii+1, jj),)

    return isolated
    
def Unpack_Fit(fit):

    non_fitted   =   ()
    accomplished =   ()
    exceded      =   ()
    fitted       =   ()
    excluded     =   ()

    for (what, (ii,jj)) in fit:

        if  what == 0:

            non_fitted  =    non_fitted +   ((ii,jj),)

        elif what == 1:
            
            accomplished      =   accomplished  +   ((ii,jj),)
            fitted            =   fitted + ((ii,jj),)

        elif what == 2:

            exceded      =   exceded  +   ((ii,jj),)
            fitted            =   fitted + ((ii,jj),)

        elif what == 3:

            excluded      =  excluded + ((ii,jj),)
    
    return (non_fitted, accomplished, exceded, fitted)

def Get_Fit_Map(n_rows, n_cols, non_fitted, exceded, excluded, fig = False, path = ''):

    """
    Funzione che ritorna e stampa matrice i cui valori sono

    0   se il fit non è stato svolto per quel (ii,jj) perchè già buona stima (solo per fit Markov)

    1   se il fit è stato fatto ed è stato portato a buon termine

    2   se il fit è stato fatto ma ha superato numero max di iterazioni senza convergere

    3   se il fit non è stato svolto perchè elemento escluso

    OSS: essendo che mi aspetto max di fit fatti, genero np.ones()


    """
    fit_map =  np.ones((n_rows, n_cols))

    for ii in range(n_rows):
        for jj in range(n_cols):

            if (ii,jj) in non_fitted:

                fit_map[ii,jj]  =   0.
            
            elif (ii, jj) in exceded:

                fit_map[ii,jj]  =   2.
            
            elif (ii,jj) in excluded:

                fit_map[ii,jj]  =   3.
    
    print('Completata Fit_Map')
    
    if fig:

        c = plt.get_cmap('Pastel1', 3)
        plt.matshow(fit_map, cmap = c)
        plt.title('Fit Map')
        plt.colorbar()
        plt.xlabel('Row Index')
        plt.ylabel('Col Idx')
        plt.savefig(path + fig+'.png')
        plt.close()
    
    return fit_map

def Verify_Fit(matrix, boni, fit, parameter, where, treshold):

    """

    Funzione che in base agli argomenti where e treshold verifica se il fit scelto per il parametro parameter ha un valore minore/maggiore di 
    treshold 
    In caso affermativo, plotta il fit

    """
    if fit == 'markov':
        fit_params = 'Markov_Fit_Params'
        y_fit = 'y_markov_fit'
    elif fit == 'tot':
        fit_params = 'Tot_Fit_Params'
        y_fit = 'y_fit'
    else:
        raise ValueError ("Specify which fit to verify: 'markov' or 'tot \n")

    for (ii, jj) in boni:

        if where == 'lower':

            condition = getattr(matrix[ii][jj], fit_params)[parameter]['Values'] < treshold

        elif where == 'upper':

            condition = getattr(matrix[ii][jj], fit_params)[parameter]['Values'] > treshold

        else: raise ValueError('Uncorrect where selection: type "upper" to verify which values are above treshold\n"lower to below\n ')
            
        if condition:

            plt.figure()
            plt.plot(matrix[ii][jj].x_freq, matrix[ii][jj].y, '+', label = data)
            plt.plot(matrix[ii][jj].x_freq, getattr(matrix[ii][jj], y_fit), label = y_fit)
            plt.title(str((ii,jj)))


def Get_Parameter_Map(fit, parameter, matrix, n_rows, n_cols, excluded, Deltas = False):


    if fit == 'markov':
        fit_params = 'Markov_Fit_Params'
    elif fit == 'tot':
        fit_params = 'Tot_Fit_Params'
    else:
        raise ValueError ("Specify which fit to verify: 'markov' or 'tot \n")
    if Deltas:
        value = 'StdErrs'
    else:
        value = 'Values'

    p_map   =   np.zeros((n_rows, n_cols))
    nans = ()
    for ii in range(n_rows):
        for jj in range (n_cols):
            if  (ii, jj) not in excluded:
                p_map[ii,jj]    =   getattr(matrix[ii][jj], fit_params)[parameter][value]
                if Deltas:
                    #rendo percentuale
                    p_map[ii,jj] = (100*p_map[ii,jj])/getattr(matrix[ii][jj], fit_params)[parameter]['Values']
            else:
                p_map[ii,jj]    = np.nan
                nans = nans +((ii,jj),)


    print('Completata Parameter_Map per '+parameter)
    print('Ho trovato {} elementi saturati'.format(len(nans)))

    return  (p_map, nans)

def Interpolate_Parameter_Map(p_map, parameter, fit, matrix, n_rows, n_cols, inf, sup):

    if fit == 'markov':
        parameters = 'Markov_Fit_Params'           
    elif fit == 'tot':
        parameters = 'Tot_Fit_Params'

    new_map = np.zeros((n_rows, n_cols))
    for ii in range(p_map.shape[0]):
        for jj in range (p_map.shape[1]): 
            
         
            #print(str((ii,jj)))
            
            neigh = Get_Neighbours2D(ii,jj,n_rows, n_cols)
            neigh_params = [getattr(matrix[kk][ll], parameters).T.Values[parameter] for kk, ll in neigh if hasattr(matrix[kk][ll], parameters)]
            #print(neigh_params)
            ave = []
            count = 0
            for (mm,nn) in neigh:

                is_nan = np.isnan(p_map[mm,nn])
                is_high = p_map[mm,nn] > sup 
                is_low = p_map[mm,nn] < inf
                if hasattr(matrix[mm][nn], parameters):
                    is_different = Is_Far_by_N_sigma_from_Mean(p_map[mm,nn], neigh_params, N = 1.5)
                    if is_different : count+=1
                else: is_different = True

                if not (is_nan | is_high | is_low | is_different):
                    ave.append(p_map[mm, nn])
 
            new_map[ii,jj] = np.average(ave,) if len(ave) != 0 else np.nanmean(neigh_params)

            #print('Su {} ho {} outlier'.format(len(neigh_params), count))


    print('Completata Interpolazione per elementi di {} {} map \n'.format(fit, parameter))


    return  new_map

def Print_Parameter_Map(p_map, inf, sup, parameter, fit, name, pix_scale, filename, path = './', **kwargs):

    if parameter == 'Omega':
        title = 'Brillouin shift $\Omega$ {} Map\n({})'.format(fit, name)
    elif parameter == 'Gamma':
        title = 'Brillouin width $\Gamma$ {} Map\n({})'.format(fit, name)
    elif parameter == 'Tau':
        title = 'Constant of relaxation tau {} Map\n({})'.format(fit, name)
    elif parameter == 'Delta':
         title = 'Relaxation width $\Delta$ {} Map\n({})'.format(fit, name)
    else: raise ValueError('Specifica parametro')

    f, ax = plt.subplots()

    img = ax.matshow(p_map, cmap = 'jet')
    cbar = plt.colorbar(img, ax = ax, label = 'GHz', orientation = 'vertical')
    img.set_clim(inf, sup)

    ax.set_title(title, fontsize = 14)
    ax.set_xlabel('Col Index')
    ax.set_ylabel('Row Idx')
    ax.tick_params('x', bottom = True, top = True, labelbottom = True, labeltop = False)
    ax.tick_params('y', left = True, right = True)

    bar = scalebar.ScaleBar(300, units = 'nm', location = kwargs['bar_loc'] if 'bar_loc' in kwargs else 'lower right', length_fraction=0.1, height_fraction=0.005, font_properties = {'size' : 10}, color = kwargs['bar_color'] if 'bar_color' in kwargs else 'white', frameon = kwargs['frameon'] if 'frameone' in kwargs else False)
    plt.gca().add_artist(bar)

    plt.tight_layout()
    f.savefig(os.path.join(path, filename+'.pdf'), format = 'pdf', bbox_to_inches = (0,0,1,1))

    plt.show()


def Escludi_a_Mano(to_add, excluded):

    for (ii, jj) in to_add:

        excluded = excluded + [(ii,jj),]

    excluded.sort()
    return excluded


def Whose_Param_Too_High(param, treshold, fit,  matrix, fitted, verbose = False):


    if fit == 'markov':
        attr = 'Markov_Fit_Params'
    elif fit == 'tot':
        attr = 'Tot_Fit_Params'
    else:raise ValueError("Select a type of fit: 'markov' or 'tot' ")
    too_high    =   ()

    for (ii,jj) in (fitted):
        if getattr(matrix[ii][jj],attr)[param]['Values'] > treshold:
            too_high    =   too_high    +   ((ii,jj),)
            if verbose:
                print(str((ii,jj))+' ha '+param+'= %3.2f'%(getattr(matrix[ii][jj],attr)[param]['Values']))

    return too_high

def Whose_Param_Too_Low(param, treshold, fit, matrix, fitted, verbose = False):

    if fit == 'markov':
        attr = 'Markov_Fit_Params'
    elif fit == 'tot':
        attr = 'Tot_Fit_Params'
    else:raise ValueError("Select a type of fit: 'markov' or 'tot' ")

    too_low    =   ()

    for (ii,jj) in (fitted):
        if getattr(matrix[ii][jj],attr)[param]['Values'] <= treshold:
            too_low    =   too_low    +   ((ii,jj),)
            if verbose:
                print(str((ii,jj))+' ha '+param+'= %3.2f'%(getattr(matrix[ii][jj],attr)[param]['Values'][param]['Values']))

    return too_low


def Save_Fitted_Info(path, n_rows, n_cols, fitted):

    with open(path+'fitted.txt', 'w') as f:
            
        f.write('(')
        for ii in range(0, n_rows,1):
            for jj in range(n_cols):
                if (ii,jj) in fitted:
                    f.write(str((ii,jj))+',')
        f.write(')')
    print('Stampato fitted info su file')

def Save_Fit_Info(fit, filename = 'fit.txt', path = './'):

    with open(path+filename, 'w') as f:
        f.write('(')
        for ft in fit:
                    f.write(str(ft)+',')
        f.write(')')
    print('Stampato fitted info su file '+path+filename)


def Verify_Initial_Conditions(matrix, ver = (), init = ()):
    
    plt.title('Verify Initial Conditions for '+str((ver[0],ver[1]))+'with '+ str((init[0],init[1]))+ 'parameters')
    plt.plot(matrix[ver[0]][ver[1]].x_freq, matrix[ver[0]][ver[1]].y, '+', label = 'dati')
    plt.plot(matrix[ver[0]][ver[1]].x_freq, matrix[init[0]][init[1]].Gauss_Convolve_Theoretical_Response_Fast(matrix[init[0]][init[1]].p0.values[0]))

    print("Il parametro è :", matrix[init[0]][init[1]].p0.values[0])

def Save_Markov_Fit_Parameters(matrix, fitted, out_filename = 'markov_fit_params.txt' , path = './'):

    """
    Salvo nella cartella di analisi un file di nome 'markov_fit_params.txt' nella cartella di analisi associata alla presa dati
    la struttura è 
    ogni riga contiene un dizionario già strutturato per file .json che DataFrame è in grado di acquisire

    FONDAMENTALE è l'ordine delle righe che corrisponde a quello di fitted()

    """

    with open(path+out_filename, 'w') as f_out:

        for (ii,jj) in fitted:
   
            f_out.write(json.dumps(matrix[ii][jj].Markov_Fit_Params.to_dict())+'\n')

    print('Stampato parametri fit su file '+path+out_filename)

def Save_Tot_Fit_Parameters(matrix, fitted, out_filename = 'tot_fit_params.txt' , path = './'):

    with open(path+out_filename, 'w') as f_out:

        for (ii,jj) in fitted:
   
            f_out.write(json.dumps(matrix[ii][jj].Tot_Fit_Params.to_dict())+'\n')

    print('Stampato parametri fit su file '+path+out_filename)

def Save_XY_VIPA(x,y, out_filename = 'xy_VIPA.txt' , path = './'):

    with open(path+out_filename, 'w') as f_out:
                f_out.write('# x,y of VIPA in couples of lines\n')
                f_out.write(np.array2string(x, max_line_width = 100000)+'\n')
                f_out.write(np.array2string(y, max_line_width = 100000)+'\n')


def Save_XY_position(matrix, n_rows, n_cols, out_filename = 'xy.txt' , path = './'):

    with open(path+out_filename, 'w') as f_out:
        f_out.write('# x,y of each spectra in couples of lines\n')
        for (ii,jj) in serpentine_range(n_rows, n_cols, 'right'):
                f_out.write(np.array2string(matrix[ii][jj].x_freq, max_line_width = 10000)+'\n')
                f_out.write(np.array2string(matrix[ii][jj].y, max_line_width = 10000)+'\n')

def Save_y_markov_fit(matrix,fitted, out_filename = 'y_markov_fit.txt', path = './'):

    with open(path+out_filename, 'w') as f_out:
        for (ii,jj) in fitted:
            f_out.write(np.array2string(matrix[ii][jj].y_markov_fit, max_line_width = 100000)+'\n')


def Save_y_fit(matrix, fitted, out_filename = 'y_tot_fit.txt', path = './'):

    with open(path+out_filename, 'w') as f_out:
        for (ii,jj) in fitted:
            f_out.write(np.array2string(matrix[ii][jj].y_fit, max_line_width = 10000)+'\n')

def Save_cost_markov(matrix, fitted, out_filename = 'cost_markov.txt', path = './'):

    with open(path+out_filename, 'w') as f_out:
        for (ii,jj) in fitted:
            f_out.write(np.array2string(matrix[ii][jj].cost_markov, max_line_width = 10000)+'\n')

def Save_cost_tot(matrix, fitted, out_filename = 'cost_tot.txt', path = './'):

    with open(path+out_filename, 'w') as f_out:
        for (ii,jj) in fitted:
            f_out.write(np.array2string(matrix[ii][jj].cost_tot, max_line_width = 10000)+'\n')

def Plot_Elements_Spectrum(matrix, elements_iterable, fit = False, pix = False, peaks = False, x_range = (), y_range = ()):

    if pix:

        attribute = 'x'
        
    else:
        attribute = 'x_freq'

    for (ii,jj) in elements_iterable:
        
        print(str((ii,jj)))
        plt.figure()
        plt.plot(getattr(matrix[ii][jj], attribute), matrix[ii][jj].y, '.', label = 'data')

        if fit:
                
            if fit == 'markov':
                plt.plot(getattr(matrix[ii][jj], attribute), matrix[ii][jj].y_markov_fit, label = 'markov fit')
                print(matrix[ii][jj].Markov_Fit_Params)

            elif fit == 'tot':
                plt.plot(getattr(matrix[ii][jj], attribute), matrix[ii][jj].y_fit, label = 'tot fit')
                print(matrix[ii][jj].Tot_Fit_Params)
            
            elif fit == 'dho':
                plt.plot(getattr(matrix[ii][jj], attribute), matrix[ii][jj].y_DHO_fit, label = 'DHO fit')
                print(matrix[ii][jj].DHO_Fit_Params)
            
        if peaks:
            #anche se non funziona con x_freq i picchi, o forse sì?
            plt.plot(getattr(matrix[ii][jj], attribute)[matrix[ii][jj].peaks['idx']], matrix[ii][jj].y[matrix[ii][jj].peaks['idx']], '+', label = 'peaks')
        
        
        if x_range:
            plt.xlim(x_range[0], x_range[1])
        if y_range:
            plt.ylim(y_range[0], y_range[1])
        plt.title(str((ii,jj)))
        plt.legend()
        plt.show()

def Get_cost_map(matrix, fit, n_rows, n_cols, fig, inf = 0, sup = 1000, cmap = 'seismic', path = './'):

    if fit == 'markov':
        attr = 'cost_markov'
    elif fit == 'tot':
        attr = 'cost_tot'

    cost_matrix = np.zeros((n_rows, n_cols))

    for ii in range(n_rows):
        for jj in range(n_cols):
            if hasattr(matrix[ii][jj], attr):
                cost_matrix[ii,jj] = getattr(matrix[ii][jj], attr)
            else:
                cost_matrix[ii,jj] = np.nan
    
    

    cm = plt.get_cmap(cmap)
    cm.set_bad(color='k')
    plt.matshow(cost_matrix, cmap = cmap)
    plt.clim(inf, sup)
    plt.title('Cost Map for '+attr)
    plt.colorbar()
    plt.xlabel('Row Index')
    plt.ylabel('Col Idx')
    plt.savefig(os.path.join(path, fig+'.pdf'), format = 'pdf')
    plt.show()

    return cost_matrix


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

def Check_Peaks_Number(matrix, iterable_elements, n, v = False, fig = False):

    n_peaks = ()

    for (ii,jj) in iterable_elements:
        if matrix[ii][jj].n_peaks == n:
            n_peaks += ((ii,jj),)
            print('Spettro {} ha {} picchi!\n'.format(str((ii,jj)), n, ))
            if v:
                print(matrix[ii][jj].peaks)
            if fig:
                plt.figure()
                plt.title(str((ii,jj)))
                plt.plot(matrix[ii][jj].x, matrix[ii][jj].y, label = 'data')
                plt.plot(matrix[ii][jj].x[matrix[ii][jj].peaks['idx']], matrix[ii][jj].y[matrix[ii][jj].peaks['idx']], '+', label = 'peaks')
                plt.legend()
                plt.show()
                plt.close()
    print('Ho trovato {} spettri con {} picchi\n'.format(len(n_peaks), n))
    return n_peaks

def Get_Bad_Elements(matrix, iterable, treshold, fit):

    if fit == 'markov':
        cost = 'cost_markov'
    elif fit == 'tot':
        cost = 'cost_tot'

    too_bad = ()
    
    for (ii,jj) in iterable:
        if getattr(matrix[ii][jj],cost) > treshold:

            too_bad += ((ii,jj),)

    print('I found {} bad elements out of {}\n'.format(len(too_bad), len(iterable)))
    return too_bad

def Get_Good_Elements(matrix, iterable, treshold, fit ):

    if fit == 'markov':
        cost = 'cost_markov'
    elif fit == 'tot':
        cost = 'cost_tot'

    too_good = ()
    for (ii,jj) in iterable:
        if getattr(matrix[ii][jj],cost) < treshold:

            too_good += ((ii,jj),)

    print('I found {} good elements out of {}\n'.format(len(too_good), len(iterable)))
    return too_good


    
def Get_p0_by_Neighbours(matrix, columns, ii_0, jj_0, n_rows, n_cols):

    p0s = []
    neigh = Get_Neighbours2D(ii_0, jj_0, n_rows, n_cols)
    for (ii,jj) in neigh:
        if hasattr(matrix[ii][jj], 'Markov_Fit_Params'):
            if True not in np.isnan(matrix[ii][jj].Markov_Fit_Params.T['StdErrs'].values):
                condition_delta = ('delta_position' in getattr(matrix[ii][jj], 'Markov_Fit_Params').keys())
                if (columns == cols_mark_nodelta) | ((columns == cols_mark) & condition_delta) :
                    p0s.append(matrix[ii][jj].Markov_Fit_Params.T['Values'][list(columns)])
                else: continue
            else: continue

    return p0s

def serpentine_range(n_rows, n_cols, start):

    new_boni = []
    
    for ii in range(n_rows):
        if start == 'right':
            for jj in np.arange(n_cols-1, -1, -1):
                new_boni.append((ii,jj))
            start = 'left'
            continue
            
        elif start == 'left':
            for jj in np.arange(0, n_cols, 1):
                new_boni.append((ii,jj))
            start = 'right'
            continue
            
            
    return new_boni



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


def Check_Settings_From_Terminal(recover_markov, skip_tot, exclude_delta, fit_algorithm):

    print('I will show the most important operatives inputs I took from config.ini, just to verify their correctness\n')
    delay = 3.

    if recover_markov:
        print('You decided to recover markov fit from {} \nIt is correct? Enter "ok" if so, any other key to change this option'.format(analysis_path))  
    else:
        print('You will perform markov fit. Enter "ok" to continue, any other key to modify this opt\n')
    
    a = input()

    if a == 'ok':
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
            print('You will skip fit tot')
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
    
    print('Inserted fit algorithm: {}, that means {}\n Enter ok to confirm, any other to switch to the other'.format(fit_algorithm, 'Trust Reflective Region' if fit_algorithm == 'trf' else 'Levenberg-Marquardart'))


    a = input()
    if a == 'ok':
        print('Ok, you confirmed {}'.format(fit_algorithm))
        time.sleep(delay)
    else:
        print('You choose to switch to {}'.format('Trust Reflective Region' if fit_algorithm != 'trf' else 'Levenberg-Marquardart'))
        fit_algorithm = 'trf' if fit_algorithm != 'trf' else 'lm'
        time.sleep(delay)

    print("Now, I am beginning to analyze data. You won't do anything from now on")
    time.sleep(delay)
    return recover_markov, skip_tot, exclude_delta, fit_algorithm

def Check_Execution_Mode(argv):

# argv[1] = -f è Jupyter su Linux, su windows ho un ipylauncher in argv[0]

    return 'terminal' if ((argv[1] != '-f') & ('ipykernel_launcher' not in argv[0])) else 'interactive'

def Find_Problems(four, which, what_is, value):

    if what_is == 'greater':
        where = np.where(np.array(which) > value)
    elif what_is == 'lesser':
        where = np.where(np.array(which) < value)
    else: raise ValueError('You can insert "greater" or "lesser"\n')

    where_is = []
    count = 0
    for (ii,jj) in four:
        if count in where[0]:
            where_is.append((ii,jj),)
        count+=1
    print('Ho trovato {} elementi che soddisfano la condizione inserita {} di {}\n'.format(len(where_is), what_is, value))
    return where_is
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

from    Models              import S_Dynamical_Form_Factor_2, S_2_Generate, S_Dynamical_Form_Factor_0, S_0_Generate        
from    Alessandria         import *
from    lmfit               import Model

free_spectral_range =   29.9702547 #GHz
cols      = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')


#       PARTE   2       Classi 


class Spectrum  :

    #   Oggetto che sarà alla base di ogni spettro sperimentale

    def __init__(self, name, n_params = 12):

        """
        Definisci il nome dell'oggetto per le stampe e il numero di parametri del tuo modello

        """
        self.name       =   name
        self.n_params   =   n_params
        self.saturation =   False
        self.Brillouin_HIgher   = False

    def __str__(self):

        return '{}'.format(self.name)
        
    def Get_Spectrum(self, y, offset = 183., cut = False, cut_range = None, fig = False):

        """

        #  Acquisisce spettro : o tramite passaggio


            'cut' = True è per gli spettri troppo lunghi e va associata a cut_range (min, max)

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

        self.peaks      =   find_peaks(self.y, **syg_kwargs)
        self.n_peaks    =   self.peaks[0].size


    def Get_Spectrum_4_Peaks_by_Height(self):

        """
        Funzione che in ogni caso mi ritorna un dict con informazioni sui 4 picchi dello spettro, disposti in ordine 
        elastico, brillouin stockes, brillouin antistockes, elastico

        """
        self.peaks      =   Find_Highest_n_peaks(self.peaks, 4)
        self.n_peaks    =   self.peaks['peaks_idx'].size
        
    def Get_Spectrum_4_Peaks_by_Order(self):

        if self.n_peaks ==  6:
            
            self.peaks  =   Find_First_n_peaks(self.peaks, 5, exclude = [3])

        elif    self.n_peaks == 5:

            self.peaks  =   Find_First_n_peaks(self.peaks, 4)

        elif    self.n_peaks == 7:

            self.peaks  =   Find_First_n_peaks(self.peaks, 6, exclude = [1, 3])

        else:

            raise ValueError("Problema: numero di picchi non previsto dal codice %d"%(self.n_peaks))


        self.n_peaks    =   self.peaks['peaks_idx'].size

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

    def Get_VIPA_tif(self, tif_filename, path='', offset = 'same', **img_kwargs):

        print ('ATTENZIONE funzione da aggiornare\n\n\n\nATTENZIONE guarda Get_VIPA_mat')

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

    def Check_Spectrum(self, saturation_height = 40000, saturation_width = 15.):

        if (self.n_peaks <= 3) | (self.n_peaks >= 7):
            
            print('Spettro invisibile')

            return  3

        elif (self.n_peaks >= 4) & (self.n_peaks < 7):


            pk_max_idx  =   np.argmax(self.peaks[1]['peak_heights'])
            condition_peaks_pos     =   ((self.y[self.peaks[0][0]] < self.y[self.peaks[0][1]]) | (self.y[self.peaks[0][3]] < self.y[self.peaks[0][2]]))
            condition_peaks_height  =   (self.peaks[1]['peak_heights'] < 1000).all()

            if (self.y.max() >= saturation_height)  &  (self.peaks[1]['widths'][pk_max_idx] > saturation_width):
                
                print('spettro saturato')
                self.saturation =   True
                return          1

            if condition_peaks_pos:
                
                if condition_peaks_height:
                        
                    print('Spettro con Brillouin più alti')
                    self.Brillouin_Higher   = True
                    return          2
                
                else:

                    print('Spettro invisibile')
                    return  3

        


    def Knowing_Spectrum_Peaks(self):

        self.peaks      =   Find_Highest_n_peaks(self.y, 4, height = 10., distance = 20, width = 1.)

    def Knowing_VIPA_Peaks(self):

        self.VIPA_peaks =   Find_Highest_n_peaks(self.y_VIPA, 5, height = 10., distance = 100, width = 1.)

    def Knowing_Peaks(self):

        peaks    =   ()

        for (n_peaks, y, dist)    in      zip((4,5),(self.y, self.y_VIPA), (20, 100)):

            peaks    =   peaks    +   (Find_Highest_n_peaks(y, n_peaks, height = 10., distance = dist, width = 1.),)

        self.peaks     =   peaks[0]
        self.VIPA_peaks         =   peaks[1]

    def How_Many_Peaks_To(self, n_peaks = 4, delta = 1., distance = 20., width = 5., treshold = 5, fig = False, verbose = False, i_know_it_is = False):

        pk      =   find_peaks(self.y, height = treshold, distance = distance, width = width)
        height  =   np.max(pk[1]['peak_heights'])/2
        
        while True:

            pk      =   find_peaks(self.y, height = height, distance = distance, width = width)

            if (height > treshold) & (pk[0].size == n_peaks):

                if verbose:
                    print("Ho trovato valore dell'altezza per avere %d picchi: %f\n"%(n_peaks, height), pk)
                    _ = Analyze_Peaks(self.x_freq, self.y, 'GHz', fig = fig, verbose = verbose, height= height, distance = distance, width = width )
    
                break
            
            elif (height <= treshold):

                print(pk)
                raise ValueError('Errore: superata altezza minima %f\nQualcosa è andato storto'%(treshold))

            else: 

                height-=delta
        
        #check che i picchi Brillouin siano dentro agli elastici, premesso che so che possano essere più alti degli elastici

        condition_peaks =   ((self.y[pk[0][0]] < self.y[pk[0][1]]) | ((self.y[pk[0][3]] < self.y[pk[0][2]])))
        
        if condition_peaks &  (i_know_it_is == False):
            raise ValueError("Picchi Brillouin non sono interni o sono più alti degli elastici, controllare")

        self.spectrum_cut_height        =   height
        self.spectrum_peaks_dist        =   distance
        
    def How_Many_Peaks_To_VIPA(self, n_GHz_peaks = 5, n_gauss_peaks = 3, delta = 1., distance = 150., width = 1, treshold = 50, fig = False, verbose = False):

        h_save = ()

        for n_peaks in (n_GHz_peaks, n_gauss_peaks):
                
            pk      =   find_peaks(self.y_VIPA, height = treshold, distance = distance, width = width)
            height  =   np.max(pk[1]['peak_heights'])/2
            
            while True:

                pk      =   find_peaks(self.y_VIPA, height = height, distance = distance, width = width)

                if (height > treshold) & (pk[0].size == n_peaks):

                    h_save  =   h_save + (height,)

                    if verbose:
                        print("Ho trovato valore dell'altezza per avere %d picchi: %f\n"%(n_peaks, height), pk)
                        _ = Analyze_Peaks(self.x_VIPA, self.y_VIPA, 'GHz', fig = fig, verbose = verbose, height= height, distance = distance, width = width )
                    break
                
                elif (height <= treshold):

                    print(pk)
                    raise ValueError('Errore: superata altezza minima %f\nQualcosa è andato storto'%(treshold))

                else: 

                    height-=delta

        self.GHz_fit_height     =   h_save[0]
        self.gauss_fit_height   =   h_save[1]
        self.VIPA_peaks_dist    =   distance

    def Fit_VIPA_Gaussian(self, fig = False, verbose = False):

        peaks_idx   =   find_peaks(self.y_VIPA, height = self.gauss_fit_height, distance = self.VIPA_peaks_dist, width = 1)[0]
    
        gmod = Model(gaussian)
        result = gmod.fit(self.y_VIPA[peaks_idx], x = self.x_VIPA_freq[peaks_idx], A = 1., mu = 1, sigma = 1)
        
        A = result.values['A']
        mu = result.values['mu']
        sigma = result.values['sigma']

        #mi salvo anche i valori della gaussiana come valori iniziali per il fit

        self.p0['A']        =   [A]
        self.p0['mu']       =   [mu]
        self.p0['sigma']    =   [sigma]   

        if verbose:

            print ('Ho stimato i parametri della gaussiana come A = %3.2f\tmu  = %3.2f\tsigma = %3.2f' %(A, mu, sigma))
            print ('E li ho aggiunti ai parametri iniziali per il fit. Ora conosco %d parametri su %d \n' %(self.p0.size, self.n_params))

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
   
    def Fit_Pixel2GHz(self,  fig = False):

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

            x_fit = np.linspace(Delta_Pixel[0], Delta_Pixel[len(Delta_Pixel)-1], 100)
            y_fit = self.Poly2GHz[0]*(x_fit**3) + self.Poly2GHz[1]*(x_fit**2) + self.Poly2GHz[2]*x_fit + self.Poly2GHz[3]

            plt.figure()
            plt.plot(Delta_Pixel, Delta_GHz, '.', label='Peaks Data')
            plt.plot(x_fit, y_fit, label='Fit Polinomiale')
            plt.legend()
            plt.title('GHz vs Pixel')
            #plt.savefig('figure/'+fig+'.png')
            plt.show()
        
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

    def Spectrum_Pix2GHz (self, fig = False):

        #modifico in GHz le ascisse dello spettro
        # --> prima devo convertire in DeltaPixel

        self.x      =   self.x  -   self.x[self.y.argmax()]
        self.x_freq     =   ((self.x**3)*self.Poly2GHz[0])+ ((self.x**2)*self.Poly2GHz[1]) + (self.x*self.Poly2GHz[2]) + (self.Poly2GHz[3])
        
        if fig:

            fig1  = plt.figure()
            plt.xlabel('GHz')
            plt.plot(self.x_freq, self.y)
            plt.title('Spettro Exp in GHz')

    def Cut_n_Estimate_Spectrum(self, cut = True, estimate = False, distanza = 2/3, verbose = False ):
        
        """

        Funzione che esegue 
        #
        # taglio      :     trovo i valori dei picchi elastici e delle relative ampiezze 
        #                   (grazie scipy) e banalmente mi sposto di distanza rispetto 
        #                   alle posizioni dei picchi sull'array delle frequenze
        #                   il taglio sulle y è consequenziale, il tutto grazie agli indici
        #
        # stima dei parametri iniziali :
        #                    riesco a stimare qualche parametro
        #  
        #   OSS: tutto funziona perchè mi aspetto due picchi elastici ad aprire e
        #        chiudere lo spettro
        #      

        """

        query_min               =   self.x_freq[self.peaks['peaks_idx'][0]] + (self.peaks['peaks_width'][0]*distanza)
        _, idx_min              =   Find_Nearest(self.x_freq, query_min)
            
        query_max               =   self.x_freq[self.peaks['peaks_idx'][3]] - (self.peaks['peaks_width'][3]*distanza)
        _, idx_max              =   Find_Nearest(self.x_freq, query_max)
        
        # STIMA PARAMETRI INIZIALI della funzione teorica
        # (quelli che posso, e devo farlo  prima di tagliare)
        
        self.p0  =   pd.DataFrame({}, columns = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width','delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset'), index = ['Initials'])

        # 1)    stima dei parametri dai dati
        #       Omega come la media della posizione dei massimi brillouin
        #       Gamma (=Delta)come la media dell'ampiezza restituita da find_peaks per i Brillouin /10 (mah..)
        #       offset come la media dei valori dello spettro nella zona tra i due picchi

        if  estimate:
                
            self.p0['Omega']            =   [np.absolute(self.x_freq[self.peaks['peaks_idx'][2]] - self.x_freq[self.peaks['peaks_idx'][1]])*0.5]
            self.p0['Gamma']            =   [(self.peaks['peaks_width'][2]+self.peaks['peaks_width'][1])/20]
            self.p0['Delta']            =   self.p0['Gamma']
            self.p0['tau'  ]            =   [100.]
            self.p0['offset']           =   np.mean(self.y[self.peaks['peaks_idx'][1]:self.peaks['peaks_idx'][2]])
            

            # 2)i parametri iniziali che dovrebbero andare bene sempre
            self.p0['Co']               =   [1.]#amplitude factor
            self.p0['shift']            =   [0.]
            self.p0['delta_amplitude']  =   [1.]
            self.p0['delta_width']      =   [0.5]


        if verbose:

            print ("\n\nHo stimato %d parametri iniziali per il fit che andrai a fare\n" %(self.p0.size))

            for p in self.p0.T.index :

                print(p, ' = %4.3f \n' %(self.p0[p][0]))


        
        # procedo con il taglio se è voluto

        if cut:           
            
            self.x_freq         =   self.x_freq[idx_min:idx_max]
            self.x          =   self.x[idx_min:idx_max]
            self.y              =   self.y[idx_min:idx_max]
            self.y_err          =   self.y_err[idx_min:idx_max]

    def Estimate_Initial_Parameters(self, p0, treshold, columns, **kwargs):

        """

        Passo a p0 i parametri ottenuti dal primo Fit
        In più aggiundo un  delta uguale a Gamma e un Tau infinito

        """
        self.cost           =   0.5*np.sum(self.Residuals(p0, self.y)**2)
        print('costo = '+str(self.cost))
        self.p0         =   pd.DataFrame({idx : value for (idx, value) in zip(cols, p0)}, index = ['Values'])

        if self.cost >= treshold:
            
            percents        =   ('positive', 0.2, 'positive', 'positive', 'positive', 0.1, 0.1, 0.1,  np.inf, np.inf)
            self.Get_Fit_Bounds(percents, columns = cols_mark)
            p0              =   np.concatenate((p0[:3], p0[5:]))

            self.Non_Linear_Least_Squares_Markov(p0, bound = (self.bounds['down'].values, self.bounds['up'].values), **kwargs)

            self.p0             = pd.DataFrame(self.Fit_Params.T['Values']).T
            self.p0['Delta']    = self.p0['Gamma']
            self.p0['tau']     = [100.]

            self.p0  =   pd.DataFrame(self.p0, columns = cols)
            self.cost           =   0.5*np.sum(self.Residuals(self.p0.values[0], self.y)**2)
            print('costo dopo fit = '+str(self.cost))
            
            if self.res_lsq['success']:
                return  1
            else:
                return  2

        else:
            return  0

    def Recover_Initial_Parameters(self, p0):

        self.p0         =   pd.DataFrame(p0, index = cols_mark).T
        self.p0['Delta']    = self.p0['Gamma']
        self.p0['tau']     = [100.]
        self.p0  =   pd.DataFrame(self.p0, columns = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width','delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset'))
        self.cost           =   0.5*np.sum(self.Residuals(self.p0.values[0], self.y)**2)
        print('costo dopo fit = '+str(self.cost))
            


    def Take_A_Look_Before_Fitting(self):
        
        self.y_Gauss_convolution = self.Gauss_Convolve_Theoretical_Response_Fast(self.p0.values[0])
        
        print("Valore stimato della cost function prima del fit:\n{}".format(self.cost))

        plt.figure()
        plt.plot(self.x_freq, self.y_Gauss_convolution, '-', label = 'Initial_Guess')
        plt.plot(self.x_freq, self.y, '*', label = 'Data')
        plt.title('Goodness of initial guess, cost = %f'%(self.cost))
        plt.xlabel('Freq(GHz)')
        plt.show()

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
        #   param delta         p[5] = delta amplitude
        #                       p[6] = delta factor
        #
        #   param gauss         p[7] = A
        #                       p[8] = mu 
        #                       p[9] = sigma 
        #
 
        #   questi ultimi       p[10] = shift
        #   per forza           p[11] = offset       

        if (p.size != 12):
            raise ValueError("Calcola stai a usa modello Real = Markov + Exp\n p deve avere 12 elementi! Ne ha {}".format(p.size))

        if fantoccio :

            conv_range = np.linspace(fantoccio[0], fantoccio[1], 200)

        else:
            
            conv_range = self.x_freq

        self.y_convolution      =       np.zeros(conv_range.size)
        w_j_VIPA                =       np.linspace(-35,35, self.x_VIPA_freq.size)
        kernel                  =       self.Interpolate_VIPA(w_j_VIPA)
        kernel                  =       kernel/(p[7]*(np.exp(-((w_j_VIPA-p[8])**2)/(2*(p[9]**2)))))
        
        for  ii in range(len(conv_range)):

            delta_w                 =   conv_range[ii] -   w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_2(delta_w-p[10], *p[0:7])
            self.y_convolution[ii]  =   np.sum(theor*kernel)

        if fig:

                plt.figure()
                plt.plot(conv_range, self.y_convolution)
                plt.xlabel('GHz')
                plt.title('convoluzione dati VIPA con funz S2 for '+self.__str__())
                plt.show
        
 

        self.y_Gauss_convolution   =   p[11] + self.y_convolution*p[7]*np.exp(-((conv_range - p[8])**2)/(2*(p[9]**2)))

        if fig:

            plt.figure()
            plt.plot(conv_range, self.y_Gauss_convolution)
            plt.title('convoluzione dati VIPA con funz S2 più moltiplicazione Gauss')

    def Convolve_Theoretical_Response_Fast (self, p, p_gauss):

        """
        Versione per non fittare parametri della gaussiana nel fit totale

        Funziona con modello Real (Markov + Rilassamento exp) --> 12 parametri
        """

        self.y_convolution      =       np.zeros(self.x_freq.size)
        _ , idx_min             =       Find_Nearest(self.x_VIPA_freq, -35.)
        _ , idx_max             =       Find_Nearest(self.x_VIPA_freq, 35.)
        w_j_VIPA                =       self.x_VIPA_freq[idx_min:idx_max]#-1 per azzeccare dim coi bins
        VIPA_w_j                =       self.y_VIPA[idx_min:idx_max-1]
        Delta_w_j_VIPA          =       Get_Delta_Between_Array_Elements(w_j_VIPA)
        w_j_VIPA                =       w_j_VIPA[:w_j_VIPA.size-1]
        kernel                  =       VIPA_w_j*Delta_w_j_VIPA/(p_gauss[0]*(np.exp(-((w_j_VIPA-p_gauss[1])**2)/(2*(p_gauss[2]**2)))))
        
        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_2(delta_w-p[7], *p[0:7])
            self.y_convolution[ii]  =   np.sum(theor*kernel)

        self.y_Gauss_convolution   =   p[8] + self.y_convolution*p_gauss[0]*np.exp(-((self.x_freq - p_gauss[1])**2)/(2*(p_gauss[2]**2)))

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
        kernel                  =       VIPA_w_j*Delta_w_j_VIPA/(p[7]*(np.exp(-((w_j_VIPA-p[8])**2)/(2*(p[9]**2)))))
        
        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_2(delta_w-p[10], *p[0:7])
            self.y_convolution[ii]  =   np.sum(theor*kernel)
 
        self.y_Gauss_convolution   =   p[11] + self.y_convolution*p[7]*np.exp(-((self.x_freq - p[8])**2)/(2*(p[9]**2)))

        return self.y_Gauss_convolution

    def Gauss_Convolve_Markovian_Response (self, p, fantoccio = False, fig = False, compare = False):
        
        #       nomenclatura:   p[0] = Co
        #                       p[1] = Omega
        #                       p[2] = Gamma

        #   param delta         p[3] = delta amplitude
        #                       p[4] = delta factor        #
        #   param gauss         p[5] = A
        #                       p[6] = mu 
        #                       p[7] = sigma 
        #   questi ultimi       p[8] = shift
        #   per forza           p[9] = offset       

        p = np.asarray(p)
        
        if (p.size != 10):

            raise ValueError("Calcola stai a usa modello Markov \n p deve avere 10 elementi! Ne ha {}".format(p.size))

        if fantoccio :

            conv_range = np.linspace(fantoccio[0], fantoccio[1], 200)

        else:
            
            conv_range = self.x_freq

        self.y_markov_convolution      =       np.zeros(conv_range.size)
        w_j_VIPA                       =       np.linspace(-35,35, self.x_VIPA_freq.size)
        kernel                         =       self.Interpolate_VIPA(w_j_VIPA)
        kernel                         =       kernel/(p[5]*(np.exp(-((w_j_VIPA-p[6])**2)/(2*(p[7]**2)))))
        
        for  ii in range(len(conv_range)):

            delta_w                         =   conv_range[ii] -   w_j_VIPA
            theor                           =   S_Dynamical_Form_Factor_0(delta_w-p[8], *p[0:5])
            self.y_markov_convolution[ii]   =   np.sum(theor*kernel)

        if fig:

                plt.figure()
                plt.plot(conv_range, self.y_markov_convolution, label = 'initial guess')
                plt.xlabel('GHz')
                plt.title('convoluzione dati VIPA con funz S0 for '+self.__str__())

        self.y_Gauss_markov_convolution   =   p[9] + self.y_markov_convolution*p[5]*np.exp(-((conv_range - p[6])**2)/(2*(p[7]**2)))

        if fig:

            plt.figure()
            plt.plot(conv_range, self.y_Gauss_markov_convolution, label = 'initial_guess')
            plt.title('convoluzione dati VIPA con funz S0 più moltiplicazione Gauss')
            if compare:
                plt.plot(self.x_freq, self.y, '+', label = 'data')
                plt.legend()
                plt.show

    def Gauss_Convolve_Markovian_Response_Fast (self, p):

        """
        Versione veloce della funzione Gauss_Convolve_Markovian
        Funziona con modello Markov --> 10 parametri

        """

        self.y_markov_convolution      =       np.zeros(self.x_freq.size)
        _ , idx_min             =       Find_Nearest(self.x_VIPA_freq, -35.)
        _ , idx_max             =       Find_Nearest(self.x_VIPA_freq, 35.)
        w_j_VIPA                =       self.x_VIPA_freq[idx_min:idx_max]#-1 per azzeccare dim coi bins
        VIPA_w_j                =       self.y_VIPA[idx_min:idx_max-1]
        Delta_w_j_VIPA          =       Get_Delta_Between_Array_Elements(w_j_VIPA)
        w_j_VIPA                =       w_j_VIPA[:w_j_VIPA.size-1]
        kernel                  =       VIPA_w_j*Delta_w_j_VIPA/(p[5]*(np.exp(-((w_j_VIPA-p[6])**2)/(2*(p[7]**2)))))
        
        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_0(delta_w-p[8], *p[0:5])
            self.y_markov_convolution[ii]  =   np.sum(theor*kernel)
 
        self.y_Gauss_markov_convolution   =   p[9] + self.y_markov_convolution*p[5]*np.exp(-((self.x_freq - p[6])**2)/(2*(p[7]**2)))

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

    def Residuals(self, p, y):
        
        return (self.Gauss_Convolve_Theoretical_Response_Fast(p) - y)/self.y_err


    def Residuals_NoGauss(self, p, y, p_gauss):
        
        return (self.Convolve_Theoretical_Response_Fast(p, p_gauss) - y)/self.y_err
    

    def Non_Linear_Least_Squares (self, p0, p_gauss, columns, bound = (-np.inf, np.inf), fig = False, **kwargs):

        """

        OSS le kwargs sono quelle di least_squares() di scipy.optimize

        """

        start            =    time.process_time()        
        self.res_lsq     =    least_squares(self.Residuals_NoGauss, p0, args= ([self.y, p_gauss]), bounds = bound, **kwargs)
        print("s impiegati a fare il fit totale ", time.process_time()-start, '\n')

        Parameters       =    self.res_lsq.x

        try:
            J                =    self.res_lsq.jac
            cov              =    np.linalg.inv(J.T.dot(J))
            Delta_Parameters =    np.sqrt(np.diagonal(cov))

        except  np.linalg.LinAlgError as err:

            if 'Singular matrix' in str(err):
                print('Ho trovato matrice singolare')
                Delta_Parameters    =   np.empty(len(p0))
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
        
        
        self.Fit_Params = pd.DataFrame((Parameters, Delta_Parameters, p0), index = ('Values', 'StdErrs', 'Initials'), columns = columns)
        if (self.res_lsq['success'] == False):
            return 2
        else: return 1



    def Get_Fit_Bounds(self, percents, columns): 

        """
        Funzione che definisce un attributo self.bounds, un dataframe con le colonne up and down
        ogni valore di self.p0 ha due elementi associati che corrispondono ai limiti che si vogliono imporre nel fit

        il calcolo dei bounds si basa sulla funzione Get_Around e sul parametro da passare percents, il quale DEVE essere di lunghezza uguale al 
        p0 di riferimento, con quell'ordine di paramentri

        se un elemento di percents è np.inf, la funzione ritorna l'asse reale come limiti per il parametro corrispondente
        se un elemento di percents è 'positive', la funzione ritorna [0, np.inf] come limiti #

        la funzione funziona sia per Markov che per Real

        """

        if len(percents) != len(columns):

            raise ValueError ("Lunghezza della tupla percents = %d errata, deve essere uguale alla dimensione delle colonne passate = %d"%(len(percents), len(columns)))
                
        self.bounds = pd.DataFrame( {}, index= columns, columns=('down', 'up'))

        for (col, frac) in zip(columns, percents):
    
            if  frac    ==  np.inf:

                self.bounds['down'][col]     =    -np.inf
                self.bounds['up'][col]       =    np.inf
            
            elif    frac == 'positive':

                self.bounds['down'][col]     =    0
                self.bounds['up'][col]       =    np.inf
            
            else:
                    
                bound   =   Get_Around(self.p0[col]['Values'], frac)
                self.bounds['down'][col]     =   bound[0]
                self.bounds['up'][col]       =   bound[1]

    def Residuals_Markov(self, p, y):
        
        return (self.Gauss_Convolve_Markovian_Response_Fast(p) - y)/self.y_err
    
    def Non_Linear_Least_Squares_Markov (self, p0, bound = (-np.inf, np.inf), max_iter = None, verbose = 0, fig = False, **kwargs):

        start            =    time.process_time()
        self.res_lsq     =    least_squares(self.Residuals_Markov, p0, args= ([self.y]), bounds = bound, max_nfev = max_iter, verbose = verbose, **kwargs)
        print("s impiegati a fare il fit ", time.process_time()-start, '\n')
        Parameters       =    self.res_lsq.x
        
        try:
            J                =    self.res_lsq.jac
            cov              =    np.linalg.inv(J.T.dot(J))
            Delta_Parameters =    np.sqrt(np.diagonal(cov))

        except  np.linalg.LinAlgError as err:

            if 'Singular matrix' in str(err):
                print('Ho trovato matrice singolare')
                Delta_Parameters        =   np.empty(len(p0))
                Delta_Parameters[:]     =   np.nan
            else :
                raise

        self.y_markov_fit  = self.y_Gauss_markov_convolution
      
        if fig:

            plt.figure()
            plt.title('Fit for '+self.__str__())
            plot(self.x_freq, self.y, '+', label='Data')
            plot(self.x_freq, self.y_markov_fit, label= 'Fit')
            plt.legend()
           
        self.Fit_Params = pd.DataFrame((Parameters, Delta_Parameters, p0), index = ('Values', 'StdErrs', 'Initials'), columns = ('Co', 'Omega', 'Gamma',  'delta_width','delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset'))
        
        if (self.res_lsq['success'] == False):
            return 2
        else: return 1
    
    def Recover_Fit_Params(self, dictio_string):

        """
        Funzione che tramite libreria json ricostruisce DataFrame Fit_Params per l'oggetto
        dictio_string contiene sotto forma di stringa le info in formato json (che sono un dizionario)

        """
        self.Fit_Params =   pd.DataFrame(json.loads(dictio_string))

def Initialize_Matrix(n_rows, n_cols):

    matrix = ()
    
    for ii in range(n_rows):

        riga = ()

        for jj in range(n_cols):

            riga  = riga + (Spectrum('Element ('+str(ii)+str(jj)+')'),)
    

        matrix = matrix + (riga,)

    print('Ho inizializzato una matrice %dx%d, per un totale di %d spettri'%(len(matrix), len(matrix[0]), len(matrix)*len(matrix[0])  ))

    return matrix

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
        else :

            isolated    =   isolated    +   ((ii+1, jj),)

    return isolated
    
    
def Unpack_Fit(fit):

    non_fitted   =   ()
    accomplished =   ()
    exceded      =   ()
    fitted       =   ()

    for (what, (ii,jj)) in fit:

        if  what == 0:

            non_fitted  =    non_fitted +   ((ii,jj),)

        elif what == 1:
            
            accomplished      =   accomplished  +   ((ii,jj),)
            fitted            =   fitted + ((ii,jj),)

        elif what == 2:

            exceded      =   exceded  +   ((ii,jj),)
            fitted            =   fitted + ((ii,jj),)
    
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

       
        plt.matshow(fit_map, cmap = 'ocean')
        plt.title('Fit Map')
        plt.colorbar()
        plt.xlabel('Row Index')
        plt.ylabel('Col Idx')
        plt.savefig(path + fig+'.png')
        plt.close()
    
    return fit_map

def Gamma_Verify_Markov_Fit(matrix, fitted):

    for (ii, jj) in fitted:

        if matrix[ii][jj].Fit_Params['Gamma']['Values'] > 1.:

            plt.figure()
            plt.plot(matrix[ii][jj].y)
            plt.plot(matrix[ii][jj].y_fit)
            plt.title(str((ii,jj)) +' successo: '+ str(matrix[ii][jj].res_lsq['success']))
            print(matrix[ii][jj].res_lsq['message'])


def Get_Parameter_Map(parameter, columns, matrix, n_rows, n_cols, fitted, excluded, fig = False, path = ''):

    if parameter not in matrix[fitted[0][0]][fitted[0][1]].Fit_Params.columns:
            
            print ('Parametro scelto non in quelli fittati: scegliere uno tra\n', columns, '\nse si sta mappando dopo fit totale, markoviano leva tau e delta\n')
            raise ValueError('COJONE')

    p_map   =   np.zeros((n_rows, n_cols))
    
    for ii in range(n_rows):
        for jj in range (n_cols):
            if (ii,jj)  in fitted:     
                if  (ii, jj) not in excluded:
                    p_map[ii,jj]    =   matrix[ii][jj].Fit_Params[parameter]['Values']

    print('Completata Paramter_Map per '+parameter)
    
    if fig:

        plt.matshow(p_map, cmap = 'ocean')
        plt.title(parameter+' Map')
        plt.colorbar()
        plt.xlabel('Row Index')
        plt.ylabel('Col Idx')
        plt.savefig(path + fig+'.png')
        plt.close()

    return  p_map


def Escludi_a_Mano(to_add, excluded):

    for (ii, jj) in to_add:

        excluded = excluded + ((ii,jj),)

    return excluded

def Whose_Gamma_Too_High(treshold, matrix, fitted):

    too_high    =   ()
    for (ii,jj) in (fitted):
        if matrix[ii][jj].Fit_Params['Gamma']['Values'] > treshold:
            too_high    =   too_high    +   ((ii,jj),)
            print(str((ii,jj))+' ha gamma = %3.2f'%(matrix[ii][jj].Fit_Params['Gamma']['Values']))

    return too_high



def Save_Params(path, n_rows, n_cols, matrix, fitted):

    f = open(path+'save_params.txt', 'w')

    for ii in range(0, n_rows,1 ):
        for jj in range(n_cols):
            if (ii,jj) in fitted:
                f.write(str(matrix[ii][jj].Fit_Params.values[0])+'\n')
    f.close()

def Save_Fitted_Info(path, n_rows, n_cols, fitted):

    f = open(path+'fitted.txt', 'w')
    f.write('(')
    for ii in range(0, n_rows,1):
        for jj in range(n_cols):
            if (ii,jj) in fitted:
                f.write(str((ii,jj))+',')
    f.write(')')
    f.close()



def Save_Fit_Info(fit, filename = 'fit.txt', path = './'):

    f = open(path+filename, 'w')
    f.write('(')
    for ft in fit:
                f.write(str(ft)+',')
    f.write(')')
    f.close()


def Parse_Parameter_Save(parameter_file = 'save_params.txt', path = ''):

    parameters  = ()
    f   =   open(path+parameter_file, 'r')
    lines   =   f.readlines()
    f.close()
    indices =   np.arange(0,len(lines), 3)
    for ii in indices:
            
        
        a   =   np.array(lines[ii].split()[1:], dtype= np.float64)
        b   =   np.array(lines[ii+1].split()[:], dtype= np.float64)
        c   =   np.array(lines[ii+2].replace(']', ' ').split()[:], dtype= np.float64)

        parameters  =   parameters  +   (np.concatenate((a,b,c)),)

    return parameters

def Verify_Initial_Conditions(matrix, ver = (), init = ()):
    
    plt.title('Verify Initial Conditions for '+str((ver[0],ver[1]))+'with '+ str((init[0],init[1]))+ 'parameters')
    plt.plot(matrix[ver[0]][ver[1]].x_freq, matrix[ver[0]][ver[1]].y, '+', label = 'dati')
    plt.plot(matrix[ver[0]][ver[1]].x_freq, matrix[init[0]][init[1]].Gauss_Convolve_Theoretical_Response_Fast(matrix[init[0]][init[1]].p0.values[0]))

    print("Il parametro è :", matrix[init[0]][init[1]].p0.values[0])

def Save_Fit_Parameters(matrix, fitted, out_filename = 'fit_params.txt' , path = './'):

    """
    Salvo nella cartella di analisi un file di nome 'fit_params.txt' nella cartella di analisi associata alla presa dati
    la struttura è 
    ogni riga contiene un dizionario già strutturato per file .json che DataFrame è in grado di acquisire

    FONDAMENTALE è l'ordine delle righe che corrisponde a quello di fitted()

    """

    with open(path+out_filename, 'w') as f_out:

        for (ii,jj) in fitted:
   
            f_out.write(json.dumps(matrix[ii][jj].Fit_Params.to_dict())+'\n')


def Plot_Elements_Spectrum(matrix, elements_iterable):

    for (ii,jj) in elements_iterable:

        plt.figure()
        plt.plot(matrix[ii][jj].x_freq, matrix[ii][jj].y)
        plt.title(str((ii,jj)))
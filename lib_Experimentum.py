import  time
import  numpy               as np
import  pandas              as pd
import  h5py
from    matplotlib.pyplot   import plot
import  matplotlib.pyplot   as plt
import  scipy.signal        as syg
from    scipy.optimize      import leastsq
from    scipy.optimize      import least_squares


from    Models              import S_Dynamical_Form_Factor_2, S_2_Generate        
from    Alessandria         import *
from    lmfit               import Model

free_spectral_range =   29.9702547 #GHz

#       PARTE   2       Classi 


class Spectrum  :

    #   Oggetto che sarà alla base di ogni spettro sperimentale

    def __init__(self, name, n_params = 12):

        """
        Definisci il nome dell'oggetto per le stampe e il numero di parametri del tuo modello

        """
        self.name       =   name
        self.n_params   =   n_params

    def __str__(self):

        return '{}'.format(self.name)
        
    def Get_Spectrum(self, y, offset = 183., **img_kwargs):

        #   y       è l'array dello spettro
        #   
        #   nu0     è la frequenza del laser con cui 
        #           lo spettro è stato acquisito
        #           (in realtà probably numero associato, v. parametro v Claudia)
        #           

        print ("Cojone  modifica le cose in modo che la funzione apra da sola il file.mat\n\n")

        self.x_pix      =   np.arange(1, len(y)+1, 1)
        self.offset     =   offset
        self.y          =   y - self.offset

        if 'fig' in img_kwargs :

            plt.figure()
            plt.plot(self.x_pix, self.y)
            plt.xlabel('Pixels')
            plt.title('Experimental spectrum for '+self.__str__())
            plt.savefig(img_kwargs['save_path']+img_kwargs['fig']+'.png')
            plt.show()


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
            self.y_VIPA         =   VIPA[:,self.nu0] - self.offset
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

    
    def Fit_VIPA_Gaussian(self, fig = False, **syg_kwargs):

        peaks_idx, _    =   Analyze_Peaks(self.x_VIPA_freq, self.y_VIPA, 'GHz', **syg_kwargs)
    
        gmod = Model(gaussian)
        result = gmod.fit(self.y_VIPA[peaks_idx], x = self.x_VIPA_freq[peaks_idx], A = 1., mu = 1, sigma = 1)
        
        A = result.values['A']
        mu = result.values['mu']
        sigma = result.values['sigma']

        #mi salvo anche i valori della gaussiana come valori iniziali per il fit

        self.p0['A']        =   [A]
        self.p0['mu']       =   [mu]
        self.p0['sigma']    =   [sigma]   

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

   

    def Fit_Pixel2GHz(self,  fig = False, **syg_kwargs):

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

        peaks               =   syg.find_peaks(self.y_VIPA, **syg_kwargs)
        peaks_pix           =   np.array(peaks[0])+1 #array con indici(anche pixel dunque) dei picchi
        peaks_counts        =   self.y_VIPA[peaks_pix]
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
            plt.xlim(-100,100)
            plt.show()

    def Spectrum_Pix2GHz (self, fig = False):

        #modifico in GHz le ascisse dello spettro
        # --> prima devo convertire in DeltaPixel

        self.x_pix      =   self.x_pix  -   self.x_pix[self.y.argmax()]
        self.x_freq     =   ((self.x_pix**3)*self.Poly2GHz[0])+ ((self.x_pix**2)*self.Poly2GHz[1]) + (self.x_pix*self.Poly2GHz[2]) + (self.Poly2GHz[3])
        
        if fig:

            fig1  = plt.figure()
            plt.xlabel('GHz')
            plt.plot(self.x_freq, self.y)
            plt.title('Spettro Exp in GHz')

    def Cut_n_Estimate_Spectrum(self, cut = True, distanza = 2/3, verbose = False, **syg_kwargs ):
        
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

        peaks_idx, peaks_width  =   Analyze_Peaks(self.x_freq, self.y, 'GHz', verbose = True, fig = False, **syg_kwargs)

        query_min               =   self.x_freq[peaks_idx[0]] + (peaks_width[0]*distanza)
        _, idx_min              =   Find_Nearest(self.x_freq, query_min)
            
        query_max               =   self.x_freq[peaks_idx[3]] - (peaks_width[3]*distanza)
        _, idx_max              =   Find_Nearest(self.x_freq, query_max)
        
        # STIMA PARAMETRI INIZIALI della funzione teorica
        # (quelli che posso, e devo farlo  prima di tagliare)
        
        
        self.p0 = pd.DataFrame({})
        # 1)    stima dei parametri dai dati
        #       Omega come la media della posizione dei massimi brillouin
        #       Gamma (=Delta)come la media dell'ampiezza restituita da find_peaks per i Brillouin /10 (mah..)
        #       offset come la media dei valori dello spettro nella zona tra i due picchi

        self.p0['Omega']            =   [np.absolute(self.x_freq[peaks_idx[2]] - self.x_freq[peaks_idx[1]])*0.5]
        self.p0['Gamma']            =   [(peaks_width[2]+peaks_width[1])/20]
        self.p0['Delta']            =   self.p0['Gamma']
        self.p0['offset']           =   np.mean(self.y[peaks_idx[1]:peaks_idx[2]])
        

        # 2)i parametri iniziali che dovrebbero andare bene sempre
        self.p0['Co']               =   [1000.]#amplitude factor
        self.p0['shift']            =   [0.]
        self.p0['tau']              =   [1.]
        self.p0['delta_amplitude']  =   self.y[[peaks_idx[self.y[peaks_idx].argmax()]]]
        self.p0['delta_width']      =   [0.5]

        if verbose:

            print ("\n\nHo stimato %d parametri iniziali per il fit che andrai a fare\n" %(self.p0.size))

            for p in self.p0.T.index :

                print(p, ' = %4.3f \n' %(self.p0[p][0]))


        
        # procedo con il taglio se è voluto

        if (cut == True):           
            
            self.x_freq         =   self.x_freq[idx_min:idx_max]
            self.x_pix          =   self.x_pix[idx_min:idx_max]
            self.y              =   self.y[idx_min:idx_max]


        
        

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
        #   questi ultimi       p[5] = shift
        #   per forza           p[6] = offset
        #
        #   param gauss         p[7] = mu 
        #                       p[8] = sigma 
        #
        #   param delta         p[9] = delta amplitude
        #                       p[10] = delta factor 
        #
        #   OSS shift è ovviamente messo dopo la moltiplicaz per gauss

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

        #return necessario per Residuals in least_squares
        return self.y_Gauss_convolution
        """
    def Sinc_Convolve_Theoretical_Response (self, p, conv_range, fig = False):

       
        Versione funzione Convolve con inviluppo Sinc invece che Gaussiano
       

        self.y_convolution      =       np.zeros(self.x_freq.size)
        w_j                     =       np.linspace(-15,15,self.x_freq.size)
        w_j_VIPA                =       np.linspace(-35,35, self.x_VIPA_freq.size)
        kernel                  =       self.Interpolate_VIPA(w_j_VIPA)
        kernel                  =       kernel*(((np.sin((self.x_freq-p[7])/p[8]))/((self.x_freq-p[7])/p[8]))**2)
        
        for  ii in range(len(self.x_freq)):

            delta_w                 =   self.x_freq[ii] -   w_j_VIPA
            theor                   =   S_Dynamical_Form_Factor_2(delta_w-p[5], p[0], p[1], p[2], p[3], p[4], p[9], p[10])
           
            self.y_convolution[ii]  =   np.sum(theor*kernel)

        if fig:

                plt.figure()
                plt.plot(self.x_freq, self.y_convolution)
                plt.title('convoluzione dati VIPA con funz S2')
        
 
        self.y_Sinc_convolution    =   p[6] + self.y_convolution/(((np.sin((self.x_freq-p[7])/p[8]))/((self.x_freq-p[7])/p[8]))**2)

        if fig:

            plt.figure()
            plt.plot(self.x_freq, self.y_Gauss_convolution)
            plt.title('convoluzione dati VIPA con funz S2 più moltiplicazione Gauss')
            plt.savefig(fig)
    """


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
        
        return (self.Gauss_Convolve_Theoretical_Response(p) - y)
    
    def Non_Linear_Least_Squares (self, p0, my_method = None, bound = (-np.inf, np.inf), max_iter = None, verbose = 0, fig = False, **kwargs):

        if (p0 == 'auto'):
            
            #riordino le colonne e le assumo nel vettore p0 da passare alla funzione

            p0  =   pd.DataFrame(self.p0, columns = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'shift', 'offset', 'delta_width','delta_amplitude', 'A', 'mu', 'sigma'))
            p0  =   self.p0.values[0]

            print("\nAssumo vettore dei parametri iniziali da quelli stimati\n", p0, '\n')
        
        else:  print("\nAssumo vettore dei parametri iniziali da quello inserito come arg\n", p0, '\n')


        start            =    time.process_time()
        
        if ( my_method == 'lsq'):
            
            self.res_lsq     =    leastsq(self.Residuals, p0, args= self.y, full_output = True, **kwargs)
            print("s impiegati a fare il fit ", time.process_time()-start, '\n')
            Parameters       =   self.res_lsq[0]
            Delta_Parameters =   np.zeros(Parameters.size)
            
        elif (my_method == 'least_squares'):
            
            self.res_lsq     =    least_squares(self.Residuals, p0, args= ([self.y]), bounds = bound, max_nfev = max_iter, verbose = verbose, **kwargs)
            print("s impiegati a fare la convoluzione ", time.process_time()-start, '\n')
            Parameters       =    self.res_lsq.x
            J                =    self.res_lsq.jac
            cov              =    np.linalg.inv(J.T.dot(J))
            Delta_Parameters =    np.sqrt(np.diagonal(cov))
            
        else : raise ValueError("Specificare se usare metodo scipy.optimize.least_squares() o scipy.optimize.lsq()\n")
           
        y_fit  = self.Gauss_Convolve_Theoretical_Response(Parameters)

        if fig:

            plt.figure()
            plt.title('Fit for '+self.__str__())
            plot(self.x_freq, self.y, '+', label='Data')
            plot(self.x_freq, y_fit, label= 'Fit')
            plt.legend()
        
        
        self.Fit_Params = pd.DataFrame((Parameters, Delta_Parameters, p0), index = ('Values', 'StdErrs', 'Initials'), columns = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'shift', 'offset', 'Mu_Sinc', 'Sigma_Sinc', 'Delta_Width', 'Delta_Factor'))
        
        


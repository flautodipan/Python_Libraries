
import  numpy               as      np
import  matplotlib.pyplot   as      plt
from    PIL                 import  Image
from    scipy.signal        import  find_peaks
from    scipy.io            import  loadmat


"""    

##############################################################################################
#############################################################################################

Funzioni generiche 
Insieme di funzioni utili per vari scopi in tutto il lavoro


"""

def Get_Delta_Between_Array_Elements(a):

    """
    Restituisce un array di lunghezza N-1 con gli intervalli tra i valori di a

    """

    a_up    =   a[1:]
    a_down  =   a[:a.size-1]

    return  a_up-a_down 

def Find_First_n_peaks(pk, n_peaks, exclude = None):

    """

    exclude deve essere sotto forma di iterabile

    """

    if exclude:
            
        idx      =   range(pk['idx'].size)

        for ii in exclude:

            idx      =   np.delete(idx, ii-1)
    
        idx      =   idx[:n_peaks]
    
    else:

        idx     =   range(n_peaks)

    return          {'idx' : pk['idx'][idx], 'heights': pk['heights'][idx], 'widths' : pk['widths'][idx]}

def Find_Highest_n_peaks(pk, n_peaks):

    from_heighest    =   np.flip(pk['heights'].argsort())
    from_heighest    =   from_heighest[:n_peaks]
    idx              =   np.sort(from_heighest)

    return          {'idx' : pk['idx'][idx], 'heights': pk['heights'][idx], 'widths' : pk['widths'][idx]}

def Analyze_Peaks(x, y, x_dim, fig = False, verbose = False, **syg_kwargs):

    """
    Grazie a scipy.signal.find_peaks(), la funzione analizza i picchi della funzione campionata da y
    passare i parametri di find_peaks voluti come kwargs

    > stampa le cose che trova se verbose == True
    > Ritorna dizionario {#picchi, indici_picchi nell'array y, ampiezza, altezza}
    
    """
    y                   =   np.asarray(y)
    peaks               =   find_peaks(y, **syg_kwargs)
    peaks_idx           =   np.array(peaks[0])
    peaks_width         =   peaks[1]['widths']
    peak_heights        =   peaks[1]['peak_heights']

    print("\n\n Ho trovato %d picchi nel tuo spettro sperimentale con le caratteristiche richieste\n Altezza > %3.2f \n Spessore > %3.2f \n\n" %(peaks_idx.size, syg_kwargs['height'], syg_kwargs['width']))
    
    if fig:

        plt.figure()
        plt.plot(x, y, label = 'Signal')
        plt.plot(x[peaks_idx], y[peaks_idx], '*', label = 'Peaks')
        plt.legend()

    if verbose:
            
        for kk in range(len(peaks_idx)):

            print("\n Il picco %d ha: \t indice = %d \t x_value (%s) = %3.2f \t ampiezza(%s) = %3.2f \t altezza = %3.2f \n" %(kk+1, peaks_idx[kk], x_dim, x[peaks_idx[kk]], x_dim, peaks_width[kk], peak_heights[kk]))
                
    return {'n_peaks' : peaks_idx.size, 'idx' : peaks_idx, 'peaks_width' : peaks_width, 'peak_heights' : peak_heights}



def Get_Around(value, factor):

    """
    Ritorna una tupla contenente (inf, sup) di un intervallo centrato su value
    la larghezza dell'intervallo è 2*factor*value --> factor tra (0,1)

    Utile per bounds nei fit

    """
    inf     =   value - np.abs(value*factor)
    sup     =   value + np.abs(value*factor)
    
    return (inf, sup)
    
def Estimate_FWHM (x,y):

    half    =   (np.max(y)+np.min(y))/2

        #suddivido array in due parti
    control =   0 #cerco di eliminare parte asimmetrica

    x_left  =   x[control:np.argmax(y)]
    y_left  =   y[control:np.argmax(y)]

    #plt.figure()
    #plt.plot(x_left,y_left) 

    x_right =   x[np.argmax(y):len(x) - control]
    y_right =   y[np.argmax(y):len(x) - control]

    #plt.plot(x_right,y_right) 

    temp_left   =   Find_Nearest(y_left, half)
    idx_left    =   np.where(y_left==temp_left)
    x_half_left =   x_left[idx_left]

    temp_right   =   Find_Nearest(y_right, half)
    idx_right    =   np.where(y_right==temp_right)
    x_half_right =   x_right[idx_right]

    return (x_half_right-x_half_left)

def Find_Nearest(array, value):
    # trova il valore più vicino a 'value'
    # all'interno di 'array'

    array   =   np.asarray(array)
    temp    =   np.abs(array-value)
       
    #print('stamo array')
    #print(array)
    #print('sto a stampa')
    idx     =   temp.argmin()

    del temp
    return (array[idx], idx)

def Find_Nearest_Array(array, value):

    # Versione della funzione precedente in cui
    # value è un array --> stessa operazione di prima ma per più valori
    
    array   =   np.asarray(array)
    idx     =   np.zeros(np.size(value), dtype=int)
    
    for ii in range (np.size(value)):

        temp        =   np.absolute(array-value[ii])
        idx[ii]     =   temp.argmin()

    return (array[idx], idx)

def Import_from_Matlab (mat_filename, path, transpose = True, **kwargs):

    """
    Funzione che importa da matlab, file 'path+matfilename'
    a seconda del tipo di dato, va strutturata con kwargs

    -   matlab cell     =       tensori, vuole una key per il dictionary
                                es.  loadmat(..)['y']

    """

    dati    =   loadmat(path+mat_filename)

    if 'var_name' in kwargs:

        dati    =   dati[kwargs['var_name']]
    
    if transpose: return dati.T
    else:     return dati


def Import_TIF(filename, path = './'):

    #   funzione che prende file TIF e lo riporta sotto forma di np.array

    image   =   Image.open(path+filename)
    imarray =   np.array(image)

    return  imarray

def Get_Neighbours2D(ii_0, jj_0, n_rows, n_cols):
    """
    Dato elemento matrice, ritorna tupla con gli 8 primi vicini
    """
    neighbours = []


    for ii in range(ii_0 - 1, ii_0+2, 1):
        for jj in range(jj_0 -1, jj_0 +2, 1):
            condition_inside = (ii >= 0) & (ii < n_rows) & (jj >= 0) & (jj < n_cols)
            if ((ii,jj) != (ii_0, jj_0)) & condition_inside:
                neighbours.append((ii,jj))
                
    return neighbours

"""
########################################################################################################
########################################################################################################

Funzioni matematiche assenti dalle librerie varie

"""

def gaussian  (x,A, mu, sigma):

    return (A*np.exp(-0.5*((x-mu)/sigma)**2))

def lorentian (x, mu, sigma, A):
    return (A/np.pi)*(sigma/((x - mu)**2 + sigma**2))

def delta_function (x, position, width, amplitude):

    return amplitude*np.exp(-((x-position)/width)**2/(np.abs(width)*np.sqrt(np.pi)))


import  numpy               as      np
import  matplotlib.pyplot   as      plt
from    PIL                 import  Image


"""                     
##############################################################################################
##############################################################################################

Funzioni generiche 
Insieme di funzioni utili per vari scopi in tutto il lavoro


"""

def Get_Around(value, factor):

    """
    Ritorna una tupla contenente (inf, sup) di un intervallo centrato su value
    la larghezza dell'intervallo è 2*factor*value --> factor tra (0,1)

    Utile per bounds nei fit

    """
    inf     =   value - (value*factor)
    sup     =   value - (value*factor)

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

def Import_from_Matlab ():
    pass

def Import_TIF(filename):

    #   funzione che prende file TIF e lo riporta sotto forma di np.array

    image   =   Image.open(filename)
    imarray =   np.array(image)

    return  imarray


"""
########################################################################################################
########################################################################################################

Funzioni matematiche assenti dalle librerie varie

"""


def lorentian (x, A, mu, sigma):
    return (A/np.pi)*(sigma/((x - mu)**2 + sigma**2))

def delta_function (x, width, amplitude):

    return amplitude*np.exp(-(x/width)**2/(np.abs(width)*np.sqrt(np.pi)))
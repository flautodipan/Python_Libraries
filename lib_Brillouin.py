###################################################################################################
###################################################################################################
#
# LIB_BRILLOUIN.py
#
# Principal library used by Brillouin.py,
# the script that performs fits of Brillouin experimental data
# 
# OSS:
# All the functions are provided of a brief explanation of their action in the 
# commented section at the beginning of their definition code
#
# 
# Contents:
#
# Part I:   class Spectrum
#           The reference class of a single data point, constituted by an element (ii,jj) 
#           of the matrix NxM that represents the whole z-stack acquisition
# 
# Part II:  External functions
#           The set of all functions that are useful to define, manage, visualize the data matrix
#           II a : 
#           II b : data visualization
#
###################################################################################################
###################################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


###################################################################################################
###################################################################################################





###################################################################################################
###################################################################################################

#######   ||||  ||||    EXTERNAL FUNCTIONS
#######    ||    ||    
#######    ||    ||
#######    ||    ||
#######   ||||  ||||



#######  II b: data visualization

def Plot_Elements_Spectrum(matrix, elements_iterable, fit = False, pix = False, peaks = False, x_range = (), y_range = ()):

    """
    This functions plots in different figures the spectrum of each matrix element contained 
    in the list elements_iterable

    Takes 
    -   the whole data matrix 
    -   the N indices of matrix elements (in the form of the list 
        elements_iterable = [(ii_0, jj_0), ..., (ii_N, jj_N)])
    -   the parameter "pix" to establish whether the data are in pixel or frequencies form
    -   the binary parameter "fit" to plot the performed fit (if so) together with the spectrum
    -   the binary parameter "peaks" to draw a star on each saved peak of the spectrum 
        (v. Spectrum.Get_Spectrum_Peaks)
    -   the x(y)_range parameter to specify a particular region for the plot

    """

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


###################################################################################################
###################################################################################################
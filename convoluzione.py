#%%
import  lib_Experimentum    as      Exp
from    Alessandria         import  *
from    Models              import  S_2_Generate, S_Dynamical_Form_Factor_2, S_Dynamical_Form_Factor_0, S_0_Generate

import pandas as pd
import  numpy               as      np
from    scipy.io            import  loadmat

from    matplotlib.pyplot   import  plot
import  matplotlib.pyplot   as      plt


import  time

#I/O 
now_path            =   '../BRILLOUIN/TDP43/ARS_13_02/'
spectra_filename    =   'ARS_13_02'
VIPA_filename       =   'NO_ARS_13_02_VIPA_not_sat.tif'


log_file            =   'log_'+spectra_filename
analysis_dir       =   'analysis_best/'

#operatives
syg_kwargs          =   {'height': 80, 'distance': 31, 'width': 3.}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
syg_kwargs_brill    =  {'height': 18, 'distance': 31, 'width': 3.}
VIPA_treshold       =   6
sat_height          =   50000
sat_width           =   13.5

#quanto mi allontano dal VIPA
pre_cut             =   False
cut                 =   False
cut_distance        =   0.25

cols_smart  =  ('Co', 'Omega', 'Gamma', 'delta_position',  'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_basic  = ('Co', 'Omega', 'Gamma', 'delta_position', 'delta_width',  'delta_amplitude')
cols        = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position',  'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_position', 'delta_width',  'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_position', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')

######
#%%
#0) importo dati e inizializzo oggetti
i   =   0
j   =   19
y    =   Import_from_Matlab(spectra_filename, now_path, var_name = 'y3')[i][j]
Data    =   Exp.Spectrum(name = str((i,j)))
Data.Get_Spectrum(y = np.resize(y, np.max(y.shape)), offset = 183., cut = pre_cut, cut_range = (10, 175))
Data.Get_Spectrum_Peaks(**syg_kwargs)
Data.Get_VIPA_tif(VIPA_filename, now_path, fig = 'VIPA_img', save_path = now_path)

plt.figure()
plt.plot(Data.x, Data.y)
plt.title(str((i,j)))

#%%
check = Data.Check_Spectrum_Saturation(saturation_height = sat_height, saturation_width = sat_width)

if      check == 1  : nature    = 'saturo'
elif    (Data.n_peaks == 2) :

    Data.Get_Spectrum_Peaks(**syg_kwargs_brill)
    if (Data.y[Data.peaks['idx'][2]] > Data.y[Data.peaks['idx'][1]]) & (Data.y[Data.peaks['idx'][2]] > Data.y[Data.peaks['idx'][3]]):
        nature    = 'brillouin_highest_dx'
    elif (Data.y[Data.peaks['idx'][1]] > Data.y[Data.peaks['idx'][2]]) & (Data.y[Data.peaks['idx'][1]] > Data.y[Data.peaks['idx'][0]]):
        nature    = 'brillouin_highest_sx'
    else:
        raise ValueError ('Non ho riconosciuto lo spettro\n')

elif (Data.n_peaks == 3) :
    Data.Get_Spectrum_Peaks(**syg_kwargs_brill) 
    nature = 'brillouin_higher'
else:                 nature    = 'bono'
print('Lo spettro è '+nature)
print('Lo spettro ha {} picchi'.format(Data.n_peaks))
plt.figure()
plt.plot(Data.x, Data.y)
plt.plot(Data.x[Data.peaks['idx']], Data.y[Data.peaks['idx']], '*')
plt.title(str((i,j)))

#%%
#1) Operazioni di modifica Spettro

Data.How_Many_Peaks_To_VIPA(treshold = VIPA_treshold, **syg_kwargs_VIPA)
Data.Fit_Pixel2GHz(fig = True)
Data.VIPA_Pix2GHz(fig = True)

if nature == 'brillouin_higher':
        Data.Get_Spectrum_Peaks(**syg_kwargs_brill)

Data.Spectrum_Pix2GHz(fig = True)
Data.Align_Spectrum()
Data.Cut_n_Estimate_Spectrum(estimate = True, cut = cut , verbose = True, distanza = cut_distance)
Data.Fit_VIPA_Gaussian(fig = True)

x_VIPA_interpolated = np.linspace(-35, 35, 2000)
y_VIPA_interpolated = Data.Interpolate_VIPA(x_VIPA_interpolated)
plt.plot(x_VIPA_interpolated, y_VIPA_interpolated)

#%%
_ = Data.Gauss_Convolve_Markovian_Response(Data.p0.values[0], fig = True, compare = True, fantoccio = (-30, 60))


# %%
from Alessandria import gaussian
p0 = np.array([1.00000000, 7.67396040e+00, 1.00000000e-01, 0.00000000e+00,
       0.010000000, 622.5462821, 1.19622722e+01, 1.67908016e+01,
       0.00000000e+00, 0.00000000e+00])
       
y_conv = Data.Gauss_Convolve_Markovian_Response_Smart_Fast(Data.p0[list(cols_smart)].values[0])#, fantoccio = (-100,100))
print(y_conv)
plt.figure()
plt.plot(Data.x_freq, Data.y, '*', label = 'data')
plt.plot(Data.x_freq, y_conv, '--', label = 'conv')
plt.plot(Data.x_freq, gaussian(Data.x_freq, p0[5], p0[6], p0[7]), label = 'gauss')
plt.legend()
# %%
from Models import S_Dynamical_Form_Factor_0_nodelta

x = np.linspace(-20, 20 , Data.x_VIPA.size)
y_VIPA = Data.Interpolate_VIPA(x)
y = y_VIPA + S_Dynamical_Form_Factor_0_nodelta(x, *p0[0:3])

plt.plot(x, y)

# %%
Data.Get_cost_markov(Data.p0[list(cols_smart)].values[0])

# %%

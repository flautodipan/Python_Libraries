#%%

from lib_Experimentum import Spectrum
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from    scipy.io            import loadmat
from    matplotlib.pyplot   import plot
import  matplotlib.pyplot   as plt
from lmfit import Model
from scipy.stats import norm

def gaussian  (x,A, mu, sigma):

    return (A*np.exp(-0.5*((x-mu)/sigma)**2))


# %%
VIPA = Spectrum()

VIPA.Get_VIPA_mat('VIPA.mat', path = '', tunable = 0, offset = 183., fig = 'VIPA_func')

#%%

VIPA.Fit_Pixel2GHz(altezza = 1, fig = 'Fit_Pix2GHz')
VIPA.VIPA_Pix2GHz(fig='Conversione_VIPA')


#%%
Picchi = find_peaks(VIPA.y_VIPA, height=9, width=0.001)

gmod = Model(gaussian)
result = gmod.fit(Picchi[1]['peak_heights'], x = VIPA.x_VIPA_freq[Picchi[0]], A = 1., mu = 1, sigma = 1)
print(result.fit_report)
A = result.values['A']
mu = result.values['mu']
sigma = result.values['sigma']

x = np.linspace(VIPA.x_VIPA_freq[Picchi[0]].min(), VIPA.x_VIPA_freq[Picchi[0]].max(), 1000)
y = gaussian(x, A, mu, sigma)

plot(VIPA.x_VIPA_freq[Picchi[0]], Picchi[1]['peak_heights'], '*')
plot(x,y)
plot(VIPA.x_VIPA_freq, VIPA.y_VIPA)
plt.xlim(-50,50)
plt.show()

print(result.values)

# %%

#%%
import numpy as np
import matplotlib.pyplot as plt
from lib_Experimentum import Print_Parameter_Map
from Alessandria import gaussian
from scipy.optimize import curve_fit

path = '../BRILLOUIN/TDP43/'
save_path = '../Scrittek/figures/'


#%%
# prendo i dati

name_NO_ARS_48 = '-ARS 48h'
pix_scale_NO_ARS_48 = 350#nm
Omegas_NO_ARS_48 = np.load(path +'markov_Omegas'+name_NO_ARS_48+'.npy', )
Gammas_NO_ARS_48 = np.load(path +'markov_Gammas'+name_NO_ARS_48+'.npy', )
omega_map_NO_ARS_48 = np.load(path +'markov_omega_map'+name_NO_ARS_48+'.npy', )
gamma_map_NO_ARS_48 = np.load(path +'markov_gamma_map'+name_NO_ARS_48+'.npy',)


name_ARS_48 = '+ARS 48h'

pix_scale_ARS_48 = 250#nm
Omegas_ARS_48 = np.load(path +'markov_Omegas'+name_ARS_48+'.npy', )
Gammas_ARS_48 = np.load(path +'markov_Gammas'+name_ARS_48+'.npy', )
omega_map_ARS_48 = np.load(path +'markov_omega_map'+name_ARS_48+'.npy', )
gamma_map_ARS_48 = np.load(path +'markov_gamma_map'+name_ARS_48+'.npy',)


name_ARS_72 = '+ARS 72h'

pix_scale_ARS_72 = 240#nm

Omegas_ARS_72 = np.load(path +'markov_Omegas'+name_ARS_72+'.npy', )
Gammas_ARS_72 = np.load(path +'markov_Gammas'+name_ARS_72+'.npy', )
omega_map_ARS_72 = np.load(path +'markov_omega_map'+name_ARS_72+'.npy', )
gamma_map_ARS_72 = np.load(path +'markov_gamma_map'+name_ARS_72+'.npy',)
#%%
# Omega distributions fit gaussiano

omega_max = 7.9
omega_min = 7.3


# NO ARS 48

y, x = np.histogram(Omegas_NO_ARS_48, bins = 50)
x = (x[1:] + x[:-1])/2
x = x[0:20]
y = y[0:20]
p_gauss_NO_ARS_48, _ = curve_fit(gaussian, x, y, p0 = [800, 7.5, .091] )
plt.figure()
plt.title(name_NO_ARS_48+'Omega distribution over sample')
bins = plt.hist(Omegas_NO_ARS_48, bins = 50, label = 'omega',histtype = 'bar', stacked = True, rwidth= 0.8,  color = 'maroon')
x_ = np.linspace(7.2, 7.5)
plt.plot(x, gaussian(x, *p_gauss_NO_ARS_48), label = 'water peak $\Omega$ = {:3.2f}'.format(p_gauss_NO_ARS_48[1]), lw = 3., c = 'yellowgreen')
omega_mean_NO_ARS_48 = np.mean(Omegas_NO_ARS_48)
omega_std_NO_ARS_48 = np.std(Omegas_NO_ARS_48)
plt.xlim(omega_min -0.2, omega_max + 0.2 )
plt.legend(title='Mean = {:3.2f}\nStd = {:3.2f}'.format(omega_mean_NO_ARS_48, omega_std_NO_ARS_48,))
plt.show()


#ARS 48h

y, x = np.histogram(Omegas_ARS_48, bins = 50)
x = (x[1:] + x[:-1])/2
x = x[3:22]
y = y[3:22]
p_gauss_ARS_48, _ = curve_fit(gaussian, x, y, p0 = [800, 7.45, .08] )
plt.figure()
plt.title(name_ARS_48+'Omega distribution over sample')
bins = plt.hist(Omegas_ARS_48, bins = 50, label = 'omega',histtype = 'bar', stacked = True, rwidth= 0.8,  color = 'maroon')
x_ = np.linspace(7.2, 7.5)
plt.plot(x, gaussian(x, *p_gauss_ARS_48), label = 'water peak $\Omega$ = {:3.2f}'.format(p_gauss_ARS_48[1]), lw = 3., c = 'yellowgreen')
omega_mean_ARS_48 = np.mean(Omegas_ARS_48)
omega_std_ARS_48 = np.std(Omegas_ARS_48)
plt.xlim(omega_min -0.2, omega_max + 0.2 )
plt.legend(title='Mean = {:3.2f}\nStd = {:3.2f}'.format(omega_mean_ARS_48, omega_std_ARS_48,))
plt.xlim(omega_min -0.2, omega_max + 0.2 )
plt.show()


#72h
y, x = np.histogram(Omegas_ARS_72, bins = 50)
x = (x[1:] + x[:-1])/2
x = x[5:22]
y = y[5:22]
p_gauss_ARS_72, _ = curve_fit(gaussian, x, y, p0 = [800, 7.45, .1] )


plt.figure()
plt.title(name_ARS_72+'Omega distribution over sample')
bins = plt.hist(Omegas_ARS_72, bins = 50, label = 'omega',histtype = 'bar', stacked = True, rwidth= 0.8,  color = 'maroon')
x_ = np.linspace(7.2, 7.5)
plt.plot(x, gaussian(x, *p_gauss_ARS_72), label = 'water peak $\Omega$ = {:3.2f}'.format(p_gauss_ARS_72[1]), lw = 3., c = 'yellowgreen')
omega_mean_ARS_72 = np.mean(Omegas_ARS_72)
omega_std_ARS_72 = np.std(Omegas_ARS_72)
plt.xlim(omega_min -0.2, omega_max + 0.2 )
plt.legend(title='Mean = {:3.2f}\nStd = {:3.2f}'.format(omega_mean_ARS_72, omega_std_ARS_72,))
plt.xlim(omega_min -0.2, omega_max + 0.2 )
plt.show()



#%%
#stampa mappe in Scrittek
ref = p_gauss_ARS_72[1]
Print_Parameter_Map(omega_map_ARS_72, omega_min, omega_max, 'Omega', 'markov', name_ARS_72, pix_scale_ARS_72, 'markov_omega_ARS_72', save_path)

shift = ref - p_gauss_ARS_48[1]
omega_map_ARS_48_new = omega_map_ARS_48 + shift
Print_Parameter_Map(omega_map_ARS_48, omega_min, omega_max, 'Omega', 'markov', name_ARS_48, pix_scale_ARS_48, 'markov_omega_ARS_48', save_path)

shift = ref - p_gauss_NO_ARS_48[1]
omega_map_NO_ARS_48_new = omega_map_NO_ARS_48 + shift
Print_Parameter_Map(omega_map_NO_ARS_48, omega_min, omega_max, 'Omega', 'markov', name_NO_ARS_48, pix_scale_NO_ARS_48, 'markov_omega_NO_ARS_48', save_path)



# %%


#%%
# Gamma distributions

gamma_max = 0.235
gamma_min = 0.1

# NO ARS 48

plt.figure()
plt.title(name_NO_ARS_48+' Gamma distribution over sample')
bins = plt.hist(Gammas_NO_ARS_48, bins = 50, label = 'gamma',histtype = 'bar', stacked = True, rwidth= 0.8,  color = 'royalblue')
x_ = np.linspace(7.2, 7.5)
gamma_mean_NO_ARS_48 = np.mean(Gammas_NO_ARS_48)
gamma_std_NO_ARS_48 = np.std(Gammas_NO_ARS_48)
plt.xlim(gamma_min -0.03, gamma_max + 0.03 )
plt.legend(title='Mean = {:3.2f}\nStd = {:3.2f}'.format(gamma_mean_NO_ARS_48, gamma_std_NO_ARS_48,))
plt.show()


# ARS 48

plt.figure()
plt.title(name_ARS_48+' Gamma distribution over sample')
bins = plt.hist(Gammas_ARS_48, bins = 50, label = 'gamma',histtype = 'bar', stacked = True, rwidth= 0.8,  color = 'royalblue')
x_ = np.linspace(7.2, 7.5)
gamma_mean_ARS_48 = np.mean(Gammas_ARS_48)
gamma_std_ARS_48 = np.std(Gammas_ARS_48)
plt.xlim(gamma_min -0.03, gamma_max + 0.03 )
plt.legend(title='Mean = {:3.2f}\nStd = {:3.2f}'.format(gamma_mean_ARS_48, gamma_std_ARS_48,))
plt.show()

# ARS 48

plt.figure()
plt.title(name_ARS_72+' Gamma distribution over sample')
bins = plt.hist(Gammas_ARS_72, bins = 50, label = 'gamma',histtype = 'bar', stacked = True, rwidth= 0.8,  color = 'royalblue')
x_ = np.linspace(7.2, 7.5)
gamma_mean_ARS_72 = np.mean(Gammas_ARS_72)
gamma_std_ARS_72 = np.std(Gammas_ARS_72)
plt.xlim(gamma_min -0.03, gamma_max + 0.03 )
plt.legend(title='Mean = {:3.2f}\nStd = {:3.2f}'.format(gamma_mean_ARS_72, gamma_std_ARS_72,))
plt.show()


# %%

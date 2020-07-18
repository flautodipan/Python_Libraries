#%%
path = '../Scrittek/figures/'
from Models import *
import numpy as np
import matplotlib.pyplot as plt

#%%
### ESEMPIO MODELLO MARKOVIANO


x = np.linspace(-15, 15, 1000)

p1 = (7.9, 0.2)
p2 = (6.5, 0.3)
p3 =  (4.3, 0.2)
p4 = (8.7, 0.3,)

y1 = S_Dynamical_Form_Factor_0_nodelta(x, 10., *p1)
y2 = S_Dynamical_Form_Factor_0_nodelta(x, 10., *p2)
y3 = S_Dynamical_Form_Factor_0_nodelta(x, 10., *p3)
y4 = S_Dynamical_Form_Factor_0_nodelta(x, 10., *p4)


fig, axs = plt.subplots(2, 2, constrained_layout=True)
fig.suptitle(r'Markovian model for $I(\omega) \propto S_{\rho\rho}(\omega)$', fontsize = 14)


axs[0,0].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(*p1), fontsize = 8.5)
axs[0,0].set_ylabel('Intensity')
axs[0, 0].plot(x, y1, ls = 'solid', color = 'darkgoldenrod')
plt.ylim(-0.5, 4.5)


axs[0,1].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(*p2), fontsize = 8.5)
plt.ylim(-0.5, 4.5)
axs[0,1].plot(x, y2, ls = 'solid', color = 'red')
axs[1,0].set_ylabel('Intensity')

axs[1,0].plot(x, y3, ls = 'solid', color = 'cornflowerblue')
axs[1,0].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(*p3), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,0].set_xlabel('$\omega$ (GHz)')

axs[1,1].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(*p4), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,1].plot(x, y4, ls = 'solid', color = 'darkolivegreen', label = '$\Omega$ = {}\n$\Gamma$ = {}'.format(7.4, 0.27))
plt.xlabel('$\omega$ (GHz)')

#fig.title('Theoretical function by Markovian model')
#fig.tight_layout()
plt.show()
plt.close()
#fig.savefig(path+'Brillouin_markov_example.pdf', format = 'pdf')

print(path)
#%%
#%%
### ESEMPIO MODELLO RELAXATION


x = np.linspace(-15, 15, 1000)
p1 = (7.9, 0.2, 0.8, 100000)
p2 = (6.5, 0.3, 0.8, 100000)
p3 =  (4.3, 0.2,0.8, 100000)
p4 = (8.7, 0.3, 0.8, 100000)
y1 = S_Dynamical_Form_Factor_2_nodelta(x, 10., *p1)
y2 = S_Dynamical_Form_Factor_2_nodelta(x, 10., *p2)
y3 = S_Dynamical_Form_Factor_2_nodelta(x, 10., *p3)
y4 = S_Dynamical_Form_Factor_2_nodelta(x, 10., *p4)


fig, axs = plt.subplots(2, 2, constrained_layout=True)
fig.suptitle(r'Relaxation model for $I(\omega) \propto S_{\rho\rho}(\omega)$', fontsize = 14)


axs[0,0].set_title(r'$\Omega = {} \;\; \Gamma = {} \;\;\Delta = {} \;\; \tau = {}$'.format(*p1), fontsize = 8.5)
axs[0,0].set_ylabel('Intensity')
axs[0, 0].plot(x, y1, ls = 'solid', color = 'darkgoldenrod')
plt.ylim(-0.5, 4.5)


axs[0,1].set_title(r'$\Omega = {} \;\; \Gamma = {} \;\;\Delta = {} \;\; \tau = {}$'.format(*p2), fontsize = 8.5)
plt.ylim(-0.5, 4.5)
axs[0,1].plot(x, y2, ls = 'solid', color = 'red')
axs[1,0].set_ylabel('Intensity')

axs[1,0].plot(x, y3, ls = 'solid', color = 'cornflowerblue')
axs[1,0].set_title(r'$\Omega = {} \;\; \Gamma = {} \;\;\Delta = {} \;\; \tau = {}$'.format(*p3), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,0].set_xlabel('$\omega$ (GHz)')

axs[1,1].set_title(r'$\Omega = {} \;\; \Gamma = {} \;\;\Delta = {} \;\; \tau = {}$'.format(*p4), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,1].plot(x, y4, ls = 'solid', color = 'darkolivegreen', label = '$\Omega$ = {}\n$\Gamma$ = {}'.format(7.4, 0.27))
plt.xlabel('$\omega$ (GHz)')

#fig.title('Theoretical function by Markovian model')
#fig.tight_layout()
plt.show()
plt.close()
#fig.savefig(path+'Brillouin_tot_example.pdf', format = 'pdf')
# %%
# MODELLO RILASSAMENTO

path = '../Scrittek/figures/'
from Models import *
import numpy as np
import matplotlib.pyplot as plt

# %%
#Spettro Rayleigh Brillouin



x = np.linspace(-15, 15, 1000)
y1 = S_Dynamical_Form_Factor_2(x, 10., 7.9, 0.2, 0.3, 0.01, 0, 8, 0.6)
f, ax = plt.subplots()
ax.plot(x, y1, c = 'yellowgreen')
ax.set_title('Example of Rayleigh Brillouin Spectrum')
ax.set_xlabel('$\omega$ - frequency shift (GHz)')
ax.set_ylabel('Intensity')
ax.set_yticklabels([])
plt.savefig(path+'spectrum_example.pdf', format = 'pdf')

#%%

#differenti valori di Omega tau

x = np.linspace(-15, 15, 1000)
Omega = 7.9
tau1 = 0.0002
tau2 = 0.02
tau3 = 20

y1 = S_Dynamical_Form_Factor_2(x, 10., Omega, 0.2, 0.3, tau1, 0, 8, 0.3)
y2 = S_Dynamical_Form_Factor_2(x, 10., Omega, 0.2, 0.3, tau2, 0, 8, 0.3)
y3 = S_Dynamical_Form_Factor_2(x, 10., Omega, 0.2, 0.3, tau3, 0, 8, 0.3)

f , ax = plt.subplots()
ax.plot(x, y1, label = r'$\Omega \tau \ll 1$')
ax.plot(x, y2, label = r'$\Omega \tau \simeq 1$')
ax.plot(x, y3, label = r'$\Omega \tau \gg 1$', c = 'yellowgreen')

ax.legend()
ax.set_title('Rayleigh Brillouin Spectrum\nfor different values of the relaxation time')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('Frequency shift')
ax.set_ylabel('Intensity')
#f.savefig(path+'brillouin_different_tau.pdf', format = 'pdf')
plt.show()

#%%

#differenti valori di Omega tau
# TRUCCATO

x = np.linspace(-15, 15, 1000)
Omega = 7.9
tau1 = 0.0002
tau2 = 0.02
tau3 = 20

y1 = S_Dynamical_Form_Factor_2(x, 10., Omega-2, 0.2, 0.2, tau1, 0, 8, 0.3)
y2 = S_Dynamical_Form_Factor_2(x, 10., Omega, 0.25, 0.3, tau2, 0, 8, 0.3)
y3 = S_Dynamical_Form_Factor_2(x, 10., Omega+2, 0.2, 0.2, tau3, 0, 8, 0.3)

"""
y1 = S_Dynamical_Form_Factor_2_nodelta(x, 10., Omega, 0.2, 0.3, tau1,)
y2 = S_Dynamical_Form_Factor_2_nodelta(x, 10., Omega, 0.2, 0.3, tau2,)
y3 = S_Dynamical_Form_Factor_2_nodelta(x, 10., Omega, 0.2, 0.3, tau3,)
"""
f , ax = plt.subplots()
ax.plot(x, y1, label = r'$\Omega \tau \ll 1$')
ax.plot(x, y2, label = r'$\Omega \tau \simeq 1$')
ax.plot(x, y3, label = r'$\Omega \tau \gg 1$', c = 'yellowgreen')

ax.legend()
ax.set_title('Rayleigh Brillouin Spectrum\nfor different values of the relaxation time')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('Frequency shift')
ax.set_ylabel('Intensity')
#f.savefig(path+'brillouin_different_tau.pdf', format = 'pdf')
# %%

take_path = '../GROMACS/'
from BioAlessandria import Parse_xvg_skip

for what in ['potential', 'temp', 'press']:

    time, y = Parse_xvg_skip(what+'.xvg', take_path, skip_lines=24)
        
    plt.figure()
    if what[1] == 'o':
        y /= 1000
    plt.plot(time,y, label = what, c = 'darkred')
    plt.legend()
    if what[1] == 'o':
            plt.title('Energy Equilibration')
            plt.ylabel('Pot Energy (MJ/mol)')

    elif what[1] == 'e':
        plt.title('Temperature Equilibration')
        plt.ylabel('Temperature (K)')
    else:
        plt.title('Pressure Equilibration')
        plt.ylabel('Pressure (bar)')

    plt.xlabel('Time (ps)')
    plt.savefig(path+what+'_eq.pdf', format = 'pdf')


# %%
# SPREADING PLOT
Kds_errs = np.array([0.9, 37, 150, 135,150, 125, 150, 450, 350,  600, 500, 400 ])
Kds = np.array([4, 105, 650, 650,  700, 750, 850, 1300, 1320, 1360, 1500, 1800])
names = [r'$RNA_{12}$', 'Apt1.2', 'Apt2.2', 'Apt3.1','Apt3.2', 'Apt2.1', 'Apt1.1', 'Apt4.2', 'NegApt1.2',  r'$NegRNA_{12}$', 'Apt5.2', 'Apt6.2']

f, ax = plt.subplots()
idx=np.arange(0,len(Kds), dtype=int)
md_idx = np.array([0, 3, 4, 8, 9, 11], dtype=int)

ax.bar(idx, Kds, yerr = Kds_errs,  color = 'yellowgreen', capsize = 3, alpha = 0.7)
ax.bar(md_idx, Kds[list(md_idx)], yerr = Kds_errs[list(md_idx)], color = 'firebrick', alpha = 0.7, label = 'runned in MD')
#ax.scatter(idx, Kds, c = colormap[np.array(colors, dtype=int)],)
ax.set_ylim(-50, 2500)
ax.set_xticks(idx)
ax.set_xticklabels(names)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
ax.set_ylabel('$K_d$ (nM)')
ax.set_title("RNA Aptamers dissociation constant $K_d$", fontsize=12.5, pad = 10)
ax.legend(loc = 'upper left')
plt.tight_layout()

f.savefig(path+'Kd_spreading_plot.pdf', format = 'pdf')

# %%

# NMR - test comparison - initial 

take_path = '../GROMACS/'
from BioAlessandria import Parse_xvg_skip

RMSD_comp = Parse_xvg_skip(take_path+'NMR_comparison.xvg', skip = 50)

RMSD_NMR = Parse_xvg_skip(take_path+'/WTC1/rmsd_wtc1.xvg', skip = 50)
RMSD_test = Parse_xvg_skip(take_path+'/WTC1_h/rmsd_wtc1_h.xvg', skip = 50)

f, ax = plt.subplots()
ax.set_title('RMSD between NMR and WTC1 structures in time', fontsize = 13.5, pad = 10)
#ax.hlines(np.mean(RMSD_eq), 0, 2000, color = 'k', linestyles = 'dashed', label = 'equilibrium RMSD mean')
#ax.fill_betweenx(np.arange(ax.get_ylim()[0], 1.55, 0.1), 1000, 2000, color = 'firebrick', alpha = 0.2, label = 'Equilibrium time range')
ax.plot(RMSD_comp[0]/1000, RMSD_comp[1], label = 'WTC1 vs NMR mean', color = 'firebrick', alpha = 0.8)
ax.plot(RMSD_NMR[0][1:]/1000, RMSD_NMR[1][1:], label = 'NMR vs NMR initial', color = 'goldenrod', alpha = 0.8)
ax.plot(RMSD_test[0][2:]/1000, RMSD_test[1][2:], label = 'WTC1 vs WTC1 initial', color = 'green', alpha = 0.8)

ax.set_xlim(0,1000)
ax.set_ylim(ax.get_ylim()[0], 1.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('RMSD (nm)')
ax.legend()
plt.tight_layout()

f.savefig(path+'NMR_comparison.pdf', format = 'pdf')


# %%
# %%
#AIRY

def FP_airy(delta, R):

    F = 4*R/(1-R)**2
    finesse = (np.pi*np.sqrt(R))/(1-R)
    return (1/(1 + F*np.sin(delta)**2), R, finesse)

import matplotlib.pyplot as plt
import numpy as np

f, ax = plt.subplots()

x = np.linspace(-10,10, 1000)
ax.set_title('FP etalon intensity profile')
ax.set_ylabel('Transmission $I/I_0$')
y = FP_airy(x, 0.6)
ax.plot(x, y[0], label = 'R = {:3.1f}\nfin = {:1d}'.format(y[1], int(y[2])), color = 'darkgoldenrod')
y = FP_airy(x, 0.2)
ax.plot(x, y[0], label = 'R = {:3.1f}\nfin = {:1d}'.format(y[1],int(y[2])), color = 'darkred', ls = '--')
y = FP_airy(x, 0.9)
ax.plot(x, y[0], label = 'R = {:3.1f}\nfin = {:1d}'.format(y[1], int(y[2])), color = 'royalblue', ls='dotted')

ax.set_xticks([])
#ax.xaxis.set_visible(False)
ax.set_xlabel('Frequency')
ax.legend(loc = 'lower left')
plt.xlim(-4.5, 4.5)
plt.savefig(path+'FP_airy.pdf', format = 'pdf', bbox_to_inches = 'tight')



# %%
#FABRY PEROT

def FP_airy(delta, R):

    F = 4*R/(1-R)**2
    finesse = (np.pi*np.sqrt(R))/(1-R)
    return (1/(1 + F*np.sin(delta)**2), R, finesse)

import matplotlib.pyplot as plt
import numpy as np

#%%
# FP risoluzione finesse

f, ax = plt.subplots(1, 2)

x = np.linspace(-10,10, 1000)
y = FP_airy(x, 0.2)
y_1 = FP_airy(x + 0.5, 0.2)
ax[0].plot(x, y[0], color = 'lightcoral', ls='dotted', label = r'freq $\nu_1$')
ax[0].plot(x, y_1[0], color = 'darkred', ls='dashed', label = r'freq $\nu_2$')
ax[0].plot(x, y[0]+y_1[0] - 0.2, label = 'Total field')


ax[0].set_xticks([])
ax[0].set_yticks([])
ax[0].set_ylabel('Intensity')
ax[0].set_xlabel('Frequency')

ax[0].set_title('Finesse = {:1d}'.format(int(y[2])))


ax[0].legend()
ax[0].set_xlim(-4.8, 1.8)
plt.savefig(path+'FP_bad_resolution.pdf', format = 'pdf', bbox_to_inches = 'tight')

x = np.linspace(-10,10, 1000)
y = FP_airy(x, 0.9)
y_1 = FP_airy(x + 0.5, 0.9)
ax[1].plot(x, y[0], color = 'darkorange', ls='dashed', label = r'freq $\nu_1$')
ax[1].plot(x, y_1[0], color = 'darkred', ls='dashed', label = r'freq $\nu_2$')
ax[1].plot(x, y[0]+y_1[0], label = 'Total field', alpha = 0.6)


ax[1].set_xticks([])
ax[1].set_yticks([])
ax[1].set_ylabel('Intensity')
ax[1].set_xlabel('Frequency')

ax[1].set_title('Finesse = {:1d}'.format(int(y[2])))


ax[1].legend()
ax[1].set_xlim(-4.8, 1.8)

#plt.tight_layout(rect=(0,0,2,1.25))
f.savefig(path+'FP_resolution.pdf', format = 'pdf', bbox_to_inches='tight')

# %%
# GAUSSIAN ENVELOPE

from Alessandria import gaussian
f, ax = plt.subplots(1, 1)

x = np.linspace(-10,10, 1000)

y = FP_airy(x, 0.9)[0]*gaussian(x, 0.991, 0, 3.55)
ax.plot(x, y, color = 'yellowgreen', ls='dashed', label = 'Airy function')
ax.plot(x, gaussian(x, 0.961, 0, 3.55), ls= 'solid', label = 'Gaussian envelope G($\omega$)', color = 'firebrick', alpha = .8)

plt.legend(loc = 'upper right')
ax.set_ylabel('Intensity')
ax.set_xlabel('Frequency shift $\omega$')


ax.set_xticks([])
ax.set_yticks([])
plt.savefig(path+'VIPA_Gauss_envelop.pdf', format = 'pdf', quality = 100)
# %%
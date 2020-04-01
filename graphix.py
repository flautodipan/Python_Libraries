#%%
path = '../Scrittek/figures/'
from Models import *
import numpy as np
import matplotlib.pyplot as plt

#%%
### ESEMPIO MODELLO MARKOVIANO


x = np.linspace(-15, 15, 1000)
y1 = S_Dynamical_Form_Factor_0_nodelta(x, 10., 7.9, 0.2)
y2 = S_Dynamical_Form_Factor_0_nodelta(x, 10., 6.5, 0.35)
y3 = S_Dynamical_Form_Factor_0_nodelta(x, 10., 1.4, 0.27)
y4 = S_Dynamical_Form_Factor_0_nodelta(x, 10., 8.7, 0.3)


fig, axs = plt.subplots(2, 2, constrained_layout=True)
fig.suptitle('Dynamical form factors in Markovian approximation', fontsize = 15)


axs[0,0].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(7.9, 0.2), fontsize = 8.5)
axs[0,0].set_ylabel('Intensity')
axs[0, 0].plot(x, y1, ls = 'solid', color = 'darkgoldenrod')
plt.ylim(-0.5, 4.5)


axs[0,1].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(6.5, 0.35), fontsize = 8.5)
plt.ylim(-0.5, 4.5)
axs[0,1].plot(x, y2, ls = 'solid', color = 'red')
axs[1,0].set_ylabel('Intensity')

axs[1,0].plot(x, y3, ls = 'solid', color = 'cornflowerblue')
axs[1,0].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(1.4, 0.27), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,0].set_xlabel('$\omega$ (GHz)')

axs[1,1].set_title('$\Omega$ = {}\t$\Gamma$ = {}'.format(8.7, 0.3), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,1].plot(x, y4, ls = 'solid', color = 'darkolivegreen', label = '$\Omega$ = {}\n$\Gamma$ = {}'.format(7.4, 0.27))
plt.xlabel('$\omega$ (GHz)')

#fig.title('Theoretical function by Markovian model')
#fig.tight_layout()
plt.show()
plt.close()
fig.savefig(path+'Brillouin_markov_example.pdf', format = 'pdf')

print(path)
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
plt.figure()
plt.plot(x, y1)
plt.title('Example of Rayleigh Brillouin Spectrum')
plt.xlabel('$\omega$ - frequency shift (GHz)')
#plt.set_yticklabels([])
plt.savefig(path+'spectrum_example.pdf', format = 'pdf')


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

Kds = np.array([4, 105, 650, 700, 750, 850, 1300, 1320, 1350, 1360, 1500, 1800])
names = ['RNA12', 'Apt1.2', 'Apt2.2', 'Apt3.2', 'Apt2.1', 'Apt1.1', 'Apt4.2', 'NegApt1.2', 'Apt3.1', 'NegRNA12', 'Apt5.2', 'Apt6.2']
md_idx = np.array([0, 3, 7, 8, 9, 11], dtype=int)
idx=np.arange(0,len(Kds), dtype=int)
other_idx = np.array([1,2,4,5,6,10], dtype = int)
colors = []
labels = []
colormap = np.array(['dodgerblue', 'orange'])
for ii in idx:
    if ii in md_idx: 
        colors.append(1)
        labels.append('Runned in MD')

    else : 
        colors.append(0)
        labels.append(None)

f, ax = plt.subplots()
ax.scatter(idx, Kds, c = colormap[np.array(colors, dtype=int)])
ax.set_xticks(idx)
ax.set_xticklabels(names)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
ax.set_ylabel('$K_d$ (nM)')
ax.set_title("Aptamer's dissociation constant $K_d$", fontsize=12.5)
ax.legend(labels)
plt.tight_layout()

f.savefig(path+'Kd_spreading_plot.pdf', format = 'pdf')

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
from Alessandria import gaussian
f, ax = plt.subplots(1, 1)

x = np.linspace(-10,10, 1000)

y = FP_airy(x, 0.9)[0]*gaussian(x, 0.991, 0, 3.55)
ax.plot(x, y, color = 'royalblue', ls='dashed', label = 'Airy function')
ax.plot(x, gaussian(x, 0.961, 0, 3.55), ls= 'solid', label = 'Gaussian envelop', color = 'darkorange')

plt.legend()
ax.set_ylabel('Intensity')
ax.set_xlabel('Frequency')


ax.set_xticks([])
ax.set_yticks([])
plt.savefig(path+'VIPA_Gauss_envelop.pdf', format = 'pdf', quality = 100)
# %%

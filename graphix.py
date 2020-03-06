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



x = np.linspace(-15, 15, 1000)
y1 = S_Dynamical_Form_Factor_2_nodelta(x, 10., 7.9, 0.2, 0.3, 0.01)
y2 = S_Dynamical_Form_Factor_2_nodelta(x, 10., 6.5, 0.35, 0.8, 0.1)
y3 = S_Dynamical_Form_Factor_2_nodelta(x, 10., 1.4, 0.27, 0.1, 10)
y4 = S_Dynamical_Form_Factor_2_nodelta(x, 10., 8.7, 0.3, 0.9, 1)


fig, axs = plt.subplots(2, 2, constrained_layout=True)
fig.suptitle('Dynamical form factors in viscoelastic approximation', fontsize = 15)


axs[0,0].set_title('$\Omega$ = {}\t$\Gamma$={}\t$\Delta$={}\t \u03C4={}'.format(7.9, 0.2,0.3, 0.1), fontsize = 8.5)
axs[0,0].set_ylabel('Intensity')
axs[0, 0].plot(x, y1, ls = 'solid', color = 'darkgoldenrod')
plt.ylim(-0.5, 4.5)


axs[0,1].set_title('$\Omega$ = {}\t$\Gamma$ = {}\t$\Delta$={}\t \u03C4={}'.format(6.5, 0.35, 0.8, 0.1), fontsize = 8.5)
plt.ylim(-0.5, 4.5)
axs[0,1].plot(x, y2, ls = 'solid', color = 'red')
axs[1,0].set_ylabel('Intensity')

axs[1,0].plot(x, y3, ls = 'solid', color = 'cornflowerblue')
axs[1,0].set_title('$\Omega$ = {}\t$\Gamma$ = {}\t$\Delta$={}\t \u03C4={}'.format(1.4, 0.27, 0.1, 10), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,0].set_xlabel('$\omega$ (GHz)')

axs[1,1].set_title('$\Omega$ = {}\t$\Gamma$ = {}\t$\Delta$={}\t \u03C4={}'.format(8.7, 0.3, 0.9, 10), fontsize = 8.5)
plt.ylim(-0.5, 4.5)

axs[1,1].plot(x, y4, ls = 'solid', color = 'darkolivegreen', label = '$\Omega$ = {}\n$\Gamma$ = {}'.format(7.4, 0.27))
plt.xlabel('$\omega$ (GHz)')

#fig.title('Theoretical function by Markovian model')
#fig.tight_layout()
plt.show()
plt.close()
fig.savefig(path+'Brillouin_real_example.pdf', format = 'pdf')

# %%

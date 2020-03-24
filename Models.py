import  numpy               as      np
import  matplotlib.pyplot   as      plt
from    Alessandria         import  *


#       PARTE   1.1     Funzioni Teoriche del Modello
#                       derivate dalla risoluzione dell'equazione di memoria
#                       con il formalismo Mori-Zwanzig
#                       S(k,w) è dynamic form factor, la I (k, w) è la corrente I(k,w)=S(k,w)*w**2
#
#       Notazione:      (omega-shift)   =   variabile indipendente frequenza in GHz
#
#       parametri       Omega   =   frequenza caratteristica        
#                       Gamma   =   fattore di larghezza parte costante della funz di memoria
#                       Delta   =   fatt di larghezza parte decadimento esponenziale della funz di memoria
#                       tau     =   costante di tempo del decadimento esponenziale                     
#                       Co      =   S(k, t=0) è il prefattore della funzione
#                                   OSS: è il valore nei tempi, non freq
# 

def S_Dynamical_Form_Factor_0_nodelta(omega, Co, Omega, Gamma):
    return (Co * Gamma * (4 * Gamma ** 2 + Omega ** 2 + omega ** 2) / (Gamma ** 2 + (Omega - omega) ** 2 / 4) / np.pi / (Gamma ** 2 + (Omega + omega) ** 2 / 4) / 8)


def S_Dynamical_Form_Factor_0(omega, Co, Omega, Gamma, delta_position, delta_width, delta_amplitude):
    return delta_function(omega, delta_position, delta_width, delta_amplitude)+(Co * Gamma * (4 * Gamma ** 2 + Omega ** 2 + omega ** 2) / (Gamma ** 2 + (Omega - omega) ** 2 / 4) / np.pi / (Gamma ** 2 + (Omega + omega) ** 2 / 4) / 8)

def I_Dynamical_Current_0(omega, Co, Omega, Gamma):
        return (omega**2)*(Co * Gamma * (4 * Gamma ** 2 + Omega ** 2 + omega ** 2) / (Gamma ** 2 + (Omega - omega) ** 2 / 4) / np.pi / (Gamma ** 2 + (Omega + omega) ** 2 / 4) / 8)

# --> caso rilassamento esponenziale    :       K(p) proporz exp(-t/tau)            'Exp'

def S_Dynamical_Form_Factor_1(omega, Co, Omega, Delta, tau):
    return (2 * Co * (omega ** 4 * tau ** 2 + (tau ** 2 * Omega ** 2 - 4 * Delta * tau + 1) * omega ** 2 + 4 * Delta ** 2 + Omega ** 2) * Delta / (omega ** 4 * tau ** 2 + 2 * tau ** 2 * Omega * omega ** 3 + (tau ** 2 * Omega ** 2 - 4 * Delta * tau + 1) * omega ** 2 + (-4 * tau * Delta * Omega + 2 * Omega) * omega + 4 * Delta ** 2 + Omega ** 2) / (omega ** 4 * tau ** 2 - 2 * tau ** 2 * Omega * omega ** 3 + (tau ** 2 * Omega ** 2 - 4 * Delta * tau + 1) * omega ** 2 + (4 * tau * Delta * Omega - 2 * Omega) * omega + 4 * Delta ** 2 + Omega ** 2) / np.pi)

def I_Dynamical_Current_1(omega, Co, Omega, Delta, tau):
    return (omega**2)*(2 * Co * (omega ** 4 * tau ** 2 + (tau ** 2 * Omega ** 2 - 4 * Delta * tau + 1) * omega ** 2 + 4 * Delta ** 2 + Omega ** 2) * Delta / (omega ** 4 * tau ** 2 + 2 * tau ** 2 * Omega * omega ** 3 + (tau ** 2 * Omega ** 2 - 4 * Delta * tau + 1) * omega ** 2 + (-4 * tau * Delta * Omega + 2 * Omega) * omega + 4 * Delta ** 2 + Omega ** 2) / (omega ** 4 * tau ** 2 - 2 * tau ** 2 * Omega * omega ** 3 + (tau ** 2 * Omega ** 2 - 4 * Delta * tau + 1) * omega ** 2 + (4 * tau * Delta * Omega - 2 * Omega) * omega + 4 * Delta ** 2 + Omega ** 2) / np.pi)

# --> caso rilassamento exp + Markov    :       K(p) = somma casi precedenti        'Real'


def S_Dynamical_Form_Factor_2_nodelta(omega, Co, Omega, Gamma, Delta, tau):
    return (Gamma * omega ** 2 * tau ** 2 + Delta + Gamma) * Co * (omega ** 4 * tau ** 2 / 4 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / (omega ** 4 * tau ** 2 / 4 - tau ** 2 * Omega * omega ** 3 / 2 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + Omega * (Delta * tau - 0.1e1 / 0.2e1) * omega + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / np.pi / (omega ** 4 * tau ** 2 / 4 + tau ** 2 * Omega * omega ** 3 / 2 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + (-tau * Delta * Omega + Omega / 2) * omega + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / 2

def S_Dynamical_Form_Factor_2(omega, Co, Omega, Gamma, Delta, tau, delta_position, delta_width, delta_amplitude):

    return lorentian(omega, delta_position, delta_amplitude, delta_width)+(Gamma * omega ** 2 * tau ** 2 + Delta + Gamma) * Co * (omega ** 4 * tau ** 2 / 4 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / (omega ** 4 * tau ** 2 / 4 - tau ** 2 * Omega * omega ** 3 / 2 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + Omega * (Delta * tau - 0.1e1 / 0.2e1) * omega + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / np.pi / (omega ** 4 * tau ** 2 / 4 + tau ** 2 * Omega * omega ** 3 / 2 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + (-tau * Delta * Omega + Omega / 2) * omega + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / 2

def I_Dynamical_Current_2 (omega, Co, Omega, Gamma, Delta, tau):
    return (omega**2)*(Gamma * omega ** 2 * tau ** 2 + Delta + Gamma) * Co * (omega ** 4 * tau ** 2 / 4 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / (omega ** 4 * tau ** 2 / 4 - tau ** 2 * Omega * omega ** 3 / 2 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + Omega * (Delta * tau - 0.1e1 / 0.2e1) * omega + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / np.pi / (omega ** 4 * tau ** 2 / 4 + tau ** 2 * Omega * omega ** 3 / 2 + (0.1e1 / 0.4e1 + (Gamma ** 2 + Omega ** 2 / 4) * tau ** 2 - Delta * tau) * omega ** 2 + (-tau * Delta * Omega + Omega / 2) * omega + Delta ** 2 + 2 * Delta * Gamma + Gamma ** 2 + Omega ** 2 / 4) / 2


#       PARTE   1.2     Generatrici delle funzioni teoriche di cui sopra
#                       generano in un intervallo dato e con densità di punti data

def S_0_Generate(Co, Omega, Gamma, delta_width, delta_factor, x_min,  x_max, density, fig = False):
    
    x   =   np.linspace(x_min,x_max, density)
    y   =    S_Dynamical_Form_Factor_0(x, Co, Omega, Gamma,delta_width, delta_factor)

    
    if fig:
        
        plt.figure()
        plt.title('S_0 generated data')
        plt.plot(x,y,'.')


    return (x,y)

def I_0_Generate(Co, Omega, Gamma, x_min, x_max, density, fig = False):
    
    x   =   np.linspace(x_min,x_max, density)
    y   =   I_Dynamical_Current_0(x, Co, Omega, Gamma)

    
    if fig:
        
        plt.figure()
        plt.title('I_0 generated data')
        plt.plot(x,y,'.')


    return (x,y)

def S_1_Generate(Co, Omega, Delta, tau, x_min, x_max, density, fig = False):
    
    x   =   np.linspace(x_min,x_max,  density)
    y   =    S_Dynamical_Form_Factor_1(x, Co, Omega, Delta, tau)

    
    if fig:
        
        plt.figure()
        plt.title('S_1 generated data')
        plt.plot(x,y,'.')


    return (x,y)

def I_1_Generate(Co, Omega, Delta, tau, x_min, x_max, density, fig = False):
    
    x   =   np.linspace(x_min,x_max, density)
    y   =   I_Dynamical_Current_1(x, Co, Omega, Delta, tau)

    
    if fig:
        
        plt.figure()
        plt.title('I_1 generated data')
        plt.plot(x,y,'.')


    return (x,y)

def S_2_Generate(Co, Omega, Gamma, Delta, tau, delta_position, delta_width, delta_amplitude, x_min, x_max, density, fig = False):

    x   =   np.linspace(x_min,x_max,  density)
    y   =   lorentian(x, delta_position, delta_width, delta_amplitude) + S_Dynamical_Form_Factor_2_nodelta(x, Co, Omega, Gamma, Delta, tau)

    
    if fig:
        
        plt.figure()
        plt.title('S_2 generated data')
        plt.plot(x,y,'.')


    return (x,y)

def I_2_Generate(Co, Omega, Gamma, Delta, tau, x_min, x_max, density, fig = False):

    x   =   np.linspace(x_min,x_max,  density)
    y   =   I_Dynamical_Current_2(x, Co, Omega, Gamma, Delta, tau)

    
    if fig:
        
        plt.figure()
        plt.title('I_2 generated data')
        plt.plot(x,y,'.')


    return (x,y)


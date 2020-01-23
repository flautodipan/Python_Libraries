import numpy as np
from    Alessandria         import * 

class Lorentian_Fit :

    # classe che fa fit di una lorentzia
    # con la libreria symfit

    def __init__(self, X_DATA, Y_DATA, X_Delta = None, Y_Delta = None):

        self.x_data =   X_DATA
        self.y_data =   Y_DATA
        self.x_Delta  =   X_Delta
        self.y_Delta  =   Y_Delta
        
    def Get_Initial_Parameters(self):

        self.a      =   np.max(self.y_data)
        self.mu     =   self.x_data[np.argmax(self.y_data)]
        self.sigma  =   Estimate_FWHM(self.x_data,self.y_data)
    
    def Fit_Data(self, fig = False):

        A       =       sft.Parameter('A', value=(self.a))#, min=0., max=self.a)
        MU      =       sft.Parameter('MU', value = self.mu)#, min=(self.mu-self.sigma), max=(self.mu+self.sigma))
        SIGMA   =       sft.Parameter('SIGMA', value= self.sigma)#, min = 0., max=2*self.sigma)
        X,Y     =       sft.variables('X,Y')
        model   =       {Y  :   (A/np.pi)*(SIGMA/((X - MU)**2 + SIGMA**2))}

        self.fit        =   sft.Fit(model, X = self.x_data, Y = self.y_data, sigma_Y = self.y_Delta)
        self.fit_result =   self.fit.execute()

        # aggiorno valori dei parametri a quelli finali

        self.a          =   self.fit_result.value(A)
        self.Delta_a    =   self.fit_result.stdev(A)

        self.mu         =   self.fit_result.value(MU)
        self.Delta_mu   =   self.fit_result.stdev(MU)

        self.sigma      =   self.fit_result.value(SIGMA)
        self.Delta_sigma=   self.fit_result.stdev(A)

        
        #print(self.fit_result.__str__())

    def Plot_Results(self, plot_title):

        x       =   self.x_data
        y       =   self.y_data

        #infittisco plot del fit

        x_fit      =   np.linspace(np.min(x), np.max(x), 100)
        y_fit      =   lorentian(x_fit, self.a, self.mu, self.sigma)
        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)

        ax1.set_title(plot_title)
        ax1.set_xlabel('Frequencies (GHz)')
        ax1.set_ylabel('Counts')

        ax1.plot(x,y,'.')
        ax1.plot(x_fit, y_fit)

        plt.close(fig1)

from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import distance as dis
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from operator import methodcaller

from os.path import join

class Protein:

    def __init__(self, pdb_filename, chain_id = 'A', from_web = False, model = 1):
        
        """
        Classe dedicata alle proteine, che vengono lette dal file pdb indicato con pdb_filename
        
        Parametri

        pdb_filen   se il file è presente in cartella, indicare estensione
                    se from_web = True, senza estensione
        
        chain_id    di default considera pdb in cui non si ha che una proteina
                    in caso contrario indicare il chain_id      
        
        model       DEVE ESSERE SETTATO SU 'None' SE NEL PDB NON SI HANNO PIÙ MODELLI
                    altrimenti deve essere un'intero, anche se per adesso non sono stato
                    in grado di farlo per n diverso da 1
                    noi prendiamo il primo..poi si vedrà

        """

        #1) carico pdb

        if from_web: 

            self.pdb = PandasPdb().fetch_pdb(pdb_filename)            
        else:
            self.pdb = PandasPdb().read_pdb(pdb_filename)

        #2A) seleziono modello

        if model:


            n= model-1 #aggiusto per combaciare indici

            if model == 1:

                models = self.pdb.df['OTHERS']

                start_index = np.array(models[self.pdb.df['OTHERS']['record_name']== 'MODEL']['line_idx'])[n]
                stop_index  = np.array(models[self.pdb.df['OTHERS']['record_name']== 'ENDMDL']['line_idx'])[n]

            #estraggo indici partenza/arresto lettura dal pdb
            # lo faccio facendo la query delle righe con MODEL e prendendo gli 
            #indici di linea --> a meno di un addendo costante è la numerosità da prendere
            # per i valori identificati dalle keys ['ATOM']

            elif (model >1):
                models = self.pdb.df['OTHERS'][self.pdb.df['OTHERS']['record_name']== 'MODEL']
                start_index = np.array(models['line_idx'])[n]
                stop_index = np.array(models['line_idx'])[n+1]
                
            # ho contato a mano che sbaglia di 4 righe considerando MODEL,ENDMDL e TER
            # che ignora perchè io gli chiedo ATOM come key e costruisco start e stop
            # in base a MODEL
            
            Delta      = stop_index-start_index-4        
        
            #3A) seleziono atomi e chain_id

            self.atoms  =   self.pdb.df['ATOM'][n:Delta]
            self.atoms  =   self.atoms[self.pdb.df['ATOM']['chain_id']==chain_id]
            self.initial=   self.atoms.residue_number[0]

        #self.terminus = self.pdb.df['OTHERS']
        
        #2B) prendo atomi e seleziono chain_id
        else:
            
            self.atoms  =   self.pdb.df['ATOM']
            self.atoms  =   self.atoms[self.atoms.chain_id == chain_id] 
            self.initial=   self.atoms.residue_number.values[0]
            #self.terminus =  self.pdb.df['OTHERS']
         
        
    def Get_Atom_Coord(self, atom_name): 

        # assegno dinamicamente il nome dell'attributo in base al nome dell'atomo

        setattr(self, atom_name, self.atoms[self.pdb.df['ATOM']['atom_name']==atom_name]) 
        setattr(self, atom_name+'_xcoord', np.array(getattr(self, atom_name)['x_coord']))
        setattr(self, atom_name+'_ycoord', np.array(getattr(self, atom_name)['y_coord']) )
        setattr(self, atom_name+'_zcoord', np.array(getattr(self, atom_name)['z_coord']) )
        setattr(self, atom_name+'_Coord', np.matrix([(getattr(self,atom_name+'_xcoord')), (getattr(self,atom_name+'_ycoord')), (getattr(self,atom_name+'_zcoord'))]).T)

    def Get_lenght(self):

        if hasattr(self, 'CA'):
            setattr(self, 'lenght', len(getattr(self, 'CA_Coord')))
        else: raise ValueError('Non so come stimare lunghezza Prot\n\n')
        
    def Get_Protein_Sequence(self):

        setattr(self, 'sequence', self.pdb.amino3to1()[self.pdb.amino3to1().chain_id == 'A'].residue_name.values)

class RNA:

    def __init__(self, pdb_filename, chain_id = 'A', from_web = False, model = 1, initial = None):
        
        """
        Classe dedicata all'RNA, che viene letto dal file pdb indicato con pdb_filename
        
        Parametri

        pdb_filen   se il file è presente in cartella, indicare estensione
                    se from_web = True, senza estensione
        
        chain_id    di default considera pdb in cui non si ha che un RNA
                    in caso contrario indicare il chain_id      
        
        model       DEVE ESSERE SETTATO SU 'None' SE NEL PDB NON SI HANNO PIÙ MODELLI
                    altrimenti deve essere un'intero, anche se per adesso non sono stato
                    in grado di farlo per n diverso da 1
                    noi prendiamo il primo..poi si vedrà

        """

        #1) carico pdb

        if from_web: 

            self.pdb = PandasPdb().fetch_pdb(pdb_filename)            
        else:
            self.pdb = PandasPdb().read_pdb(pdb_filename)

        #2A) seleziono modello

        if model:

            n= model-1
            
            if model == 1:

                models = self.pdb.df['OTHERS']

                start_index = np.array(models[self.pdb.df['OTHERS']['record_name']== 'MODEL']['line_idx'])[n]
                stop_index  = np.array(models[self.pdb.df['OTHERS']['record_name']== 'ENDMDL']['line_idx'])[n]


            elif model > 1:
                models = self.pdb.df['OTHERS'][self.pdb.df['OTHERS']['record_name']== 'MODEL']
                start_index = np.array(models['line_idx'])[n]
                stop_index = np.array(models['line_idx'])[n+1]
                
            # ho contato a mano che sbaglia di 4 righe considerando MODEL,ENDMDL e TER
            # che ignora perchè io gli chiedo ATOM come key e costruisco start e stop
            # in base a MODEL
            
            Delta      = stop_index-start_index-4        
        
            #3A) seleziono atomi e chain_id

            self.atoms  =   self.pdb.df['ATOM'][n:Delta]
            self.atoms  =   self.atoms[self.pdb.df['ATOM']['chain_id']==chain_id]

            #self.terminus = self.pdb.df['OTHERS']
        
        #2B) prendo atomi e seleziono chain_id

        else:
            
            self.atoms  =   self.pdb.df['ATOM']
            self.atoms  =   self.atoms[self.atoms.chain_id ==chain_id] 
            #self.terminus =  self.pdb.df['OTHERS']

        if initial:

            self.initial = initial
        else:
            self.initial = self.atoms.residue_number.values[0]

        
    def Get_Atom_Coord(self, atom_name): 

        # assegno dinamicamente il nome dell'attributo in base al nome dell'atomo
        real_atom_name = atom_name
        if atom_name == "O5'":
            atom_name = 'O5'
        setattr(self, atom_name, self.atoms[self.pdb.df['ATOM']['atom_name']==real_atom_name]) 

        setattr(self, atom_name+'_xcoord', np.array(getattr(self, atom_name)['x_coord']))
        setattr(self, atom_name+'_ycoord', np.array(getattr(self, atom_name)['y_coord']) )
        setattr(self, atom_name+'_zcoord', np.array(getattr(self, atom_name)['z_coord']) )
        setattr(self, atom_name+'_Coord', np.matrix([(getattr(self,atom_name+'_xcoord')), (getattr(self,atom_name+'_ycoord')), (getattr(self,atom_name+'_zcoord'))]).T)

    def Get_lenght(self):

        if hasattr(self, 'O5'):
            setattr(self, 'lenght', len(getattr(self, 'O5_Coord')))
        elif (not hasattr(self, 'O5')) & hasattr(self, 'P'):
            setattr(self, 'lenght', len(getattr(self, 'P_Coord')+1))
        else: raise ValueError('Non so come stimare lunghezza RNA\n\n')

    def Get_RNA_Sequence(self, ):

        if hasattr(self, 'O5'):
            setattr(self, 'sequence', self.atoms.residue_name[self.atoms.atom_name == "O5'"].values)
        else: raise ValueError('Non ho atomi O5 per capire, sequenza puoi prenderla solo dopo')

    
class Trajectory():

    """

    Classe per analisi degli autovalori e autovettori della traiettoria effettuata in simulazione
    MD con il Software GROMACS
    --> indicare il path dove si trovano i file eigenval.xvg e eigenvec.trr ù

    """
    def __init__(self, bio_name):

        self.bio_name   =   bio_name

    def __str__(self):

        return '{}'.format(self.bio_name)

    def Get_EigenValues(self, xvg_filename = 'eigenval.xvg', path = './', **kwargs):

        _ , self.EigenValues    = Parse_xvg(path+xvg_filename)

    def Print_EigenValues(self, path = './', **kwargs):

        
        plt.figure()
        plt.plot(self.EigenValues[2:], '.', c = kwargs['color'])
        plt.plot(self.EigenValues[0:2], '*', c = kwargs['contrastcolor'], label = 'First two eigenvalues\nExplained Variance Ratio = {:3.2f}'.format(self.perceigen) )
        plt.title('Eigenvalues of MD trajetory for '+self.__str__())
        plt.xlabel('Eig index (ordered by explained variance)') 
        plt.legend()          
        plt.savefig(path+kwargs['fig']+'.pdf', format = 'pdf')
        plt.show()

    def Set_Time_Info(self, n_frames, time_range, timestep = 100):

        self.n_frames       =   n_frames
        self.initial_time   =  time_range[0]
        self.final_time     =  time_range[1]
        self.timestep       =  timestep
        self.time_range     = time_range

        time_step = (self.final_time-self.initial_time)/(self.n_frames-1)
        if time_step != self.timestep:
            raise ValueError("Errore in acquisizione: timestep inserito = {}\ntimestep calcolato = {}".format(self.timestep, time_step))

        print('Acquisito correttamente informazioni temporali:\n')
        print('Tempo iniziale =\t{} ps \nTempo finale =\t{} ps \ntime step =\t{} ps\nN frame =\t{}\n'.format(self.initial_time, self.final_time, self.timestep, self.n_frames))


    def Get_2D_Traj(self, xvg_filename = '2dproj.xvg', path = './',  **kwargs):


        self.x_proj, self.y_proj    =   Parse_xvg(xvg_filename, path = path)
        if (self.n_frames != self.x_proj.size):

           raise ValueError("Dimensione spazio essenziale differente da #frame \n #frames = %d \t len(essential) = %d \n Sicuro hai scelto il file giusto?\n" %(self.n_frames, self.RMSD.size))


        print   ('Acquisito proiezioni su PC1 e PC2 della traiettoria.\n')
        print   ('Numero di punti nello spazio essenziale : {}\n'.format(self.n_frames))

        if ('fig' in kwargs) :

            fig, ax = plt.subplots()

            c = range(self.n_frames)
            scatt = ax.scatter(self.x_proj,self.y_proj, c=c, cmap= 'viridis', s = 10)
            cbar = plt.colorbar(scatt, ax = ax, label = 'Time (ns)')
            cbar.ax.set_label('Time (ns)',)
            ax.set_xlabel(' PC 1 ')
            ax.set_ylabel(' PC 2 ')
            ax.set_title('PCA of MD traj for '+self.__str__())
            plt.tight_layout()
            fig.savefig(path+kwargs['fig']+'.pdf', format = 'pdf', bbox_inches = 'tight')  
            plt.show()

    
    def Get_RMSD(self, xvg_filename, equilibrium = False, mode = 'tot', scale = 'ns', path='./', skip = False, **kwargs):

        if equilibrium:
            
            if mode == 'RNA':
                RMSD = 'RMSD_eq_RNA'
            elif mode == 'BS':
                RMSD = 'RMSD_eq_BS'

            time_range = 'time_range_eq'
            n_frames = 'n_frames_eq'

        else:
            n_frames = 'n_frames'
            time_range = 'time_range'
            if mode == 'RNA':
                RMSD = 'RMSD_RNA'
            elif mode == 'BS':
                RMSD = 'RMSD_BS'
            else:
                RMSD = 'RMSD'

        setattr(self, RMSD, Parse_xvg_skip(xvg_filename, skip = skip, path = path)[1])
        print('lunghezza iniz = {}'.format(len(getattr(self, RMSD))))
        if equilibrium:
            setattr(self, RMSD, getattr(self, RMSD)[self.idx_eq_left:self.idx_eq_right+1])
        print('lunghezza fin = {}'.format(len(getattr(self, RMSD))))

        if getattr(self, n_frames) != getattr(self, RMSD).size:
            raise ValueError("Incompatibilità tra numero frame = {} e len({}) = {}".format(getattr(self, n_frames), RMSD, len(getattr(self, RMSD))))
  
        
        if ('fig' in kwargs):

            plt.figure()
            if scale == 'ns':
                plt.plot(np.arange(getattr(self, time_range)[0], getattr(self, time_range)[1]+self.timestep, self.timestep)/1000, getattr(self, RMSD), ls = '--', lw = 0.3, color=kwargs['color'], alpha = 0.8)
                plt.xlabel('Time (ns)')
            else:
                plt.plot(np.arange(getattr(self, time_range)[0], getattr(self, time_range)[1]+self.timestep, self.timestep), getattr(self, RMSD), ls = '--', lw = 0.3, color=kwargs['color'], alpha = 0.8)
                plt.xlabel('Time (ps)')
            if 'ylim' in kwargs:
                plt.ylim(kwargs['ylim'][0], kwargs['ylim'][1])
            
            plt.ylabel('{} (nm)'.format(RMSD))
            plt.title('{} for {}'.format(RMSD, self.__str__()))
            plt.savefig(join(path, kwargs['fig']+'.pdf'), format = 'pdf')
            

        if ('histo' in kwargs):

            plt.figure()
            _ =   plt.hist(getattr(self, RMSD), bins = kwargs['bins'], histtype='bar', rwidth=0.8, color=kwargs['color'], alpha = 0.9)
            plt.xlabel('{} (nm)'.format(RMSD))
            plt.title("{} distribution for {}".format(RMSD, self.__str__()))
            plt.savefig(join(path, kwargs['histo']+'.pdf'), format = 'pdf')
            plt.show()

    def Define_Equilibrium_by_RMSD(self, time_range_eq, path = './', scale = 'ns', **kwargs):

        self.time_range_eq = time_range_eq
        self.idx_eq_left = int((time_range_eq[0] - self.initial_time)/self.timestep)
        self.idx_eq_right = int((time_range_eq[1] - self.initial_time)/self.timestep)

        if (time_range_eq[0] - self.initial_time)%self.timestep != 0:
            raise ValueError("tempo inserito  {} come iniziale non è multiplo del timetep {}".format(time_range_eq[0], self.timestep))
        
        self.RMSD_eq = self.RMSD[self.idx_eq_left:self.idx_eq_right+1]
        self.n_frames_eq = len(self.RMSD_eq)

        if 'fig' in kwargs:

            f = plt.figure()
            ax = plt.subplot(111)

            if scale == 'ns':
                ax.plot(np.arange(self.time_range[0], self.time_range[1]+self.timestep, self.timestep)/1000, self.RMSD, ls = '--', lw = 0.3, color=kwargs['color'], alpha = 0.8, label = 'RMSD')
                plt.xlabel('Time (ns)')
            else:
                ax.plot(np.arange(self.time_range[0], self.time_range[1]+self.timestep, self.timestep), self.RMSD, ls = '--', lw = 0.3, color=kwargs['color'], alpha = 0.8, label = 'RMSD')
                plt.xlabel('Time (ps)')   

            if 'ylim' in kwargs:
                plt.ylim(kwargs['ylim'][0], kwargs['ylim'][1])

            plt.ylabel('RMSD (nm)')
            plt.fill_betweenx(np.arange(ax.get_ylim()[0], 10, 0.1), time_range_eq[0]/1000, time_range_eq[1]/1000, alpha = kwargs['alpha'], color = kwargs['darkcolor'], label = 'Equilibrium time range')
            plt.title('Equilibrium RMSD for {}'.format(self.__str__()))
            
            plt.legend()
            plt.savefig(join(path, kwargs['fig']+'.pdf'), format = 'pdf') 

        print("Selezionata zona di equilibrio da {} ps  a {} ps".format(time_range_eq[0], time_range_eq[1]))
        print('Pari a un numero di frame = {}'.format(self.n_frames_eq))
        print("l'indice di riferimento da cui partire rispoetto all'array RMSD completo è  {}".format(self.idx_eq_left))
        print("l'indice di riferimento a cui finire rispoetto all'array RMSD completo è  {}".format(self.idx_eq_right+1))


    def Acquire_Atoms_List(self, xvg_filename, what, path = './' , **kwargs):

            if what == 'RNA':
                attr = 'RNA_index'
            elif what == 'BS':
                attr = 'BS_index'

            setattr(self, attr, np.array(Parse_xvg_skip(xvg_filename, path, **kwargs)[0], dtype = int))
    
    """
    def Get_RMSF_res(self, xvg_filename, now_name, text_size = 7, path = './'):
    
        self.res, self.RMSF_res = Parse_xvg_skip(xvg_filename, path)

        #scambio RNA e PROT per coerenza con altre figure 
        idx_RNA_start = np.where( self.res == 1.)[0][0]
        self.res[idx_RNA_start:] += self.res[idx_RNA_start-1]

        text = []   
        size = text_size
        start = 0

        for ii in range(int(BS.size/size)):
            text.append([str(v) for v in BS[start:start+size]])
            start+=size
        if BS.size%size != 0:
            final = [str(v) for v in BS[start:]]
            if len(final) != size:
                while True:
                    final.append(' ')
                    if len(final) == size:
                        break
                text.append(final)

        f, ax = plt.subplots(1,1)
        ax.stem(self.res[:idx_RNA_start], self.RMSF_res[:idx_RNA_start],  markerfmt='orangered', basefmt='orangered', linefmt='orange', label = 'Protein')
        ax.stem(self.res[idx_RNA_start:], self.RMSF_res[idx_RNA_start:], markerfmt=color, basefmt=color, linefmt='yellowgreen', label = 'RNA')
        ax.legend(title = 'RNA starts at self.res {}'.format(int(self.res[idx_RNA_start])))
        ax.set_title('RMSF for {} residues'.format(now_name), pad = 1.3)
        ax.table(text)
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.set_xlabel('Residue number')
        ax.set_ylabel('RMSF (nm)')

        f.savefig(path+'RMSF_res_'+now_name+'.pdf', format = 'pdf')
    """

    def Get_RMSF(self, xvg_filename, equilibrium = False, path='./', skip_lines = 17, **kwargs):

        n_tot, self.RMSF = Parse_xvg_skip(xvg_filename, path, skip_lines = skip_lines)

        self.remaining_index = list(n_tot)

        for ii in n_tot:
            if (ii in self.RNA_index) | (ii in self.BS_index):
                self.remaining_index.remove(ii)
        
        self.remaining_index = np.array(self.remaining_index, dtype = int)

        if (len(self.remaining_index)+ self.RNA_index.size + self.BS_index.size) != n_tot.size:
            raise ValueError("Dimensioni degli indici non tornano\nBS = {}\tRNA={}\tremaining={}\nTot={}".format(self.BS_index.size, self.RNA_index.size, len(self.remaining_index), n_tot.size))

        if 'fig' in kwargs:

            plt.figure()
            plt.title('Equilibrium RMSF for {}'.format(self.__str__()))
            plt.plot(self.remaining_index, self.RMSF[self.remaining_index], color = kwargs['thirdcolor'], alpha=0.65)       
            plt.plot(self.RNA_index, self.RMSF[self.RNA_index-1], '--', color = kwargs['color'], label = 'RNA')
            plt.plot(self.BS_index, self.RMSF[self.BS_index-1], '.',color = kwargs['darkcolor'], label = 'BS')
            if 'ylim' in kwargs:
                plt.ylim(kwargs['ylim'][0], kwargs['ylim'][1])
            plt.legend()
            plt.xlabel('Number of atom')
            plt.ylabel('RMSF (nm)')
            plt.savefig(join(path, kwargs['fig']+'.pdf'), format = 'pdf')

    def Get_Gyradium(self, xvg_filename, path  = './', skip_lines = 28, **kwargs):

        self.Gyradium = Parse_xvg_skip(xvg_filename, path, skip_lines=skip_lines)[1]

        if 'fig' in kwargs:

            plt.figure()
            plt.title('Radium of Gyration BS+RNA for {}'.format(self.__str__()))
            plt.plot(np.arange(self.time_range[0], self.time_range[1]+self.timestep, self.timestep)/1000, self.Gyradium, '-', c = kwargs['color'], label = 'Gyration radium')
            plt.fill_betweenx(np.arange(0, 10, 0.1), self.time_range_eq[0]/1000, self.time_range_eq[1]/1000, alpha = kwargs['alpha'], color = kwargs['darkcolor'], label = 'Equilibrium time range')
            plt.legend()
            plt.xlabel('Time (ns)')
            plt.ylabel('Gyration Radium (nm)')
            if 'ylim' in kwargs:
                plt.ylim(kwargs['ylim'][0], kwargs['ylim'][1])
            plt.savefig(join(path, kwargs['fig']+'.pdf'), format = 'pdf')


    def Get_Terminals_Dist(self, xvg_filename, path='./', skip = False, **kwargs):

        """
        Prendo distanza tra terminali da file xvg associata alla traiettoria

        - skip: fattore di cui si vuole ridurre numerosità

        - fig solita kwarg

        - histo come fig ma stampa istogramma distribuzione distanza, da accompagnare
          con n_bins come kwarg

        """

        _ ,  self.ter_dist =  Parse_xvg_skip(xvg_filename, skip_lines = kwargs['skip_lines'], skip = skip, path = path)

        if (self.ter_dist.size != self.x_proj.size):

           raise ValueError("Dimensione ter_dist differente da #frame \n #frames = %d \t len(ter_dist) = %d \n Sicuro hai scelto il file giusto?\n" %(self.n_frames, self.ter_dist.size))
        
        if ('fig' in kwargs):

            plt.figure()
            plt.plot(np.arange(0,self.x_proj.size*self.timestep, self.timestep)/1000, self.ter_dist,  ls = '--', lw = 0.3, color=kwargs['color'], alpha = 0.8)
            plt.xlabel('Time (ns)')
            plt.ylabel('terminals distance(nm)')
            plt.title('Terminals distance for '+self.__str__())
            plt.savefig(path+kwargs['fig']+'.pdf', format = 'pdf')
        

        if ('histo' in kwargs):

            plt.figure()
            _ =   plt.hist(self.ter_dist, bins = kwargs['bins'], histtype='bar', rwidth=0.8, color=kwargs['color'], alpha = 0.9, density = True)
            plt.xlabel('ter_dist (nm)')
            plt.title("Distribution of terminals distance "+self.__str__())
            plt.grid(axis = 'y', alpha = 0.45)
            plt.savefig(path+kwargs['histo']+'.pdf', format = 'pdf')
            plt.show()

    def Analyze_Traj_by_Descriptor (self, descriptor, N_gauss, p0, bins = 50, **image_kwargs):

        self.N_gauss = N_gauss
        """
        Funzione che
        -   fitta la distribuzione RMSD (con P0 iniziali) con N_gauss gaussiane (2 o 3)
            --> permette di ripetere 
        
        OSS:
        p0 deve essere di dimensione 3*N_gauss, nell'ordine (mu, sigma, A)
        kwargs importanti fig, fig_fit, path

        """
        descrittore  = descriptor #nome per fig

        if (descriptor == 'RMSD'):
            descriptor = self.RMSD
        elif (descriptor ==  'terminals distance'):
            descriptor = self.ter_dist
        else:
            raise ValueError (" Scegliere descrittore per analisi traiettoria : \n 'RMSD' o  'terminal distance' ")
        
        if (N_gauss == 3 ):
            fit_mode = trimodal
            i = 6
            idx = 'mu3'
            index = ('mu1', 'sigma1', 'A1', 'mu2', 'sigma2', 'A2', 'mu3', 'sigma3', 'A3')
        
        elif (N_gauss == 2):
            i = 3
            fit_mode = bimodal
            idx = 'mu2'
            index = ('mu1', 'sigma1', 'A1', 'mu2', 'sigma2', 'A2')
        
        #fit e plot
        y, x,  = np.histogram(self.ter_dist, bins = bins)
        x = (x[1:] + x[:-1])/2 
        self.gauss_params, cov = curve_fit(fit_mode,x,y, p0=p0)
        
        print("Parametri del fit \n")
        self.df  =   pd.DataFrame({'Values':self.gauss_params, 'StdErr':np.sqrt(np.diag(cov))}, index = index) 
        print(self.df)
        
        # divido perchè distribuzione è indipendente da ordine

        fold = self.ter_dist[self.ter_dist < self.gauss_params[i]]
        unfold = self.ter_dist[self.ter_dist >= self.gauss_params[i]]

        # plotto primo histo
        plt.figure()
        plt.hist([fold, unfold], histtype = 'bar', stacked = True, color=[image_kwargs['color_fold'], image_kwargs['color_unfold']], bins = bins, rwidth=0.8, label = ['Folded states', 'unfolded states'], alpha = 0.9)
    
        #fit e plot
        x = np.linspace(descriptor.min(), descriptor.max(), 1000)
        if N_gauss == 3:
            plt.plot(x, trimodal(x, *self.df['Values'].values), color = image_kwargs['fit_color'], linestyle = 'dashed', lw = 2, label = 'Multimodal Fit')    
        elif N_gauss == 2:
            plt.plot(x, bimodal(x, *self.df['Values'].values), color = image_kwargs['fit_color'], linestyle = 'dashed', lw = 2, label = 'Multimodal Fit')    

        plt.title(' Distribution of '+descrittore+' for '+self.__str__())
        plt.xlabel(descrittore +' (nm)')
        plt.legend(title= 'mu = {:3.2f}\u00B1{:3.2f}'.format(self.gauss_params[i], self.df['StdErr'][idx]))
        plt.savefig(image_kwargs['path']+image_kwargs['fig_fit']+'.pdf', format = 'pdf')
        plt.show()


    def Divide_Traj_by_Descriptor (self, descriptor, sigma_factor, **image_kwargs):

        """
        -   suddivide lo spazio essenziale in base al treshold = mu - sigma*factor della 
            gaussiana con RMSD più alto (corrispondente a stato UNFOLD)

        
        """
        descrittore  = descriptor #nome per fig

        if (descriptor == 'RMSD'):
            descriptor = self.RMSD
        elif (descriptor ==  'ter_dist'):
            descriptor = self.ter_dist
        else:
            raise ValueError (" Scegliere descrittore per analisi traiettoria : \n 'RMSD' o  'ter_dist' ")

        if (self.N_gauss == 3):

            treshold = self.gauss_params[6] - self.gauss_params[7]*sigma_factor

        elif (self.N_gauss == 2):

            treshold = self.gauss_params[3] - self.gauss_params[4]*sigma_factor

        print("Treshold data dalla gaussiana con "+descrittore + " più alto:\nmu - simga*%3.2f = %3.2f\n\n" %(sigma_factor, treshold))  
        
        self.x_fold     =   self.x_proj[descriptor <= treshold]
        self.y_fold     =   self.y_proj[descriptor <= treshold]

        self.x_unfold   =   self.x_proj[descriptor > treshold]
        self.y_unfold   =   self.y_proj[descriptor > treshold]
        
        self.unfold_percentage = 100*self.x_unfold.size/(self.x_fold.size+self.x_unfold.size)

        print("Percentuale di popolazione stimata unfold = {} ".format(self.unfold_percentage))

        if ('fig' in image_kwargs) :

            fig     =   plt.figure()
            mask    =   descriptor > treshold
            plt.scatter(self.x_fold, self.y_fold, s = 7, color = image_kwargs['color'])
            plt.scatter(self.x_unfold, self.y_unfold, s = 7, color = image_kwargs['darkcolor'])
            plt.xlabel(' PC 1 ')
            plt.ylabel(' PC 2 ')
            title = 'Folded/Unfolded states division for {}'.format( self.__str__())
            plt.title(title)
            plt.legend(['fold', 'unfold = {:3.2f} of population'.format(self.unfold_percentage)]) #, title = 'percentage of unfolded population = {:3.2f}'.format(self.unfold_percentage))

            fig.savefig(image_kwargs['path']+image_kwargs['fig']+'.pdf', format = 'pdf') 
            plt.show()
      
    def Which_Part_of_Traj(self, points):

        """

        Funzione che stampa i valori di frame, time, descrittore associati i punti nello spazio essenziale points
        Pensata per i risultati della clustering analyis, quindi su valori numerici esattamente uguali (no bisogno find_nearest)

        Funzionamento : compara la x di ogni riga di points, tanto l'y è allo stesso indice

        """
        for jj in range(len(points)):
            
            idx =   int(np.nonzero(self.x_proj == points[jj][0])[0])

            if (idx != int(np.nonzero(self.y_proj == points[jj][1])[0])):
                raise ValueError("Qualcosa di strano: indice x e y non corrispondono.\n x: %d \t y: %d\n\n"%(idx, np.nonzero(self.y_proj == points[jj][1])))
            
            x   =   self.x_proj[idx]
            y   =   self.y_proj[idx]

            print("Il punto nello spazio essenziale corrispondente al clusterchief %d è (%3.2f, %3.2f)\n"%(jj+1, x, y))
            print("\nCorrisponde al frame {}\nOssia al tempo {} ps\ned ha un valore di RMSD pari a {}\nun valore di distanza dei terminali di {}\n".format(idx+1, (idx+1)*self.timestep, self.RMSD[idx], self.ter_dist[idx]))   

    def Analyze_Variance_Explained_Ratio(self, percent_sigma):

        """
        
        Funzione che salva e stampa gli autovalori che spiegano
        la varianza sigma richiesta da utente 
    
        
        """

        self.Norma           = np.sum(self.EigenValues)
        self.Sigma_ERatio    = self.EigenValues/self.Norma

        n               = 0
        test            = 0.
        Sigma_Eigens    = []

        while True:

            Sigma_Eigens.append(self.EigenValues[n])
            test        += self.Sigma_ERatio[n]         
            n           += 1
            if test     > percent_sigma: break

        df              = pd.DataFrame(Sigma_Eigens)
        print('# Autovalori che spiegano la varianza richiesta = ', n, '\n\n')
        df
        print(Sigma_Eigens)
    
    def Analyze_Eig_Variance_Ratio (self, n_eig = 2):

        """
        Funzione che ritorna la percentuale di varianza del moto spiegata
        da n_eig autovalori della matrice delle covarianze
        default = 2 perchè lavoriamo in spazio essenziale di 2d

        """
        
        self.Norma          = np.sum(self.EigenValues)
        self.Sigma_ERatio   = self.EigenValues/self.Norma

        mask                = np.zeros(len(self.Sigma_ERatio), dtype=bool)
        mask[:n_eig]        = True
        self.perceigen      = np.sum(self.Sigma_ERatio, where=mask)       


        print('# Percentuale di varianza spiegata da n_eig = %d autovalori \n %f\n\n' %(n_eig, self.perceigen*100))
        print('Autovalori:\n', self.EigenValues[:n_eig])

class Cluster_2DAnalysis():

    def __init__(self, clustype, name):
        
        self.clustype   =   clustype
        self.name       =   name        
          
    def __str__(self):

        return "{}".format(self.name)

    def Get_by_Traj(self, traj):

        self.xy         =   np.array([traj.x_proj,traj.y_proj]).T
        self.timestep   =   traj.timestep

    def Get_by_Passing(self, xy, timestep):

            self.xy         =   xy
            self.timestep   =   timestep
        
    def Elbow_KMeans(self, kmax):

        """ Funzione che implementa il metodo "Elbow" sull'algoritmo k_means

            calcola e mostra la Within-Cluster Sum of Squared errors 
            per     differenti valori di k, al fine di scegliere il migliore.

        -----

            Parameters

            points sono i dati su cui clusterizzare
            variando k tra 1 e kmax

        """
    
        self.Elbow = []

        for k in range(1, kmax+1):

            kmeans = KMeans(n_clusters = k).fit(self.xy)
            centroids = kmeans.cluster_centers_
            pred_clusters = kmeans.predict(self.xy)
            curr_sse = 0
    
            # calculate square of Euclidean distance of each point from its cluster center and add to current WSS
            for i in range(len(self.xy)):

                curr_center = centroids[pred_clusters[i]]
                curr_sse += (self.xy[i, 0] - curr_center[0]) ** 2 + (self.xy[i, 1] - curr_center[1]) ** 2
      
            self.Elbow.append(curr_sse)
        

        plt.figure()
        plt.plot(range(kmax), self.Elbow)
        plt.title('Elbow check for K-means')
        plt.xlabel('k value')
        plt.savefig('Elbow_test')

    def Silhouette_KMeans(self, kmax, **kwargs):

        self.Silh = []

        # dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
        for k in range(2, kmax+1):
            kmeans = KMeans(n_clusters = k).fit(self.xy)
            labels = kmeans.labels_
            self.Silh.append(silhouette_score(self.xy, labels, metric = 'euclidean'))

        self.silhouette_max =  np.array(self.Silh).argmax() + 2 
        
        if ('fig' in kwargs):
                
            plt.figure()
            plt.plot(np.arange(0,kmax-1,1)+2, self.Silh)
            plt.title('Silohuette check for k-means algorithm for '+self.__str__())
            plt.xlabel('k value')
            plt.savefig(kwargs['path']+kwargs['fig']+'.png')
            plt.show()

    def Clusterize(self, by_silhouette = True, **kwargs):
        
        """
        - Effettua la clusterizzazione, di default prendendo il massimo dall'analisi di silhouette
        ma può essere usato specificando il keyword arg n_clust e mettendo False per by_silhouette
        - Calcola centroidi e trova elementi in xy più vicini ai centroidi,
        salvandoli in self.clusterchiefs con i corrispettivi indici in ""_idx

        - verbose = True per stampe


        """
        if (self.clustype == 'KMeans'):

            if not by_silhouette:

                self.n_clusters =   kwargs['n_clust']

            else:

                self.n_clusters = self.silhouette_max

            self.Cluster = KMeans(n_clusters=self.n_clusters).fit(self.xy)

            #definisco clusterchiefs : elementi (x,y) in xy più vicini
            # ai centroidi di ogni cluster
            self.clusterchiefs      =   ()
            self.clusterchiefs_idx  =   ()

            for jj in range(self.n_clusters):

                self.clusterchiefs      =  self.clusterchiefs + (Find_Nearest_2D(self.xy, (self.Cluster.cluster_centers_[jj,:]))[0],)
                self.clusterchiefs_idx  =  self.clusterchiefs_idx + (Find_Nearest_2D(self.xy, (self.Cluster.cluster_centers_[jj,:]))[1],)
                          
                if kwargs['verbose']:

                    print("Il punto più vicino al centroide del cluster " + str(jj+1) + " è \n" +str(self.clusterchiefs[jj])+'\n')

                    
            if ('fig' in kwargs):

                fig = plt.figure()

                plt.scatter(self.xy[:,0], self.xy[:,1], c=self.Cluster.labels_, cmap= 'viridis', s = 10)
                plt.scatter(self.Cluster.cluster_centers_[:,0],self.Cluster.cluster_centers_[:,1], marker='*', edgecolors='black', label = 'Cluster centroids')
                plt.xlabel(' PC 1 ')
                plt.ylabel(' PC 2 ')
                plt.title(str(self.n_clusters)+ ' Cluster representation for '+self.__str__())
                plt.legend()
                fig.savefig(kwargs['path']+kwargs['fig']+'.png')
                plt.show()

        else:

            raise ValueError ("Non è stato selezionato un algoritmo di clustering corretto\
                Scegliere:\n 'KMeans' per algoritmo k-means\n")    
        
    def Cluster_RMSD_Division (self, RMSD, histo = False, **kwargs):
        
        """
        Funzione che individua il cluster il cui centroide ha RMSD più alto, 
        assumendolo rappresentativo dello stato unfold

        histo se True va accompagnato dai kwargs bins, etc della funzione plt.histo()
        essendo una fig, va passato anche il path di salvataggio

        verbose stampa un po' di cose


        """

                
        if (RMSD.size != self.xy.shape[0]):

           raise ValueError("Dimensione RMSD differente da #frame \n #frames = %d \t len(RMSD) = %d \n Sicuro hai scelto il file giusto?\n" %(self.xy.shape[0], RMSD.size))
        
        # get cluster idx with the higher RMSD

        self.unfoldest_idx  =   RMSD[list(self.clusterchiefs_idx)].argmax()
        
        if ('verbose' in kwargs):

            print("Cluster con RMSD maggiore è il %d \n con RMSD = %f\n" %(self.unfoldest_idx +1, RMSD[list(self.clusterchiefs_idx)].max()))

def BioStructure (pdb_filename, n_structures, structures_names, structures_type, structures_chain_ids, structures_models):

    """
        if len(structures_type) is not n_structures:
            raise ValueError("Size of structure type must be n_structures\n")
        if (structures_type.type is not str):
            raise ValueError("structures type must be a tuple of strings \n Strings accepted are 'Protein' or 'RNA' for now\n\n")
    """

    """ IDEA 

    fare una clasese che si prende in input oggetti di classe e come attributi ha automaticamente
    - matrice coordinate
    - matrice contatti
    - matrice legami (ossia matrice contatti tagliata nella parte asimmetrica)
    
    
    
    """
    for (st_type, st_name, st_chain_id) in zip(structures_type, structures_names, structures_chain_ids):
        command = st_name+' = '+st_type+'('+pdb_filename+', chain_id ='+st_chain_id+')'
        print(command)
        exec(command)

def Dist_Matrix (Coord):  
# prende matrice coordinate x,y,z di N elementi, ritorna matrice delle distanze relative
# capire se come in C è un problema ritornare una matrice, 
# è meglio fare con i puntatori

    if (np.shape(Coord)[1] != 3):
        raise ValueError ("Dist_Matrix is defined with 3 columns matrix\n")

    N = np.shape(Coord)[0]
    Distances = np.zeros((N, N))

    i = 0
    k = 0
    for i in range (N):
        for k in range (N):
            Distances[i,k] = dis.euclidean(Coord[i,], Coord[k,]) 
    return Distances


def R_Coord (Coord):
# prende matrice coordinate x,y,z di N elementi, ritorna vettore coordinate radiali
    if (np.shape(Coord)[1] != 3):
        raise ValueError ("Dist_Matrix is defined with 3 columns matrix\n")

    N = np.shape(Coord)[0]
    RCoord = np.zeros(N)

    for l in range(N):
        RCoord[l] = dis.euclidean(np.zeros(3), Coord[l,:])
        
    return RCoord


def Contacts_Matrix (Dist, treshold, fig = False):

    #prende matrice distanze e ritorna matrice binaria
    # dei contatti definita sul limite treshold
    # NOVITA' : fig funge da nome, tanto l'importante è che non sia False !

    if np.shape(Dist)[0] != np.shape(Dist)[1]:
        raise ValueError("Matrix should be square matrix of relative distances")
    N = np.shape(Dist)[0]

    Cont_Matrix = np.zeros((N,N))

    for i in range (N):
        for k in range (N):
            if Dist[i,k] <= treshold:
                Cont_Matrix[i,k] = 1


    if fig:

        fig1 = plt.figure()
    
        ax1 = fig1.add_subplot(111)
        ax1.set_title(str(fig)+' Contacts Matrix'+'- '+str(treshold)+' angstrom')
    
        ax1.set_xlabel('Element index', labelpad = 10)
        ax1.set_ylabel('Element index')

    
        cax1 = ax1.matshow(Cont_Matrix, cmap = 'ocean', origin = 'lower')
    
        fig1.savefig(str(fig)+' Contact_matrix_'+str(treshold).replace('.', ','))

        plt.show()

    return Cont_Matrix

def Get_Covariance_Matrix(N, filename, path = './',):
    
    """

    Leggo matrice covarianza atomica in formato .dat 
    che esce da cova di GROMACS opzione -ascii

    """

    cov_matrix = np.zeros((3*N,3*N), dtype = float)

    with open (join(path,filename+'.dat'), 'r') as fin:
        lines = fin.readlines()
    

    if N != (np.sqrt(len(lines)/3)):
        raise ValueError ('Incompatibile lunghezza WTC {} con N da file di covarianza {}'.format(N, (np.sqrt(len(lines)/3))))

    init = 0
    for ii in range(3*N):
        print('riga {}'.format(ii)) 

        jj = 0
        line_scroll = lines[init:init+N]

        for line in line_scroll:
            array = np.array(line.split(), dtype= 'float')        
            cov_matrix[ii, jj] = array[0]
            cov_matrix[ii,jj+1] = array[1]
            cov_matrix[ii,jj+2] = array[2]
            jj+=3

        init+= N
    
    np.save(join(path,filename+'.npy'), cov_matrix)
    return cov_matrix

def Pearson_Matrix_from_Cov(cov_matrix):

    N = cov_matrix.shape[0]
    pearson = np.zeros((N,N))
    for ii in range(N):
        for jj in range(N):
            
            pearson[ii,jj] = cov_matrix[ii,jj]/(np.sqrt(cov_matrix[ii,ii]*cov_matrix[jj,jj]))

    return pearson

def Print_Cov_Matrix_BS(c_matrix, name, c_type, BS, idx, pearson = False, text = False, text_size = 5, path = './', **kwargs):
        
    if c_type == 'Atoms':
        ylabel = 'Atom number'
    elif (c_type == 'CAP') | (c_type == 'CAO'):
        ylabel = 'Residue Number'
    else:
        raise ValueError ('Select type of cov matrix\n Atoms or CAO\CAP \n')
    
    if text:
        text = []   
        size = text_size
        start = 0

        for ii in range(int(BS.size/size)):
            text.append([str(v) for v in BS[start:start+size]])
            start+=size
        if BS.size%size != 0:
            final = [str(v) for v in BS[start:]]
            if len(final) != size:
                while True:
                    final.append(' ')
                    if len(final) == size:
                        break
            text.append(final)


    f, ax = plt.subplots(1,1)
    
    cm = ax.imshow(c_matrix, cmap = 'jet')  

    if ((c_type == 'CAP') | (c_type == 'CAO')):
        ticks =  np.linspace(0, len(idx)-1, 8, dtype = int)
        tickslabels = idx[ticks]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(tickslabels)
        ax.set_yticklabels(tickslabels)

    if text:
        ax.table(text)
    if not pearson:
        ax.set_title('{} Covariance Matrix for {}'.format(c_type, name), pad = 1.)
    else:
        ax.set_title('{} Pearson Covariance Matrix for {}'.format(c_type, name), pad = 1.)
    
    plt.colorbar(cm) 
    ax.set_ylabel(ylabel)

    if 'clim' in kwargs:
        cm.set_clim(kwargs['clim'])
    ax.xaxis.set_ticks_position('top')

    plt.tight_layout()
    f.autolayout=True
    if not pearson:
        f.savefig(join(path, name+'_'+c_type+'_cov_matrix.pdf'), format = 'pdf', bbox_inches = 'tight')
    else:
        f.savefig(join(path, name+'_'+c_type+'_pear_matrix.pdf'), format = 'pdf', bbox_inches = 'tight')

def CAP_Cov_Matrix(cov_matrix, idx, save_filename, save_path = './'):

    """
    Calcolo covarianza degli atomi di CA per proteina e P per RNA
    a partire dalla matrice di covarianza atomica in cui CA e P hanno n_atom in idx
    """

    cov_matrix_CAP = np.zeros((idx.size, idx.size))

    for (ii, first_idx) in zip(range(idx.size), idx):

        jj_0 = first_idx*3

        for (jj, second_idx) in zip(range(idx.size), idx):

            ii_0 = second_idx*3
            cov_matrix_CAP[ii,jj] = np.mean(cov_matrix[ii_0:ii_0+3, jj_0:jj_0+3])
            
    np.savetxt(join(save_path,save_filename), cov_matrix_CAP)

    return cov_matrix_CAP




def Analyze_Bond_Residues (Cont_Matrix, structure_sizes, structure_names, first = ('Proteina', 0), second = ('RNA', 0), fig = False, verbose = False):

    """
    Funzione che analizza i legami tra i residui di recettore (proteina) e ligando (RNA) nel dettaglio a 
    partire dalla matrice dei contatti.

    Prende: 
        - Cont_Matrix, matrice (N+M*)x(N+M*) dove N è il numero di res della proteina, M* = M -1 dove
          M è il numero di residui dell'RNA (-1 uno per il problema del P)
        - Structure_sizes è una tupla dim = 2 che contiene il numero totale di residui di Proteina e RNA (conmpreso il +1)
        - structure_names è una tupla dim = 2 con due stringhe identificativive, ma attualmente inutile sarebbe 
          utile se poi si riuscisse a generalizzare il codice per altre biomolecole
        - first e second sono due tuple (<nome>, initial), dove initial è il residuo di partenza della struttura <nome>
          nel pdb. In base a questo l'output è diverso:
          first = 'RNA' stampa per ogni res di RNA quali res di 'proteina' a distanza data e viceversa

    Oss sul funzionamento:(assumendo first = 'rna')
        1) mi creo una matrice Asym_Cont dim = M* x N con i contatti tra residui RNA (righe) e residui Proteina (colonne)
           e mi creo una tupla Bonds che sarà della dim del numero di nucleotidi RNA
        2) NonZero è una tupla che contiene gli indici degli elementi zero
            es.

            NonZero = Asym_Cont.nonzero()
            (array([0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5, 6, 6, 6, 6,
                    6, 8, 8, 8]),
            array([ 82,  83,  84,  85, 162, 163, 164,  47,  48, 100, 133,  47,  48,
                    49,  14,  14,  76,  77,  14,  16,  75,  76,  77,  16,  17,  18]))

            vuol dire che per la riga 0 della matrice Asym ho un elemento di contatto all'indice 82, 83, etc.
            per questo nel secondo array della tupla ricomincia, perché scende di riga
        3) a questo punto grazie a cont_prev, cont_next e la funzione count_nonzero di numpy che mi permettono di scorrere la complessa struttura
           di NonZero, faccio un ciclo sul numero dei nucleotidi RNA (range(Asym_Cont.shape[0])) e creo una tupla Bond_i
           con il numero identificativo del residuo in contatto con l'i-esimo residuo di RNA
           OSS: questo è il VERO numero del residuo in quanto tramite operazione Bond_i = Bond_i + (int(NonZero[1][j]+second[1]),)
                sto aggiungendo la lunghezza totale della proteina (second) includendo l'1 di scarto dell'indice
                (ho controllato su pymol l'effettiva vicinanza maggiore del residuo che esce fuori)
            Alla fine del ciclo, questa tupla Bond_i diventa un elemento della tupla finale Bonds che sarà poi ritornata al main

    """


    Asym_Cont = Cont_Matrix[structure_sizes[0]:, :structure_sizes[0]]

    if fig:

        plt.figure()
        plt.matshow(Asym_Cont)
        plt.xlabel(structure_names[0]+" residue index")
        plt.ylabel(structure_names[1]+" residue index")
        plt.title("Asym_Cont Matrix of "+structure_names[0]+" + "+structure_names[1])
        plt.savefig(structure_names[0]+"+"+structure_names[1]+"_Asym_Cont_Matrix")
        plt.show()

    if (first == 'Proteina') & (second == 'RNA'):

        Asym_Cont = Asym_Cont.T
        
    NonZero =   Asym_Cont.nonzero()
    cont_prev = 0
    Bonds = ()
        
    for i in range(Asym_Cont.shape[0]):
        Bond_i = ()
        cont_next = np.count_nonzero(Asym_Cont[i, :])
        j = 0
        for j in range (cont_prev, cont_next + cont_prev):
                Bond_i = Bond_i + (int(NonZero[1][j]+second[1]),)

        if (cont_next!=0):
            if verbose:
                print ("Il "+str(i+1+first[1])+" residuo di "+first[0]+" lega con "+str(cont_next)+" residui di "+second[0]+ "\n I residui \n", np.array(Bond_i)+1,'\n')
            cont_prev = cont_prev + cont_next #aggiorno
        Bonds = Bonds + (Bond_i,)
        if (cont_next != len(Bond_i)):
              
            raise ValueError("Qualcosa è andato storto: taglia della tupla dei residui di contatto non coincide con il numero dei cont_nonzero\n\n")
        
    return Bonds

def Print_Bonds_HDOCK (Bonds, filename, distance, initial):

    f = open(filename, 'w')

    for (jj,Bond) in zip(range(len(Bonds)), Bonds):

        if (len(Bond) == 1 ):
            f.write("%d:A %d:B %d, " %(jj+initial+1, np.array(Bond[0]) +2 , int(distance)))
        elif (len(Bond) > 1):
            f.write("%d:A %d-%d:B %d, " %(jj+initial+1,np.array(Bond[0]) +2, np.array(Bond[len(Bond)-1])+2,int(distance)))
        else:
            pass
    
    f.close()

def Print_Protein_BS_old (res, Bonds, size, filename = 'BS.txt', prot_initial = 0, RNA_initial = 0):

    """
    stampa residui proteina coinvolti nel legame

    Dato che potrebbero esserci buchi nel file, è necessario avere un indice dei residui, così da fare match

    Il problema col residuo di RNA che manca non si verifica perché siamo prima 
    """
    Binding_Site = ()
    Binding_Site_RNA = ()

    count = 0
    idxs = []

    for kk in range(size):
        for Bond in Bonds:
            if (((kk+1+prot_initial) in Bond) & ((kk+1+prot_initial) not in Binding_Site)):

                Binding_Site = Binding_Site + (kk+1+prot_initial,)
                idxs = idxs + [kk+1,]
                count += 1
        
    Binding_Site = res[idxs]

    for ii in range(len(Bonds)):
          if len(Bonds[ii]) != 0:
              Binding_Site_RNA += (RNA_initial+ii,)

    Binding_Site_RNA = np.sort(Binding_Site_RNA)

    print ("La proteina utilizza nel legame (definito dalla distanza inserita) i residui\n", Binding_Site)
    print('per un totale di %d residui'%len(Binding_Site))

    BS = {'Prot' : Binding_Site, 'RNA' : Binding_Site_RNA}

    return BS

def Print_Protein_BS(Bonds, filename, path='./'):

    pass









def Parse_xvg(filename, path='./', skip = False):

    """Funzione che legge file xvg indicato da filename
    e ritorna i conseguenti  array
    """

    with open(path+filename) as file:

        temp    =  np.array([row.split() for row in filter(methodcaller('startswith', ''), file.readlines()[17:])], dtype=np.float64)
    

    if skip:
        return temp[np.arange(0,temp.shape[0], skip),0],  temp[np.arange(0,temp.shape[0], skip),1]
    
    else:
        return     temp[:,0],  temp[:,1]

def Parse_xvg_skip(filename, path='./', skip_lines = 18, skip = False):

    print('Sto aprendo file da {}'.format(join(path,filename)))
    """Funzione che legge file xvg indicato da filename
    e ritorna i conseguenti  array

    se sono più di due colonne, legge solo le prime 2

    """

    with open(join(path,filename), 'r') as f:
        parser    =  f.readlines()[skip_lines:]
    
    temp      =  np.zeros((len(parser),2))

    for row in range(len(parser)):

        temp[row]   =   np.array(str(parser[row]).split()[:2], dtype=np.float64)


    if skip:

        return temp[np.arange(0,temp.shape[0], skip),0],  temp[np.arange(0,temp.shape[0], skip),1]
    
    else:

        return     temp[:,0],  temp[:,1]

def Find_Nearest_2D(insieme, point):

    """

    Trova il punto 2d nell'insieme più vicino a 'point'
    l'insieme deve essere passato sotto forma di matrice [N x 2]
    i cui elementi  [:, 0] sono le x, [:,1] sono le y

    Ritorna una tupla contenente il punto dell'insieme più
    vicino a point e l'indice di riga dell'insieme

    """

    insieme =   np.asarray(insieme)
    dist    =   np.zeros(insieme.shape[0])

    for ii in range (insieme.shape[0]):

        dist[ii]   =   dis.euclidean([insieme[ii,0], insieme[ii,1]], [point[0], point[1]])
       
    
    idx     =   dist.argmin()


    return (insieme[idx, :], idx)

def Get_Protein_BS(protein, RNA, res, treshold, RNA_atom = 'P'):

    if RNA_atom == 'P':

        protein.Get_Atom_Coord(atom_name = 'CA')
        protein.Get_lenght()
        #protein.Get_Protein_Sequence()
        RNA.Get_Atom_Coord(atom_name='P')
        #RNA.Get_RNA_Sequence()
        RNA.Get_lenght()
        #sequence = np.concatenate([protein.sequence, RNA.sequence])
        Coord_P   = np.concatenate((protein.CA_Coord, RNA.P_Coord))
        Dist_P    = Dist_Matrix(Coord_P)
        Cont_P   = Contacts_Matrix(Dist_P, treshold)
        Bonds_P   = Analyze_Bond_Residues(Cont_P, (protein.lenght, RNA.lenght), ("protein", "RNA"), first=  ('RNA', 1), second = ('Proteina', protein.initial))
        BS_P      = Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['Prot']
        BS_RNA_P = Print_Protein_BS_old(res, Bonds_P, protein.lenght, prot_initial=protein.initial, RNA_initial=RNA.initial)['RNA']

        return BS_P, BS_RNA_P

    else: raise ValueError('Non hai ancora scritto questo codice')

def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)
    
def trimodal(x,mu1,sigma1,A1,mu2,sigma2,A2, mu3, sigma3, A3):
    
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x, mu3, sigma3, A3)
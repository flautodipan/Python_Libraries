from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import distance as dis
from mpl_toolkits.mplot3d import Axes3D

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from operator import methodcaller

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

            #estraggo indici partenza/arresto lettura dal pdb
            # lo faccio facendo la query delle righe con MODEL e prendendo gli 
            #indici di linea --> a meno di un addendo costante è la numerosità da prendere
            # per i valori identificati dalle keys ['ATOM']
    
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
            self.atoms  =   self.atoms[self.pdb.df['ATOM']['chain_id']==chain_id] 
            self.initial=   self.atoms.residue_number[0]
            #self.terminus =  self.pdb.df['OTHERS']
         
            
    def Get_All_Coord (self):
        
        pass
        
    def Get_CA_Coord(self):      
       
        self.CA = self.atoms[self.pdb.df['ATOM']['atom_name']=='CA']
        
        self.CA_xcoord = np.array(self.CA['x_coord'])
        self.CA_ycoord = np.array(self.CA['y_coord'])
        self.CA_zcoord = np.array(self.CA['z_coord'])
        

        self.CA_Coord = np.matrix([(self.CA_xcoord), (self.CA_ycoord), (self.CA_zcoord)])
        self.CA_Coord = np.matrix.transpose(self.CA_Coord)
        
        self.lenght   = self.CA_Coord.shape[0]

        if (self.lenght == len(self.CA_Coord)):

            print('Tutto ok \n Selezionati %s atomi di Carbonio pari al numero di residui \n Salvate le coordinate cartesiane di ogni atomo CA\n\n' %len(self.CA_Coord))
        
        else:
            raise ValueError ("Something went wrong in determining protein size\n\n")

class RNA:

    def __init__(self, pdb_filename, chain_id = 'A', from_web = False, model = 1):
        
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
            self.atoms  =   self.atoms[self.pdb.df['ATOM']['chain_id']==chain_id] 
            #self.terminus =  self.pdb.df['OTHERS']
    
        
    def Get_P_Coord(self): 

       
        self.P = self.atoms[self.pdb.df['ATOM']['atom_name']=='P']
        
        self.P_xcoord = np.array(self.P['x_coord'])
        self.P_ycoord = np.array(self.P['y_coord'])
        self.P_zcoord = np.array(self.P['z_coord'])
        

        self.P_Coord = np.matrix([(self.P_xcoord), (self.P_ycoord), (self.P_zcoord)])
        self.P_Coord = np.matrix.transpose(self.P_Coord)

        self.lenght   = self.P_Coord.shape[0] + 1 
        #il primo residuo 5' dell'RNA non ha il gruppo fosfato

        if (self.lenght == (len(self.P_Coord)+1)):
            print('Tutto ok \n Selezionati %s atomi di Fosforo pari al numero di residui -1\n Salvate le coordinate cartesiane di ogni atomo P\n\n' %len(self.P_Coord))
        else:
            raise ValueError("Something went wrong in determining size of RNA\n\n")
    
class Eigen_Trj():

    """
    Classe per analisi degli autovalori e autovettori della traiettoria effettuata in simulazione
    MD con il Software GROMACS
    --> indicare il path dove si trovano i file eigenval.xvg e eigenvec.trr (solitamente cartella
    supeiore a quella di questo script)
    """
    def __init__(self, path = './', fig = False):

        _ , self.EigenValues = Parse_xvg(path+'eigenval.xvg')

        if fig :
            
            plt.figure()
            plt.plot(self.EigenValues, '.')
            plt.title('Eigenvalues of MD trajetory')
            plt.xlabel('Eig index (ordered by variance significance)')
            plt.show()

    def Analyze_Variance_Explained_Ratio(self, percent_sigma):

        """
        
        Funzione che salva e stampa gli autovalori che spiegano
        la varianza sigma richiesta da utente 
    
        
        """

        Norma           = np.sum(self.EigenValues)
        Sigma_ERatio    = self.EigenValues/Norma

        n               = 0
        test            = 0.
        Sigma_Eigens    = []

        while True:

            Sigma_Eigens.append(self.EigenValues[n])
            test        += Sigma_ERatio[n]         
            n           += 1
            if test     > percent_sigma: break

        df              = pd.DataFrame(Sigma_Eigens)
        print('# Autovalori che spiegano la varianza richiesta = ', n, '\n\n')
        df
        print(Sigma_Eigens)

class Cluster_2DAnalysis():

    def __init__(self, xy, final_time, initial_time = 0, clustype = 'KMeans', fig = False):
        #time units = ps

        self.clustype       =  clustype
        self.xy             =  xy
        self.final_time     =  final_time
        self.initial_time   =  initial_time
        self.timestep       =  (self.final_time-self.initial_time)/self.xy.shape[0]
       

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

    def Silhouette_KMeans(self, kmax):

        self.Silh = []

        # dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
        for k in range(2, kmax+1):
            kmeans = KMeans(n_clusters = k).fit(self.xy)
            labels = kmeans.labels_
            self.Silh.append(silhouette_score(self.xy, labels, metric = 'euclidean'))

        self.n_clusters =  np.array(self.Silh).argmax() +2 
        
        plt.figure()
        plt.plot(np.arange(0,kmax-1,1)+2, self.Silh)
        plt.title('Silohuette check for k-means')
        plt.xlabel('k value')
        plt.savefig('Silh_test.png')

    def Clusterize(self, verbose = False):

        if (self.clustype == 'KMeans'):

            self.Cluster = KMeans(n_clusters=self.n_clusters).fit(self.xy)

            #definisco clusterchiefs : elementi (x,y) in xy più vicini
            # ai centroidi di ogni cluster
            self.clusterchiefs      =   ()
            self.clusterchiefs_idx  =   ()
            for jj in range(self.n_clusters):

                self.clusterchiefs      =  self.clusterchiefs + (Find_Nearest_2D(self.xy, (self.Cluster.cluster_centers_[jj,:]))[0],)
                self.clusterchiefs_idx  =  self.clusterchiefs_idx + (Find_Nearest_2D(self.xy, (self.Cluster.cluster_centers_[jj,:]))[1],)
                          
                if verbose:

                    print("Il punto più vicino al centroide del cluster " + str(jj+1) + " è \n" +str(self.clusterchiefs[jj])+"\nCorrisponde al frame %d\nOssia al tempo %f ps\n" %((self.clusterchiefs_idx[jj]+1), (self.clusterchiefs_idx[jj]+1)*self.timestep))

        else:
            raise ValueError ("Non è stato selezionato un algoritmo di clustering corretto\
                Scegliere:\n 'KMeans' per algoritmo k-means\n")    
        


    def Figure_Clusters(self):

        #color map data dall'appartenenza a un certo cluster
        c   =   self.Cluster.labels_

        fig = plt.figure()

        plt.scatter(self.xy[:,0], self.xy[:,1], c=c, cmap= 'viridis', s = 10)
        plt.scatter(self.Cluster.cluster_centers_[:,0],self.Cluster.cluster_centers_[:,1], marker='*', edgecolors='black')
        plt.xlabel(' PC 1 ')
        plt.ylabel(' PC 2 ')
        plt.title(str(self.n_clusters)+' Cluster representation by '+self.clustype+' algorithm')
        plt.show()
        fig.savefig('3 micros/Clusters.png')

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

def Analyze_Bond_Residues (Cont_Matrix, structure_sizes, structure_names, first = ('Proteina', 0), second = ('RNA', 0)):

    #1) estrarre matrice dei contatti nella parte che interessa 
    # per ora a due poi servirà giocare meglio su sizes e numero strutture
    # o forse fare un ciclo nel main
    # first = 'RNA' stampa per ogni res di RNA quali res di proteina a distanza data
    # first = 'Proteina' viceversa
    # i secondi argomenti delle tuple first e secondo sono gli initial : es nel pdb 4bs2 la proteina è risolta
    # solo per i domini centrail -> initial = 96


    Asym_Cont = Cont_Matrix[structure_sizes[0]:, :structure_sizes[0]]


    plt.figure()
    plt.matshow(Asym_Cont)
    plt.xlabel(structure_names[0]+" residue index")
    plt.ylabel(structure_names[1]+" residue index")
    plt.title("Asym_Cont Matrix of "+structure_names[0]+" + "+structure_names[1])
    plt.savefig(structure_names[0]+"+"+structure_names[1]+"_Asym_Cont_Matrix")
    plt.show()

    #2) estraggo informazioni dalla tupla di array ritornata da np.nonzero()
    
    #numero da aggiungere a indice per ottenere numero residuo reale, ma 
    # per adesso deve essere 2 in quanto non ho atomi di fosforo in primo res RNA

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
                Bond_i = Bond_i + (int(NonZero[1][j]+1+second[1]),)

        if (cont_next!=0):

            print ("Il "+str(i+1+first[1])+" residuo di "+first[0]+" lega con "+str(cont_next)+" residui di "+second[0]+ "\n I residui \n", np.array(Bond_i)+1,'\n')
            cont_prev = cont_prev + cont_next #aggiorno

        if (cont_next != len(Bond_i)):
                
            raise ValueError("Qualcosa è andato storto: taglia della tupla dei residui di contatto non coincide con il numero dei cont_nonzero\n\n")
        
        Bonds = Bonds + (Bond_i,)  
 
    
        
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

def Print_Protein_BS (Bonds, size, filename = 'BS.txt', initial = 0):

    """
    stampa residui proteina coinvolti nel legame

    """
    Binding_Site = np.zeros(size)

    count = 0

    for kk in range(size):
        for Bond in Bonds:
            
            if (((kk+1+initial) in Bond) & ((kk+1+initial) not in Binding_Site)):

                Binding_Site[count] = kk+1+initial
                count += 1
        
    Binding_Site = Binding_Site[Binding_Site != 0]
    Binding_Site = np.sort(Binding_Site)

    #print ("La proteina utilizza nel legame (definito dalla distanza inserita) i residui\n")
    
    f = open(filename, 'w')

    for Bond in Binding_Site:

        f.write("%d A\n"%Bond)
    
    f.close()
        
    return Binding_Site





def Parse_xvg(filename, path='./'):

    """Funzione che legge file xvg indicato da filename
    e ritorna i conseguenti  array
    """

    with open(path+filename) as file:

        temp    =  np.array([row.split() for row in filter(methodcaller('startswith', ' '), file.readlines())], dtype=np.float64)
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


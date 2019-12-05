# #libreria utile per programmi sulle proteine
# 
# OSS 30/09/2019 :
#                 PER ORA A PARTE UNA FUNZIONE DICHIARO METODI STATICI
#                 CIOè NON PROGAMMANDO OBJ ORIENTED 
#                 per ora solamente le CA_Coord diventano un attributo

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

            #self.terminus = self.pdb.df['OTHERS']
        
        #2B) prendo atomi e seleziono chain_id
        else:
            
            self.atoms  =   self.pdb.df['ATOM']
            self.atoms  =   self.atoms[self.pdb.df['ATOM']['chain_id']==chain_id] 
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
         
        print('Tutto ok \n Selezionati %s atomi di Carbonio pari al numero di residui \n Salvate le coordinate cartesiane di ogni atomo CA\n\n' %len(self.CA_Coord))
      

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
         
        print('Tutto ok \n Selezionati %s atomi di Fosforo pari al numero di residui -1\n Salvate le coordinate cartesiane di ogni atomo CA\n\n' %len(self.P_Coord))
     
def BioStructure (pdb_filename, n_structures, structures_names, structures_type, structures_chain_ids, structures_models):

    """
        if len(structures_type) is not n_structures:
            raise ValueError("Size of structure type must be n_structures\n")
        if (structures_type.type is not str):
            raise ValueError("structures type must be a tuple of strings \n Strings accepted are 'Protein' or 'RNA' for now\n\n")
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


def Contacts_Matrix (Dist, treshold):
    #prende matrice distanze e ritorna matrice binaria
    # dei contatti definita sul limite treshold

    if np.shape(Dist)[0] != np.shape(Dist)[1]:
        raise ValueError("Matrix should be square matrix of relative distances")
    N = np.shape(Dist)[0]

    Cont_Matrix = np.zeros((N,N))

    for i in range (N):
        for k in range (N):
            if Dist[i,k] <= treshold:
                Cont_Matrix[i,k] = 1

    return Cont_Matrix

    
class Eigen_Trj():

    """
    Classe per analisi degli autovalori e autovettori della traiettoria effettuata in simulazione
    MD con il Software GROMACS
    --> indicare il path dove si trovano i file eigenval.xvg e eigenvec.trr (solitamente cartella
    supeiore a quella di questo script)
    """
    def __init__(self, path = '../', fig = False):

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

    def __init__(self, xy, clustype = 'k-means', fig = False):
        
        self.clustype   =  clustype
        self.xy         =  xy

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

        plt.figure()
        plt.plot(range(kmax-1), self.Silh)
        plt.title('Silohuette check for k-means')
        plt.xlabel('k value')
        plt.savefig('Silh_test')


from operator import methodcaller

def Parse_xvg(filepath):

    """Funzione che legge file xvg indicato da filepath
    e ritorna i conseguenti  array
    """

    with open(filepath) as file:

        temp    =  np.array([row.split() for row in filter(methodcaller('startswith', ' '), file.readlines())], dtype=np.float64)
        return     temp[:,0],  temp[:,1]

"Sto facendo una modifica a cazzo de merda"

"Ne sto facendo un'altra vivaddio"
"Mannaggia alla madonna"
def Knowing_Spectrum_Peaks(self):

        self.peaks      =   Find_Highest_n_peaks(self.y, 4, height = 10., distance = 20, width = 1.)

    def Knowing_VIPA_Peaks(self):

        self.VIPA_peaks =   Find_Highest_n_peaks(self.y_VIPA, 5, height = 10., distance = 100, width = 1.)

    def Knowing_Peaks(self):

        peaks    =   ()

        for (n_peaks, y, dist)    in      zip((4,5),(self.y, self.y_VIPA), (20, 100)):

            peaks    =   peaks    +   (Find_Highest_n_peaks(y, n_peaks, height = 10., distance = dist, width = 1.),)

        self.peaks     =   peaks[0]
        self.VIPA_peaks         =   peaks[1]

    def How_Many_Peaks_To(self, n_peaks = 4, delta = 1., distance = 20., width = 5., treshold = 5, fig = False, verbose = False, i_know_it_is = False):

        pk      =   find_peaks(self.y, height = treshold, distance = distance, width = width)
        height  =   np.max(pk[1]['peak_heights'])/2
        
        while True:

            pk      =   find_peaks(self.y, height = height, distance = distance, width = width)

            if (height > treshold) & (pk[0].size == n_peaks):

                if verbose:
                    print("Ho trovato valore dell'altezza per avere %d picchi: %f\n"%(n_peaks, height), pk)
                    _ = Analyze_Peaks(self.x_freq, self.y, 'GHz', fig = fig, verbose = verbose, height= height, distance = distance, width = width )
    
                break
            
            elif (height <= treshold):

                print(pk)
                raise ValueError('Errore: superata altezza minima %f\nQualcosa è andato storto'%(treshold))

            else: 

                height-=delta
        
        #check che i picchi Brillouin siano dentro agli elastici, premesso che so che possano essere più alti degli elastici

        condition_peaks =   ((self.y[pk[0][0]] < self.y[pk[0][1]]) | ((self.y[pk[0][3]] < self.y[pk[0][2]])))
        
        if condition_peaks &  (i_know_it_is == False):
            raise ValueError("Picchi Brillouin non sono interni o sono più alti degli elastici, controllare")

        self.spectrum_cut_height        =   height
        self.spectrum_peaks_dist        =   distance
        

    def Interpolate_VIPA (self, freq):

        # funzione che interpola il valore della funzione di trasf
        # richiesta nel valore freq a partire dai dati sperim
        # funziona per valore singolo e per array di valori


        freq = np.array(freq)

        if (freq.size == 1):
        
            _ , idx             =   Find_Nearest(self.x_VIPA_freq, freq)
            x_fittino           =   np.array([self.x_VIPA_freq[idx-1], self.x_VIPA_freq[idx], self.x_VIPA_freq[idx+1]])
            y_fittino           =   np.array([self.y_VIPA[idx-1], self.y_VIPA[idx], self.y_VIPA[idx+1]])
            Parameters          =   np.polyfit(x_fittino, y_fittino, 2)

            interpolate         =   ((freq**2) *(Parameters[0])) + (freq * (Parameters[1])) + (Parameters[2])

        else:

            interpolate         =   np.zeros(np.size(freq))
            _ , idx             =   Find_Nearest_Array(self.x_VIPA_freq, freq)

            for ii in range(np.size(freq)): 

                x_fittino           =   np.array([self.x_VIPA_freq[idx[ii]-1], self.x_VIPA_freq[idx[ii]], self.x_VIPA_freq[idx[ii]+1]])
                y_fittino           =   np.array([self.y_VIPA[idx[ii]-1], self.y_VIPA[idx[ii]], self.y_VIPA[idx[ii]+1]])
                Parameters          =   np.polyfit(x_fittino, y_fittino, 2)
                interpolate[ii]     =   ((freq[ii]**2) *(Parameters[0])) + (freq[ii] * (Parameters[1])) + (Parameters[2])

        return              interpolate

        ################
        # ##############  funzioni fuori da classe spectrum



def Parse_Parameter_Save(parameter_file = 'save_params.txt', path = ''):

    parameters  = ()
    f   =   open(path+parameter_file, 'r')
    lines   =   f.readlines()
    f.close()
    indices =   np.arange(0,len(lines), 3)
    for ii in indices:
            
        
        a   =   np.array(lines[ii].split()[1:], dtype= np.float64)
        b   =   np.array(lines[ii+1].split()[:], dtype= np.float64)
        c   =   np.array(lines[ii+2].replace(']', ' ').split()[:], dtype= np.float64)

        parameters  =   parameters  +   (np.concatenate((a,b,c)),)

    return parameters


def Save_Gaussia_Info(matrix, fitted, out_filename, path = './'):

    with open(path+out_filename, 'w') as f_out:
        
        for (ii,jj) in fitted:

            f_out.write(json.dumps(matrix[ii][jj].Fit_Params[list(cols_gauss)].to_dict())+'\n')

    print('Stampato parametri fit su file '+path+out_filename)


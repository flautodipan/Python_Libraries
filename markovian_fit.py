######################################################################################################################


#######   ||||  ||||  ||||   MARKOVIAN FIT: - opero fit markoviano con tutti i parametri
#######    ||    ||    ||                     tranne Delta e tau (quindi anche Gauss)
#######    ||    ||    ||
#######    ||    ||    ||
#######   ||||  ||||  ||||


import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Experimentum    import  *
from        Alessandria         import  *

##########
#variabili globali
cols      = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')

############

def Exec(matrix, boni, excluded, inputs, analysis_path):

    n_rows  =   len(matrix)
    n_cols  =   len(matrix[0])
    dim     =   n_cols*n_rows
    fit = ()

    if inputs['DEFAULT']['markov_recover'] == False:

        ####### faccio er fit

        isolated = Get_Isolated_Elements(excluded)
        percents        =   eval(inputs['DEFAULT']['percents_bound_markov'])
        for (ii,jj) in boni:
            print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))

            if (ii,jj) in isolated:

                matrix[ii][jj].Get_p0(matrix[0][0].p0.values[0], cols_mark)
            
            else:

                if ii == 0:
                    matrix[ii][jj].Get_p0(matrix[ii][jj-1].p0.values[0], cols_mark)
                else:
                    matrix[ii][jj].Get_p0(matrix[ii-1][jj].p0.values[0], cols_mark)

            matrix[ii][jj].Get_Fit_Bounds(percents, cols_mark)
            fit = fit + ((matrix[ii][jj].Non_Linear_Least_Squares_Markov(bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values), max_nfev = int(inputs['DEFAULT']['max_nfev_markov'])),(ii,jj)),)

            del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].res_lsq, matrix[ii][jj].bounds


        non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)
        Save_Fit_Info(fit, filename = 'markov_fit.txt', path=analysis_path)
        Save_Fit_Parameters(matrix, fitted, out_filename = 'markov_fit_params.txt', path = analysis_path)

    else:

        ############## faccio er ricovery

        with open(analysis_path+'markov_fit.txt', 'r') as fin:
            fit     =   eval(fin.read())

        non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

        with open(analysis_path+'markov_fit_params.txt', 'r') as fin:
            lines   =   fin.readlines()

        if (len(fitted) != len(lines)):
            raise ValueError("Incompatibilità tra lunghezza file parametri ({}) e informazioni fit ({})".format(len(fitted), len(lines)))

        for (line, (ii,jj)) in zip(lines, fitted) :
            matrix[ii][jj].Recover_Fit_Params(line)

    return matrix

######################################################################################################################


#######   ||||   ||      ||   TOT FIT: - opero fit totale con tutti i parametri
#######    ||     ||    ||               viscoelastici ma senza più i gaussiani
#######    ||      ||  ||
#######    ||       ||||   
#######   ||||       ||

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
    fit_tot = ()
    percents = eval(inputs['DEFAULT']['percents_bound_tot'])
    p_gauss = matrix[0][0].p0[list(cols_gauss)].values[0]

    for (ii,jj) in boni:

        print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
        matrix[ii][jj].Initials_Parameters_from_Markov(matrix[ii][jj].Fit_Params.T['Values'].values)
        matrix[ii][jj].Get_Fit_Bounds(percents, columns = cols_real)
        fit_tot =   fit_tot + (((matrix[ii][jj].Non_Linear_Least_Squares(p_gauss, cols_real, bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values), max_nfev = int(inputs['DEFAULT']['max_nfev_tot']))), (ii,jj)),)

    #after fit
    non_fitted_tot, accomplished_tot, exceded_tot, fitted_tot = Unpack_Fit(fit_tot)
    Save_Fit_Info(fit_tot, filename = 'tot_fit.txt', path=analysis_path)
    Save_Fit_Parameters(matrix, fitted_tot, out_filename = 'tot_fit_params.txt', path = analysis_path)

    


#########################                                                   #######################
#########################           VERSIONE 1 FIT CON MODELLO MARKOVIANO   #######################

#%%


import      os
now_path        =   '../BRILLOUIN/TDP43/ARS_13_02/'
spectra_filename    =   'ARS_13_02'
VIPA_filename       =   'NO_ARS_13_02_VIPA_not_sat.tif'

os.system('cd ' + now_path +' & mkdir ' + now_path+'analysis/')
analysis_path            =   now_path +'analysis/'


cols      = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_mark   = ('Co', 'Omega', 'Gamma', 'delta_width', 'delta_amplitude', 'A', 'mu', 'sigma', 'shift', 'offset')
cols_real   = ('Co', 'Omega', 'Gamma', 'Delta', 'tau', 'delta_width', 'delta_amplitude','shift', 'offset')
cols_gauss  = ( 'A', 'mu', 'sigma')


recover_markov = False

#%%
import      numpy               as      np
import      matplotlib.pyplot   as      plt
from        lib_Experimentum    import  *
from        Alessandria         import  *
import      time
super_start         =   time.process_time()
tempo               =   ()

syg_kwargs          =   {'height': 20, 'distance': 20, 'width': 3.}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
syg_kwargs_brill    =  {'height': 15.8, 'distance': 25, 'width': 3.}
VIPA_treshold       =   6
sat_height          =   50000
sat_width           =   13.5

# %%
#0) Acquisisco dati e inizializzo oggetti Spectrum per ognuno su una matrice (n_rows, n_cols)
#   e compio alcune operazioni di sistema utili per salvataggio dati

#import dati spettro
dati    =   Import_from_Matlab(spectra_filename, now_path, var_name = 'y3')
n_rows  =   len(dati)
n_cols  =   len(dati[0])
#matrix, rows, cols = Initialize_Matrix(1,8, 2+1, 10)
matrix, rows, cols = Initialize_Matrix(0,0, n_rows, n_cols)
dim     =   len(rows)*len(cols)
#definisco quantità di interesse

invisible           =   () 
saturated           =   () 
brillouin_higher    =   ()
brillouin_highest   =   ()
boni                =   ()
excluded            =   ()

# %%
#1) Acquisisco VIPA e Spettri
start = time.process_time()
matrix[0][0].Get_VIPA_tif(VIPA_filename, now_path, offset = 183.)


for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
        
        matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = 183., cut = False, cut_range = (200, 600))
        matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs)
        matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
        matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA

        check   =   matrix[ii][jj].Check_Spectrum(saturation_height = sat_height, saturation_width = sat_width)
        
        if (check == 1):

            saturated   =   saturated   +   ((ii,jj), )
            excluded    =   excluded    +   ((ii,jj), )

        elif (check == 2):

            brillouin_higher    =   brillouin_higher    +   ((ii,jj),)
            boni                =   boni                +   ((ii,jj),)

        elif (check == 3):

            invisible           =   invisible           +   ((ii,jj),)
            excluded    =   excluded    +   ((ii,jj), )
        
        else:

            boni                =   boni                +   ((ii,jj),)


acq_time    =   time.process_time()-start
tempo       =   tempo + (('acquisizione', acq_time),)
print('tempo impiegato per acquisizione spettri: %f s'%(acq_time))
print('Totale spettri saturati : %d\n'%(len(saturated)), saturated)


print('Totale spettri con Brillouin più alti : %d\n'%(len(brillouin_higher)), brillouin_higher)
print('Totale spettri con Brillouin invisibili: %d\n'%(len(invisible)), invisible)
print('Totale spettri boni :   %d'%(len(boni)))
print('Totale   spettri : %d\ndi cui %d inutilizzabili'%(dim, len(invisible)+len(saturated)))
print('ossia il %3.2f percento'%(float(len(invisible)+len(saturated))*100/dim))
print('Di cui il {} invisibili e il {} saturati'.format(len(invisible)*100/dim, len(saturated)*100/dim))

# %%
#2) Faccio operazioni di modifica spettro

start = time.process_time()

matrix[0][0].How_Many_Peaks_To_VIPA(treshold = VIPA_treshold, **syg_kwargs_VIPA)
matrix[0][0].Fit_Pixel2GHz()
matrix[0][0].VIPA_Pix2GHz()
matrix[0][0].Spectrum_Pix2GHz()
matrix[0][0].Get_Spectrum_4_Peaks_by_Height()
matrix[0][0].Cut_n_Estimate_Spectrum(estimate = True, distanza = 0.12)
matrix[0][0].Fit_VIPA_Gaussian()

####recupero gli invisibili buoni
##DA INSERIRE


for ii in range(len(rows)):
    for jj in range(len(cols)):
        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
        
        
        matrix[ii][jj].x_VIPA_freq   =   matrix[0][0].x_VIPA_freq
        matrix[ii][jj].y_VIPA        =   matrix[0][0].y_VIPA
        
        matrix[ii][jj].Poly2GHz      =   matrix[0][0].Poly2GHz
        
        
        if (ii,jj) in brillouin_higher:

            matrix[ii][jj].Spectrum_Pix2GHz(align = False)
        else:
            matrix[ii][jj].Spectrum_Pix2GHz()

        if ((ii,jj) != (0,0)) & ((ii,jj) not in excluded):
                    

            if (ii,jj) in brillouin_higher:
            
                matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs_brill)
                matrix[ii][jj].Get_Spectrum_4_Peaks_by_Order()

                # se è Brillouin highest dx
                if (matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][2]] > matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][1]]) & (matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][2]] > matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][3]]):
                    print('ok è brillouin highest dx')
                    brillouin_highest = brillouin_highest + ((ii,jj),)
                    matrix[ii][jj].Align_Brillouin_Highest('dx')

                #se è brillouin highest sx
                elif (matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][1]] > matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][2]]) & (matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][1]] > matrix[ii][jj].y[matrix[ii][jj].peaks['peaks_idx'][0]]):
                    print('ok è brillouin highest sx')
                    brillouin_highest = brillouin_highest + ((ii,jj),)
                    matrix[ii][jj].Align_Brillouin_Highest('sx')
                            
            else :
                #boni
                matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                

            matrix[ii][jj].Cut_n_Estimate_Spectrum(distanza = 0.25)
            #del matrix[ii][jj].x, matrix[ii][jj].x_VIPA, matrix[ii][jj].Poly2GHz, matrix[ii][jj].peaks

mod_time    =   time.process_time()-start
tempo       =   tempo + (('modifica',mod_time), )
print('tempo impiegato per modifica spettri: %f s'%(mod_time))


# salvo info spettri e VIPA
Save_XY_position(matrix, len(rows), len(cols), path = analysis_path)
Save_XY_VIPA(matrix[0][0].x_VIPA_freq, matrix[0][0].y_VIPA, path = analysis_path)
print('\n I saved xy info on xy.txt and xy_VIPA.txt in your analysis directory\n\n')


#%%
#3) faccio il fit markoviano

if recover_markov == False:
        
    print('\n\n You chose to do the markovian fit\n\n')
    fit                 =   ()
    start = time.process_time()
    isolated = Get_Isolated_Elements(excluded)
    percents        =   ('positive', 0.2, 'positive', 'positive', 'positive', 0.1, 0.1, 0.1,  np.inf, np.inf)
    
    for (ii,jj) in boni:

        print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))

        if (ii,jj) in isolated:

            matrix[ii][jj].Get_p0(matrix[0][0].p0.values[0], cols_mark)
        
        else:

            if ii == 0:
                matrix[ii][jj].Get_p0(matrix[ii][jj-1].p0.values[0], cols_mark)
            else:
                matrix[ii][jj].Get_p0(matrix[ii-1][jj].p0.values[0], cols_mark)

        matrix[ii][jj].Get_cost_markov(matrix[ii][jj].p0.values[0])
        print('Cost before fitting = {}'.format(matrix[ii][jj].cost_markov))
        matrix[ii][jj].Get_Fit_Bounds(percents, cols_mark)
        fit = fit + ((matrix[ii][jj].Non_Linear_Least_Squares_Markov(bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values),  max_nfev = 100),(ii,jj)),)
        matrix[ii][jj].Get_cost_markov(matrix[ii][jj].Markov_Fit_Params.values[0])
        print('Cost after fitting = {}\n'.format(matrix[ii][jj].cost_markov))

        del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].y_markov_convolution
        #print(fit)
else:

    print('\n\n You chose to SKIP the markovian fit and recover info \n\n')

    ############## faccio er ricovery

    with open(analysis_path+'markov_fit.txt', 'r') as fin:
        fit     =   eval(fin.read())

    non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

    with open(analysis_path+'markov_fit_params.txt', 'r') as fin:
        lines   =   fin.readlines()

    if (len(fitted) != len(lines)):
        raise ValueError("Incompatibilità tra lunghezza file parametri ({}) e informazioni fit ({})".format(len(fitted), len(lines)))

    for (line, (ii,jj)) in zip(lines, fitted) :
        matrix[ii][jj].Recover_Markov_Fit_Params(line)

markov_time     =   time.process_time()-start
tempo           =   tempo + (('fit markoviano', markov_time),)

print('tempo impiegato per fit markoviani: %f s'%(markov_time))
print('tempo impiegato ore = %3.2f'%(markov_time/3600))


# 4) after - fit markoviano

non_fitted, accomplished, exceded, fitted = Unpack_Fit(fit)

#too_markov         =   Whose_Param_Too_High('Gamma', 2., matrix, fitted)
Save_Fit_Info(fit, filename = 'markov_fit.txt', path=analysis_path)
Save_Markov_Fit_Parameters(matrix, fitted, out_filename = 'markov_fit_params.txt', path = analysis_path)
Save_y_markov_fit(matrix, boni, path = analysis_path)
Save_cost_markov(matrix, boni, path = analysis_path)

#%%
##################################################################################################################################
# fit tot

fit_tot = ()
percents = (0.2, 0.1, 0.15, 'positive', 'positive', 0.15, 0.15, np.inf, np.inf)

start = time.process_time()
for (ii,jj) in boni:

    print('Passo row = %d/%d col = %d/%d'%(ii,len(rows)-1, jj, len(cols)-1))
    p_gauss = matrix[ii][jj].Markov_Fit_Params[list(cols_gauss)].values[0]
    matrix[ii][jj].Initials_Parameters_from_Markov(matrix[ii][jj].Markov_Fit_Params.T['Values'].values)
    matrix[ii][jj].Get_Fit_Bounds(percents, columns = cols_real)
    matrix[ii][jj].Get_cost_tot(matrix[ii][jj].p0.values[0], p_gauss)
    print('\nCost before fitting = {}\n'.format(matrix[ii][jj].cost_tot))
    fit_tot =   fit_tot + (((matrix[ii][jj].Non_Linear_Least_Squares(p_gauss, cols_real, bound = (matrix[ii][jj].bounds['down'].values, matrix[ii][jj].bounds['up'].values), max_nfev = 35)), (ii,jj)),)
    matrix[ii][jj].Get_cost_tot(matrix[ii][jj].Tot_Fit_Params.values[0], p_gauss)
    print('\nCost after fitting = {}\n'.format(matrix[ii][jj].cost_tot))
    #del matrix[ii][jj].y_Gauss_markov_convolution, matrix[ii][jj].res_lsq, matrix[ii][jj].bounds



tot_time     =   time.process_time()-start
tempo           =   tempo + (('fit tot', tot_time),)

#after fit
non_fitted_tot, accomplished_tot, exceded_tot, fitted_tot = Unpack_Fit(fit_tot)

Save_Fit_Info(fit_tot, filename = 'tot_fit.txt', path=analysis_path)
Save_Tot_Fit_Parameters(matrix, fitted_tot, out_filename = 'tot_fit_params.txt', path = analysis_path)
Save_y_fit(matrix, boni, path = analysis_path)
Save_cost_tot(matrix, boni, path = analysis_path)


super_time = super_start - time.process_time()

print('tempo impiegato per esecuzione dello script ore = %3.2f\n '%(super_time/3600))

for (what,t) in tempo:

    print('di cui %f secondi =  %f  ore in %s \n' %(t, t/3600, what))


# %%


#######   ||||  ||||    DATA ACQUISITION AND TREATMENT : - acquisiamo i file di INPUT
#######    ||    ||                                  - eseguiamo operazioni su cartelle
#######    ||    ||
#######    ||    ||
#######   ||||  ||||

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

def Exec(analysis_path, inputs):

        
    #0) import dati spettro da file 

    dati    =   Import_from_Matlab(inputs['I/O']['spectra_filename'], inputs['I/O']['now_path'], var_name = 'y')
    if inputs['DEFAULT']['dim'] == 'tot':
        n_rows  =   len(dati)
        n_cols  =   len(dati[0])
        dim     =   n_cols*n_rows
    elif inputs['DEFAULT']['dim'] == 'partial':
        n_rows = int(inputs['DEFAULT']['n_rows'])
        n_cols = int(inputs['DEFAULT']['n_cols'])
        dim     =   n_cols*n_rows

    else:
        raise ValueError('Value of dim inserted in config file is wrong\n Choose either "tot" or "partial", indicating in this case n_rows and n_cols \n')
    
    matrix = Initialize_Matrix(n_rows,n_cols)
    syg_kwargs   =   {'height': float(inputs['DEFAULT']['norm_height']), 'distance':float(inputs['DEFAULT']['norm_distance']), 'width': float(inputs['DEFAULT']['norm_width'])}
    syg_kwargs_brill = {'height': float(inputs['DEFAULT']['brill_height']), 'distance':float(inputs['DEFAULT']['brill_distance']), 'width': float(inputs['DEFAULT']['brill_width'])}
    cut_range =  (int(inputs['DEFAULT']['cut_left']), int(inputs['DEFAULT']['cut_right']))

    ###########################################################################################################################

    #definisco quantità di interesse

    invisible           =   () 
    saturated           =   () 
    brillouin_higher    =   ()
    boni                =   ()
    excluded            =   ()

    #################### 1) Acquisisco VIPA e Spettri

    matrix[0][0].Get_VIPA_tif(inputs['I/O']['VIPA_filename'], inputs['I/O']['now_path'], offset = float(inputs['DEFAULT']['offset']))
    print(cut_range)
    print(type(cut_range))
    
    for ii in range(n_rows):
        for jj in range(n_cols):
            print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))
            
            matrix[ii][jj].Get_Spectrum(y = np.resize(dati[ii][jj],np.max(dati[ii][jj].shape)) , offset = float(inputs['DEFAULT']['offset']), cut = inputs['DEFAULT']['cut'], cut_range = cut_range)
            matrix[ii][jj].Get_Spectrum_Peaks(**syg_kwargs)
            matrix[ii][jj].x_VIPA   =   matrix[0][0].x_VIPA
            matrix[ii][jj].y_VIPA   =   matrix[0][0].y_VIPA
            check   =   matrix[ii][jj].Check_Spectrum(saturation_width = float(inputs['DEFAULT']['saturation_width']))
            if (check == 1):
                saturated   =   saturated   +   ((ii,jj), )
                excluded    =   excluded    +   ((ii,jj), )
            elif (check == 2):
                brillouin_higher    =   brillouin_higher    +   ((ii,jj),)
                boni                =   boni                +   ((ii,jj),)
            elif (check == 3):
                invisible           =   invisible   +   ((ii,jj),)
                excluded            =   excluded    +   ((ii,jj),)
            else:
                boni                =   boni                +   ((ii,jj),)


    ###################### 2) Modifico gli spettri


    matrix[0][0].How_Many_Peaks_To_VIPA(treshold = float(inputs['DEFAULT']['VIPA_height_treshold']))
    matrix[0][0].Fit_Pixel2GHz(fig = True)
    matrix[0][0].VIPA_Pix2GHz(fig=True)

    matrix[0][0].Spectrum_Pix2GHz(fig=True)
    matrix[0][0].Get_Spectrum_4_Peaks_by_Height()
    matrix[0][0].Cut_n_Estimate_Spectrum(estimate = True, columns = cols, distanza = float(inputs['DEFAULT']['width_factor']))
    matrix[0][0].Fit_VIPA_Gaussian()


    for ii in range(n_rows):
        for jj in range(n_cols):

            print('Passo row = %d/%d col = %d/%d'%(ii,n_rows, jj, n_cols))

            if ((ii,jj) != (0,0)) & ((ii,jj) not in excluded):
                        
                matrix[ii][jj].x_VIPA_freq   =   matrix[0][0].x_VIPA_freq
                matrix[ii][jj].y_VIPA        =   matrix[0][0].y_VIPA
                
                matrix[ii][jj].Poly2GHz      =   matrix[0][0].Poly2GHz
                matrix[ii][jj].Spectrum_Pix2GHz()

                if (ii,jj) in brillouin_higher:
                
                    matrix[ii][jj].Get_Spectrum_Peaks(height = 5., distance = 50, width = 5)
                    matrix[ii][jj].Get_Spectrum_4_Peaks_by_Order()
                
                else :
                    #boni
                    matrix[ii][jj].Get_Spectrum_4_Peaks_by_Height()
                    

                matrix[ii][jj].Cut_n_Estimate_Spectrum(distanza = 0.25)

                del matrix[ii][jj].x, matrix[ii][jj].x_VIPA, matrix[ii][jj].Poly2GHz, matrix[ii][jj].peaks, matrix[ii][jj].offset



    # salvo info spettri e VIPA
    Save_XY_position(matrix, n_rows, n_cols, path = analysis_path)
    Save_XY_VIPA(matrix[0][0].x_VIPA_freq, matrix[0][0].y_VIPA, path = analysis_path)
    print('\n I saved xy info on xy.txt and xy_VIPA.txt in your analysis directory\n\n')

    with open(analysis_path+'log.txt', 'a') as f_log:
        f_log.write(('\n I saved xy info on xy.txt and xy_VIPA.txt in your analysis directory\n\n'))
        f_log.write('\nTotale spettri saturati : %d\n %s \n\n'%(len(saturated), str(saturated)))
        f_log.write('\nTotale spettri con Brillouin più alti : %d\n %s \n\n'%(len(brillouin_higher), str(brillouin_higher)))
        f_log.write('\nTotale spettri con Brillouin invisibili: %d\n %s \n\n'%(len(invisible), str(invisible)))
        f_log.write('\nTotale spettri boni :   %d'%(len(boni)))
        f_log.write('\nTotale   spettri : %d\ndi cui %d inutilizzabili'%(dim, len(invisible)+len(saturated)))
        f_log.write('\nossia il %3.2f percento'%(float(len(invisible)+len(saturated))*100/dim))


    return (matrix, boni, excluded)


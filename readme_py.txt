#I/O 
now_path        =   '../BRILLOUIN/TDP43/NO_ARS_12_02/'
spectra_filename    =   'NO_ARS_12_02'
VIPA_filename       =   'NO_ARS_12_02_VIPA_quasisat.tif'
log_file            =   'log_'+spectra_filename
analysis_dir        =   'analysis_new_nocut/'

#operatives

#esclusi a mano
to_add              =   []

syg_kwargs          =   {'height': 119, 'distance': 31, 'width': 3.}
syg_kwargs_VIPA     =   {'distance':70, 'width': 1}
syg_kwargs_brill    =  {'height': 23, 'distance': 31, 'width': 3.}
VIPA_treshold       =   6
sat_height          =   50000
sat_width           =   13.5
almost_treshold     =   15000

#quanto mi allontano dal VIPA
pre_cut             =   False
cut                 =   False

mean_dist_01 = 37
mean_dist_23 = 34
#markov_fit
p0_normal = np.array([ 1.07378474e-01,  7.57148558e+00,  1.49128813e-01,  1.19015861e-01,
        1.448930518e-01,  8.34614271,  4.79747192e+03, -1.00904973e+01,
        1.58007162e+01,  2.11019859e-01, -3.10388495e-01])
p0_brillouin = np.array([ 1.07378474e-01,  7.57148558e+00,  1.49128813e-01,  1.19015861e-01,
        1.48930518e-01,  2.34614271e-01,  4.79747192e+03, -1.00904973e+01,
        1.58007162e+01,  2.11019859e-01, -3.10388495e-01])
p0_almost = np.array([ 1.08633225e-01,  7.70983143e+00,  1.58967633e-01,  1.70455195e+00,
        6.40427573e-01,  2.20351667e+00,  5.23638443e+03, -9.18245455e+00,
        1.43788115e+01,  2.73907418e-01,  8.73821212e+00])

recover_markov = False
percents_markov     =   ('positive', 0.2, 'positive', np.inf, 'positive', 'positive', 0.2, 0.01, 0.01,  np.inf, np.inf)
#tot fit
skip_tot = True
percents_tot        = (0.1, 0.1, 0.5, 'positive', 'positive', 0.15,  0.15, 0.15, np.inf, np.inf)

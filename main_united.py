#%%
######################################################################################################################
#######   ||||    PREAMBOLO: - acquisiamo i file di INPUT
#######    ||                - eseguiamo operazioni su cartelle
#######    ||   
#######    ||    
#######   |||| 

import  Spectra_Acquisition_Modify
import markovian_fit
import tot_fit
import configparser
import os
import time

super_start         =   time.process_time()
tempo               =   ()

### ACQUISISCO INPUT 

inputs = configparser.ConfigParser()
with open('K27M.ini', 'r') as f:
    inputs.read_file(f)


### OPERAZIONI SU CARTELLE

while True:

    if os.path.exists(inputs['I/O']['now_path'] +'analysis/'):

        print("\nDirectory named analysis already exists.\nPress 'o' to overwrite it, or 'n' to generate a new directory for analysis\n")
        ans = input()
        if (ans == 'n'):

            print("Insert name of new directory\n")
            new_name = input()
            os.system('cd '+inputs['I/O']['now_path'] +' && mkdir '+new_name)
            analysis_path = inputs['I/O']['now_path']  + new_name +'/'
            break

        elif (ans == 'o'):
            os.system('cd '+inputs['I/O']['now_path'] +' && rm -r analysis_backup/')
            os.system('cd '+inputs['I/O']['now_path'] +' && cp -r analysis/ analysis_backup/')
            print('I backed up your directory, for any inconvenience...\n')
            #os.system('cd '+inputs['I/O']['now_path'] +' && mkdir analysis/')
            analysis_path = inputs['I/O']['now_path']  + 'analysis/'
            break
        
        else:
            print('\nValue inserted not correct\n Try again motherfucker\n')
    else:
        os.system('cd '+inputs['I/O']['now_path'] +' && mkdir analysis/')
        analysis_path = inputs['I/O']['now_path']  + 'analysis/'
        break

with open(analysis_path+'log.txt', 'w') as f_log:
    f_log.write('#This is a log file: info about functionment on the script\n')

print('\nOn file "log" in your working directory, you will find some informations on script running\n')

#%%
######################################################################################################################

#######   ||||  ||||    DATA ACQUISITION AND TREATMENT : - acquisiamo i file di INPUT
#######    ||    ||                                  - eseguiamo operazioni su cartelle
#######    ||    ||
#######    ||    ||
#######   ||||  ||||


print("\n\nI'm beginning acquiring and modifying spectra\n\n")
start = time.process_time()
matrix, boni, excluded = Spectra_Acquisition_Modify.Exec(analysis_path, inputs)
acq_time    =   time.process_time()-start
tempo       =   tempo + (('acquisizione', acq_time),)
with open(analysis_path+'log.txt', 'a') as f_log:
    f_log.write('\ntempo totale impiegato per acquisizione e modifica spettri: %f s\n'%(acq_time))


#%%
#######################################################################################################################

#######   ||||  ||||  ||||   MARKOVIAN FIT: - opero fit markoviano con tutti i parametri
#######    ||    ||    ||                     tranne Delta e tau (quindi anche Gauss)
#######    ||    ||    ||
#######    ||    ||    ||
#######   ||||  ||||  ||||

print("\n\nI'm beginning markovian fit\n\n")
start = time.process_time()
matrix = markovian_fit.Exec(matrix, boni, excluded, inputs, analysis_path)
markov_time     =   time.process_time()-start
tempo           =   tempo + (('fit markoviano', markov_time),)

with open(analysis_path+'log.txt', 'a') as f_log:
    f_log.write('\ntempo impiegato per fit markonviani: %f s\n'%(markov_time))
    f_log.write('\ntempo impiegato ore = %3.2f\n'%(markov_time/3600))

#%%
######################################################################################################################


#######   ||||   ||      ||   TOT FIT: - opero fit totale con tutti i parametri
#######    ||     ||    ||               viscoelastici ma senza più i gaussiani
#######    ||      ||  ||
#######    ||       ||||   
#######   ||||       ||

print("\n\nI'm beginning total fit\n\n")
start = time.process_time()
tot_fit.Exec(matrix, boni, excluded, inputs, analysis_path)
tot_time    = time.process_time()-start
tempo       = tempo +(('fit totale',tot_time),)

with open(analysis_path+'log.txt', 'a') as f_log:
    f_log.write('\ntempo impiegato per fit totali: %f s\n'%(tot_time))
    f_log.write('\ntempo impiegato ore = %3.2f\n'%(tot_time/3600))

##### stampe finali

super_time    = time.process_time()-super_start
with open(analysis_path+'log.txt', 'a') as f_log:

    f_log.write('tempo impiegato per esecuzione dello script ore = %3.2f\n '%(super_time/3600))

    for (what,t) in tempo:

        f_log.write('di cui %f secondi =  %f  ore in %s \n' %(t, t/3600, what))

"""
Mass processing script to run everything.


#NOTE THAT THERE IS NO CALL TO PYTHON IN subprocess.run() BECAUSE ITS INSTALLED VIA SETUPTOOLS

"""

import os
import subprocess
import time
import sys

### need to fix this section later to handle plots, etc
#import odf_bottle_combine as btl_combine
#import odf_codes_plots as plots

def merge_files(file1, file2, file_out_path):
    with open(file_out_path, 'w') as file_out:
        with open(file1, 'r') as f1:
            for line in f1:
                file_out.write(line)
        with open(file2, 'r') as f2:
            for line in f2:
                file_out.write(line)
    return True

def process_all():
    prefix = 'nbp1802_'
    btl = '.btl'
    ros = '.ros'
    cnv = '.cnv'
    ext_hex = '.hex'
    xmlcon = '.XMLCON'
    csv = '.csv'
    #directories and index files
    ssscc_file = 'data/ssscc.csv'
    dir_logs = 'data/logs/'
    qual_dir_temp_primary = 'quality_code/temp_primary/'
    qual_dir_temp_secondary = 'quality_code/temp_secondary/'
    qual_dir_cond_primary = 'quality_code/cond_primary/'
    qual_dir_cond_secondary = 'quality_code/cond_secondary/'
    qual_flag_temp = 'quality_flag_temp'
    qual_flag_cond = 'quality_flag_cond'
    fit_t1 = 'fitting_t1'
    fit_t2 = 'fitting_t2'
    fit_c1 = 'fitting_c1'
    fit_c2 = 'fitting_c2'

    #load ssscc from file
    ssscc = []
    with open(ssscc_file, 'r') as filename:
        ssscc = [line.strip() for line in filename]

    ### TODO ADD IN REFT PROCESSOR AND MOVE FILES TO CORRECT DIRECTORY

    #check for already converted files to skip later
    time_start = time.perf_counter()
    cnv_dir_list = os.listdir('data/converted/')
    time_dir_list = os.listdir('data/time/')

    for x in ssscc:
        if '{}.pkl'.format(x) in cnv_dir_list:
            continue
        #convert hex to ctd
        subprocess.run(['odf_convert_sbe.py', 'data/raw/' + x + '.hex', 'data/raw/' + x + '.XMLCON', '-o', 'data/converted'], stdout=subprocess.PIPE)
        print('odf_convert_sbe.py SSSCC: ' + x + ' done')

    time_convert = time.perf_counter()

    for x in ssscc:
        if '{}_time.pkl'.format(x) in time_dir_list:
            continue
        subprocess.run(['odf_sbe_metadata.py', 'data/converted/' + x + '.pkl'], stdout=subprocess.PIPE)
        print('odf_sbe_metadata.py SSSCC: ' + x + ' done')

    btl_dir_list = os.listdir('data/bottle/')
    for x in ssscc:
        if '{}_btl_mean.pkl'.format(x) in btl_dir_list:
            continue
        #process bottle file
        subprocess.run(['odf_process_bottle.py', 'data/converted/' + x + '.pkl', '-o', 'data/bottle/'], stdout=subprocess.PIPE)
        print('odf_process_bottle.py SSSCC: ' + x + ' done')
    time_bottle = time.perf_counter()

    #btl_combine.main(None)
    ################################################################################

        #process pressure offset
    subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-press'], stdout=subprocess.PIPE)
    print('odf_calibrate_ctd.py pressure averaging: done')
    time_pressure_calibrate = time.perf_counter()

    for x in ssscc:
        #apply pressure offset to selected files
        subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-pres'], stdout=subprocess.PIPE)
        print('odf_fit_ctd.py pressure fit SSSCC: ' + x + ' done')
    time_pressure_fit = time.perf_counter()

    ################################################################################

    for x in range(2):
        #temperature fit against reftemp data
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'P', '-primary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}{qual_flag_temp}.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_pressure{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_pressure{csv}'], stdout=subprocess.PIPE)
        print('Primary temperature wrt P calibrated')


        #one pass on secondary sensors
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}2_pressure{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_pressure{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt P calibrated')

        #apply temperature fits to cast data (time)
        for x in ssscc:
            subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-temp'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py temp coefficients (pressure) appplied to SSSCC: ' + x + ' done')

        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'T', '-primary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_temperature{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_temperature_1{csv}'], stdout=subprocess.PIPE)
        print('Primary temperature wrt T calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'T', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}2_temperature{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_temperature_2{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt T calibrated')
        time_temperature_calibrate = time.perf_counter()

        time_temperature_calibrate = time.perf_counter()
        for x in ssscc:
            subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-temp'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py temp coefficients (temperature) appplied to SSSCC: ' + x + ' done')
        time_temperature_fit = time.perf_counter()

    ################################################################################

    for x in range(2):
        #conductivity fit against salt data
        #using 6000 because nominal bottom depth - change if not to not bias data

        #conductivity fit against salt data
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}_pressure_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_1{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')

        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}_pressure{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')

        #apply conductivity fits to cast data (time)
        for x in ssscc:
            subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-cond'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py cond coefficients (pressure) appplied to SSSCC: ' + x + ' done')
        time_conductivity_fit = time.perf_counter()

        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}_conductivity_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_1{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'C', '-secondary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}_conductivity{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')

        #apply conductivity fits to cast data (time)
        for x in ssscc:
            subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-cond'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py cond coefficients (conductivity) appplied to SSSCC: ' + x + ' done')
        time_conductivity_fit = time.perf_counter()

    ################################################################################

    #apply oxygen fits to cast data (time)
#    for x in ssscc:
#        subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-oxy', 'data/oxygen/' + x], stdout=subprocess.PIPE)
#        print('odf_fit_ctd.py oxy coefficients appplied to SSSCC: ' + x + ' done')

    subprocess.run(['oxy_fit_script.py'], stdout=subprocess.PIPE)
    time_oxygen_fit = time.perf_counter()

    subprocess.run(['ctd_to_bottle.py'], stdout=subprocess.PIPE)
    time_bottle_file = time.perf_counter()

    #plots.main(None)

    time_end = time.perf_counter()
    print(str(time_end - time_start) + ' total seconds elapsed via perf_counter(), ' + str((time_end - time_start)/60) + ' minutes elapsed')
    print(str(time_convert - time_start) + ' seconds elapsed converting raw files')
    print(str(time_bottle - time_convert) + ' seconds elapsed bottle crunching')
    print(str(time_pressure_calibrate - time_bottle) + ' seconds elapsed pressure calib')
    print(str(time_pressure_fit - time_pressure_calibrate) + ' seconds elapsed pressure fit')
    print(str(time_temperature_calibrate - time_pressure_fit) + ' seconds elapsed temperature calib')
    print(str(time_temperature_fit - time_temperature_calibrate) + ' seconds elapsed temperature fit')
    print(str(time_conductivity_calibrate - time_temperature_fit) + ' seconds elapsed conductivity calib')
    print(str(time_conductivity_fit - time_conductivity_calibrate) + ' seconds elapsed conductivity fit')
    print(str(time_oxygen_fit - time_conductivity_fit) + ' seconds elapsed oxygen fit')
    print(str(time_bottle_file - time_oxygen_fit) + ' seconds elapsed bottle file creation')

def main(argv):
    '''Run everything.
    '''
    process_all()

if __name__ == '__main__':
    main(sys.argv[1:])

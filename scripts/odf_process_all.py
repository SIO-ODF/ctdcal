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

def merge_files_v2(files, file_out_path):
    '''merge all files that are passed in as a list'''
    memory_buffer = f''
    for x in files:
        with open(x, 'r') as f_in:
            memory_buffer = memory_buffer + f_in.read()
            # for line in f_in:
            #     memory_buffer = memory_buffer + line
    with open(file_out_path,'w') as f_out:
        f_out.write(memory_buffer)
    return True

def fit_merging(side, param, max_cuts, file_out_path):
    path_logs = f'data/logs/quality_code/'
    c1_dir = f'cond_primary/'
    c2_dir = f'cond_secondary/'
    c1_file_cond = f'fitting_c1_conductivity_'
    c2_file_cond = f'fitting_c2_conductivity_'
    c1_file_pres = f'fitting_c1_pressure_'
    c2_file_pres = f'fitting_c2_pressure_'
    csv = f'.csv'

    c_dir = ''
    c_file = ''

    # determine which strings to use in paths
    if side == 1:
        c_dir = c1_dir
        if param == 'pressure':
            c_file = c1_file_pres
        if param == 'conductivity':
            c_file = c1_file_cond
    if side == 2:
        c_dir = c2_dir
        if param == 'pressure':
            c_file = c2_file_pres
        if param == 'conductivity':
            c_file = c2_file_cond

    #compile list of path names
    l = []
    for x in range(1, max_cuts + 1):
        l.append(f'{path_logs}{c_dir}{c_file}{x}{csv}')

    return merge_files_v2(l, file_out_path)

def code_merging(side, param, max_cuts, file_out_path):
    path_logs = f'data/logs/quality_code/'
    c1_dir = f'cond_primary/'
    c2_dir = f'cond_secondary/'
    c1_file_cond = f'quality_flag_cond_conductivity_'
    c2_file_cond = f'quality_flag_cond2_conductivity_'
    c1_file_pres = f'quality_flag_cond_pressure_'
    c2_file_pres = f'quality_flag_cond2_pressure_'
    csv = f'.csv'

    c_dir = ''
    c_file = ''

    # determine which strings to use in paths
    if side == 1:
        c_dir = c1_dir
        if param == 'pressure':
            c_file = c1_file_pres
        if param == 'conductivity':
            c_file = c1_file_cond
    if side == 2:
        c_dir = c2_dir
        if param == 'pressure':
            c_file = c2_file_pres
        if param == 'conductivity':
            c_file = c2_file_cond

    #compile list of path names
    l = []
    for x in range(1, max_cuts + 1):
        l.append(f'{path_logs}{c_dir}{c_file}{x}{csv}')

    return merge_files_v2(l, file_out_path)

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

    ssscc_t1 = f'data/ssscc/ssscc_t1.csv'
    ssscc_t2 = f'data/ssscc/ssscc_t2.csv'
    ssscc_t3 = f'data/ssscc/ssscc_t3.csv'
    ssscc_c1 = f'data/ssscc/ssscc_c1.csv'
    ssscc_c2 = f'data/ssscc/ssscc_c2.csv'
    ssscc_c3 = f'data/ssscc/ssscc_c3.csv'
    ssscc_c4 = f'data/ssscc/ssscc_c4.csv'
    ssscc_c5 = f'data/ssscc/ssscc_c5.csv'
    ssscc_c6 = f'data/ssscc/ssscc_c6.csv'
    ssscc_c7 = f'data/ssscc/ssscc_c7.csv'
    ssscc_c8 = f'data/ssscc/ssscc_c8.csv'



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

    #     #process pressure offset
    # subprocess.run(['odf_calibrate_ctd.py', ssscc_file, '-press'], stdout=subprocess.PIPE)
    # print('odf_calibrate_ctd.py pressure averaging: done')
    time_pressure_calibrate = time.perf_counter()
    #
    # for x in ssscc:
    #     #apply pressure offset to selected files
    #     subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-pres'], stdout=subprocess.PIPE)
    #     print('odf_fit_ctd.py pressure fit SSSCC: ' + x + ' done')
    time_pressure_fit = time.perf_counter()

    ################################################################################

    for x in range(2):
        #temperature fit against reftemp data
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t1, '-temp', '-calib', 'P', '-primary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}{qual_flag_temp}.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_pressure_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_pressure_1_1{csv}'], stdout=subprocess.PIPE)
        print('Primary temperature wrt P calibrated')
        # using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t2, '-temp', '-calib', 'P', '-primary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_pressure_2{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_pressure_1_2{csv}'], stdout=subprocess.PIPE)
        print('Primary temperature wrt P calibrated')
        # using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t3, '-temp', '-calib', 'P', '-primary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_pressure_3{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_pressure_1_3{csv}'], stdout=subprocess.PIPE)
        print('Primary temperature wrt P calibrated')

        merge_files(f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_pressure_1_1{csv}',f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_pressure_1_2{csv}', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_swap.csv')
        merge_files(f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_swap.csv',f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_pressure_1_3{csv}', f'{dir_logs}{fit_t1}.csv')


        #one pass on secondary sensors
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t1, '-temp', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}2_pressure_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_pressure_2_1{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt P calibrated')
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t2, '-temp', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}2_pressure_2{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_pressure_2_2{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t3, '-temp', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '800:5000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}2_pressure_3{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_pressure_2_3{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt P calibrated')

        merge_files(f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_pressure_2_1{csv}',f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_pressure_2_2{csv}', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}swap2.csv')
        merge_files(f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}swap2.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_pressure_2_3{csv}',f'{dir_logs}{fit_t2}.csv')

        #apply temperature fits to cast data (time)
        for x in ssscc:
            subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-temp'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py temp coefficients (pressure) appplied to SSSCC: ' + x + ' done')

        # subprocess.run(['odf_calibrate_ctd.py', ssscc_t1, '-temp', '-calib', 'T', '-primary', '-order', '1'], stdout=subprocess.PIPE)
        # subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_temperature_1{csv}'], stdout=subprocess.PIPE)
        # subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_temperature_1_1{csv}'], stdout=subprocess.PIPE)
        # print('Primary temperature wrt T calibrated')
        # subprocess.run(['odf_calibrate_ctd.py', ssscc_t2, '-temp', '-calib', 'T', '-primary', '-order', '1'], stdout=subprocess.PIPE)
        # subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_temperature_2{csv}'], stdout=subprocess.PIPE)
        # subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_temperature_1_2{csv}'], stdout=subprocess.PIPE)
        # print('Primary temperature wrt T calibrated')
        # subprocess.run(['odf_calibrate_ctd.py', ssscc_t3, '-temp', '-calib', 'T', '-primary', '-order', '1'], stdout=subprocess.PIPE)
        # subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_primary}{qual_flag_temp}1_temperature_3{csv}'], stdout=subprocess.PIPE)
        # subprocess.run(['cp', f'{dir_logs}{fit_t1}.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_temperature_1_3{csv}'], stdout=subprocess.PIPE)
        # print('Primary temperature wrt T calibrated')
        #
        # merge_files(f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_temperature_1_1{csv}', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_temperature_1_2{csv}',f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_swap.csv')
        # merge_files(f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_swap.csv', f'{dir_logs}{qual_dir_temp_primary}{fit_t1}_temperature_1_3{csv}',f'{dir_logs}{fit_t1}.csv')
        #
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t1, '-temp', '-calib', 'T', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}_temperature_2_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_temperature_2_1{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt T calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t2, '-temp', '-calib', 'T', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}_temperature_2_2{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_temperature_2_2{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt T calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_t3, '-temp', '-calib', 'T', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_temp.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}_temperature_2_3{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_t2}.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_temperature_2_3{csv}'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt T calibrated')

        merge_files(f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_temperature_2_1{csv}', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_temperature_2_2{csv}',f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}swap2.csv')
        merge_files(f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}swap2.csv', f'{dir_logs}{qual_dir_temp_secondary}{fit_t2}_temperature_2_3{csv}',f'{dir_logs}{fit_t2}.csv')

        merge_files(f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}_temperature_2_1{csv}', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}_temperature_2_2{csv}',f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}swap2.csv')
        merge_files(f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}swap2.csv', f'{dir_logs}{qual_dir_temp_secondary}{qual_flag_temp}_temperature_2_3{csv}',f'{dir_logs}{qual_flag_temp}.csv')


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
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c1, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_1{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c2, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_2{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_2{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c3, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_3{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_3{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c4, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_4{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_4{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c5, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_5{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_5{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c6, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_6{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_6{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c7, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_7{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_7{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c8, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_pressure_8{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_pressure_8{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')
        fit_merging(1, 'pressure', 8, f'{dir_logs}{fit_c1}.csv')

        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c1, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_1{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c2, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_2{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_2{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c3, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_3{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_3{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c4, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_4{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_4{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        #using 6000 because nominal bottom depth - change if not to not bias data
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c5, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_5{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_5{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c6, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_6{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_6{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c7, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_7{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_7{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c8, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:5000'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_pressure_8{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_pressure_8{csv}'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')
        fit_merging(2, 'pressure', 8, f'{dir_logs}{fit_c2}.csv')

        #apply conductivity fits to cast data (time)
        for x in ssscc:
            subprocess.run(['odf_fit_ctd.py', 'data/time/' + x + '_time.pkl', '-cond'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py cond coefficients (pressure) appplied to SSSCC: ' + x + ' done')
        time_conductivity_fit = time.perf_counter()

        subprocess.run(['odf_calibrate_ctd.py', ssscc_c1, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_1{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c2, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_2{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_2{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c3, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_3{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_3{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c4, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_4{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_4{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c5, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_5{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_5{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c6, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_6{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_6{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c7, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_7{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_7{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c8, '-cond', '-calib', 'C', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_primary}{qual_flag_cond}1_conductivity_8{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c1}.csv', f'{dir_logs}{qual_dir_cond_primary}{fit_c1}_conductivity_8{csv}'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt C calibrated')

        fit_merging(1, 'conductivity', 8, f'{dir_logs}{fit_c1}.csv')

        subprocess.run(['odf_calibrate_ctd.py', ssscc_c1, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_1{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_1{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c2, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_2{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_2{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c3, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_3{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_3{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c4, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_4{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_4{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c5, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_5{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_5{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c6, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_6{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_6{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c7, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_7{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_7{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')
        subprocess.run(['odf_calibrate_ctd.py', ssscc_c8, '-cond', '-calib', 'C', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}quality_flag_cond.csv', f'{dir_logs}{qual_dir_cond_secondary}{qual_flag_cond}2_conductivity_8{csv}'], stdout=subprocess.PIPE)
        subprocess.run(['cp', f'{dir_logs}{fit_c2}.csv', f'{dir_logs}{qual_dir_cond_secondary}{fit_c2}_conductivity_8{csv}'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt C calibrated')

        fit_merging(2, 'conductivity', 8, f'{dir_logs}{fit_c2}.csv')
        code_merging(2, 'conductivity', 8, f'{dir_logs}quality_flag_cond.csv')

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

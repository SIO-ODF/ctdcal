"""
Mass processing script to run everything.

"""

import os
import subprocess
import time
import sys

def process_all():
    prefix = 'nbp1706_'
    btl = '.btl'
    ros = '.ros'
    cnv = '.cnv'
    ext_hex = '.hex'
    xmlcon = '.XMLCON'
    ssscc_file = '../ssscc.csv'

    #load ssscc from file
    ssscc = []
    with open(ssscc_file, 'r') as filename:
        ssscc = [line.strip() for line in filename]

    #check for already converted files to skip later
    time_start = time.perf_counter()
    cnv_dir_list = os.listdir('data/converted/')

    for x in ssscc:
        if '{}_cnv.pkl'.format(x) in cnv_dir_list:
            continue
        #convert hex to ctd
        subprocess.run(['python3', './odf_convert_sbe.py', 'data/raw/' + x + '.hex', 'data/raw/' + x + '.XMLCON', '-o', 'data/converted'], stdout=subprocess.PIPE)
        print('odf_convert_sbe.py SSSCC: ' + x + ' done')

    time_convert = time.perf_counter()

    btl_dir_list = os.listdir('data/bottle/')
    for x in ssscc:
        if '{}_btl.pkl'.format(x) in btl_dir_list:
            continue
        #process bottle file
        subprocess.run(['python3', './odf_process_bottle.py', 'data/converted/' + x + '_cnv.pkl', '-o', 'data/bottle/'], stdout=subprocess.PIPE)
        print('odf_process_bottle.py SSSCC: ' + x + ' done')
    time_bottle = time.perf_counter()

    ################################################################################

        #process pressure offset
    subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-press'], stdout=subprocess.PIPE)
    print('odf_calibrate_ctd.py pressure averaging: done')
    time_pressure_calibrate = time.perf_counter()

    for x in ssscc:
        #apply pressure offset to selected files
        subprocess.run(['python3', './odf_fit_ctd.py', 'data/time/' + x + '_time.csv', '-pres'], stdout=subprocess.PIPE)
        print('odf_fit_ctd.py pressure fit SSSCC: ' + x + ' done')
    time_pressure_fit = time.perf_counter()

    ################################################################################

    for x in range(2):
        #temperature fit against reftemp data
        #using 5000 because no casts went past 5000
        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'P', '-primary', '-order', '2', '-xRange', '1000:6000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()
        print('Primary temperature wrt P calibrated')

        #one pass on secondary sensors
        #using 5000 because no casts went past 5000
        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '1000:6000'], stdout=subprocess.PIPE)
        time_temperature_calibrate = time.perf_counter()

        print('Secondary temperature wrt P calibrated')

        #apply temperature fits to cast data (time)
        for x in ssscc:
            subprocess.run(['python3', './odf_fit_ctd.py', 'data/time/' + x + '_time.csv', '-temp'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py temp coefficients appplied to SSSCC: ' + x + ' done')

        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'T', '-primary', '-order', '1'], stdout=subprocess.PIPE)
        print('Primary temperature wrt T calibrated')
        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-temp', '-calib', 'T', '-secondary', '-order', '1'], stdout=subprocess.PIPE)
        print('Secondary temperature wrt T calibrated')

        time_temperature_calibrate = time.perf_counter()
        for x in ssscc:
            subprocess.run(['python3', './odf_fit_ctd.py', 'data/time/' + x + '_time.csv', '-temp'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py temp coefficients appplied to SSSCC: ' + x + ' done')
        time_temperature_fit = time.perf_counter()

    ################################################################################

    for x in range(2):
        #conductivity fit against salt data
        #using 5000 because no casts went past 5000
        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'P', '-primary', '-order', '2', '-xRange', '200:6000'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Primary conductivity wrt P calibrated')

        #using 5000 because no casts went past 5000
        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'P', '-secondary', '-order', '2', '-xRange', '200:6000'], stdout=subprocess.PIPE)
        time_conductivity_calibrate = time.perf_counter()
        print('Secondary conductivity wrt P calibrated')

        #apply conductivity fits to cast data (time)
        for x in ssscc:
            subprocess.run(['python3', './odf_fit_ctd.py', 'data/time/' + x + '_time.csv', '-cond'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py cond coefficients appplied to SSSCC: ' + x + ' done')
        time_conductivity_fit = time.perf_counter()

        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'T', '-primary', '-order', '2'], stdout=subprocess.PIPE)
        print('Primary conductivity wrt T calibrated')
        subprocess.run(['python3', './odf_calibrate_ctd.py', ssscc_file, '-cond', '-calib', 'T', '-secondary', '-order', '2'], stdout=subprocess.PIPE)
        print('Secondary conductivity wrt T calibrated')

        #apply conductivity fits to cast data (time)
        for x in ssscc:
            subprocess.run(['python3', './odf_fit_ctd.py', 'data/time/' + x + '_time.csv', '-cond'], stdout=subprocess.PIPE)
            print('odf_fit_ctd.py cond coefficients appplied to SSSCC: ' + x + ' done')
        time_conductivity_fit = time.perf_counter()

    ################################################################################

    #apply oxygen fits to cast data (time)
    for x in ssscc:
        subprocess.run(['python3', './odf_fit_ctd.py', 'data/time/' + x + '_time.csv', '-oxy', 'data/oxygen/' + x], stdout=subprocess.PIPE)
        print('odf_fit_ctd.py oxy coefficients appplied to SSSCC: ' + x + ' done')
    time_oxygen_fit = time.perf_counter()

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

def main(argv):
    '''Run everything.
    '''
    process_all()

if __name__ == '__main__':
    main(sys.argv[1:])

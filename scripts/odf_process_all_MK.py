'''
Attempt to write a cleaner processing script from scratch.
'''

# import necessary packages
import os
import time
import sys
import subprocess
import pandas as pd
import numpy as np
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import ctdcal.rinko as rinko
import odf_salt_parser as salt_parser
import odf_reft_parser as reft_parser

def process_all():

    #####
    # Step 0: Define necessary variables
    #####

    # define cruise and file information/extensions
    prefix = 'nbp1802_'
    btl = '.btl'
    ros = '.ros'
    cnv = '.cnv'
    ext_hex = '.hex'
    xmlcon = '.XMLCON'
    csv = '.csv'

    # directories and index files
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

    # segmented station/cast files
    ssscc_t1 = f'data/ssscc/ssscc_t1.csv'
    ssscc_t2 = f'data/ssscc/ssscc_t2.csv'
    ssscc_c1 = f'data/ssscc/ssscc_c1.csv'
    ssscc_c2 = f'data/ssscc/ssscc_c2.csv'
    ssscc_c3 = f'data/ssscc/ssscc_c3.csv'
    ssscc_c4 = f'data/ssscc/ssscc_c4.csv'
    ssscc_c5 = f'data/ssscc/ssscc_c5.csv'
    ssscc_c6 = f'data/ssscc/ssscc_c6.csv'

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    ssscc_list = []
    with open(ssscc_file, 'r') as filename:
        ssscc_list = [line.strip() for line in filename]

    # make list of already converted files to skip later
    cnv_dir_list = os.listdir('data/converted/')
    time_dir_list = os.listdir('data/time/')
    btl_dir_list = os.listdir('data/bottle/')
    salt_dir_list = os.listdir('data/salt/')
    reft_dir_list = os.listdir('data/reft/')

    # convert hex to ctd (TODO: convert this to function form)
    for ssscc in ssscc_list:
        if '{}.pkl'.format(ssscc) not in cnv_dir_list:
            subprocess.run(['odf_convert_sbe.py', 'data/raw/' + ssscc + '.hex', 'data/raw/' + ssscc + '.XMLCON', '-o', 'data/converted'], stdout=subprocess.PIPE)
            print('odf_convert_sbe.py SSSCC: ' + ssscc + ' done')

    # ??? (TODO: convert this to function form)
    for ssscc in ssscc_list:
        if '{}_time.pkl'.format(ssscc) not in time_dir_list:
            subprocess.run(['odf_sbe_metadata.py', 'data/converted/' + ssscc + '.pkl'], stdout=subprocess.PIPE)
            print('odf_sbe_metadata.py SSSCC: ' + ssscc + ' done')

    # process bottle file (TODO: convert this to function form)
    # does this generate "ondeck_pressure.csv"?
    for ssscc in ssscc_list:
        if '{}_btl_mean.pkl'.format(ssscc) not in btl_dir_list:
            subprocess.run(['odf_process_bottle.py', 'data/converted/' + ssscc + '.pkl', '-o', 'data/bottle/'], stdout=subprocess.PIPE)
            print('odf_process_bottle.py SSSCC: ' + ssscc + ' done')

    # generate salt files
    for ssscc in ssscc_list:
        if '{}_salts.csv'.format(ssscc) not in salt_dir_list:
            salt_parser.process_salts(ssscc,'data/salt/')

    # generate ref temp files
    for ssscc in ssscc_list:
        if '{}_reft.csv'.format(ssscc) not in reft_dir_list:
            reft_parser.process_reft(ssscc,'data/reft/')

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    btl_cols = ['index', 'CTDTMP1', 'CTDTMP2', 'CTDPRS', 'CTDCOND1', 'CTDCOND2',
       'CTDSAL', 'CTDOXY1', 'CTDOXYVOLTS', 'CTDXMISS', 'ALT', 'REF_PAR',
       'GPSLAT', 'GPSLON','new_fix', 'pressure_temp_int', 'pump_on', 'btl_fire', 
       'scan_datetime', 'btl_fire_num']
    btl_data_all = process_ctd.load_all_ctd_files(ssscc_list, 'bottle', btl_cols)
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list, 'time', None)

    # process pressure offset
    pressure_log = process_ctd.load_pressure_logs('data/logs/ondeck_pressure.csv')
    p_offset = process_ctd.get_pressure_offset(pressure_log.ondeck_start_p, 
                                                pressure_log.ondeck_end_p)

    btl_data_all = fit_ctd.apply_pressure_offset(btl_data_all, 'CTDPRS', p_offset)
    time_data_all = fit_ctd.apply_pressure_offset(time_data_all, 'CTDPRS', p_offset)


def main(argv):
    '''Run everything.
    '''
    process_all()

if __name__ == '__main__':
    main(sys.argv[1:])

'''
Attempt to write a cleaner processing script from scratch.
'''

# import necessary packages
import os
import glob
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
    # Step 0: Load and define necessary variables
    #####

    import config

    # define cruise and file information/extensions
    prefix = 'nbp1802_'

    # directories and index files
    qual_flag_temp = 'quality_flag_temp'
    qual_flag_cond = 'quality_flag_cond'
    fit_t1 = 'fitting_t1'
    fit_t2 = 'fitting_t2'
    fit_c1 = 'fitting_c1'
    fit_c2 = 'fitting_c2'

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    ssscc_list = []
    with open(config.directory['ssscc_file'], 'r') as filename:
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

    btl_data_all = fit_ctd.apply_pressure_offset(btl_data_all, config.column['p'], p_offset)
    time_data_all = fit_ctd.apply_pressure_offset(time_data_all, config.column['p'], p_offset)

    # temperature calibration
    # steps: 1) remove non-finite data
    #        2) flag points w/ large deviations
    #        3) calculate fit parameters (on data w/ flag 2) -> save them too!
    #        4) apply fit
    #        5) qualify flag file

    ssscc_files_t = sorted(glob.glob('data/ssscc/ssscc_*t*.csv'))
    qual_flag_t1 = pd.DataFrame()
    qual_flag_t2 = pd.DataFrame()

    for f in ssscc_files_t:
        # 0) grab ssscc chunk to fit
        ssscc_list_t = pd.read_csv(f, header=None, dtype='str', squeeze=True).to_list()

        # 1) remove non-finite data
        df_temp_good = process_ctd.prepare_fit_data(btl_data_all[btl_data_all['SSSCC'].isin(ssscc_list_t)], config.column['reft'])

        # 2 & 3) calculate fit params 
        coef_temp_prim,df_ques_t1 = process_ctd.calibrate_param(df_temp_good[config.column['t1_btl']], df_temp_good[config.column['reft']],
                                                                df_temp_good[config.column['p_btl']], 'T', 1, df_temp_good['SSSCC'],
                                                                df_temp_good['btl_fire_num'], xRange='1000:5000')
        coef_temp_sec,df_ques_t2 = process_ctd.calibrate_param(df_temp_good[config.column['t2_btl']], df_temp_good[config.column['reft']],
                                                                df_temp_good[config.column['p_btl']], 'T', 1, df_temp_good['SSSCC'],
                                                                df_temp_good['btl_fire_num'], xRange='1000:5000')
                                                                
        # 4) apply fit
        btl_data_all.loc[btl_data_all['SSSCC'].isin(ssscc_list_t).values][config.column['t1_btl']] = fit_ctd.temperature_polyfit(btl_data_all.loc[btl_data_all['SSSCC'].isin(ssscc_list_t)][config.column['t1_btl']],
                                                                                                                    btl_data_all.loc[btl_data_all['SSSCC'].isin(ssscc_list_t)][config.column['p_btl']], coef_temp_prim)
        btl_data_all.loc[btl_data_all['SSSCC'].isin(ssscc_list_t).values][config.column['t2_btl']] = fit_ctd.temperature_polyfit(btl_data_all[btl_data_all['SSSCC'].isin(ssscc_list_t)][config.column['t2_btl']],
                                                                                                                    btl_data_all.loc[btl_data_all['SSSCC'].isin(ssscc_list_t)][config.column['p_btl']], coef_temp_sec)
        time_data_all.loc[time_data_all['SSSCC'].isin(ssscc_list_t).values][config.column['t1']] = fit_ctd.temperature_polyfit(time_data_all.loc[time_data_all['SSSCC'].isin(ssscc_list_t)][config.column['t1']],
                                                                                                                    time_data_all.loc[time_data_all['SSSCC'].isin(ssscc_list_t)][config.column['p']], coef_temp_prim)
        time_data_all.loc[time_data_all['SSSCC'].isin(ssscc_list_t).values][config.column['t2']] = fit_ctd.temperature_polyfit(time_data_all[time_data_all['SSSCC'].isin(ssscc_list_t)][config.column['t2']],
                                                                                                                    time_data_all.loc[time_data_all['SSSCC'].isin(ssscc_list_t)][config.column['p']], coef_temp_sec)                                                                                                            

        # 5) handle quality flags
        qual_flag_t1 = pd.concat([qual_flag_t1,df_ques_t1])
        qual_flag_t2 = pd.concat([qual_flag_t2,df_ques_t2])

    # export temp quality flags
    qual_flag_t1.to_csv(config.directory['logs'] + 'qual_flag_t1.csv', index=False)
    qual_flag_t2.to_csv(config.directory['logs'] + 'qual_flag_t2.csv', index=False)


def main(argv):
    '''Run everything.
    '''
    process_all()

if __name__ == '__main__':
    main(sys.argv[1:])

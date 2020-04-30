#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 03:55:51 2018

@author: k3jackson
"""

import os
import subprocess
import time
import sys
import configparser
#sys.path.append('ctdcal/')
import pandas as pd
import numpy as np
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import ctdcal.rinko as rinko

# MK: deprecated 04/29/20, use odf_process_all_MK.py instead

def process_all_new():

    ssscc_file = 'data/ssscc.csv'
    iniFile = 'data/ini-files/configuration.ini'
    config = configparser.RawConfigParser()
    config.read(iniFile)

    # maybe reformat config file?
    # fold this code into config module
    p_col = config['analytical_inputs']['p']
    p_btl_col = config['inputs']['p']
    t1_col = config['analytical_inputs']['t1']
    t1_btl_col = config['inputs']['t1']
    t2_col = config['analytical_inputs']['t2']
    t2_btl_col = config['inputs']['t2']
    c1_col = config['analytical_inputs']['c1']
    c1_btl_col = config['inputs']['c1']
    c2_col = config['analytical_inputs']['c2']
    c2_btl_col = config['inputs']['c2']
    sal_col = config['analytical_inputs']['salt']
    reft_col = config['inputs']['reft']
    sample_rate = config['ctd_processing']['sample_rate']
    search_time = config['ctd_processing']['roll_filter_time']

    expocode = config['cruise']['expocode']
    sectionID = config['cruise']['sectionid']
    raw_directory = config['ctd_processing']['raw_data_directory']
    time_directory = config['ctd_processing']['time_data_directory']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    oxygen_directory = config['ctd_processing']['oxygen_directory']
    btl_directory = config['ctd_processing']['bottle_directory']
    o2flask_file = config['ctd_processing']['o2flask_file']
    log_directory = config['ctd_processing']['log_directory']
    sample_rate = config['ctd_processing']['sample_rate']
    search_time = config['ctd_processing']['roll_filter_time']
    ctd = config['ctd_processing']['ctd_serial']

    # MK: import from settings/config instead of hard code?
    btl_data_prefix = 'data/bottle/'
    btl_data_postfix = '_btl_mean.pkl'
    time_data_prefix = 'data/time/'
    time_data_postfix = '_time.pkl'
    p_log_file = 'data/logs/ondeck_pressure.csv'
    t1_btl_col = 'CTDTMP1'
    t2_btl_col = 'CTDTMP2'
    t1_col = 'CTDTMP1'
    t2_col = 'CTDTMP2'
    p_btl_col = 'CTDPRS'
    p_col = 'CTDPRS'
    c1_btl_col = 'CTDCOND1'
    c2_btl_col = 'CTDCOND2'
    c1_col = 'CTDCOND1'
    c2_col = 'CTDCOND2'
    reft_col = 'T90'
    refc_col = 'BTLCOND'
    sal_col = 'CTDSAL'

    #time_column_data = config['time_series_output']['data_names'].split(',')
    time_column_data = config['time_series_output']['data_output']
    time_column_names = config['time_series_output']['column_name'].split(',')
    time_column_units = config['time_series_output']['column_units'].split(',')
    time_column_format = config['time_series_output']['format']

    #pressure_column_data = config['time_series_output']['data_names'].split(',')
    p_column_data = config['pressure_series_output']['data'].split(',')
    p_column_names = config['pressure_series_output']['column_name'].split(',')
    p_column_units = config['pressure_series_output']['column_units'].split(',')
    p_column_format = config['pressure_series_output']['format']
    p_column_qual = config['pressure_series_output']['qual_columns'].split(',')
    p_column_one = list(config['pressure_series_output']['q1_columns'].split(','))

    # Columns from btl file to be read:
    cols = ['index', 'CTDTMP1', 'CTDTMP2', 'CTDPRS', 'CTDCOND1', 'CTDCOND2',
       'CTDSAL', 'CTDOXY1', 'CTDOXYVOLTS', 'CTDXMISS', 'ALT', 'REF_PAR', 'GPSLAT', 'GPSLON','new_fix', 'pressure_temp_int', 'pump_on', 'btl_fire', 'scan_datetime',
       'btl_fire_num']


    # Load ssscc from file
    ssscc = []
    with open(ssscc_file, 'r') as filename:
        ssscc = [line.strip() for line in filename]

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

    # generate salts file here: (salt parser)
    salt_dir = './data/salt/'
    import scripts.odf_salt_parser as salt_parser

    file_list = os.listdir(salt_dir)
    files = []
    for file in file_list:
        if '.' not in file:
            files.append(file)

    for file in files:
        print(file)
        salt_path = salt_dir + file
        saltDF = salt_parser.salt_loader(saltpath=salt_path)
        salt_parser.salt_df_parser(saltDF, salt_dir)

    # generate reft file here
    reft_path = './data/reft/'
    import scripts.odf_reft_parser as reft_parser
    for station in ssscc:
        with open(reft_path + station + '.cap', 'r') as ssscc_reftemp:
            # create csv files
            df_part = reft_parser.parse(ssscc_reftemp, station)

        df_part.to_csv(reft_path + station + '_reft.csv',index=False)

    # generate oxy file here (separate out from calib)

    ###########################################################################

    ### File I/O

    # Load in all bottle, time, ref_data files into DataFrame

    btl_data_all = process_ctd.load_all_ctd_files(ssscc,btl_data_prefix,
                                                  btl_data_postfix,'bottle',cols)
    time_data_all = process_ctd.load_all_ctd_files(ssscc,time_data_prefix,
                                                   time_data_postfix,'time',None)
                                                   ##
   # end of data collection, begin calibration
    ################################################################################

    ### Pressure Calibration

    # Determine Pressure offset from logs

    pressure_log = process_ctd.load_pressure_logs(p_log_file)
    p_off = process_ctd.get_pressure_offset(pressure_log.ondeck_start_p,pressure_log.ondeck_end_p) # insert logic to check if good value?

    btl_data_all = fit_ctd.apply_pressure_offset(btl_data_all,p_btl_col, p_off) # not actually bottle data
    time_data_all = fit_ctd.apply_pressure_offset(time_data_all,p_col,p_off) # continuous data
    # btl_data_all[p_btl_col] = fit_ctd.apply_pressure_offset(btl_data_all[p_btl_col], p_off)
    # time_data_all[p_col] = fit_ctd.apply_pressure_offset(time_data_all[p_col],p_off)


    ###########################################################################

    ### Temperature Calibration
    for x in range(2): # 2 passes checks for outliers twice (i.e. before/after fit)

     # Second order calibration

        df_temp_good = process_ctd.prepare_fit_data(btl_data_all, 'T90')

        df_ques_reft = process_ctd.quality_check(df_temp_good[t2_btl_col], df_temp_good[t1_btl_col], df_temp_good[p_btl_col], df_temp_good['SSSCC'], df_temp_good['btl_fire_num'], 'quest')
        df_ques_reft['Parameter'] = 'REF_TEMP'

        # sort out fit parameters, what order P, what order T
        # only fit `linear' regions, i.e. not mixed/benthic layers
        # built better way to handle xRange?
        coef_temp_1,df_ques_t1 = process_ctd.calibrate_param(df_temp_good[t1_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'TP', 2, df_temp_good.SSSCC, df_temp_good.btl_fire_num, xRange='800:6000')
        coef_temp_2,df_ques_t2 = process_ctd.calibrate_param(df_temp_good[t2_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'TP', 2, df_temp_good.SSSCC, df_temp_good.btl_fire_num, xRange='1500:6000')

    # Apply fitting coef to data

        btl_data_all[t1_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t1_btl_col], btl_data_all[p_btl_col], coef_temp_1)
        btl_data_all[t2_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t2_btl_col], btl_data_all[p_btl_col], coef_temp_2)


        time_data_all[t1_col] = fit_ctd.temperature_polyfit(time_data_all[t1_col], time_data_all[p_col], coef_temp_1)
        time_data_all[t2_col] = fit_ctd.temperature_polyfit(time_data_all[t2_col], time_data_all[p_col], coef_temp_2)


    # Construct Quality Flag file

        qual_flag_temp = process_ctd.combine_quality_flags([df_ques_reft,df_ques_t1,df_ques_t2])

    ## First order calibtation

        df_temp_good = process_ctd.prepare_fit_data(btl_data_all, 'T90')

#        df_ques_reft = process_ctd.quality_check(df_temp_good[t2_btl_col], df_temp_good[t1_btl_col], df_temp_good[p_btl_col], df_temp_good['SSSCC'], df_temp_good['btl_fire_num'], 'quest')
#        df_ques_reft['Parameter'] = 'REF_TEMP'

        coef_temp_prim,df_ques_t1 = process_ctd.calibrate_param(df_temp_good[t1_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'T', 1, df_temp_good.SSSCC, df_temp_good.btl_fire_num)
        coef_temp_sec,df_ques_t2 = process_ctd.calibrate_param(df_temp_good[t2_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'T', 1, df_temp_good.SSSCC, df_temp_good.btl_fire_num)

        btl_data_all[t1_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t1_btl_col], btl_data_all[p_btl_col], coef_temp_prim)
        btl_data_all[t2_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t2_btl_col], btl_data_all[p_btl_col], coef_temp_sec)


    # Apply fitting coef to data

        btl_data_all[t1_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t1_btl_col], btl_data_all[p_btl_col], coef_temp_prim)
        btl_data_all[t2_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t2_btl_col], btl_data_all[p_btl_col], coef_temp_sec)

        time_data_all[t1_col] = fit_ctd.temperature_polyfit(time_data_all[t1_col], time_data_all[p_col], coef_temp_prim)
        time_data_all[t2_col] = fit_ctd.temperature_polyfit(time_data_all[t2_col], time_data_all[p_col], coef_temp_sec)

        qual_flag_temp = process_ctd.combine_quality_flags([df_ques_reft,df_ques_t1,df_ques_t2])

    ###########################################################################

#
    ### Conductivity Calibration
    for x in range(2):
    # Update ref conductivity from salinometer

        btl_data_all['BTLCOND'] = fit_ctd.CR_to_cond(btl_data_all.CRavg,btl_data_all.BathTEMP,btl_data_all.CTDTMP1,btl_data_all.CTDPRS)
        df_cond_good = process_ctd.prepare_fit_data(btl_data_all,'BTLCOND')

        df_ques_refc = process_ctd.quality_check(df_cond_good['CTDCOND2'], df_temp_good['CTDCOND1'],df_temp_good['CTDPRS'], df_temp_good['SSSCC'], df_temp_good['btl_fire_num'], 'quest')
        df_ques_refc['Parameter'] = 'REF_COND'

    # Second Order Calibration

        # can probably use same xRange as temp, unless C sensor has some quirk
        coef_cond_1,df_ques_c1 = process_ctd.calibrate_param(df_cond_good.CTDCOND1,df_cond_good.BTLCOND,df_cond_good.CTDPRS,'CP',2,df_cond_good.SSSCC,df_cond_good.btl_fire_num,xRange='800:6000')
        coef_cond_2,df_ques_c2 = process_ctd.calibrate_param(df_cond_good.CTDCOND2,df_cond_good.BTLCOND,df_cond_good.CTDPRS,'CP',2,df_cond_good.SSSCC,df_cond_good.btl_fire_num,xRange='1500:6000')


    # Apply fitting coef to data

        # fit_ctd.conductivity_polyfit has salinity commented out? (fit_ctd.py, ln 460) MK
        btl_data_all['CTDCOND1'] = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND1'],btl_data_all['CTDTMP1'],btl_data_all['CTDPRS'],coef_cond_1)
        btl_data_all['CTDCOND2'] = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND2'],btl_data_all['CTDTMP2'],btl_data_all['CTDPRS'],coef_cond_2)
        # btl_data_all['CTDCOND1'],btl_data_all['CTDSAL'] = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND1'],btl_data_all['CTDTMP1'],btl_data_all['CTDPRS'],coef_cond_1)
        # btl_data_all['CTDCOND2'],sal_2 = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND2'],btl_data_all['CTDTMP2'],btl_data_all['CTDPRS'],coef_cond_2)

        # fit_ctd.conductivity_polyfit has salinity commented out? (fit_ctd.py, ln 460) MK
        time_data_all['CTDCOND1'] = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND1'], time_data_all['CTDTMP1'], time_data_all['CTDPRS'], coef_cond_1)
        time_data_all['CTDCOND2'] = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND2'], time_data_all['CTDTMP2'], time_data_all['CTDPRS'], coef_cond_2)
        # time_data_all['CTDCOND1'],time_data_all['CTDSAL'] = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND1'], time_data_all['CTDTMP1'], time_data_all['CTDPRS'], coef_cond_1)
        # time_data_all['CTDCOND2'],sal2 = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND2'], time_data_all['CTDTMP2'], time_data_all['CTDPRS'], coef_cond_2)


        qual_flag_cond = process_ctd.combine_quality_flags([df_ques_c1,df_ques_c2,df_ques_refc])

    # Finally WRT to Cond

        btl_data_all['BTLCOND'] = fit_ctd.CR_to_cond(btl_data_all.CRavg,btl_data_all.BathTEMP,btl_data_all.CTDTMP1,btl_data_all.CTDPRS)
        df_cond_good = process_ctd.prepare_fit_data(btl_data_all,'BTLCOND')

        coef_cond_prim,df_ques_c1 = process_ctd.calibrate_param(df_cond_good.CTDCOND1,df_cond_good.BTLCOND,df_cond_good.CTDPRS,'C',2,df_cond_good.SSSCC,df_cond_good.btl_fire_num)
        coef_cond_sec,df_ques_c2 = process_ctd.calibrate_param(df_cond_good.CTDCOND2,df_cond_good.BTLCOND,df_cond_good.CTDPRS,'C',2,df_cond_good.SSSCC,df_cond_good.btl_fire_num)

    # Apply fitting coef to data

        # fit_ctd.conductivity_polyfit has salinity commented out? (fit_ctd.py, ln 460) MK
        btl_data_all['CTDCOND1'] = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND1'],btl_data_all['CTDTMP1'],btl_data_all['CTDPRS'],coef_cond_prim)
        btl_data_all['CTDCOND2'] = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND2'],btl_data_all['CTDTMP2'],btl_data_all['CTDPRS'],coef_cond_sec)
        # btl_data_all['CTDCOND1'],sal_1 = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND1'],btl_data_all['CTDTMP1'],btl_data_all['CTDPRS'],coef_cond_prim)
        # btl_data_all['CTDCOND2'],sal_2 = fit_ctd.conductivity_polyfit(btl_data_all['CTDCOND2'],btl_data_all['CTDTMP2'],btl_data_all['CTDPRS'],coef_cond_sec)

        # fit_ctd.conductivity_polyfit has salinity commented out? (fit_ctd.py, ln 460) MK
        time_data_all['CTDCOND1'] = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND1'], time_data_all['CTDTMP1'], time_data_all['CTDPRS'], coef_cond_prim)
        time_data_all['CTDCOND2'] = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND2'], time_data_all['CTDTMP2'], time_data_all['CTDPRS'], coef_cond_sec)
        # time_data_all['CTDCOND1'],time_data_all['CTDSAL'] = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND1'], time_data_all['CTDTMP1'], time_data_all['CTDPRS'], coef_cond_prim)
        # time_data_all['CTDCOND2'],sal2 = fit_ctd.conductivity_polyfit(time_data_all['CTDCOND2'], time_data_all['CTDTMP2'], time_data_all['CTDPRS'], coef_cond_sec)

        qual_flag_cond = process_ctd.combine_quality_flags([df_ques_c1,df_ques_c2,df_ques_refc])
    ###########################################################################
#
#    ## Oxygen Calibration
#    btl_data_all,time_data_all = oxy_fitting.calibrate_oxygen(time_data_all,btl_data_all,ssscc)

#    subprocess.run(['oxy_fit_script.py'], stdout=subprocess.PIPE)
#
#
#    subprocess.run(['ctd_to_bottle.py'], stdout=subprocess.PIPE)
#
#
#   ##########################################################################

    #####
    # Oxygen calibration (MK)
    # look over this code again, maybe refence matlab scripts from PMEL

    import ctdcal.oxy_fitting as oxy_fitting
    import gsw

    # definitions
    import ctdcal.settings as settings
    hex_prefix = settings.ctd_processing_dir['hex_prefix']
    hex_postfix = settings.ctd_processing_dir['hex_postfix']
    xml_prefix = settings.ctd_processing_dir['xml_prefix']
    xml_postfix = settings.ctd_processing_dir['xml_postfix']
    oxy_btl_col = settings.bottle_inputs['btl_oxy']
    rinko_volts = settings.ctd_inputs['rinko_oxy']
    rinko_btl_volts = settings.bottle_inputs['rinko_oxy']
    dov_col = settings.ctd_inputs['dov']

    # unsure about these
    lat_col = config['analytical_inputs']['lat']
    lon_col = config['analytical_inputs']['lon']
    lat_btl_col = config['inputs']['lat']
    lon_btl_col = config['inputs']['lon']
    sal_btl_col = config['analytical_inputs']['btl_salt']
    t_btl_col = config['inputs']['t']
    t_col = config['analytical_inputs']['t']

    # calculate sigma
    btl_data_all['sigma_btl'] = oxy_fitting.sigma_from_CTD(btl_data_all[sal_btl_col], btl_data_all[t_btl_col], btl_data_all[p_btl_col], btl_data_all[lon_btl_col], btl_data_all[lat_btl_col])
    time_data_all['sigma_ctd'] = oxy_fitting.sigma_from_CTD(time_data_all[sal_col], time_data_all[t_col], time_data_all[p_col], time_data_all[lon_col], time_data_all[lat_col])

    btl_data_all[oxy_btl_col] = oxy_fitting.calculate_bottle_oxygen(ssscc, btl_data_all['SSSCC'], btl_data_all['TITR_VOL'], btl_data_all['TITR_TEMP'], btl_data_all['FLASKNO'])
    btl_data_all[oxy_btl_col] = oxy_fitting.oxy_ml_to_umolkg(btl_data_all[oxy_btl_col], btl_data_all['sigma_btl'])
    btl_data_all['OXYGEN_FLAG_W'] = oxy_fitting.flag_winkler_oxygen(btl_data_all[oxy_btl_col])

    # Calculate SA and PT
    btl_data_all['SA'] = gsw.SA_from_SP(btl_data_all[sal_btl_col], btl_data_all[p_btl_col], btl_data_all[lon_btl_col], btl_data_all[lat_btl_col])
    btl_data_all['PT'] = gsw.pt0_from_t(btl_data_all['SA'], btl_data_all[t_btl_col], btl_data_all[p_btl_col])

    time_data_all['SA'] = gsw.SA_from_SP(time_data_all[sal_col], time_data_all[p_col], time_data_all[lon_col], time_data_all[lat_col])
    time_data_all['PT'] = gsw.pt0_from_t(time_data_all['SA'], time_data_all[t_col], time_data_all[p_col])

    # Calculate OS in µmol/kg

    btl_data_all['OS_btl'] = oxy_fitting.os_umol_kg(btl_data_all['SA'], btl_data_all['PT'])
    time_data_all['OS_ctd'] = oxy_fitting.os_umol_kg(time_data_all['SA'], time_data_all['PT'])

    # collect data by stations

    btl_data_all['oxy_stn_group'] = btl_data_all['SSSCC'].str[0:3]
    time_data_all['oxy_stn_group'] = time_data_all['SSSCC'].str[0:3]

    station_list = time_data_all['oxy_stn_group'].unique()
    station_list = station_list.tolist()
    station_list.sort()

    btl_data_oxy = btl_data_all.copy()#loc[btl_data_all['OXYGEN'].notnull()]

    rinko_coef0 = rinko.rinko_o2_cal_parameters()

    all_rinko_df = pd.DataFrame()
    all_sbe43_df = pd.DataFrame()
    rinko_dict = {}
    sbe43_dict = {}
    for station in station_list:
        
        #time_data = time_data_all[time_data_all['SSSCC'].str[0:3] == station].copy()
        #btl_data = btl_data_oxy[btl_data_oxy['SSSCC'].str[0:3] == station].copy()
        
        time_data = time_data_all[time_data_all['oxy_stn_group'] == station].copy()
        btl_data = btl_data_oxy[btl_data_oxy['oxy_stn_group'] == station].copy()
        
        if (btl_data['OXYGEN_FLAG_W'] == 9).all():
            rinko_dict[station] = np.full(8,np.nan)
            sbe43_dict[station] = np.full(7,np.nan)
            print(station + ' skipped, all oxy data is NaN')
            continue

        rinko_coef, rinko_oxy_df = rinko.rinko_oxygen_fit(btl_data[p_btl_col],btl_data[oxy_btl_col],btl_data['sigma_btl'],time_data['sigma_ctd'],
                            time_data['OS_ctd'],time_data[p_col],time_data[t_col],time_data[rinko_volts],rinko_coef0, btl_data['SSSCC'])

        station_ssscc = time_data['SSSCC'].values[0]
        hex_file = hex_prefix + station_ssscc + hex_postfix
        xml_file = xml_prefix + station_ssscc + xml_postfix
        sbe_coef0 = oxy_fitting.get_SB_coef(hex_file, xml_file)
        sbe_coef, sbe_oxy_df = oxy_fitting.sbe43_oxy_fit(btl_data[p_btl_col], btl_data[oxy_btl_col], btl_data['sigma_btl'],time_data['sigma_ctd'],
                                                        time_data['OS_ctd'],time_data[p_col],time_data[t_col],time_data[dov_col],
                                                        time_data['scan_datetime'],sbe_coef0,btl_data['SSSCC'])
        
        rinko_dict[station] = rinko_coef
        sbe43_dict[station] = sbe_coef
        all_rinko_df = pd.concat([all_rinko_df,rinko_oxy_df])
        all_sbe43_df = pd.concat([all_sbe43_df,sbe_oxy_df])
        
        print(station + ' Done!')

    sbe43_coef_df = oxy_fitting.create_coef_df(sbe43_dict)
    rinko_coef_df = oxy_fitting.create_coef_df(rinko_dict)

    btl_data_all = btl_data_all.merge(all_rinko_df, left_on=['SSSCC',p_btl_col], right_on=['SSSCC_rinko','CTDPRS_rinko_btl'],how='left')

    btl_data_all = btl_data_all.merge(all_sbe43_df, left_on=['SSSCC',p_btl_col], right_on=['SSSCC_sbe43','CTDPRS_sbe43_btl'],how='left')

    btl_data_all.drop(list(btl_data_all.filter(regex = 'rinko')), axis = 1, inplace = True)

    btl_data_all.drop(list(btl_data_all.filter(regex = 'sbe43')), axis = 1, inplace = True)

    ### Handle Missing Values

    X = btl_data_all['CTDOXYVOLTS_x'], btl_data_all['CTDPRS'], btl_data_all[t_btl_col], btl_data_all['dv_dt_x'], btl_data_all['OS_btl']
    rinko_X = btl_data_all[p_btl_col], btl_data_all[t_btl_col], btl_data_all['CTDRINKOVOLTS'], btl_data_all['OS_btl']

    for station in btl_data_all['SSSCC']:
        coef_43 = sbe43_coef_df.loc[station[0:3]]
        coef_rinko = rinko_coef_df.loc[station[0:3]]
        btl_data_all['CTDOXY_fill'] = oxy_fitting.oxy_equation(X, coef_43[0], coef_43[1], coef_43[2], coef_43[3], coef_43[4], coef_43[5], coef_43[6])
        btl_data_all['CTDRINKO_fill'] = rinko.rinko_curve_fit_eq(rinko_X, coef_rinko[0], coef_rinko[1], coef_rinko[2], coef_rinko[3], coef_rinko[4], coef_rinko[5], coef_rinko[6], coef_rinko[7])

    btl_data_all.loc[btl_data_all['CTDOXY'].isnull(),'CTDOXY_FLAG_W'] = 2
    btl_data_all.loc[btl_data_all['CTDOXY'].isnull(),'CTDOXY'] = btl_data_all.loc[btl_data_all['CTDOXY'].isnull(),'CTDOXY_fill']
    btl_data_all.loc[btl_data_all['CTDRINKO'].isnull(),'CTDRINKO_FLAG_W'] = 2
    btl_data_all.loc[btl_data_all['CTDRINKO'].isnull(),'CTDRINKO'] = btl_data_all.loc[btl_data_all['CTDRINKO'].isnull(),'CTDRINKO_fill']

    btl_data_all = oxy_fitting.flag_oxy_data(btl_data_all)
    btl_data_all = oxy_fitting.flag_oxy_data(btl_data_all,ctd_oxy_col='CTDRINKO',flag_col='CTDRINKO_FLAG_W')

    btl_data_all['res_rinko'] = btl_data_all['OXYGEN'] - btl_data_all['CTDRINKO']
    btl_data_all['res_sbe43'] = btl_data_all['OXYGEN'] - btl_data_all['CTDOXY']
    btl_data_all.loc[np.abs(btl_data_all['res_rinko']) >=6 , 'CTDRINKO_FLAG_W'] = 3 # what's the benefit of np.abs() vs abs()?
    btl_data_all.loc[np.abs(btl_data_all['res_sbe43']) >=6 , 'CTDOXY_FLAG_W'] = 3

    time_data_all['CTDOXY'] = '-999'
    time_data_all['CTDRINKO'] = '-999'
    btl_data_all.sort_values(by='sigma_btl',inplace=True)
    time_data_all.sort_values(by='sigma_ctd',inplace=True)
    for station in station_list:
        rinko_coef = rinko_coef_df.loc[station].values
        if np.isnan(rinko_coef).all():
            print(station + ' data bad, leaving -999 values and flagging as 9')
            time_data_all.loc[time_data_all['SSSCC'].str[0:3] == station,'CTDRINKO_FLAG_W'] = 9
            time_data_all.loc[time_data_all['SSSCC'].str[0:3] == station,'CTDOXY_FLAG_W'] = 9
            continue

        time_data = time_data_all[time_data_all['oxy_stn_group'] == station].copy()
        time_data['CTDOXY'] = oxy_fitting.SB_oxy_eq(sbe43_coef_df.loc[station],time_data[dov_col],time_data[p_col],time_data[t_col],time_data['dv_dt'],time_data['OS_ctd'])
        time_data['CTDRINKO'] = rinko.rinko_curve_fit_eq((time_data[p_col],time_data[t_col],time_data[rinko_volts],time_data['OS_ctd']),rinko_coef[0],rinko_coef[1],
                                                        rinko_coef[2],rinko_coef[3],rinko_coef[4],rinko_coef[5],rinko_coef[6],rinko_coef[7])
        
        time_data_all.loc[time_data_all['SSSCC'].str[0:3] == station,'CTDOXY'] = time_data['CTDOXY']
        time_data_all.loc[time_data_all['SSSCC'].str[0:3] == station,'CTDRINKO'] = time_data['CTDRINKO']
        time_data_all.loc[time_data_all['SSSCC'].str[0:3] == station,'CTDRINKO_FLAG_W'] = 2
        time_data_all.loc[time_data_all['SSSCC'].str[0:3] == station,'CTDOXY_FLAG_W'] = 2
        print(station + ' time data done')

    #####

    #######################
    # turn correction instantaneous CTD measurements into points @ bottle firings
    # subprocess.run(['ctd_to_bottle.py'], stdout=subprocess.PIPE)

################ Clean and export data #######################
    # btl_data_all = process_ctd.merge_cond_flags(btl_data_all,qual_flag_cond, c_btl_col) # make code not require specific column input
    btl_data_all = process_ctd.merge_refcond_flags(btl_data_all,qual_flag_cond)
    # btl_data_all = process_ctd.merge_temp_flags(btl_data_all, qual_flag_temp, t_btl_col)
    btl_data_all = process_ctd.merged_reftemp_flags(btl_data_all,qual_flag_temp)

### Export Quality Flags

    qual_flag_temp.to_csv('data/logs/qual_flag_temp_new.csv',index=False)
    qual_flag_cond.to_csv('data/logs/qual_flag_cond_new.csv',index=False)

### Clean up Bottle Data by removing rows with no ctd data

    btl_data_all = oxy_fitting.clean_oxygen_df(btl_data_all) # clean up column names
    btl_data_all =  btl_data_all.dropna(subset=cols)

### Create CT Files(oxygen done by hand)

    # needed to run when there's no data available yet (i.e. at sea)
    # btl_data_all['OXYGEN'] = -999
    # btl_data_all['OXYGEN'] = -999
    # time_data_all['CTDOXY'] = -999

    # depreciated according to process_ctd.py (MK)
    # process_ctd.export_btl_data(btl_data_all,expocode,sectionID,expocode)
    # process_ctd.export_time_data(time_data_all,ssscc,int(sample_rate),int(search_time),expocode,sectionID,ctd,p_column_names,p_column_units)

    ### MK
    # create depth_log.csv
    station_list = time_data_all['SSSCC'].str[0:3].unique()
    depth_dict = {}
    for station in station_list:
        print(station)
        time_data = time_data_all[time_data_all['SSSCC'].str[0:3] == station].copy()
        max_depth = process_ctd.find_cast_depth(time_data['CTDPRS'],time_data['GPSLAT'],time_data['ALT'])
        depth_dict[station] = max_depth

    depth_df = pd.DataFrame.from_dict(depth_dict,orient='index')
    depth_df.reset_index(inplace=True)
    depth_df.rename(columns={0:'DEPTH', 'index':'STNNBR'}, inplace=True)

    depth_df.to_csv('data/logs/depth_log.csv',index=False)

    # needs depth_log.csv, manual_depth_log.csv
    process_ctd.export_ct1(time_data_all,ssscc,expocode,sectionID,ctd,p_column_names,p_column_units)

    # section for making plots/etc.

    # section to output stdev/etc. (maybe in calib to use for flagging data)


def main(argv):
    '''Run everything.
    '''
    process_all_new()

if __name__ == '__main__':
    main(sys.argv[1:])

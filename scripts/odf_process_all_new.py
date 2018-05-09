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
#import ctdcal
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import sys
sys.path.append('ctdcal/')
import ctdcal.oxy_fitting as oxy_fitting
import configparser
#import pandas as pd
#import numpy as np

def process_all_new():

    prefix = 'nbp1707_'
    btl = '.btl'
    ros = '.ros'
    cnv = '.cnv'
    ext_hex = '.hex'
    xmlcon = '.XMLCON'
    csv = '.csv'
    #directories and index files
    ssscc_file = 'data/ssscc.csv'
    dir_logs = 'data/logs/'
    qual_dir_temp_primary = 'data/quality_code/temp_primary/'
    qual_dir_temp_secondary = 'data/quality_code/temp_secondary/'
    qual_dir_cond_primary = 'data/quality_code/cond_primary/'
    qual_dir_cond_secondary = 'data/quality_code/cond_secondary/'
    qual_flag_temp = 'quality_flag_temp'
    qual_flag_cond = 'quality_flag_cond'
    fit_t1 = 'fitting_t1'
    fit_t2 = 'fitting_t2'
    fit_c1 = 'fitting_c1'
    fit_c2 = 'fitting_c2'

    iniFile = 'data/ini-files/configuration.ini'
    config = configparser.RawConfigParser()
    config.read(iniFile)

###################################
        #Initialise Configuration Parameters
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

    time_zone = config['inputs']['time_zone']
    p_col = config['analytical_inputs']['p']
    p_btl_col = config['inputs']['p']
    t_col = config['analytical_inputs']['t']
    t_btl_col = config['inputs']['t']
    t1_col = config['analytical_inputs']['t1']
    t1_btl_col = config['inputs']['t1']
    t2_col = config['analytical_inputs']['t2']
    t2_btl_col = config['inputs']['t2']
    c_col = config['analytical_inputs']['c']
    c1_col = config['analytical_inputs']['c1']
    c1_btl_col = config['inputs']['c1']
    c2_col = config['analytical_inputs']['c2']
    c2_btl_col = config['inputs']['c2']
    sal_col = config['analytical_inputs']['salt']
    sal_btl_col = config['inputs']['salt']
    btl_sal_col = config['analytical_inputs']['btl_salt']
    dov_col = config['analytical_inputs']['dov']
    dov_btl_col = config['inputs']['dov']
    dopl_col = config['analytical_inputs']['dopl']
    dopl_btl_col = config['inputs']['dopl']
    dopkg_col = config['analytical_inputs']['dopkg']
    btl_oxy_col = config['analytical_inputs']['btl_oxy']
    xmis_col = config['analytical_inputs']['xmis']
    fluor_col = config['analytical_inputs']['fluor']
    timedate = config['analytical_inputs']['datetime']
    lat_col = config['analytical_inputs']['lat']
    lat_btl_col = config['inputs']['lat']
    lon_col = config['analytical_inputs']['lon']
    lon_btl_col = config['inputs']['lon']
    reft_col = config['inputs']['reft']
    btl_num_col = config['inputs']['btl_num']

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
#####ADDED BY KJ ##################
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
    
    ###########################################################################
    
    ### File I/O
    
    # Load in all bottle, time, ref_data files into DataFrame
    
    btl_data_all = process_ctd.load_all_ctd_files(ssscc,btl_data_prefix,
                                                  btl_data_postfix,'bottle')
    time_data_all = process_ctd.load_all_ctd_files(ssscc,time_data_prefix,
                                                   time_data_postfix,'time')
    
    ################################################################################
    
    ### Pressure Calibration
    
    # Determine Pressure offset from logs
    
    pressure_log = process_ctd.load_pressure_logs(p_log_file)
    p_off = process_ctd.get_pressure_offset(pressure_log)
#    p_off = np.around(p_off,decimals=4)

        
    # Apply Pressure offset to data
    
    btl_data_all = fit_ctd.apply_pressure_offset(btl_data_all,p_off)
    time_data_all = fit_ctd.apply_pressure_offset(time_data_all,p_off)
    
    
    ###########################################################################
    
    
    ### Temperature Calibration
    
    
    # Remove unmeasured datapoints
    for x in range(2):
        
     ## Second order calibration
             
        df_temp_good = process_ctd.prepare_fit_data(btl_data_all,reft_col)
        
        coef_temp_1,df_ques_t1,df_reft = process_ctd.calibrate_temperature(df_temp_good,2,'P',1,xRange='800:6000')
        coef_temp_2,df_ques_t2,df_reft2 = process_ctd.calibrate_temperature(df_temp_good,2,'P',2,xRange='1500:6000')
    
    # Apply fitting coef to data
#    
        btl_data_all[t1_btl_col] = fit_ctd.temperature_polyfit(btl_data_all,coef_temp_1,t1_btl_col,p_btl_col)
        time_data_all[t1_col] = fit_ctd.temperature_polyfit(time_data_all,coef_temp_1,t1_col,p_col)        
    
    
        btl_data_all[t2_btl_col] = fit_ctd.temperature_polyfit(btl_data_all,coef_temp_2,t2_btl_col,p_btl_col)
        time_data_all[t2_col] = fit_ctd.temperature_polyfit(time_data_all,coef_temp_2,t2_col,p_col)   
        
    # Construct Quality Flag file
    
        qual_flag_temp = process_ctd.combine_quality_flags(df_ques_t1,df_ques_t2,df_reft)
   
    ## First order calibtation
    
        df_temp_good = process_ctd.prepare_fit_data(btl_data_all,reft_col)
        
        coef_prim_temp,df_ques_t1,df_reft = process_ctd.calibrate_temperature(df_temp_good,1,'T',1,coef=[])
        coef_sec_temp,df_ques_t2,df_quest_reft2 = process_ctd.calibrate_temperature(df_temp_good,1,'T',2,coef=[])
        
    # Apply fitting coef to data
    
        btl_data_all[t1_btl_col] = fit_ctd.temperature_polyfit(btl_data_all,coef_prim_temp,t1_btl_col,p_btl_col)
        time_data_all[t1_col] = fit_ctd.temperature_polyfit(time_data_all,coef_prim_temp,t1_col,p_col)
    
        btl_data_all[t2_btl_col] = fit_ctd.temperature_polyfit(btl_data_all,coef_sec_temp,t2_btl_col,p_btl_col)
        time_data_all[t2_col] = fit_ctd.temperature_polyfit(time_data_all,coef_sec_temp,t2_col,p_col)
            
        qual_flag_temp = process_ctd.combine_quality_flags(df_ques_t1,df_ques_t2,df_reft)
    
    ###########################################################################
    
#    
    ### Conductivity Calibration
    
    df_cond_good = process_ctd.prepare_fit_data(btl_data_all,refc_col)
    #df_cond_good,refc_data_clean = process_ctd.prepare_conductivity_data(ssscc,btl_data_all,refc_data_all)

    
    #coef_cond_1,df_ques_c1,df_refc = process_ctd.calibrate_conductivity(df_cond_good,refc_data_clean,2,'P',1,xRange='800:6000')
    coef_cond_1,df_ques_c1,df_refc = process_ctd.calibrate_conductivity(df_cond_good,2,'P',1,xRange='800:6000')
    
    #coef_cond_2,df_ques_c2,df_refc2 = process_ctd.calibrate_conductivity(df_cond_good,refc_data_clean,2,'P',2,xRange='1500:6000')
    coef_cond_2,df_ques_c2,df_refc2 = process_ctd.calibrate_conductivity(df_cond_good,2,'P',2,xRange='1500:6000')
    
    
    # Apply fitting coef to data
    
    btl_data_all[c1_btl_col], btl_data_all[sal_col] = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_1,t1_btl_col,c1_btl_col,p_btl_col)
    time_data_all[c1_col], time_data_all[sal_col] = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_1,t1_col,c1_col,p_col)
    
#    btl_data_all[c1_btl_col] = fit_ctd.conductivity_polyfit(coef_cond_1,btl_data_all[p_btl_col],btl_data_all[t1_btl_col],btl_data_all[c1_btl_col])
#    time_data_all[c1_col] = fit_ctd.conductivity_polyfit(coef_cond_1,time_data_all[p_col],time_data_all[t1_col],time_data_all[c1_col])
    
    btl_data_all[c2_btl_col], sensor_2_sal = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_2,t2_btl_col,c2_btl_col,p_btl_col)
    time_data_all[c2_col], sensor_2_sal = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_2,t2_col,c2_col,p_col)

#    btl_data_all[c2_btl_col] = fit_ctd.conductivity_polyfit(coef_cond_2,btl_data_all[p_btl_col],btl_data_all[t2_btl_col],btl_data_all[c2_btl_col])
#    time_data_all[c2_col] = fit_ctd.conductivity_polyfit(coef_cond_2,time_data_all[p_col],time_data_all[t2_col],time_data_all[c2_col])
#    
    qual_flag_cond = process_ctd.combine_quality_flags(df_ques_c1,df_ques_c2,df_refc)

#    # Repeat WRT Temp
#    
    df_cond_good = process_ctd.prepare_fit_data(btl_data_all,refc_col)
    
    #coef_cond_1,df_ques_c1,df_refc = process_ctd.calibrate_conductivity(df_cond_good,refc_data_clean,2,'T',1)
    coef_cond_1,df_ques_c1,df_refc = process_ctd.calibrate_conductivity(df_cond_good,2,'T',1)
    
    #coef_cond_2,df_ques_c2,df_refc2 = process_ctd.calibrate_conductivity(df_cond_good,refc_data_clean,2,'T',2)
    coef_cond_2,df_ques_c2,df_refc2 = process_ctd.calibrate_conductivity(df_cond_good,2,'T',2)
    
    
    # Apply fitting coef to data
    
    btl_data_all[c1_btl_col], btl_data_all[sal_col] = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_1,t1_btl_col,c1_btl_col,p_btl_col)
    time_data_all[c1_col], time_data_all[sal_col] = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_1,t1_col,c1_col,p_col)
    
    btl_data_all[c2_btl_col], sensor_2_sal = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_2,t2_btl_col,c2_btl_col,p_btl_col)
    time_data_all[c2_col], sensor_2_sal = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_2,t2_col,c2_col,p_col)
    
    qual_flag_cond = process_ctd.combine_quality_flags(df_ques_c1,df_ques_c2,df_refc)
    
#    btl_data_all[c1_btl_col] = fit_ctd.conductivity_polyfit(coef_cond_1,btl_data_all[p_btl_col],btl_data_all[t1_btl_col],btl_data_all[c1_btl_col])
#    time_data_all[c1_col] = fit_ctd.conductivity_polyfit(coef_cond_1,time_data_all[p_col],time_data_all[t1_col],time_data_all[c1_col])
#    
#    btl_data_all[c2_btl_col] = fit_ctd.conductivity_polyfit(coef_cond_2,btl_data_all[p_btl_col],btl_data_all[t2_btl_col],btl_data_all[c2_btl_col])
#    time_data_all[c2_col] = fit_ctd.conductivity_polyfit(coef_cond_2,time_data_all[p_col],time_data_all[t2_col],time_data_all[c2_col])
#    
    # Finally WRT to Cond
    
    df_cond_good = process_ctd.prepare_fit_data(btl_data_all,refc_col)
    
    #coef_cond_1,df_ques_c1,df_refc = process_ctd.calibrate_conductivity(df_cond_good,refc_data_clean,2,'C',1)
    coef_cond_1,df_ques_c1,df_refc = process_ctd.calibrate_conductivity(df_cond_good,2,'C',1)
    
    #coef_cond_2,df_ques_c2,df_refc2 = process_ctd.calibrate_conductivity(df_cond_good,refc_data_clean,2,'C',2)
    coef_cond_2,df_ques_c2,df_refc2 = process_ctd.calibrate_conductivity(df_cond_good,2,'C',2)
    
    # Apply fitting coef to data
    
    btl_data_all[c1_btl_col], btl_data_all[sal_col] = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_1,t1_btl_col,c1_btl_col,p_btl_col)
    time_data_all[c1_col], time_data_all[sal_col] = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_1,t1_col,c1_col,p_col)
    
    btl_data_all[c2_btl_col], sensor_2_sal = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_2,t2_btl_col,c2_btl_col,p_btl_col)
    time_data_all[c2_col], sensor_2_sal = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_2,t2_col,c2_col,p_col)
    
    qual_flag_cond = process_ctd.combine_quality_flags(df_ques_c1,df_ques_c2,df_refc)
    
#    btl_data_all[c1_btl_col] = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_1,t1_btl_col,c1_btl_col,p_btl_col)
#    time_data_all[c1_col] = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_1,t1_col,c1_col,p_col)
#    
#    btl_data_all[c2_btl_col] = fit_ctd.conductivity_polyfit(btl_data_all,coef_cond_2,t2_btl_col,c2_btl_col,p_btl_col)
#    time_data_all[c2_col] = fit_ctd.conductivity_polyfit(time_data_all,coef_cond_2,t2_col,c2_col,p_col)
##
#    qual_flag_cond = process_ctd.combine_quality_flags(df_ques_c1,df_ques_c2,df_refc)
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

    # Clean and export data
#     ### THIS STILL NEEDS TO BE TESTED
#    process_ctd.export_time_data(time_data_all,ssscc,int(sample_rate),int(search_time),expocode,sectionID,ctd,p_column_names,p_column_units)
#    cruise_line = 'P06'
#    process_ctd.export_btl_data(btl_data_all,expocode,sectionID,cruise_line)


def main(argv):
    '''Run everything.
    '''
    process_all_new()

if __name__ == '__main__':
    main(sys.argv[1:])
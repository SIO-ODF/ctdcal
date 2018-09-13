#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 11:42:27 2018

@author: k3jackson
"""

import os
import time
import sys
sys.path.append('ctdcal/')
import settings
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import pandas as pd
import gsw
import ctdcal.oxy_fitting as oxy_fitting
import sbe_convert

def process_all_new():

    
    # Directory and file information
    expocode = settings.cruise['expocode']
    sectionID = settings.cruise['sectionid']
    raw_directory = settings.ctd_processing_dir['raw_data_directory']
    time_directory = settings.ctd_processing_dir['time_data_directory']
    converted_directory = settings.ctd_processing_dir['converted_directory']
    pressure_directory = settings.ctd_processing_dir['pressure_data_directory']
    oxygen_directory = settings.ctd_processing_dir['oxygen_directory']
    btl_directory = settings.ctd_processing_dir['bottle_directory']
    o2flask_file = settings.ctd_processing_dir['o2flask_file']
    log_directory = settings.ctd_processing_dir['log_directory']
    p_log_file = settings.ctd_processing_dir['pressure_log']
    hex_prefix = settings.ctd_processing_dir['hex_prefix']
    hex_postfix = settings.ctd_processing_dir['hex_postfix']
    xml_prefix = settings.ctd_processing_dir['xml_prefix']
    xml_postfix = settings.ctd_processing_dir['xml_postfix']
    
    # CTD Data Inputs
    p_col = settings.ctd_inputs['p']
    t_col = settings.ctd_inputs['t']
    t1_col = settings.ctd_inputs['t1']
    t2_col = settings.ctd_inputs['t2']
    c_col = settings.ctd_inputs['c']
    c1_col = settings.ctd_inputs['c1']
    c2_col = settings.ctd_inputs['c2']
    sal_col = settings.ctd_inputs['salt']
    dov_col = settings.ctd_inputs['dov']
    lat_col = settings.ctd_inputs['lat']
    lon_col = settings.ctd_inputs['lon']
    time_col = settings.ctd_inputs['scan_datetime']

        # Bottle Data Inputs
    p_btl_col = settings.bottle_inputs['p']
    t_btl_col = settings.bottle_inputs['t']
    t1_btl_col = settings.bottle_inputs['t1']
    t2_btl_col = settings.bottle_inputs['t2']
    c_btl_col = settings.bottle_inputs['c']
    c1_btl_col = settings.bottle_inputs['c1']
    c2_btl_col = settings.bottle_inputs['c2']
    reft_col = settings.bottle_inputs['reft']
    cond_col = settings.bottle_inputs['btl_cond']
    cr_avg = settings.bottle_inputs['cond_ratio']
    bath_temp = settings.bottle_inputs['bath_temp']
    sal_btl_col = settings.bottle_inputs['salt']
    dov_btl_col = settings.bottle_inputs['dov']
    lat_btl_col = settings.bottle_inputs['lat']
    lon_btl_col = settings.bottle_inputs['lon']
    oxy_btl_col = settings.bottle_inputs['btl_oxy']
    time_btl_col = settings.bottle_inputs['scan_datetime']
    
    # CTD Information    
    sample_rate = settings.ctd_processing_constants['sample_rate']
    search_time = settings.ctd_processing_constants['roll_filter_time']
    ctd = settings.ctd_processing_constants['ctd_serial']
    
    p_column_names = settings.pressure_series_output['column_name']
    p_column_units = settings.pressure_series_output['column_units']
    
    btl_data_prefix = 'data/bottle/'
    btl_data_postfix = '_btl_mean.pkl'
    time_data_prefix = 'data/time/'
    time_data_postfix = '_time.pkl'
    p_log_file = 'data/logs/ondeck_pressure.csv'
 

    
    # Columns from btl and ctd file to be read:
    btl_cols = settings.btl_input_array
    ctd_cols = settings.ctd_input_array
    
    ssscc = settings.ssscc

#    time_start = time.perf_counter()
    cnv_dir_list = os.listdir(converted_directory)
    time_dir_list = os.listdir(time_directory)
    btl_dir_list = os.listdir(btl_directory)

    for station in ssscc:
        if '{}.pkl'.format(station) in cnv_dir_list:
            continue
        #convert hex to ctd
        hex_file = hex_prefix + station + hex_postfix
        xml_file = xml_prefix + station + xml_postfix
        
        sbe_convert.convert_sbe(station, hex_file, xml_file, converted_directory)
        print('Converted_sbe SSSCC: ' + station + ' done')

#    time_convert = time.perf_counter()

    for station in ssscc:
        if '{}_time.pkl'.format(station) in time_dir_list:
            continue
        sbe_convert.sbe_metadata(station)
        print('sbe_metadata SSSCC: ' + station + ' done')

    for station in ssscc:
        if '{}_btl_mean.pkl'.format(station) in btl_dir_list:
            continue
        #process bottle file
        sbe_convert.process_bottle(station)
        print('process_bottle SSSCC: ' + station + ' done')
#    time_bottle = time.perf_counter()
    
    ###########################################################################
    
    ### File I/O
    
    # Load in all bottle, time, ref_data files into DataFrame
    
    btl_data_all = process_ctd.load_all_ctd_files(ssscc,btl_data_prefix,
                                                  btl_data_postfix,'bottle',btl_cols)
    time_data_all = process_ctd.load_all_ctd_files(ssscc,time_data_prefix,
                                                   time_data_postfix,'time',ctd_cols)
    
    ################################################################################
    
    ### Pressure Calibration
    
    # Determine Pressure offset from logs
    
    pressure_log = process_ctd.load_pressure_logs(p_log_file)
    p_off = process_ctd.get_pressure_offset(pressure_log.ondeck_start_p,pressure_log.ondeck_end_p)
    
    btl_data_all[p_btl_col] = fit_ctd.apply_pressure_offset(btl_data_all[p_btl_col], p_off)
    time_data_all[p_col] = fit_ctd.apply_pressure_offset(time_data_all[p_col],p_off)
    
    
    ###########################################################################

    df_ques_t1 = pd.DataFrame()
    df_ques_t2 = pd.DataFrame()
    
    df_ques_c1 = pd.DataFrame()
    df_ques_c2 = pd.DataFrame()

    ### Temperature Calibration   
    for x in range(2):
        
     # Second order calibration
             
        df_temp_good = process_ctd.prepare_fit_data(btl_data_all, reft_col)
        
        df_ques_reft = process_ctd.quality_check(df_temp_good[t2_btl_col], df_temp_good[t1_btl_col], df_temp_good[p_btl_col], df_temp_good['SSSCC'], df_temp_good['btl_fire_num'], 'quest')
        df_ques_reft['Parameter'] = 'REF_TEMP'
        


        if settings.do_primary == 1:
            coef_temp_1,df_ques_t1 = process_ctd.calibrate_param(df_temp_good[t1_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'TP', 2, df_temp_good.SSSCC, df_temp_good.btl_fire_num, xRange='800:6000')
            btl_data_all[t1_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t1_btl_col], btl_data_all[p_btl_col], coef_temp_1)
            time_data_all[t1_col] = fit_ctd.temperature_polyfit(time_data_all[t1_col], time_data_all[p_col], coef_temp_1)
        
        elif settings.do_secondary == 1:
            coef_temp_2,df_ques_t2 = process_ctd.calibrate_param(df_temp_good[t2_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'TP', 2, df_temp_good.SSSCC, df_temp_good.btl_fire_num, xRange='1500:6000')
            btl_data_all[t2_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t2_btl_col], btl_data_all[p_btl_col], coef_temp_2)
            time_data_all[t2_col] = fit_ctd.temperature_polyfit(time_data_all[t2_col], time_data_all[p_col], coef_temp_2)
    
    # Apply fitting coef to data
  
    # Construct Quality Flag file
    
        qual_flag_temp = process_ctd.combine_quality_flags([df_ques_reft,df_ques_t1,df_ques_t2])
   
    ## First order calibtation
    
        df_temp_good = process_ctd.prepare_fit_data(btl_data_all, reft_col)
        
#        df_ques_reft = process_ctd.quality_check(df_temp_good[t2_btl_col], df_temp_good[t1_btl_col], df_temp_good[p_btl_col], df_temp_good['SSSCC'], df_temp_good['btl_fire_num'], 'quest')
#        df_ques_reft['Parameter'] = 'REF_TEMP'
        if settings.do_primary == 1:
            coef_temp_prim,df_ques_t1 = process_ctd.calibrate_param(df_temp_good[t1_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'T', 1, df_temp_good.SSSCC, df_temp_good.btl_fire_num)
            btl_data_all[t1_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t1_btl_col], btl_data_all[p_btl_col], coef_temp_prim)
            time_data_all[t1_col] = fit_ctd.temperature_polyfit(time_data_all[t1_col], time_data_all[p_col], coef_temp_prim)
        
        elif settings.do_secondary == 1:
            coef_temp_sec,df_ques_t2 = process_ctd.calibrate_param(df_temp_good[t2_btl_col], df_temp_good[reft_col], df_temp_good[p_btl_col], 'T', 1, df_temp_good.SSSCC, df_temp_good.btl_fire_num)
            btl_data_all[t2_btl_col] = fit_ctd.temperature_polyfit(btl_data_all[t2_btl_col], btl_data_all[p_btl_col], coef_temp_sec)
            time_data_all[t2_col] = fit_ctd.temperature_polyfit(time_data_all[t2_col], time_data_all[p_col], coef_temp_sec)
        
        
        
    # Apply fitting coef to data
            
        qual_flag_temp = process_ctd.combine_quality_flags([df_ques_reft,df_ques_t1,df_ques_t2])
    
    ###########################################################################
    
#    
    ### Conductivity Calibration
    for x in range(2):
        
        btl_data_all[cond_col] = fit_ctd.CR_to_cond(btl_data_all[cr_avg], btl_data_all[bath_temp], btl_data_all[t1_btl_col], btl_data_all[p_btl_col])
        df_cond_good = process_ctd.prepare_fit_data(btl_data_all, cond_col)

        df_ques_refc = process_ctd.quality_check(df_cond_good[c2_btl_col], df_temp_good[c1_btl_col], df_temp_good[p_btl_col], df_temp_good['SSSCC'], df_temp_good['btl_fire_num'], 'quest')
        df_ques_refc['Parameter'] = 'REF_COND'

    # Second Order Calibration
        if settings.do_primary == 1:
            coef_cond_1,df_ques_c1 = process_ctd.calibrate_param(df_cond_good[c1_btl_col], df_cond_good[cond_col], df_cond_good[p_btl_col], 'CP', 2, df_cond_good['SSSCC'], df_cond_good['btl_fire_num'], xRange='800:6000')
            btl_data_all[c1_btl_col], btl_data_all[sal_btl_col] = fit_ctd.conductivity_polyfit(btl_data_all[c1_btl_col], btl_data_all[t1_btl_col], btl_data_all[p_btl_col], coef_cond_1)
            time_data_all[c1_col],time_data_all[sal_col] = fit_ctd.conductivity_polyfit(time_data_all[c1_col], time_data_all[t1_col], time_data_all[p_col], coef_cond_1)
        
        elif settings.do_secondary == 1:
            coef_cond_2,df_ques_c2 = process_ctd.calibrate_param(df_cond_good[c2_btl_col], df_cond_good[cond_col], df_cond_good[p_btl_col], 'CP', 2, df_cond_good['SSSCC'], df_cond_good['btl_fire_num'], xRange='1500:6000')
            btl_data_all[c2_btl_col], sal_2 = fit_ctd.conductivity_polyfit(btl_data_all[c2_btl_col], btl_data_all[t2_btl_col], btl_data_all[p_btl_col] ,coef_cond_2)
            time_data_all[c2_btl_col],sal2 = fit_ctd.conductivity_polyfit(time_data_all[c2_col], time_data_all[t2_col], time_data_all[p_col], coef_cond_2)

        qual_flag_cond = process_ctd.combine_quality_flags([df_ques_c1,df_ques_c2,df_ques_refc])
    
        
        btl_data_all[cond_col] = fit_ctd.CR_to_cond(btl_data_all[cr_avg], btl_data_all[bath_temp], btl_data_all[t1_btl_col], btl_data_all[p_btl_col])
        df_cond_good = process_ctd.prepare_fit_data(btl_data_all,cond_col)
    
        if settings.do_primary == 1:
            coef_cond_prim,df_ques_c1 = process_ctd.calibrate_param(df_cond_good[c1_btl_col], df_cond_good[cond_col], df_cond_good[p_btl_col], 'C', 2 , df_cond_good['SSSCC'], df_cond_good['btl_fire_num'])
            btl_data_all[c1_btl_col], btl_data_all[sal_btl_col] = fit_ctd.conductivity_polyfit(btl_data_all[c1_btl_col], btl_data_all[t1_btl_col], btl_data_all[p_btl_col], coef_cond_prim)        
            time_data_all[c1_col],time_data_all[sal_col] = fit_ctd.conductivity_polyfit(time_data_all[c1_col], time_data_all[t1_col], time_data_all[p_col], coef_cond_prim)
        
        elif settings.do_secondary == 1:
            coef_cond_sec,df_ques_c2 = process_ctd.calibrate_param(df_cond_good.CTDCOND2,df_cond_good.BTLCOND,df_cond_good.CTDPRS,'C',2,df_cond_good.SSSCC,df_cond_good.btl_fire_num)
            btl_data_all[c2_btl_col], sal_2 = fit_ctd.conductivity_polyfit(btl_data_all[c2_btl_col], btl_data_all[t2_btl_col], btl_data_all[p_btl_col], coef_cond_sec)
            time_data_all[c2_col],sal2 = fit_ctd.conductivity_polyfit(time_data_all[c2_col], time_data_all[t2_col], time_data_all[p_col], coef_cond_sec)
        
    
        qual_flag_cond = process_ctd.combine_quality_flags([df_ques_c1, df_ques_c2, df_ques_refc])
    ###########################################################################
#
#    ## Oxygen Calibration 
    
    # Calculate Sigma
    btl_data_all['sigma_btl'] = oxy_fitting.sigma_from_CTD(btl_data_all[sal_btl_col], btl_data_all[t_btl_col], btl_data_all[p_btl_col], btl_data_all[lon_btl_col], btl_data_all[lat_btl_col])
    time_data_all['sigma_ctd'] = oxy_fitting.sigma_from_CTD(time_data_all[sal_col], time_data_all[t_col], time_data_all[p_col], time_data_all[lon_col], time_data_all[lat_col])

    btl_data_all[oxy_btl_col] = oxy_fitting.calculate_bottle_oxygen(ssscc, btl_data_all['SSSCC'], btl_data_all['TITR_VOL'], btl_data_all['TITR_TEMP'], btl_data_all['FLASKNO'])
    btl_data_all[oxy_btl_col] = oxy_fitting.oxy_ml_to_umolkg(btl_data_all[oxy_btl_col], btl_data_all['sigma_btl'])

    # Calculate SA and PT
    btl_data_all['SA'] = gsw.SA_from_SP(btl_data_all[sal_btl_col], btl_data_all[p_btl_col], btl_data_all[lon_btl_col], btl_data_all[lat_btl_col])
    btl_data_all['PT'] = gsw.pt0_from_t(btl_data_all['SA'], btl_data_all[t_btl_col], btl_data_all[p_btl_col])

    time_data_all['SA'] = gsw.SA_from_SP(time_data_all[sal_col], time_data_all[p_col], time_data_all[lon_col], time_data_all[lat_col])
    time_data_all['PT'] = gsw.pt0_from_t(time_data_all['SA'], time_data_all[t_col], time_data_all[p_col])

    # Calculate OS in µmol/kg

    btl_data_all['OS_btl'] = oxy_fitting.os_umol_kg(btl_data_all['SA'], btl_data_all['PT'])
    time_data_all['OS_ctd'] = oxy_fitting.os_umol_kg(time_data_all['SA'], time_data_all['PT'])   

    oxy_df = pd.DataFrame()
    coef_dict = {}
    for station in ssscc:
        btl_data = btl_data_all[btl_data_all['SSSCC'] == station]
        time_data = time_data_all[time_data_all['SSSCC'] == station]
    
        hex_file = hex_prefix + station + hex_postfix
        xml_file = xml_prefix + station + xml_postfix
        coef0 = oxy_fitting.get_SB_coef(hex_file, xml_file)
        cfw_coef, df = oxy_fitting.oxy_fit(btl_data[p_btl_col], btl_data[oxy_btl_col], 
                                           btl_data['sigma_btl'], time_data['sigma_ctd'], 
                                           time_data['OS_ctd'], time_data[p_col], 
                                           time_data[t_col], time_data[dov_col], 
                                           time_data[time_col], coef0)
        df['SSSCC'] = station
        coef_dict[station] = cfw_coef
        oxy_df = pd.concat([oxy_df,df])
        
        print(station, ' Completed')
        
    coef_df = oxy_fitting.create_coef_df(coef_dict)
    oxy_df = oxy_fitting.flag_oxy_data(oxy_df)
    
    # Merge oxygen fitting DF to btl_data_all
    
    btl_data_all = oxy_fitting.merge_oxy_df(btl_data_all,oxy_df)
    
    # Apply coef to Time Data
    
    time_data_all = oxy_fitting.apply_oxygen_coef_ctd(time_data_all, coef_df, ssscc)



################ Clean and export data #######################
    btl_data_all = process_ctd.merge_cond_flags(btl_data_all,qual_flag_cond)
    btl_data_all = process_ctd.merge_refcond_flags(btl_data_all,qual_flag_cond)
    btl_data_all = process_ctd.merged_reftemp_flags(btl_data_all,qual_flag_temp)
### Export Quality Flags
    
    qual_flag_temp.to_csv('data/logs/qual_flag_temp_new.csv',index=False)
    qual_flag_cond.to_csv('data/logs/qual_flag_cond_new.csv',index=False)
    
### Clean up Bottle Data by removing rows with no ctd data
    
    btl_data_all =  btl_data_all.dropna(subset=btl_cols)

### Add DATE and TIME
    
    btl_data_all['DATE'] = ''
    btl_data_all['TIME'] = ''
    for station in ssscc:
        df = btl_data_all.loc[btl_data_all['SSSCC'] == station].copy()
        btl_data_all.loc[btl_data_all['SSSCC'] == station] = process_ctd.get_btl_time(df,'btl_fire_num',time_col)

### Create CT Files and HY files
    
    
    process_ctd.export_btl_data(btl_data_all,expocode,sectionID,expocode)
    process_ctd.export_time_data(time_data_all,ssscc,int(sample_rate),int(search_time),expocode,sectionID,ctd,p_column_names,p_column_units)
    
    
def main(argv):
    '''Run everything.
    '''
    process_all_new()

if __name__ == '__main__':
    main(sys.argv[1:])
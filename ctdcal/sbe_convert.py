#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 10:53:36 2018

@author: k3jackson
"""

import os
import sys
sys.path.append('ctdcal/')
import settings
import ctdcal.process_ctd as process_ctd
import pandas as pd
import ctdcal.report_ctd as report_ctd
import ctdcal.convert as cnv
import ctdcal.process_bottle as btl
import ctdcal.sbe_reader as sbe_reader
import warnings
#import ctd
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Get information from settings.py file

expocode = settings.cruise['expocode']
sectionID = settings.cruise['sectionid']
raw_directory = settings.ctd_processing_dir['raw_data_directory']
time_directory = settings.ctd_processing_dir['time_data_directory']
pressure_directory = settings.ctd_processing_dir['pressure_data_directory']
converted_directory = settings.ctd_processing_dir['converted_directory']
oxygen_directory = settings.ctd_processing_dir['oxygen_directory']
btl_directory = settings.ctd_processing_dir['bottle_directory']
o2flask_file = settings.ctd_processing_dir['o2flask_file']
log_directory = settings.ctd_processing_dir['log_directory']
p_log_file = settings.ctd_processing_dir['pressure_log']
hex_prefix = settings.ctd_processing_dir['hex_prefix']
hex_postfix = settings.ctd_processing_dir['hex_postfix']
xml_prefix = settings.ctd_processing_dir['xml_prefix']
xml_postfix = settings.ctd_processing_dir['xml_postfix']
tc1_align = settings.ctd_processing_constants['TC_primary_align']
tc2_align = settings.ctd_processing_constants['TC_secondary_align']
do_align = settings.ctd_processing_constants['DO_align']
input_parameters = settings.input_array
time_column_names = settings.time_series_output['column_name']
time_column_units = settings.time_series_output['column_units']
time_column_data = settings.time_series_output['data_output']
time_column_format = settings.time_series_output['format']
conductivity_startup = settings.ctd_processing_constants['conductivity_start']

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
dopl_col = settings.ctd_inputs['dopl']
lat_col = settings.ctd_inputs['lat']
lon_col = settings.ctd_inputs['lon']
time_col = settings.ctd_inputs['scan_datetime']
alt_col = settings.ctd_inputs['alt']

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
btl_fire_num = settings.bottle_inputs['btl_num']

# CTD Information    
sample_rate = settings.ctd_processing_constants['sample_rate']
search_time = settings.ctd_processing_constants['roll_filter_time']
ctd = settings.ctd_processing_constants['ctd_serial']


btl_data_prefix = 'data/bottle/'
btl_data_postfix = '_btl_mean.pkl'
time_data_prefix = 'data/time/'
time_data_postfix = '_time.pkl'
p_log_file = 'data/logs/ondeck_pressure.csv'
cnv_dir = 'data/cnv/'
def sbe_metadata(ssscc):
    ssscc_file = converted_directory + ssscc + '.pkl'
    raw_data = pd.read_pickle(ssscc_file)
    raw_data = process_ctd.ondeck_pressure(ssscc, p_col, c1_col, c2_col, time_col, raw_data, float(conductivity_startup), log_directory+'ondeck_pressure.csv')
    
    if not c1_col in raw_data.keys():
        print('c1_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c1_col, float(tc1_align))

    if not c2_col in raw_data.keys():
        print('c2_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c2_col, float(tc2_align))

    if not dopl_col in raw_data.keys():
        print('do_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, dopl_col, float(do_align))
    
    # Filter data
    filter_data = process_ctd.raw_ctd_filter(raw_data, 'triangle', 24, input_parameters)
    
    stime, etime, btime, startP, maxP, btm_lat, btm_lon, btm_alt, cast_data = process_ctd.cast_details(ssscc, log_directory+'cast_details.csv', p_col, time_col, lat_col, lon_col, alt_col, filter_data)
    report_ctd.report_time_series_data(ssscc, time_directory, expocode, time_column_names, time_column_units, time_column_data, time_column_format, cast_data)
    return

def sbe_metadata_2(ssscc):
    input_parameters = ['CTDOXYVOLTS']
    ssscc_file = converted_directory + ssscc + '.pkl'
    raw_data = pd.read_pickle(ssscc_file)
    raw_data = process_ctd.ondeck_pressure(ssscc, p_col, c1_col, c2_col, time_col, raw_data, float(conductivity_startup), log_directory+'ondeck_pressure.csv')
    
    if not c1_col in raw_data.keys():
        print('c1_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c1_col, float(tc1_align))

    if not c2_col in raw_data.keys():
        print('c2_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c2_col, float(tc2_align))

    if not dopl_col in raw_data.keys():
        print('do_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, dopl_col, float(do_align))
    
    # Filter data
    filter_data = process_ctd.raw_ctd_filter(raw_data, 'triangle', 24, input_parameters)
    
    stime, etime, btime, startP, maxP, btm_lat, btm_lon, btm_alt, cast_data = process_ctd.cast_details(ssscc, log_directory+'cast_details.csv', p_col, time_col, lat_col, lon_col, alt_col, filter_data)
    report_ctd.report_time_series_data(ssscc, time_directory, expocode, time_column_names, time_column_units, time_column_data, time_column_format, cast_data)
    return

def process_bottle(ssscc):
    
    cnv_file = converted_directory + ssscc + '.pkl'
    imported_df = cnv.importConvertedFile(cnv_file, False)
    bottle_df = btl.retrieveBottleData(imported_df, False)
    if bottle_df.empty:
        #errPrint("Bottle fire data not found in converted file")
        sys.exit(1)
    else:
        total_bottles_fired = bottle_df[btl_fire_num].iloc[-1]
        bottle_num = 1

        while bottle_num <= total_bottles_fired:
            bottle_num += 1

    mean_df = btl.bottle_mean(bottle_df)
    
    meanfilePath = btl_data_prefix + ssscc + btl_data_postfix
    cnv.saveConvertedDataToFile(mean_df, meanfilePath)
    
    return 

def convert_sbe(ssscc, hexFile, xmlFile, outdir=converted_directory):
    import ctd
    sbeReader = sbe_reader.SBEReader.from_paths(hexFile, xmlFile)
    # Build Output Directory exists
    if outdir:
        if os.path.isdir(outdir):
            outputDir = outdir
        else:
            
            try:
                os.mkdir(outdir)
            except:
                print('ERROR: Could not create output directory:', outdir)
                sys.exit(1)
            else:
                outputDir = outdir
                
    converted_df = cnv.convertFromSBEReader(sbeReader, False)
    ### Test pickle file conversion
    pickle_file_name = ssscc + '.pkl'
    pickle_file_path = os.path.join(outputDir, pickle_file_name)
    ### Added to handle I06 sensor issues, use this if you are using seabird postprocessing PO values with ODF postprocessing
#     file = 'data/cnv/' + ssscc + '_WE_corr_tri24.cnv'
#     sbe_df = ctd.from_cnv('data/cnv/' + ssscc + '_WE_corr_tri24.cnv')
#     converted_df['CTDTMP1'] = sbe_df['t090C']
#     converted_df['CTDTMP2'] = sbe_df['t190C']
#     converted_df['CTDCOND1'] = sbe_df['c0mS/cm']
#     converted_df['CTDCOND2'] = sbe_df['c1mS/cm']
#     converted_df['CTDPRS'] = sbe_df['Pressure [dbar]']
    
    converted_df.to_pickle(pickle_file_path) 
    return
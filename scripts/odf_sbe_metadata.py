#! /usr/bin/env python
import sys
import os

import ctdcal.sbe_reader as sbe_reader
import numpy as np
import pandas as pd
import ctdcal.convert as cnv
import ctdcal.process_ctd as process_ctd
import ctdcal.report_ctd as report_ctd
import pickle
from pathlib import PurePosixPath
import xarray as xr

#remove and streamline imports below later
import argparse
import configparser
#import matplotlib.pyplot as plt

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

#File extension to use for raw output
RAW_SUFFIX = '_raw'

#File extension to use for converted output
CONVERTED_SUFFIX = '_cnv'

PKL_EXT = 'pkl'


# these two can be removed in favor of logger
def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# def main(argv):
def main(pkl_path):
    # parser = argparse.ArgumentParser(description='Convert SBE raw data to a converted, csv-formatted text file')
    # parser.add_argument('cnv_pkl', metavar='pkl_file', help='the converted .pkl data file to process')

    # args = parser.parse_args()

    #Import Cruise Configuration File
    iniFile = 'data/ini-files/configuration.ini'
    config = configparser.RawConfigParser()
    config.read(iniFile)

    #Initialise Configuration Parameters
    expocode = config['cruise']['expocode']
    sectionID = config['cruise']['sectionid']
    time_directory = config['ctd_processing']['time_data_directory']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    log_directory = config['ctd_processing']['log_directory']
    conductivity_startup = config['ctd_processing']['conductivity_start']
    tc1_align = config['ctd_processing']['TC_primary_align']
    tc2_align = config['ctd_processing']['TC_secondary_align']
    do_align = config['ctd_processing']['DO_align']
    sample_rate = config['ctd_processing']['sample_rate']
    search_time = config['ctd_processing']['roll_filter_time']
    H1 = config['ctd_processing']['hysteresis_1']
    H2 = config['ctd_processing']['hysteresis_2']
    H3 = config['ctd_processing']['hysteresis_3']
    ctd = config['ctd_processing']['ctd_serial']

    lat_col = config['inputs']['lat']
    lon_col = config['inputs']['lon']
    alt_col = config['inputs']['alt']
    input_parameters = config['inputs']['input_array'].split("\n")
    p_col = config['inputs']['p']
    t1_col = config['inputs']['t1']
    t2_col = config['inputs']['t2']
    c1_col = config['inputs']['c1']
    c2_col = config['inputs']['c2']
    sal_col = config['inputs']['salt']
    dopl_col = config['inputs']['dopl']
    dopkg_col = config['analytical_inputs']['dopkg']
#    xmis_col = config['inputs']['xmis']
#    fluor_col = config['inputs']['fluor']
#    v2_col = config['inputs']['backscatter']
#    v3_col = config['inputs']['rinko_oxy']
#    v4_col = config['inputs']['rinko_tmp']
#    time_zone = config['inputs']['time_zone']
    nmea_time_col = config['inputs']['nmea_datetime']
    scan_time_col = config['inputs']['scan_datetime']

    #time_column_data = config['time_series_output']['data_names'].split(',')
    time_column_data = config['time_series_output']['data_output'].split(',')
    time_column_names = config['time_series_output']['column_name'].split(',')
    time_column_units = config['time_series_output']['column_units'].split(',')
    time_column_format = config['time_series_output']['format']

    #pressure_column_data = config['time_series_output']['data_names'].split(',')
    p_column_data = config['pressure_series_output']['data_output'].split(',')
    p_column_names = config['pressure_series_output']['column_name'].split(',')
    p_column_units = config['pressure_series_output']['column_units'].split(',')
    p_column_format = config['pressure_series_output']['format']
    p_column_qual = config['pressure_series_output']['qual_columns'].split(',')
    p_column_one = list(config['pressure_series_output']['q1_columns'].split(','))

    #
    # filename_ext = os.path.basename(args.cnv_pkl) # original filename with ext
    # filename_base = os.path.splitext(filename_ext)[0] # original filename w/o ext

    # converted_df = pd.read_pickle(args.cnv_pkl)

    filename_ext = PurePosixPath(pkl_path).name  # returns ssscc.pkl
    filename_base = PurePosixPath(pkl_path).stem  # returns ssscc

    converted_df = pd.read_pickle(pkl_path)
    ds_path = pkl_path.split(".pkl")[0] + ".nc"
    try:
        with xr.open_dataset(ds_path, group="raw") as ds:
            raw_ds = ds
    except OSError:
        raw_ds = xr.open_dataset(ds_path)

    # Construct NDarray - fix this serialization asap
    #raw_data = process_ctd.dataToNDarray(convertedfilePath,None,list(converted_df.columns.insert(0,'index')),',',2)
    raw_data = converted_df.to_records()
    #import pdb; pdb.set_trace()

    if nmea_time_col in raw_ds.variables:
        time_col = nmea_time_col
    else:
        time_col = scan_time_col

    #this gets the ondeck pressure for future pressure calibration
    raw_data = process_ctd.ondeck_pressure(filename_base, p_col, c1_col, c2_col, time_col, raw_ds, float(conductivity_startup), log_directory+'ondeck_pressure.csv')

    if c1_col not in raw_data.variables:
        errPrint('c1_col data not found, skipping ctd_align')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c1_col, float(tc1_align))

    if c2_col not in raw_data.variables:
        errPrint('c2_col data not found, skipping ctd_align')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c2_col, float(tc2_align))

    if dopl_col not in raw_data.variables:
        errPrint('do_col data not found, skipping ctd_align')
    else:
        raw_data = process_ctd.ctd_align(raw_data, dopl_col, float(do_align))
        #hysteresis_matrix = process_ctd.hysteresis_correction(float(H1),float(H2), float(H3), raw_matrix)

    # Filter data
    filter_data = process_ctd.raw_ctd_filter(raw_data, 'triangle', 24, input_parameters)

    # Cast Details
    cast_data = process_ctd.cast_details(filename_base, log_directory+'cast_details.csv', p_col, time_col, lat_col, lon_col, alt_col, filter_data)

    # Write out time data
    # TODO: appending to data/converted/SSSCC.nc with group="time" doesn't work...
    cast_data.to_netcdf(time_directory + filename_base + "_time.nc")

    ###
    # MK 05/29/20: Is it necessary to do this step?
    # Raw data are being sequenced, flagged, and exported to _ct1 files
    # They will be overwritten when the code runs all the way through
    # and outputs calibrated data...
    ###
    # Pressure Sequence
    # pressure_seq_data = process_ctd.pressure_sequence(cast_data,p_col,2.0,stime,startP,'down',int(sample_rate),int(search_time))
    # Convert dissolved oxygen from ml/l to umol/kg
    # dopkg = process_ctd.o2pl2pkg(p_col, t1_col, sal_col, dopl_col, dopkg_col, lat_col, lon_col, pressure_seq_data)

    # Add quality codes to data
    # qual_pseq_data = process_ctd.ctd_quality_codes(None, None, None, False, p_column_qual, p_column_one, pressure_seq_data)

    # Write time data to file
    # depth = -999
    # report_ctd.report_pressure_series_data(filename_base, expocode, sectionID, btime, btm_lat, btm_lon, depth, btm_alt, ctd, pressure_directory, p_column_names, p_column_units, p_column_data, qual_pseq_data, dopkg, pressure_seq_data)

    debugPrint('Done!')


# if __name__ == '__main__':
#     main(sys.argv[1:])

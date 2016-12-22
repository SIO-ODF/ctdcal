#! /usr/bin/env python
import sys
import os
import argparse
import libODF_sbe_reader as sbe_reader
import numpy as np
import pandas as pd
import libODF_convert as cnv
import libODF_process_ctd as process_ctd
import libODF_report_ctd as report_ctd
import configparser
import matplotlib.pyplot as plt

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

#File extension to use for raw output
RAW_SUFFIX = '_raw'
 
#File extension to use for converted output
CONVERTED_SUFFIX = '_cnv'


def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# -------------------------------------------------------------------------------------
# Main function of the script should it be run as a stand-alone utility.
# -------------------------------------------------------------------------------------
def main(argv):

    parser = argparse.ArgumentParser(description='Convert SBE raw data to a converted, csv-formatted text file')
    parser.add_argument('hexFile', metavar='hex_file', help='the .hex data file to process')
    parser.add_argument('xmlconFile', metavar='XMLCON_file', help='the .XMLCON data file to process')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # raw output
    parser.add_argument('-r', '--raw', action='store_true', help='return the raw data values')

    # output directory
    parser.add_argument('-o', metavar='dest_dir', dest='outDir', help='location to save output files')

    # Process Command-line args
    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify hex file exists
    if not os.path.isfile(args.hexFile):
        errPrint('ERROR: Input hex file:', args.hexFile, 'not found\n')
        sys.exit(1)

    # Verify xmlcon file exists
    if not os.path.isfile(args.xmlconFile):
        errPrint('ERROR: Input xmlcon file:', args.xmlconFile, 'not found\n')
        sys.exit(1)

    # Set the default output directory to be the same directory as the hex file
    outputDir = os.path.dirname(args.hexFile)

    # Used later for building output file names
    filename_ext = os.path.basename(args.hexFile) # original filename with ext
    filename_base = os.path.splitext(filename_ext)[0] # original filename w/o ext

    # Parse the input files
    debugPrint("Parsing", args.hexFile, "and", args.xmlconFile + '... ', end='')
    sbeReader = sbe_reader.SBEReader.from_paths(args.hexFile, args.xmlconFile)
    debugPrint("Success!")

    # Build Output Directory exists
    if args.outDir:
        if os.path.isdir(args.outDir):
            outputDir = args.outDir
        else:
            debugPrint("Creating output directory...", end='')
            try:
                os.mkdir(args.outDir)
            except:
                errPrint('ERROR: Could not create output directory:', args.outDir)
                sys.exit(1)
            else:
                outputDir = args.outDir
                debugPrint('Success!')

    # Save the raw scans as csv
    if args.raw:

        debugPrint('Building raw dataset... ', end='')

        # Retrieve parsed scans
        rawData = sbeReader.parsed_scans()

        # Convert raw data to dataframe
        raw_df = pd.DataFrame(rawData)
        raw_df.index.name = 'index'
        raw_df = raw_df.apply(pd.to_numeric, errors="ignore")

        debugPrint('Success!')

        rawfileName = str(filename_base + RAW_SUFFIX + '.' + FILE_EXT)
        rawfilePath = os.path.join(outputDir, rawfileName)

        debugPrint('Saving raw data to:', rawfilePath + '... ', end='')
        try:
            raw_df.to_csv(rawfilePath)
        except:
            errPrint('ERROR: Could not save raw data to file')
        else:
            debugPrint('Success!')


    debugPrint("Converting raw scans to scientific units... ")
    converted_df = cnv.convertFromSBEReader(sbeReader, False)
   
    convertedfileName  = filename_base + CONVERTED_SUFFIX + '.' + FILE_EXT
    convertedfilePath = os.path.join(outputDir, convertedfileName)

    debugPrint('Saving converted data to:', convertedfilePath + '... ', end='')
    if cnv.saveConvertedDataToFile(converted_df, convertedfilePath, False):
        debugPrint('Success!')
    else:
        errPrint('ERROR: Could not save converted data to file')
    
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
    input_parameters = config['analytical_inputs']['input_array'].split("\n")
    p_col = config['analytical_inputs']['p']
    t1_col = config['analytical_inputs']['t1']
    t2_col = config['analytical_inputs']['t2']
    c1_col = config['analytical_inputs']['c1']
    c2_col = config['analytical_inputs']['c2']
    do_col = config['analytical_inputs']['do']
    xmis_col = config['analytical_inputs']['xmis']
    fluor_col = config['analytical_inputs']['fluor']
    time_zone = config['inputs']['time_zone']
    nmea_time_col = config['inputs']['nmea_datetime']
    scan_time_col = config['inputs']['scan_datetime']
    
    #time_column_data = config['time_series_output']['data_names'].split(',')
    time_column_data = config['time_series_output']['data_output']
    time_column_names = config['time_series_output']['column_name'].split(',')
    time_column_units = config['time_series_output']['column_units'].split(',')
    time_column_format = config['time_series_output']['format']

    #pressure_column_data = config['time_series_output']['data_names'].split(',')
    p_column_data = config['pressure_series_output']['data_output']
    p_column_names = config['pressure_series_output']['column_name'].split(',')
    p_column_units = config['pressure_series_output']['column_units'].split(',')
    p_column_format = config['pressure_series_output']['format']
    p_column_qual = config['pressure_series_output']['qual_columns'].split(',')
    p_column_one = list(config['pressure_series_output']['q1_columns'].split(','))

    if nmea_time_col in converted_df.columns:
        time_col = nmea_time_col
    else:
        time_col = scan_time_col

    # Construct NDarray
    raw_data = process_ctd.dataToNDarray(convertedfilePath,None,list(converted_df.columns.insert(0,'index')),',')
    raw_data = process_ctd.ondeck_pressure(filename_base, p_col, c1_col, c2_col, time_col, raw_data, float(conductivity_startup), log_directory+'ondeck_pressure.csv')

    if not c1_col in raw_data.dtype.names:
        errPrint('c1_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c1_col, float(tc1_align))

    if not c2_col in raw_data.dtype.names:
        errPrint('c2_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, c2_col, float(tc2_align))

    if not do_col in raw_data.dtype.names:
        errPrint('do_col data not found, skipping')
    else:
        raw_data = process_ctd.ctd_align(raw_data, do_col, float(do_align))
        #hysteresis_matrix = process_ctd.hysteresis_correction(float(H1),float(H2), float(H3), raw_matrix) 

    # Filter data
    filter_data = process_ctd.raw_ctd_filter(raw_data, 'triangle', 24, input_parameters)

    # Cast Details
    stime, etime, btime, startP, maxP, btm_lat, btm_lon, btm_alt, cast_data = process_ctd.cast_details(filename_base, log_directory+'cast_details.csv', p_col, time_col, lat_col, lon_col, alt_col, filter_data)

    # Write time data to file
    report_ctd.report_time_series_data(filename_base, time_directory, expocode, time_column_names, time_column_units, time_column_data, time_column_format, cast_data)

    # Pressure Sequence
    pressure_seq_data = process_ctd.pressure_sequence(filename_base, p_col, time_col, 2.0, stime, startP, 'down', int(sample_rate), int(search_time), cast_data)

    # Add quality codes to data
    qual_pseq_data = process_ctd.ctd_quality_codes(None, None, None, False, p_column_qual, p_column_one, pressure_seq_data)

    # Write time data to file
    depth = -999
    report_ctd.report_pressure_series_data(filename_base, expocode, sectionID, btime, btm_lat, btm_lon, depth, btm_alt, ctd, pressure_directory, p_column_names, p_column_units, qual_pseq_data, pressure_seq_data)

    debugPrint('Done!')

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

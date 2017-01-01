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
import libODF_fit_ctd as fit_ctd
import configparser
import matplotlib.pyplot as plt

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

#File extension to use for raw output
BTL_SUFFIX = '_btl'
 
#File extension to use for raw output
FIT_SUFFIX = '_fit'
 
#File extension to use for raw output
MEAN_SUFFIX = '_mean'
 
#File extension to use for raw output
TIME_SUFFIX = '_time'
 
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
    #parser.add_argument('hexFile', metavar='hex_file', help='the .hex data file to process')
    parser.add_argument('timeFile', metavar='time_file', help='the .csv data file to fit by bottle data')
    #parser.add_argument('xmlconFile', metavar='XMLCON_file', help='the .XMLCON data file to process')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # raw output
    #parser.add_argument('-r', '--raw', action='store_true', help='return the raw data values')
    parser.add_argument('-oxy', '--oxygen', type=argparse.FileType('r'), help='return the oxygen data file')

    # output directory
    #parser.add_argument('-o', metavar='dest_dir', dest='outDir', help='location to save output files')

    # Process Command-line args
    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify hex file exists
    if not os.path.isfile(args.timeFile):
        errPrint('ERROR: Input time dependent .csv file:', args.timeFile, 'not found\n')
        sys.exit(1)

    # Verify xmlcon file exists
    #if not os.path.isfile(args.xmlconFile):
    #    errPrint('ERROR: Input xmlcon file:', args.xmlconFile, 'not found\n')
    #    sys.exit(1)

    # Used later for building output file names
    filename_ext = os.path.basename(args.timeFile) # original filename with ext
    filename_base = os.path.splitext(filename_ext)[0] # original filename w/o ext
   
    if '_' in filename_base:
        filename_base = filename_base.split('_')[0]

    # Parse the input files
    #debugPrint("Parsing", args.timeFile, "and", args.xmlconFile + '... ', end='')
    #sbeReader = sbe_reader.SBEReader.from_paths(args.timeFile, args.xmlconFile)
    #debugPrint("Success!")

    # Build Output Directory exists
    #if args.outDir:
    #    if os.path.isdir(args.outDir):
    #        outputDir = args.outDir
    #    else:
    #        debugPrint("Creating output directory...", end='')
    #        try:
    #            os.mkdir(args.outDir)
    #        except:
    #            errPrint('ERROR: Could not create output directory:', args.outDir)
    #            sys.exit(1)
    #        else:
    #            outputDir = args.outDir
    #            debugPrint('Success!')

    # Save the raw scans as csv
    #if args.raw:

    #    debugPrint('Building raw dataset... ', end='')

        # Retrieve parsed scans
    #    rawData = sbeReader.parsed_scans()
        # Retrieve parsed scans

        # Convert raw data to dataframe
    #    raw_df = pd.DataFrame(rawData)
    #    raw_df.index.name = 'index'
    #    raw_df = raw_df.apply(pd.to_numeric, errors="ignore")

    #    debugPrint('Success!')

    #    rawfileName = str(filename_base + FIT_SUFFIX + '.' + FILE_EXT)
    #    rawfilePath = os.path.join(outputDir, rawfileName)

    #    debugPrint('Saving raw data to:', rawfilePath + '... ', end='')
    #    try:
    #        raw_df.to_csv(rawfilePath)
    #    except:
    #        errPrint('ERROR: Could not save raw data to file')
    #    else:
    #        debugPrint('Success!')

    #debugPrint("Converting raw scans to scientific units... ")
    #converted_df = cnv.convertFromSBEReader(sbeReader, False)
   
    #convertedfileName  = filename_base + CONVERTED_SUFFIX + '.' + FILE_EXT
    #convertedfilePath = os.path.join(outputDir, convertedfileName)

    #debugPrint('Saving converted data to:', convertedfilePath + '... ', end='')
    #if cnv.saveConvertedDataToFile(converted_df, convertedfilePath, False):
    #    debugPrint('Success!')
    #else:
    #    errPrint('ERROR: Could not save converted data to file')
    
    #Import Cruise Configuration File 
    iniFile = 'data/ini-files/configuration.ini' 
    config = configparser.RawConfigParser()
    config.read(iniFile)

    #Initialise Configuration Parameters
    expocode = config['cruise']['expocode']
    sectionID = config['cruise']['sectionid']
    time_directory = config['ctd_processing']['time_data_directory']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    oxygen_directory = config['ctd_processing']['oxygen_directory']
    btl_directory = config['ctd_processing']['bottle_directory']
    o2flask_file = config['ctd_processing']['o2flask_file']
    log_directory = config['ctd_processing']['log_directory']
    sample_rate = config['ctd_processing']['sample_rate']
    search_time = config['ctd_processing']['roll_filter_time']
    ctd = config['ctd_processing']['ctd_serial']

    #alt_col = config['inputs']['alt']
    #input_parameters = config['analytical_inputs']['input_array'].split("\n")
    time_zone = config['inputs']['time_zone']
    p_col = config['analytical_inputs']['p']
    t1_col = config['analytical_inputs']['t1']
    t2_col = config['analytical_inputs']['t2']
    c1_col = config['analytical_inputs']['c1']
    c2_col = config['analytical_inputs']['c2']
    sal_btl_col = config['inputs']['salt']
    sal_col = config['analytical_inputs']['salt']
    dopl_col = config['analytical_inputs']['dopl']
    dopkg_col = config['analytical_inputs']['dopkg']
    xmis_col = config['analytical_inputs']['xmis']
    fluor_col = config['analytical_inputs']['fluor']
    timedate = config['analytical_inputs']['datetime']
    lat_col = config['analytical_inputs']['lat']
    lon_col = config['analytical_inputs']['lon']
    
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

    timefileName = str(filename_base + TIME_SUFFIX + '.' + FILE_EXT)
    timefilePath = os.path.join(time_directory, timefileName)

    outtimefileName = str(filename_base + FIT_SUFFIX + '.' + FILE_EXT)
    outtimefilePath = os.path.join(time_directory, outtimefileName)

    pressfileName = str(filename_base + FIT_SUFFIX + '.' + FILE_EXT)
    pressfilePath = os.path.join(pressure_directory, pressfileName)

    btlfileName = str(filename_base + BTL_SUFFIX + MEAN_SUFFIX + '.' + FILE_EXT)
    btlfilePath = os.path.join(btl_directory, btlfileName)

    if args.oxygen:
        o2pkg_btl, o2pl_btl = fit_ctd.o2_calc(o2flask_file,args.oxygen.name,sal_btl_col,btlfilePath)

    # Construct NDarray
    time_data = process_ctd.dataToNDarray(timefilePath,None,True,',',1)
    print(time_data.dtype.names)

    # Get analytical Oxygen data
     

    # Write time data to file
    #report_ctd.report_time_series_data(filename_base, time_directory, expocode, time_column_names, time_column_units, time_column_data, time_column_format, time_data)

    # Pressure Sequence
    #pressure_seq_data = process_ctd.pressure_sequence(filename_base, p_col, time_col, 2.0, stime, startP, 'down', int(sample_rate), int(search_time), cast_data)

    # Convert dissolved oxygen from ml/l to umol/kg
    #dopkg = process_ctd.o2pl2pkg(p_col, t1_col, sal_col, dopl_col, dopkg_col, lat_col, lon_col, pressure_seq_data)

    # Add quality codes to data
    #qual_pseq_data = process_ctd.ctd_quality_codes(None, None, None, False, p_column_qual, p_column_one, pressure_seq_data)

    # Write time data to file
    #depth = -999
    #report_ctd.report_pressure_series_data(filename_base, expocode, sectionID, btime, btm_lat, btm_lon, depth, btm_alt, ctd, pressure_directory, p_column_names, p_column_units, qual_pseq_data, dopkg, pressure_seq_data)

    #plt.plot(pressure_seq_data[t1_col], pressure_seq_data[p_col], color='r', label='raw')
    #plt.plot(pressure_seq_data[c1_col], pressure_seq_data[p_col], color='g', label='raw')
    #plt.plot(pressure_seq_data[do_col], pressure_seq_data[p_col], color='b', label='raw')
    #plt.gca().invert_yaxis()
    #plt.axis()
    #plt.show()

    debugPrint('Done!')

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

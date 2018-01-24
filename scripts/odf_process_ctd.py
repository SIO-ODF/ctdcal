#! /usr/bin/env python
import sys
import os
import argparse
import numpy as np
import pandas as pd
import configparser
import ctdcal.process_ctd as process_ctd
import ctdcal.sbe_reader as sbe_reader
import ctdcal.convert as cnv

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

#File extension to use for processed output
PROC_SUFFIX = '_proc'


def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# -------------------------------------------------------------------------------------
# Main function of the script should it be run as a stand-alone utility.
# -------------------------------------------------------------------------------------
def main(argv):

    parser = argparse.ArgumentParser(description='General Utility for processing CTD sensors data from converted, csv-formatted text files')
    parser.add_argument('iniFile', metavar='ini_file', help='the .ini file to use for processing')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # input file
    parser.add_argument('-i', metavar='cnv_file', dest='inFile', help='the converted ctd data to process, overrides the input file defined in the .ini file')

    # output directory
    parser.add_argument('-o', metavar='dest_dir', dest='outDir', help='location to save output files')

    # Process Command-line args
    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify ini file exists
    if not os.path.isfile(args.iniFile):
        errPrint('ERROR: INI file:', args.iniFile, 'not found\n')
        sys.exit(1)
    else:
        #Import Cruise Configuration File
        config = configparser.ConfigParser()
        config.read(args.iniFile)

    # Verify input file exists
    if args.inFile:
        if not os.path.isfile(args.inFile):
            errPrint('ERROR: Input file:', args.inFile, 'not found\n')
            sys.exit(1)


    # Setting the input file
    inputFile = ''

    # Change the output directory based on other command-line arguments
    if args.inFile:
        inputFile = args.inFile
    else:
        try:
            config['ctd_processing']['input_data_file']
        except:
            debugPrint('No input file defined in ini file')
            errPrint('ERROR: No Input file defined\n')
            sys.exit(1)
        else:
            inputFile = config['ctd_processing']['input_data_file']
            if not os.path.isfile(inputFile):
                errPrint('ERROR: Input file:', inputFile, 'not found\n')
                sys.exit(1)

    debugPrint('Input data file set to:', inputFile)

    # Setting the output directory
    outputDir = ''

    # Change the output directory based on other command-line arguments
    if args.outDir:
        debugPrint('Setting output directory to:', args.outDir)
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
    elif args.inFile:
        debugPrint('Setting output directory to:', os.path.dirname(args.inFile))
        outputDir = os.path.dirname(args.inFile) # original filename with ext
    elif config['ctd_processing']['output_data_directory'] != None:
        debugPrint('Setting output directory to:', os.path.dirname(args.inFile))
        outputDir = config['ctd_processing']['output_data_directory']
    else:
        debugPrint('Setting output to current directory')
        outputDir = '.'

    debugPrint("Import converted data to dataframe... ", end='')
    imported_df = cnv.importConvertedFile(inputFile, False)
    debugPrint("Success!")

    # CTD Processing
    debugPrint('Processing data')

    #Import Cruise Configuration File
    config = configparser.ConfigParser()
    config.read(args.iniFile)

    #Initialise Configuration Parameters
    raw_directory = config['ctd_processing']['raw_data_directory']
    time_directory = config['ctd_processing']['time_data_directory']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    tc1_align = config['ctd_processing']['TC_primary_align']
    tc2_align = config['ctd_processing']['TC_secondary_align']
    do_align = config['ctd_processing']['DO_align']
    sample_rate = config['ctd_processing']['sample_rate']
    search_time = config['ctd_processing']['roll_filter_time']
    H1 = config['ctd_processing']['hysteresis_1']
    H2 = config['ctd_processing']['hysteresis_2']
    H3 = config['ctd_processing']['hysteresis_3']

    raw_matrix = process_ctd.dataToMatrix(inputFile,None,list(imported_df.columns.insert(0,'index')),',')
    debugPrint("Martix DTypes:",raw_matrix.dtype)

    data_matrix = process_ctd.ondeck_pressure(raw_matrix, float(config['ctd_processing']['conductivity_start']))

    if not 'c1_Sm' in raw_matrix.dtype.names:
        errPrint('c1_Sm data not found, skipping')
    else:
        align_matrix = process_ctd.ctd_align(raw_matrix,'c1_Sm', float(tc1_align))

    if not 'c2_Sm' in raw_matrix.dtype.names:
        errPrint('c2_Sm data not found, skipping')
    else:
        align_matrix = process_ctd.ctd_align(raw_matrix,'c2_Sm', float(tc2_align))

    if not 'o1_mll' in raw_matrix.dtype.names:
        errPrint('o1_mll data not found, skipping')
    else:
        hysteresis_matrix = process_ctd.hysteresis_correction(float(H1),float(H2), float(H3), raw_matrix)
        align_matrix = process_ctd.ctd_align(hysteresis_matrix,'o1_mll', float(do_align))


    debugPrint('Done!')

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

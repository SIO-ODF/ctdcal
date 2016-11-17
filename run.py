#! /usr/bin/env python
import sys
import os
import argparse
import sbe_reader
import numpy as np
import pandas as pd
import web_viewer as wv
import converter_scaffolding as cnv
import configparser
import process_ctd 
#import converter
#import raw_converter

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

#File extension to use for raw output
RAW_SUFFIX = '_raw'

#File extension to use for converted output
CONVERTED_SUFFIX = '_converted'

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

    parser = argparse.ArgumentParser(description='General Utility for processing SBE Data')
    parser.add_argument('hexFile', metavar='hex File', help='the .hex data file to process')
    parser.add_argument('xmlconFile', metavar='XMLCON File', help='the .XMLCON data file to process')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # raw output
    parser.add_argument('-r', '--raw', action='store_true', help='return the raw data values')

    # processed output
    parser.add_argument('-p', metavar='iniFile', dest='iniFile', help='process the ctd data, requires and addition ini file')

    # output directory
    parser.add_argument('-o', metavar='destDir', dest='outDir', help='location to save output files')

    # web viewer
    parser.add_argument('-w', metavar='destDir', dest='wvDir', nargs='?', const='default', help='build web viewer, optionally specify the destination directory for the webviewer files')

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

    # Verify ini file exists
    if args.iniFile:
        if not os.path.isfile(args.iniFile):
            errPrint('ERROR: Input ini file:', args.iniFile, 'not found\n')
            sys.exit(1)


    # Set the default output directory to be the same directory as the hex file
    outputDir = os.path.dirname(args.hexFile)

    # Used later for building output file names
    filename_ext = os.path.basename(args.hexFile) # original filename with ext
    filename_base = os.path.splitext(filename_ext)[0] # original filename w/o ext

    # Verify Output Directory exists
    if args.outDir:
        if os.path.isdir(args.outDir):
            outputDir = args.outDir
        else:
            errPrint('ERROR: Output directory:', args.outDir, 'is unreachable.\n')
            sys.exit(1)

    # Parse the input files
    debugPrint("Parsing", args.hexFile, "and", args.xmlconFile + '... ', end='')
    sbeReader = sbe_reader.SBEReader.from_paths(args.hexFile, args.xmlconFile)
    debugPrint("Success!")

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

        #debugPrint("Raw Data Types:", raw_df.dtypes)
        #debugPrint("Raw Data:", raw_df.head)

        # Retrieve Config data
        #rawConfig = sbeReader.parsed_config()

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
    converted_df = cnv.convertFromSBEReader(sbeReader, DEBUG)

    #debugPrint('converted_df:\n', converted_df.head())

    convertedfileName  = filename_base + CONVERTED_SUFFIX + '.' + FILE_EXT
    convertedfilePath = os.path.join(outputDir, convertedfileName)

    debugPrint('Saving converted data to:', convertedfilePath + '... ', end='')
    if cnv.saveConvertedDataToFile(converted_df, convertedfilePath, DEBUG):
        debugPrint('Success!')
    else:
        errPrint('ERROR: Could not save converted data to file')

    #Import Cruise Configuration File 
    config = configparser.ConfigParser()
    config.read('configuration.ini')

    #Initialise Configuration Parameters
    raw_directory = config['ctd_processing']['raw_data_directory']
    time_directory = config['ctd_processing']['time_data_directory']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    infile = convertedfilePath 
    tc1_align = config['ctd_processing']['TC_primary_align']
    tc2_align = config['ctd_processing']['TC_secondary_align']
    do_align = config['ctd_processing']['DO_align']
    sample_rate = config['ctd_processing']['sample_rate']
    search_time = config['ctd_processing']['roll_filter_time']
    H1 = config['ctd_processing']['hysteresis_1']
    H2 = config['ctd_processing']['hysteresis_2']
    H3 = config['ctd_processing']['hysteresis_3']

    # CTD Processing
    debugPrint('Processing data')
    raw_matrix = process_ctd.dataToMatrix(infile,None,list(converted_df.columns.insert(0,'index')),',') 
    print(raw_matrix.dtype)
    
    align_matrix = process_ctd.ctd_align(raw_matrix,'c1_Sm', float(tc1_align))
    align_matrix = process_ctd.ctd_align(raw_matrix,'c2_Sm', float(tc2_align))
    align_matrix = process_ctd.ctd_align(raw_matrix,'o1_mll', float(do_align))

    #data_matrix = process_ctd.ondeck_pressure(align_matrix, float(config['ctd_processing']['conductivity_start']))
    #hysteresis_matrix = process_ctd.hysteresis_correction(float(H1),float(H2), float(H3), raw_matrix) 
    
    

   
    # ---------------------------------------------------------------- #
    #
    #  Here's where you can make a call to any CTD processing classes
    #  of functions.
    #
    #  The raw scans are available at: raw_df, (pandas dataframe)
    #  The sensor config is available as: rawConfig, (object)
    #  The cooked data is available at: converted_df, (pandas dataframe)
    #  The ini file needed for processing is available at: args.iniFile
    #
    # -----------------------------------------------------------------#

    if args.wvDir:
        debugPrint('Building webviewer')
        webView = wv.WebViewer()
        if DEBUG:
            webView._setDebug()

        debugPrint('\tSetting webviewer directory to: ', end='')
        if args.wvDir == 'default':
            wvDir = os.path.join(os.path.dirname(args.hexFile),'webviewer')
            debugPrint(wvDir + '... ', end='')
        else:
            wvDir = args.wvDir
            debugPrint(wvDir + '... ', end='')

        try:
            os.mkdir(wvDir, 0o755);
        except FileExistsError as exception:
            pass
        except PermissionError as exception:
            errPrint('ERROR: do not have write permission at destination directory')
            sys.exit(1)

        debugPrint('Success!')

        webView._setWebviewerFolder(wvDir)

        debugPrint('\tCreating webviewer directory structure... ', end='')

        if webView._buildScaffolding():
            debugPrint('Success!')
        else:
            debugPrint('ERROR!')
            errPrint('Unable to created webviewer directories')
            sys.exit(1)

        debugPrint("Building data file for webviewer... ", end='')

        output = webView._buildDataFromDF(converted_df, args.hexFile)
        debugPrint('Success!')

        debugPrint('Saving webviewer data to file... ', end='')
        webView._saveData(output)
        debugPrint('Success!')
    debugPrint('Done!')

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

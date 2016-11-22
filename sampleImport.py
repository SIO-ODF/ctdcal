import sys
import os
import argparse
import pandas as pd
import libODF_convert as cnv

DEBUG = False

def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# -------------------------------------------------------------------------------------
# Main function of the script should it be run as a stand-alone utility.
# -------------------------------------------------------------------------------------
def main(argv):

    parser = argparse.ArgumentParser(description='Sample Script for importing converted SBE Data')
    parser.add_argument('convertedFile', metavar='converted_File', help='the converted data file to process')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify hex file exists
    if not os.path.isfile(args.convertedFile):
        errPrint('ERROR: Input converted file:', args.convertedFile, 'not found\n')
        sys.exit(1)

    debugPrint("Import converted data to dataframe... ", end='')
    imported_df = cnv.importConvertedFile(args.convertedFile, DEBUG)
    debugPrint("Success!")

    errPrint(imported_df.head())
    errPrint(imported_df.dtypes)

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

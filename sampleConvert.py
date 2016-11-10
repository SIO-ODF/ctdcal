import sys
import os
import argparse
import pandas as pd
import converter_scaffolding as cnv

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

    parser = argparse.ArgumentParser(description='Sample Script for converting raw SBE Data')
    parser.add_argument('hexFile', metavar='hex File', help='the .hex data file to process')
    parser.add_argument('xmlconFile', metavar='XMLCON File', help='the .XMLCON data file to process')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

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

    debugPrint("Converting raw scans to scientific units... ", end='')
    converted_df = cnv.convertFromFiles(args.hexFile, args.xmlconFile, DEBUG)
    debugPrint("Success!")

    errPrint(converted_df.head())
    errPrint(converted_df.dtypes)

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

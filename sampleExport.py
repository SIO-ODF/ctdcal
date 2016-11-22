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

    parser = argparse.ArgumentParser(description='Sample Script for importing converted, csv-formatted data to a pandas dataframe and then exporting the dataframe as a new csv-formatted file')
    parser.add_argument('convertedFile', metavar='cnv_file', help='the converted data file to import to a dataframe')
    parser.add_argument('outputFile', metavar='output_file', help='the filename to export the dataframe to')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # overwrite existing output file
    parser.add_argument('-w', '--overwrite', action='store_true', help='overwrite the pre-existing output file if found')

    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify input file exists
    if not os.path.isfile(args.convertedFile):
        errPrint('ERROR: Input file:', args.convertedFile, 'not found\n')
        sys.exit(1)

    # Verify input file exists
    if not args.overwrite and os.path.isfile(args.outputFile):
        errPrint('ERROR: Output file:', args.outputFile, 'already exists.  Use -w flag to overwrite the pre-existing file\n')
        sys.exit(1)

    debugPrint("Import converted data to dataframe... ", end='')
    imported_df = cnv.importConvertedFile(args.convertedFile, False)
    debugPrint("Success!")

    #debugPrint(imported_df.head())
    #debugPrint(imported_df.dtypes)

    debugPrint('Saving output data to:', args.outputFile + '... ', end='')
    if cnv.saveConvertedDataToFile(imported_df, args.outputFile, False):
        debugPrint('Success!')
    else:
        errPrint('ERROR: Could not save converted data to file')

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

import sys
import os
import argparse
import pandas as pd
import converter_scaffolding as cnv
import bottle_lib as btl

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

#File extension to use for converted output
CONVERTED_SUFFIX = '_converted'

#File extension to use for processed output
BOTTLE_SUFFIX = '_btl'

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

    parser = argparse.ArgumentParser(description='Build bottle file from converted file')
    parser.add_argument('convertedFile', metavar='converted File', help='the converted data file to process')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # output directory
    parser.add_argument('-o', metavar='destDir', dest='outDir', help='location to save output files')

    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify converted exists
    if not os.path.isfile(args.convertedFile):
        errPrint('ERROR: Input converted file:', args.convertedFile, 'not found\n')
        sys.exit(1)

    # Set the default output directory to be the same directory as the hex file
    outputDir = os.path.dirname(args.convertedFile)

    # Used later for building output file names
    filename_ext = os.path.basename(args.convertedFile) # original filename with ext
    filename_base = os.path.splitext(filename_ext)[0] # original filename w/o ext

    # Verify Output Directory exists
    if args.outDir:
        if os.path.isdir(args.outDir):
            outputDir = args.outDir
        else:
            errPrint('ERROR: Output directory:', args.outDir, 'is unreachable.\n')
            sys.exit(1)

    debugPrint("Import converted data to dataframe... ", end='')
    imported_df = cnv.importConvertedFile(args.convertedFile, False)
    debugPrint("Success!")

    debugPrint(imported_df.head())
    #debugPrint(imported_df.dtypes)

    # Retrieve the rows from the imported dataframe where the btl_fire_bool column == True
    # Returns a new dataframe
    bottle_df = btl.retrieveBottleData(imported_df, DEBUG)

    if bottle_df.empty:
        errPrint("Bottle fire data not found in converted file")
        sys.exit(1)
    else:
        total_bottles_fired = bottle_df["bottle_fire_num"].iloc[-1]
        debugPrint(total_bottles_fired, "bottle fire(s) detected")
        bottle_num = 1

        while bottle_num <= total_bottles_fired:
            debugPrint('Bottle:', bottle_num)
            debugPrint(bottle_df.loc[bottle_df['bottle_fire_num'] == bottle_num].head())
            bottle_num += 1

    # Build the filename for the bottle fire data
    bottlefileName  = filename_base.replace(CONVERTED_SUFFIX,'') + BOTTLE_SUFFIX + '.' + FILE_EXT
    bottlefilePath = os.path.join(outputDir, bottlefileName)

    # Save the bottle fire dataframe to file.
    debugPrint('Saving bottle data to:', bottlefilePath + '... ', end='')
    if not cnv.saveConvertedDataToFile(bottle_df, bottlefilePath, DEBUG):
        errPrint('ERROR: Could not save bottle fire data to file')
    else:
        debugPrint('Success!')


# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

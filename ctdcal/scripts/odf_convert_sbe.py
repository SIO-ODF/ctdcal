#! /usr/bin/env python
# remove and streamline imports below later
import argparse
import os
import sys

from . import convert as cnv
from . import sbe_reader as sbe_reader

DEBUG = False

# File extension to use for output files (csv-formatted)
FILE_EXT = "csv"

# File extension to use for raw output
RAW_SUFFIX = "_raw"

# File extension to use for converted output
CONVERTED_SUFFIX = "_cnv"

PKL_EXT = "pkl"


def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)


def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# -------------------------------------------------------------------------------------
# Main function of the script should it be run as a stand-alone utility.
# -------------------------------------------------------------------------------------
def main(argv):
    # MK: depreciated 04/27/20, use ctdcal.convert.hex_to_ctd instead

    parser = argparse.ArgumentParser(
        description="Convert SBE raw data to a converted, csv-formatted text file"
    )
    parser.add_argument(
        "hexFile", metavar="hex_file", help="the .hex data file to process"
    )
    parser.add_argument(
        "xmlconFile", metavar="XMLCON_file", help="the .XMLCON data file to process"
    )

    # debug messages
    # parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # raw output
    # parser.add_argument('-r', '--raw', action='store_true', help='return the raw data values')

    # output directory
    parser.add_argument(
        "-o", metavar="dest_dir", dest="outDir", help="location to save output files"
    )

    # Process Command-line args
    args = parser.parse_args()
    # if args.debug:
    #     global DEBUG
    #     DEBUG = True
    #     debugPrint("Running in debug mode")

    # Verify hex file exists
    if not os.path.isfile(args.hexFile):
        errPrint("ERROR: Input hex file:", args.hexFile, "not found\n")
        sys.exit(1)

    # Verify xmlcon file exists
    if not os.path.isfile(args.xmlconFile):
        errPrint("ERROR: Input xmlcon file:", args.xmlconFile, "not found\n")
        sys.exit(1)

    # Set the default output directory to be the same directory as the hex file
    outputDir = os.path.dirname(args.hexFile)

    # Used later for building output file names
    filename_ext = os.path.basename(args.hexFile)  # original filename with ext
    filename_base = os.path.splitext(filename_ext)[0]  # original filename w/o ext

    # Parse the input files
    debugPrint("Parsing", args.hexFile, "and", args.xmlconFile + "... ", end="")
    sbeReader = sbe_reader.SBEReader.from_paths(args.hexFile, args.xmlconFile)
    debugPrint("Success!")

    # Build Output Directory exists
    if args.outDir:
        if os.path.isdir(args.outDir):
            outputDir = args.outDir
        else:
            debugPrint("Creating output directory...", end="")
            try:
                os.mkdir(args.outDir)
            except OSError:
                errPrint("ERROR: Could not create output directory:", args.outDir)
                sys.exit(1)
            else:
                outputDir = args.outDir
                debugPrint("Success!")

    # # Save the raw scans as csv
    # if args.raw:
    #
    #     debugPrint('Building raw dataset... ', end='')
    #
    #     # Retrieve parsed scans
    #     rawData = sbeReader.parsed_scans()
    #
    #     # Convert raw data to dataframe
    #     raw_df = pd.DataFrame(rawData)
    #     raw_df.index.name = 'index'
    #     raw_df = raw_df.apply(pd.to_numeric, errors="ignore")
    #
    #     debugPrint('Success!')
    #
    #     rawfileName = str(filename_base + RAW_SUFFIX + '.' + FILE_EXT)
    #     rawfilePath = os.path.join(outputDir, rawfileName)
    #
    #     debugPrint('Saving raw data to:', rawfilePath + '... ', end='')
    #     try:
    #         raw_df.to_csv(rawfilePath)
    #     except:
    #         errPrint('ERROR: Could not save raw data to file')
    #     else:
    #         debugPrint('Success!')

    # debugPrint("Converting raw scans to scientific units... ")
    converted_df = cnv.convertFromSBEReader(sbeReader, False)

    # convertedfileName  = filename_base + CONVERTED_SUFFIX + '.' + FILE_EXT
    # convertedfilePath = os.path.join(outputDir, convertedfileName)

    ### Test pickle file conversion
    pickle_file_name = filename_base + "." + PKL_EXT
    pickle_file_path = os.path.join(outputDir, pickle_file_name)
    converted_df.to_pickle(pickle_file_path)

    # debugPrint('Saving converted data to:', convertedfilePath + '... ', end='')
    # if cnv.saveConvertedDataToFile(converted_df, convertedfilePath, False):
    #     debugPrint('Success!')
    # else:
    #     errPrint('ERROR: Could not save converted data to file')


# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])

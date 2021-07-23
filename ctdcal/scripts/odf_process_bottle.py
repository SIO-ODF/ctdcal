#! /usr/bin/env python
import argparse
import os
import sys

import pandas as pd

from . import process_bottle as btl

# File extension to use for output files (csv-formatted)
FILE_EXT = "csv"
EXT_PKL = "pkl"

# File extension to use for converted output
CONVERTED_SUFFIX = "_cnv"

# File extension to use for processed output
BOTTLE_SUFFIX = "_btl"
MEAN_SUFFIX = "_mean"
MEDIAN_SUFFIX = "_median"

# Column name for the bottle fire number
BOTTLE_FIRE_NUM_COL = "btl_fire_num"

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
    # MK: depreciated 04/27/20, use ctdcal.convert.make_btl_mean instead

    parser = argparse.ArgumentParser(
        description="Process bottle fire data from a converted, csv-formatted text file"
    )
    parser.add_argument(
        "cnv_file", metavar="cnv_file", help="the converted data file to process"
    )

    # debug messages
    parser.add_argument(
        "-d", "--debug", action="store_true", help="display debug messages"
    )

    # output directory
    parser.add_argument(
        "-o", metavar="dest_dir", dest="outDir", help="location to save output files"
    )

    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify converted exists
    if not os.path.isfile(args.cnv_file):
        errPrint("ERROR: Input converted file:", args.cnv_file, "not found\n")
        sys.exit(1)

    # Set the default output directory to be the same directory as the hex file
    outputDir = os.path.dirname(args.cnv_file)

    # Used later for building output file names
    filename_ext = os.path.basename(args.cnv_file)  # original filename with ext
    filename_base = os.path.splitext(filename_ext)[0]  # original filename w/o ext

    # Verify Output Directory exists
    if args.outDir:
        if os.path.isdir(args.outDir):
            outputDir = args.outDir
        else:
            errPrint("ERROR: Output directory:", args.outDir, "is unreachable.\n")
            sys.exit(1)

    debugPrint("Import converted data to dataframe... ", end="")
    imported_df = pd.read_pickle(args.cnv_file)
    debugPrint("Success!")

    debugPrint(imported_df.head())
    # debugPrint(imported_df.dtypes)

    # Retrieve the rows from the imported dataframe where the btl_fire_bool column == True
    # Returns a new dataframe
    bottle_df = btl.retrieveBottleData(imported_df)

    if bottle_df.empty:
        errPrint("Bottle fire data not found in converted file")
        sys.exit(1)
    else:
        total_bottles_fired = bottle_df[BOTTLE_FIRE_NUM_COL].iloc[-1]
        debugPrint(total_bottles_fired, "bottle fire(s) detected")
        bottle_num = 1

        while bottle_num <= total_bottles_fired:
            debugPrint("Bottle:", bottle_num)
            debugPrint(
                bottle_df.loc[bottle_df[BOTTLE_FIRE_NUM_COL] == bottle_num].head()
            )
            bottle_num += 1

    # # Build the filename for the bottle fire data
    # bottlefileName  = filename_base.replace(CONVERTED_SUFFIX,'') + BOTTLE_SUFFIX + '.' + FILE_EXT
    # bottlefilePath = os.path.join(outputDir, bottlefileName)
    #
    # # Save the bottle fire dataframe to file.
    # debugPrint('Saving bottle data to:', bottlefilePath + '... ', end='')
    # if not cnv.saveConvertedDataToFile(bottle_df, bottlefilePath, DEBUG):
    #     errPrint('ERROR: Could not save bottle fire data to file')
    # else:
    #     debugPrint('Success!')
    #
    mean_df = btl.bottle_mean(bottle_df)

    # Build the filename for the bottle fire mean data
    # meanfileName  = filename_base.replace(CONVERTED_SUFFIX,'') + BOTTLE_SUFFIX + MEAN_SUFFIX + '.' + FILE_EXT
    meanfileName = (
        filename_base.replace(CONVERTED_SUFFIX, "")
        + BOTTLE_SUFFIX
        + MEAN_SUFFIX
        + "."
        + EXT_PKL
    )
    meanfilePath = os.path.join(outputDir, meanfileName)

    # Save the bottle fire mean dataframe to file.
    debugPrint("Saving mean data to:", meanfilePath + "... ", end="")
    if not mean_df.to_pickle(meanfilePath):
        errPrint("ERROR: Could not save mean fire data to file")
    else:
        debugPrint("Success!")

    #### Statistical analysis later to see if mean or median is better for values

    # median_df = btl.bottle_median(bottle_df)
    #
    # # Build the filename for the bottle fire median data
    # medianfileName  = filename_base.replace(CONVERTED_SUFFIX,'') + BOTTLE_SUFFIX + MEDIAN_SUFFIX + '.' + FILE_EXT
    # medianfilePath = os.path.join(outputDir, medianfileName)
    #
    # # Save the bottle fire median dataframe to file.
    # debugPrint('Saving median data to:', medianfilePath + '... ', end='')
    # if not cnv.saveConvertedDataToFile(median_df, medianfilePath, DEBUG):
    #     errPrint('ERROR: Could not save median fire data to file')
    # else:
    #     debugPrint('Success!')


# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])

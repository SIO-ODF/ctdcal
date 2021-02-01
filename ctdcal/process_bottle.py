"""Library to create SBE .btl equivalent files.
TODO: allow for variable bottle fire scans instead of SBE standard 36
    ex: user doesn't know how to change the config for the cast to add more scans,
    instead does it post-cast?

Joseph Gum SIO/ODF
Nov 7, 2016
"""

import csv
import statistics
import sys

import pandas as pd

BOTTLE_FIRE_COL = "btl_fire"
BOTTLE_FIRE_NUM_COL = "btl_fire_num"

DEBUG = False


def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)


def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# Retrieve the bottle data from a converted file.
def retrieveBottleDataFromFile(converted_file, debug=False):

    converted_df = pd.read_pickle(converted_file)

    return retrieveBottleData(converted_df, debug)


# Retrieve the bottle data from a dataframe created from a converted file.
def retrieveBottleData(converted_df, debug=False):
    if BOTTLE_FIRE_COL in converted_df.columns:
        converted_df[BOTTLE_FIRE_NUM_COL] = (
            (
                (converted_df[BOTTLE_FIRE_COL] is True)
                & (
                    converted_df[BOTTLE_FIRE_COL]
                    != converted_df[BOTTLE_FIRE_COL].shift(1)
                )
            )
            .astype(int)
            .cumsum()
        )
        # converted_df['bottle_fire_num'] = ((converted_df[BOTTLE_FIRE_COL] == False)).astype(int).cumsum()

        return converted_df.loc[converted_df[BOTTLE_FIRE_COL] is True]
        # return converted_df
    else:
        debugPrint("Bottle fire column:", BOTTLE_FIRE_COL, "not found")

    return pd.DataFrame()  # empty dataframe


def bottle_mean(btl_df):
    """Compute the mean for each bottle from a dataframe."""
    btl_max = int(btl_df[BOTTLE_FIRE_NUM_COL].tail(n=1))
    i = 1
    output = pd.DataFrame()
    while i <= btl_max:
        output = pd.concat(
            (
                output,
                btl_df[btl_df[BOTTLE_FIRE_NUM_COL] == i]
                .mean()
                .to_frame(name=i)
                .transpose(),
            )
        )
        i += 1
    return output


def bottle_median(btl_df):
    """Compute the median for each bottle from a dataframe."""
    btl_max = int(btl_df[BOTTLE_FIRE_NUM_COL].tail(n=1))
    i = 1
    output = pd.DataFrame()
    while i <= btl_max:
        output = pd.concat(
            (
                output,
                btl_df[btl_df[BOTTLE_FIRE_NUM_COL] == i]
                .median()
                .to_frame(name=i)
                .transpose(),
            )
        )
        i += 1
    return output

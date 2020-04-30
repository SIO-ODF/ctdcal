import csv
import pandas as pd
import numpy as np
from collections import OrderedDict

# MK: depreciated 04/22/20
# functions have been exported to ctdcal.process_ctd

def reft_loader(ssscc, reft_dir):
    reft_path = reft_dir + ssscc + ".cap"
    with open(reft_path, "r", newline="") as f:
        reftF = csv.reader(
            f, delimiter=" ", quoting=csv.QUOTE_NONE, skipinitialspace="True"
        )
        reftArray = []
        for row in reftF:
            if len(row) != 17:  # skip over 'bad' rows (empty lines, comments, etc.)
                continue
            reftArray.append(row)

    reftDF = pd.DataFrame.from_records(reftArray)
    reftDF = reftDF.replace(  # remove text columns, only need numbers and dates
        to_replace=["bn", "diff", "val", "t90", "="], value=np.nan
    )
    reftDF = reftDF.dropna(axis=1)
    # combine date/time columns (day/month/year/time are read in separately)
    reftDF[1] = reftDF[1] + " " + reftDF[2] + " " + reftDF[3] + " " + reftDF[4]
    reftDF.drop(columns=[2, 3, 4], inplace=True)
    # rename columns and recast datatypes
    columns = OrderedDict(
        [
            ("index_memory", int),
            ("datetime", object),
            ("btl_fire_num", int),
            ("diff", int),
            ("raw_value", float),
            ("T90", float),
        ]
    )
    reftDF.columns = list(columns.keys())  # name columns
    reftDF = reftDF.astype(columns)  # force dtypes
    # assign initial qality flags
    reftDF.loc[:, "REFTMP_FLAG_W"] = 2
    reftDF.loc[abs(reftDF["diff"]) >= 3000, "REFTMP_FLAG_W"] = 3
    # add in STNNBR, CASTNO columns
    # TODO: should these be objects or floats? be consistent!
    # string prob better for other sta/cast formats (names, letters, etc.)
    reftDF["STNNBR"] = ssscc[0:3]
    reftDF["CASTNO"] = ssscc[3:5]
    return reftDF


def process_reft(ssscc, reft_dir):
    try:
        reftDF = reft_loader(ssscc, reft_dir)
        reftDF.to_csv(reft_dir + ssscc + "_reft.csv", index=False)
    except FileNotFoundError:
        print("refT file for cast " + ssscc + " does not exist... skipping")
        return

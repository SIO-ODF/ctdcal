import csv
import os
from collections import OrderedDict

import gsw
import pandas as pd

# MK: depreciated 04/22/20
# functions have been exported to ctdcal.odf_io


def salt_loader(ssscc, salt_dir):
    saltpath = salt_dir + ssscc  # salt files have no extension
    with open(saltpath, newline="") as f:
        saltF = csv.reader(
            f, delimiter=" ", quoting=csv.QUOTE_NONE, skipinitialspace="True"
        )
        saltArray = []
        for row in saltF:
            saltArray.append(row)
        del saltArray[0]  # remove header

    header = OrderedDict(  # having this as a dict streamlines next steps
        [
            ("STNNBR", int),
            ("CASTNO", int),
            ("SAMPNO", int),
            ("BathTEMP", int),
            ("CRavg", float),
            ("autosalSAMPNO", int),
            ("Unknown", int),
            ("StartTime", object),
            ("EndTime", object),
            ("Attempts", int),
        ]
    )
    saltDF = pd.DataFrame.from_records(saltArray)
    # add as many "Reading#"s as needed
    for ii in range(0, len(saltDF.columns) - len(header)):
        header["Reading{}".format(ii + 1)] = float
    saltDF.columns = list(header.keys())  # name columns
    saltDF = saltDF[saltDF["autosalSAMPNO"] != "worm"]
    saltDF = saltDF.astype(header)  # force dtypes
    saltDF["SALNTY"] = gsw.SP_salinometer(
        (saltDF["CRavg"] / 2.0), saltDF["BathTEMP"]
    ).round(4)
    return saltDF


# TODO: cleanup and improve
def salt_df_parser(saltDF, outdir, stn_col="STNNBR", cast_col="CASTNO"):
    stations = saltDF[stn_col].unique()
    for station in stations:
        saltStation = saltDF[saltDF[stn_col] == station]
        casts = saltStation[cast_col].unique()
        for cast in casts:
            stn_cast_salts = saltStation[saltStation[cast_col] == cast]
            # format to SSSCC_salts.csv
            outfile = (
                outdir + "{0:03}".format(station) + "{0:02}".format(cast) + "_salts.csv"
            )
            if not os.path.exists(outfile):
                stn_cast_salts.to_csv(outfile, index=False)
            else:
                print(outfile + " already exists...skipping")


def process_salts(ssscc, salt_dir):
    try:
        saltDF = salt_loader(ssscc, salt_dir)
    except FileNotFoundError:
        print("Salt file for cast " + ssscc + " does not exist... skipping")
        return
    salt_df_parser(saltDF, salt_dir)

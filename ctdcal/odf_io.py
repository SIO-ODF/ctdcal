"""
A module for handling Autosalinomter files, in the ODF format (Carl Mattson).
"""

import csv
import io
import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from ctdcal import get_ctdcal_config

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def _salt_loader(filename, flag_file="tools/salt_flags_handcoded.csv"):
    """
    Load raw file into salt and reference DataFrames.
    """

    csv_opts = dict(delimiter=" ", quoting=csv.QUOTE_NONE, skipinitialspace="True")
    if isinstance(filename, (str, Path)):
        with open(filename, newline="") as f:
            saltF = csv.reader(f, **csv_opts)
            saltArray = [row for row in saltF]
            ssscc = Path(filename).stem
    elif isinstance(filename, io.StringIO):
        saltF = csv.reader(filename, **csv_opts)
        saltArray = [row for row in saltF]
        ssscc = "test_odf_io"
    else:
        raise NotImplementedError(
            "Salt loader only able to read in str, Path, or StringIO classes"
        )

    del saltArray[0]  # remove file header
    saltDF = pd.DataFrame.from_records(saltArray)

    cols = dict(  # having this as a dict streamlines next steps
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

    # add as many "Reading#"s as needed
    for ii in range(0, len(saltDF.columns) - len(cols)):
        cols["Reading{}".format(ii + 1)] = float
    saltDF.columns = list(cols.keys())  # name columns

    # TODO: check autosalSAMPNO against SAMPNO for mismatches?
    # TODO: handling for re-samples?

    # check for commented out lines
    commented = saltDF["STNNBR"].str.startswith(("#", "x"))
    if commented.any():
        log.debug(f"Found comment character (#, x) in {ssscc} salt file, ignoring line")
        saltDF = saltDF[~commented]

    # check end time for * and code questionable
    # (unconfirmed but * appears to indicate a lot of things from LabView code:
    # large spread in values, long time between samples, manual override, etc.)
    flagged = saltDF["EndTime"].str.contains("*", regex=False)
    if flagged.any():
        # remove asterisks from EndTime and flag samples
        log.debug(f"Found * in {ssscc} salt file, flagging value(s) as questionable")
        saltDF["EndTime"] = saltDF["EndTime"].str.strip("*")
        questionable = pd.DataFrame()
        questionable["SAMPNO"] = saltDF.loc[flagged, "SAMPNO"].astype(int)
        questionable.insert(0, "SSSCC", ssscc)
        questionable["diff"] = np.nan
        questionable["salinity_flag"] = 3
        questionable["comments"] = "Auto-flagged by processing function (had * in row)"
        questionable.to_csv(flag_file, mode="a+", index=False, header=None)

    # add time (in seconds) needed for autosal drift removal step
    saltDF["IndexTime"] = pd.to_datetime(saltDF["EndTime"], format="%H:%M:%S")
    saltDF["IndexTime"] = (saltDF["IndexTime"] - saltDF["IndexTime"].iloc[0]).dt.seconds
    saltDF["IndexTime"] += (saltDF["IndexTime"] < 0) * (3600 * 24)  # fix overnight runs

    refDF = saltDF.loc[
        saltDF["autosalSAMPNO"] == "worm", ["IndexTime", "CRavg"]
    ].astype(float)
    saltDF = saltDF[saltDF["autosalSAMPNO"] != "worm"].astype(cols)  # force dtypes

    return saltDF, refDF


def remove_autosal_drift(saltDF, refDF):
    """Calculate linear CR drift between reference values"""
    if refDF.shape != (2, 2):
        ssscc = f"{saltDF['STNNBR'].unique()[0]:03d}{saltDF['CASTNO'].unique()[0]:02d}"
        log.warning(
            f"Failed to find start/end reference readings for {ssscc}, check salt file"
        )
    else:
        # find rate of drift
        diff = refDF.diff(axis="index").dropna()
        time_coef = (diff["CRavg"] / diff["IndexTime"]).iloc[0]

        # apply offset as a linear function of time
        saltDF = saltDF.copy(deep=True)  # avoid modifying input dataframe
        saltDF["CRavg"] += saltDF["IndexTime"] * time_coef
        saltDF["CRavg"] = saltDF["CRavg"].round(5)  # match initial precision

    return saltDF.drop(labels="IndexTime", axis="columns")


def _salt_exporter(
    saltDF, outdir=cfg.dirs["salt"], stn_col="STNNBR", cast_col="CASTNO"
):
    """
    Export salt DataFrame to .csv file. Extra logic is included in the event that
    multiple stations and/or casts are included in a single raw salt file.
    """
    stations = saltDF[stn_col].unique()
    for station in stations:
        stn_salts = saltDF[saltDF[stn_col] == station]
        casts = stn_salts[cast_col].unique()
        for cast in casts:
            stn_cast_salts = stn_salts[stn_salts[cast_col] == cast].copy()
            stn_cast_salts.dropna(axis=1, how="all", inplace=True)  # drop empty columns
            outfile = Path(outdir) / f"{station:03.0f}{cast:02.0f}_salts.csv"  # SSSCC_*
            if outfile.exists():
                log.info(str(outfile) + " already exists...skipping")
                continue
            stn_cast_salts.to_csv(outfile, index=False)


def process_salts(ssscc_list, salt_dir=cfg.dirs["salt"]):
    """
    Master salt processing function. Load in salt files for given station/cast list,
    calculate salinity, and export to .csv files.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to process
    salt_dir : str, optional
        Path to folder containing raw salt files (defaults to data/salt/)

    """
    for ssscc in ssscc_list:
        if (Path(salt_dir) / f"{ssscc}_salts.csv").exists():
            log.info(f"{ssscc}_salts.csv already exists in {salt_dir}... skipping")
            continue
        else:
            try:
                saltDF, refDF = _salt_loader(Path(salt_dir) / ssscc)
            except FileNotFoundError:
                log.warning(f"Salt file for cast {ssscc} does not exist... skipping")
                continue
            saltDF = remove_autosal_drift(saltDF, refDF)
            saltDF["SALNTY"] = gsw.SP_salinometer(
                (saltDF["CRavg"] / 2.0), saltDF["BathTEMP"]
            )  # .round(4)
            _salt_exporter(saltDF, salt_dir)

def print_progress_bar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    printEnd="\r",
):
    """
    A progress bar, helpful for implementing into loops or highlighting progression through processing.
    
    https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters/13685020
    credit: u/Greenstick
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)  #   Potential to add to log
    # Print New Line on Complete
    if iteration == total:
        print()
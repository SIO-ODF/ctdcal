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
from ctdcal.common import validate_file
from ctdcal.fitting.common import (
    NodeNotFoundError,
    df_node_to_BottleFlags,
    get_node,
    save_node,
)

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def _salt_loader(filename):
    """
    Load raw ODF salt file into salt, reference, and flagged DataFrames.

    ODF salt files have a header of configuration information, including the K15 and run date.
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

    #   Check on the header
    date = pd.to_datetime(saltArray[0][-1])
    k15  = float(saltArray[0][8])
    if ((date < pd.Timestamp('1985-01-01')) | (date > pd.Timestamp.today())):
        #   Check on the run date - should this be ± a few weeks to confirm that the analyst's computer is set up right?
        log.warning(f"Salt file for {ssscc} has an erroneous timestamp in the header.")
    if ((int(saltArray[0][-2]) < 0) | int(saltArray[0][-2]) > 1000):
        #   Standardization dial ranges from 0 - 1000
        log.warning(f"Salt file for {ssscc} has an erroneous standardization dial in the header.")

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
    if "Reading6" in cols:
        log.warning("More than 5 readings found in salts file - check columns of source.")
    saltDF.columns = list(cols.keys())  # name columns

    # check for commented out lines
    commented = saltDF["STNNBR"].str.startswith(("#", "x"))
    if commented.any():
        log.debug(f"Found comment character (#, x) in {ssscc} salt file, ignoring line")
        saltDF = saltDF[~commented]
    if saltDF["autosalSAMPNO"].value_counts().get('worm',0) != 2:
        #   Confirm 2-point standardization
        log.warning(f"Suspect number of standards found in salts file {ssscc}. Confirm standardizations.")    
    checkvals = saltDF[saltDF.autosalSAMPNO == "worm"].CRavg.astype(float)
    if ((checkvals >= 1.0001 * 2 * k15) & (checkvals <= 0.9999 * 2 * k15)).any():
        #   Check the worm values rel. to batch K15 - within 1/100 %
        log.warning(f"Standards in raw salt file {ssscc} are >0.01 % off from expected K15 value.")
    if any((saltDF['BathTEMP'].astype(float) < 17.5) | (saltDF['BathTEMP'].astype(float) > 33.5)):
        #   Check bath temperature
        log.warning(f"Raw salt file for {ssscc} contains bath temperatures outside of operational ranges.")

    # check end time for * and code questionable
    # (unconfirmed but * appears to indicate a lot of things from LabView code:
    # large spread in values, long time between samples, manual override, etc.)
    flagged = saltDF["EndTime"].str.contains("*", regex=False)
    questionable = None
    if flagged.any():
        # remove asterisks from EndTime and flag samples
        log.debug(f"Found * in {ssscc} salt file, flagging value(s) as questionable")
        saltDF["EndTime"] = saltDF["EndTime"].str.strip("*")
        # construct flagged dataframe
        questionable = pd.DataFrame()
        questionable["bottle_num"] = saltDF.loc[flagged, "SAMPNO"].astype(int)
        questionable["cast_id"] = [str(s).zfill(3)[-3:] + str(c).zfill(2) for s, c in zip(saltDF.loc[flagged, "STNNBR"], saltDF.loc[flagged, "CASTNO"])]
        questionable["value"] = 3
        questionable["notes"] = "Auto-flagged by processing function (had * in row)"

    # add time (in seconds) needed for autosal drift removal step
    saltDF["IndexTime"] = pd.to_datetime(saltDF["EndTime"], format="%H:%M:%S")
    saltDF["IndexTime"] = (saltDF["IndexTime"] - saltDF["IndexTime"].iloc[0]).dt.seconds
    saltDF["IndexTime"] += (saltDF["IndexTime"] < 0) * (3600 * 24)  # fix overnight runs

    refDF = saltDF.loc[
        saltDF["autosalSAMPNO"] == "worm", ["IndexTime", "CRavg"]
    ].astype(float)
    saltDF = saltDF[saltDF["autosalSAMPNO"] != "worm"].astype(cols)  # force dtypes

    return saltDF, refDF, questionable


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

        if any(abs(diff["CRavg"] > 0.0001)):
            #   Warn analyst that the standard drifted by a lot (precise to within 5 decimal places)
            log.warning(f"Salt run at station {saltDF['STNNBR'].iloc[0]} had a CR drift in excess of 0.0001.")
        if any(diff["IndexTime"] > 43200):
            #   Warn analyst that over 12 hours passed between standardizations
            log.warning(f"Salt run at station {saltDF['STNNBR'].iloc[0]} had >12 hours between standardizations.")
        if any(saltDF["IndexTime"].diff().sum() > (diff["IndexTime"])):
            #   Warn analyst that standardization and sample timestamps may not be aligned
            #   Standards should always bookend the samples
            log.warning(f"Salt run at station {saltDF['STNNBR'].iloc[0]} has misaligned standards and samples.")

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
    if len(stations) > 1:
        log.info("Multiple stations found in salt source file. Writing out .CSVs seperately.")

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


def process_salts(ssscc_list, user_cfg=None, salt_dir=cfg.dirs["salt"]):
    """
    Master salt processing function. Load in salt files for given station/cast list,
    calculate salinity, and export to .csv files.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to process
    user_cfg : Munch object
        Dictionary of user configuration parameters
    salt_dir : str, optional
        Path to folder containing raw salt files (defaults to data/salt/)

    """
    flags_df = None
    for ssscc in ssscc_list:
        if (Path(salt_dir) / f"{ssscc}_salts.csv").exists():
            log.info(f"{ssscc}_salts.csv already exists in {salt_dir}... skipping")
            continue
        else:
            try:
                saltDF, refDF, questionable = _salt_loader(Path(salt_dir) / ssscc)
            except FileNotFoundError:
                log.warning(f"Salt file for cast {ssscc} does not exist... skipping")
                continue
            saltDF = remove_autosal_drift(saltDF, refDF)
            saltDF["SALNTY"] = gsw.SP_salinometer(
                (saltDF["CRavg"] / 2.0), saltDF["BathTEMP"]
            )  # .round(4)
            _salt_exporter(saltDF, salt_dir)

            # compile flags
            if questionable is not None:
                if flags_df is None:
                    flags_df = questionable
                else:
                    flags_df = pd.concat([flags_df, questionable], ignore_index=True)



    # save flags
    if flags_df is not None:
        if user_cfg is not None:
            flag_path = Path(user_cfg.datadir, 'flag', user_cfg.bottleflags_man)
        else:   #   No user_cfg, use get_ctdcal_config TODO: Remember to remove all the old cfg calls
            flag_path = Path(cfg.dirs('flags'), "bottleflags_man.csv")
        flag_file = validate_file(flag_path, create=True)
        try:
            salt = pd.DataFrame.from_dict(get_node(flag_file, 'salt'))
        except NodeNotFoundError:
            salt = None
        flags_df = pd.concat([flags_df, salt], ignore_index=True)
        new_flags = df_node_to_BottleFlags(flags_df)
        save_node(flag_file, new_flags, 'salt', create_new=True)


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

import csv
import io
import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from . import get_ctdcal_config

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
    saltDF["IndexTime"] = pd.to_datetime(saltDF["EndTime"])
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


def sdl_loader(filename):
    """
    For loading sdl.dat (Salinity data logger) files and rewriting as ODF formatted DataFrames.
    """
    sdl = pd.read_csv(filename, sep="\t")  #   Read tab deliminated file
    sdl = sdl.loc[sdl.ReadingNumber.isnull()]  #   Remove non-averaged readings
    try:
        cut_std = sdl[
            sdl["BottleLabel"].str.contains("P")
        ]  #   Remove rows associated with standard (P coming from IAPSO batch)
        sdl = sdl[sdl["BottleLabel"].str.contains("P") == False]
    except AttributeError:
        #   If standard was run previously
        cut_std = pd.DataFrame()

    if not sdl.empty:

        #   Create SSSCC columns. For OSNAP, CC = 01
        sdl[["ID2", "Num"]] = sdl.SampleID.str.split("#", expand=True)
        sdl["SSS"] = sdl.ID2.str[0:3]
        sdl["CC"] = sdl.ID2.str[-2:]

        #   Extract microcat samples and write Autosal lines for microcat group
        #   (Pass the microcat.csv file into sdl_std for salinity)
        commented = sdl[sdl.Comments.isnull() == False]
        if not commented.empty:
            micro = commented[commented.Comments.str.contains("cat")]
            if not micro.empty:
                #   Cut out the columns they don't want
                micro = pd.DataFrame.from_dict(
                    {
                        "STNNBR": micro.SSS,
                        "BOTTLE": micro.BottleLabel.astype(int),
                        "CRavg": 2
                        * micro.AdjustedRatio,  #   SDL writes ratio out as half of what ODF routine does
                        "CR_unadjusted": 2 * micro.UncorrectedRatio,
                        "EndTime": micro.DateTime,
                        "Comment": micro.Comments,
                    }
                )
                if cut_std.empty:
                    #   Need to pass through sdl_std
                    micro.to_csv(
                        Path(cfg.dirs["salt"])
                        / f"microcat_salt_to_correct_{micro.STNNBR.iloc[0]}.csv",
                        index=False,
                    )
                else:
                    micro["SALNTY"] = gsw.SP_salinometer((micro["CRavg"] / 2), 24)
                    micro.to_csv(
                        Path(cfg.dirs["salt"])
                        / f"microcat_salt_{micro.STNNBR.iloc[0]}.csv",
                        index=False,
                    )

        sdl = sdl[
            sdl["Comments"].isnull() == True
        ]  #   Cut out everything  that has a comment (Aaron indicating as needing a rerun or microcat, which are not fit)

        #   Make saltDF in same format as odf_io
        to_write = {
            "STNNBR": sdl.SSS,
            "CASTNO": sdl.CC,
            "SAMPNO": sdl.BottleLabel.astype(int),
            "BathTEMP": 24,  #   There is no option to export bath temp with decimals and is occasionally very wrong
            "CRavg": 2
            * sdl.AdjustedRatio,  #   SDL writes ratio out as half of what ODF routine does
            "CRraw": 2 * sdl.UncorrectedRatio,
            "autosalSAMPNO": sdl.SampleNumber,
            "Unknown": np.nan,  # Standardize dial
            "StartTime": np.nan,
            "EndTime": sdl.DateTime,
            "Attempts": np.nan,
        }
        saltDF = pd.DataFrame.from_dict(to_write)
        saltDF.reset_index(inplace=True)  #   Ignore the index given to it before
    else:
        print(
            "Warning: No entries discovered in saltDF and results in empty dataframe:",
            filename,
        )

    return saltDF, cut_std


def sdl_std(saltDF, cut_std, salt_dir=cfg.dirs["salt"], infile="standards.csv"):
    """For applying standard correction to saltDF/updating standards list in Salinity data logger"""
    outfile = Path(salt_dir) / infile
    if not outfile.exists():
        #   If the standards file does not exist, create it
        print("Salt standards file does not exist. Creating...")
        header = [
            "SampleID",
            "BottleLabel",
            "DateTime",
            "BathTemperature",
            "UncorrectedRatio",
            "UncorrectedRatioStandDev",
            "Correction",
            "AdjustedRatio",
        ]
        cut_std.to_csv(Path(outfile), columns=header, index=False)

    std_list = pd.read_csv(outfile)
    if cut_std.empty:
        #   If no standards were run in the last iteration, the CRavg = CRraw and needs a flat adjustment
        #   If a standard is run, then the CRavg already has a flat standardization adjustment

        # print(
        #     saltDF.STNNBR.iloc[0],
        #     "is missing a standard. Pulling from the appropriate standard run...",
        # )

        try:
            times = pd.to_datetime(std_list.DateTime)
            samp_time = pd.to_datetime(saltDF.EndTime.iloc[0])
            corr = std_list.Correction.loc[
                times == min(times, key=lambda d: abs(d - samp_time))
            ].item()
            # print("Pulled correction value:", corr)
        except:
            #   Read last value in as Correction
            corr = std_list.Correction.iloc[-1]
            print(
                "Latest standard run was:",
                std_list["DateTime"].iloc[-1],
                "Correction =",
                corr,
            )

        #   apply offset from last standard (CRavg = CRavg + Correction)
        saltDF["CRavg"] = saltDF["CRavg"] + corr * 2
        # print("New CRavg:\n", saltDF["CRavg"])

    else:  #   Append the standard lines to the standards file. Use current adjusted CR.
        print("Appending new standards to list.")
        std_list = pd.concat([std_list, cut_std], ignore_index=True)
        std_list = std_list.drop_duplicates(
            subset="DateTime", keep="first"
        )  #   In case of overlap
        std_list.to_csv(Path(outfile), index=False)

    return saltDF


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

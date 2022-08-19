"""
A series of functions for reading in, adjusting, and exporting salinity data from many file types.

Loads the raw .dat files from the AutoSal "Salinometer Data Logger" software

Developed for OSNAP32, 2022 with the intention of merging into odf_io
"""

import gsw
import csv
import io
import logging
import numpy as np
import pandas as pd
from pathlib import Path

from . import get_ctdcal_config

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def sdl_loader(filename):
    """
    For loading sdl.dat files and rewriting as ODF formatted DataFrames.

    TODO: Flagging?
    """
    sdl = pd.read_csv(filename, sep="\t")  #   Read tab deliminated file
    sdl = sdl.loc[sdl.ReadingNumber.isnull()]  #   Remove non-averaged readings
    cut_std = sdl[
        sdl["BottleLabel"].str.contains("P")
    ]  #   Remove rows associated with standard (P coming from IAPSO batch)
    sdl = sdl[sdl["BottleLabel"].str.contains("P") == False]
    sdl[["ID2", "Num"]] = sdl.SampleID.str.split("#", expand=True)
    sdl["SSS"] = sdl.ID2.str[0:3]
    sdl["CC"] = sdl.ID2.str[3:5]  #   Create SSSCC columns
    to_write = {
        "STNNBR": sdl.SSS,
        "CASTNO": sdl.CC,
        "SAMPNO": sdl.BottleLabel,
        "BathTEMP": sdl.BathTemperature,
        "CRavg": sdl.AdjustedRatio,
        "CRraw": sdl.UncorrectedRatio,
        "autosalSAMPNO": sdl.SampleNumber,
        "Unknown": np.nan,
        "StartTime": np.nan,
        "EndTime": sdl.DateTime,
        "Attempts": np.nan,
    }
    saltDF = pd.DataFrame.from_dict(to_write)
    saltDF.reset_index(inplace=True)  #   Ignore the index given to it before

    return saltDF, cut_std


def sdl_std(saltDF, cut_std, salt_dir=cfg.dirs["salt"], infile="standards.csv"):
    """For applying standard correction to saltDF/updating standards list"""
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
        print(
            "No standards were in the last run. Pulling from the previous standard run..."
        )
        #   Read last value in as Correction
        corr = std_list["Correction"].iloc[-1]
        print(
            "Latest standard run was:",
            std_list["DateTime"].iloc[-1],
            "Correction =",
            corr,
        )
        #   apply offset from last standard (CRavg = CRavg + Correction)
        saltDF["CRavg"] = saltDF["CRavg"] + corr
        print("New CRavg:\n", saltDF["CRavg"])

    else:  #   Append the standard lines to the standards file. Use current adjusted CR.
        print("Appending new standards to list.")
        std_list = pd.concat([std_list, cut_std], ignore_index=True)
        std_list = std_list.drop_duplicated(
            subset="DateTime", keep="first"
        )  #   In case of overlap
        std_list.to_csv(Path(outfile), index=False)

    return saltDF


def sdl_exporter(saltDF, outdir=cfg.dirs["salt"], stn_col="STNNBR", cast_col="CASTNO"):
    """
    Basically a copy of _salt_exporter. For each file only one cast is expected, as SSSCC in the Salinometer Data Logger cannot be changed during a run.

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


def osnap_salts(ssscc_list, salt_dir=cfg.dirs["salt"]):
    """
    Basically a copy of odf_io.proces_salts using modified loader and standardization functions.

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
                saltDF, cut_std = sdl_loader(Path(salt_dir) / ssscc)
            except FileNotFoundError:
                log.warning(f"Salt file for cast {ssscc} does not exist... skipping")
                continue
            saltDF = sdl_std(saltDF, cut_std)
            saltDF["SALNTY"] = gsw.SP_salinometer(
                (saltDF["CRavg"] / 2.0), saltDF["BathTEMP"]
            )  # .round(4)
            sdl_exporter(saltDF, salt_dir)

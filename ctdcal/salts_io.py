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

from ctdcal import get_ctdcal_config

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def sdl_loader_dat(filename):
    """
    For loading sdl.dat files and rewriting as ODF formatted DataFrames, in accordance
    to routines aboard OSNAP32.

    TODO: Flagging?
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

def sdl_loader_xlsx(filename):
    """
    For loading sdl.xlsx files and rewriting as ODF formatted DataFrames, in accorance
    to routines aboard NBP2401.

    Originally intended as a subset of sdl_loader_dat, though there were sufficient changes
    in routine xlsx generation that warranted new columns and to_write parameters.

    Requires openpyxl
    """
    sdl = pd.read_excel(filename)
    sdl = sdl.loc[~sdl.DateTime.isnull()]   #   Cut out the redundant "to average" lines
    sdl = sdl[sdl["Reading #"].isnull()]    #   Exctract the final values from however many readings

    sdl["IndexTime"] = pd.to_datetime(sdl["DateTime"])
    sdl["IndexTime"] = (sdl["IndexTime"] - sdl["IndexTime"].iloc[0]).dt.seconds
    sdl["IndexTime"] += (sdl["IndexTime"] < 0) * (3600 * 24)

    #   Start cutting out the standards
    try:
        cut_std = sdl[
            sdl["Sample ID"].str.contains("Cal")
        ]  #   Remove rows associated with standard (P coming from IAPSO batch)
        sdl = sdl.loc[~sdl["Calculated Salinity"].isnull()]
    except AttributeError:
        #   If standard was run previously
        cut_std = pd.DataFrame()

    if not sdl.empty:
        #   Create SSSCC columns. For NBP2401, CC = 01
        #   SSSCC GTC (dates, times, whatever) #samplenum
        sdl[["SSSCC", "Num"]] = sdl["Sample ID"].str.split("#", expand=True)
        sdl["SSS"] = sdl.SSSCC.str[0:3]
        sdl["CC"] = sdl.SSSCC.str[3:5]
        sdl["SSSCC"] = sdl["SSS"] + sdl["CC"]

        to_write = {
            "STNNBR": sdl.SSS,
            "CASTNO": sdl.CC,
            "SAMPNO": sdl.Num.astype(int),
            "BathTEMP": sdl["Bath Temperature"].str.strip("C").astype(int),  #   There is no option to export bath temp with decimals and is occasionally very wrong
            "CRavg": 2
            * sdl["Adjusted Ratio"],  #   SDL writes ratio out as half of what ODF routine does
            "CRraw": 2 * sdl["Uncorrected Ratio"],
            "autosalSAMPNO": sdl.Num,
            "Unknown": np.nan,  # Standardize dial
            "StartTime": np.nan,
            "EndTime": sdl.DateTime,
            "Attempts": np.nan,
            "IndexTime": sdl.IndexTime
        }
        saltDF = pd.DataFrame.from_dict(to_write)
        saltDF.reset_index(inplace=True)
    else:
        print(
            "Warning: No entries discovered in saltDF and results in empty dataframe:",
            filename,
        )

    return saltDF, cut_std

def portasal_remove_drift(saltDF, cut_std):
    """
    For removing the difference between the first and last samples.
    Portasal values are automatically adjusted for the first standard.

    This is basically a copy from odf_io, assuming the run was
    standardized with a bottle at the beginning and the end.
    """
    if len(cut_std) > 1:
        diff = cut_std[["Adjusted Ratio", "IndexTime"]].diff().dropna()
        time_coef = (diff["Adjusted Ratio"]*2 / diff["IndexTime"]).iloc[0]
        saltDF = saltDF.copy(deep=True)  # avoid modifying input dataframe
        saltDF["CRavg"] += saltDF["IndexTime"] * time_coef
        saltDF["CRavg"] = saltDF["CRavg"].round(5)
    elif len(cut_std) == 1:
        log.warning(f"Only 1 standardization found for run {saltDF.SSSCC.iloc[0]}")
    else:
        log.warning(f"Standardization of {saltDF.SSSCC.iloc[0]} was an unusual shape.")

    return saltDF.drop(labels="IndexTime", axis="columns")

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


def sdl_exporter(saltDF, outdir=cfg.dirs["salt"], stn_col="STNNBR", cast_col="CASTNO"):
    """
    Basically a copy of _salt_exporter. For each file only one cast is expected, as SSSCC in the Salinometer Data Logger cannot be changed during a run.

    Export salt DataFrame to .csv file. Extra logic is included in the event that
    multiple stations and/or casts are included in a single raw salt file.
    """
    saltDF = saltDF.drop("index", axis=1)  #   Not needed in written file
    stations = saltDF[stn_col].unique()
    for station in stations:
        stn_salts = saltDF[saltDF[stn_col] == station]
        casts = stn_salts[cast_col].unique()
        for cast in casts:
            stn_cast_salts = stn_salts[stn_salts[cast_col] == cast].copy()
            stn_cast_salts.dropna(axis=1, how="all", inplace=True)  # drop empty columns
            outfile = Path(outdir) / f"{station}{cast}_salts.csv"  # SSS
            if outfile.exists():
                log.info(str(outfile) + " already exists...skipping")
                continue
            stn_cast_salts.to_csv(outfile, index=False)


def osnap_microcat(salt_dir=cfg.dirs["salt"]):
    import glob

    microcat_list = glob.glob(salt_dir + "microcat_salt_to_correct*")
    if len(microcat_list) != 0:
        for microcat in microcat_list:
            saltDF = pd.read_csv(microcat)
            no_std = pd.DataFrame()
            saltDF = sdl_std(saltDF, no_std)  #   Apply standard offset
            saltDF["SALNTY"] = gsw.SP_salinometer((saltDF["CRavg"] / 2), 24)
            outfile = str(Path(salt_dir)) + "/microcat_salt_" + microcat[-7:-4] + ".csv"
            saltDF.to_csv(outfile, index=False)


def osnap_salts(ssscc_list, salt_dir=cfg.dirs["salt"]):
    """
    Basically a copy of odf_io.proces_salts using modified loader and standardization functions.
    * Read in file with .dat extension
    * Adjust the conductivity ratio by a flat correction
        * If no conductivity ratio is in the file, read the last one (many casts per day)
    * Calculate SALNTY with 2002 software exports
    * Write the file, dropping any unnessesary columns

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
                ext = ".dat"
                saltDF, cut_std = sdl_loader_dat(Path(str(Path(salt_dir) / ssscc) + ext))
            except FileNotFoundError:
                log.warning(f"Salt file for cast {ssscc} does not exist... skipping")
                continue
            saltDF = sdl_std(saltDF, cut_std)
            #   The ODF autosal writeout doubles the CRavg. Here we don't have to divide by 2.
            saltDF["SALNTY"] = gsw.SP_salinometer(
                (saltDF["CRavg"] / 2), saltDF["BathTEMP"]
            )  # .round(4)
            sdl_exporter(saltDF, salt_dir)
    osnap_microcat()  #   Correct CRavg and add SALNTY for any microcats

def portasal_salts(ssscc_list, salt_dir=cfg.dirs["salt"]):
    """
    Basically a copy of osnap_salts using modified loader and standardization functions.
    * Read in file with .xlsx extension
    * Adjust the conductivity ratio by a linear correction
        * If no conductivity ratio is in the file, read the last one (many casts per day)
    * Calculate SALNTY with 2002 software exports
    * Write the file, dropping any unnessesary columns

    Master salt processing function. Load in salt files for given station/cast list,
    calculate salinity, and export to .csv files.

    Developed for use on NBP2401: GEOTRACES

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
                ext = ".xlsx"
                saltDF, cut_std = sdl_loader_xlsx(Path(str(Path(salt_dir) / ssscc) + ext))
            except FileNotFoundError:
                log.warning(f"Salt file for cast {ssscc} does not exist... skipping")
                continue
            saltDF = portasal_remove_drift(saltDF, cut_std)   # Linear correction
            #   The ODF autosal writeout doubles the CRavg. Here we don't have to divide by 2.
            saltDF["SALNTY"] = gsw.SP_salinometer(
                (saltDF["CRavg"] / 2), saltDF["BathTEMP"]
            )  # .round(4)
            sdl_exporter(saltDF, salt_dir)
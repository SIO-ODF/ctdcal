"""Library to create SBE .btl equivalent files.
TODO: allow for variable bottle fire scans instead of SBE standard 36
    ex: user doesn't know how to change the config for the cast to add more scans,
    instead does it post-cast?

Joseph Gum SIO/ODF
Nov 7, 2016
"""

import csv
import logging
import statistics
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from . import flagging as flagging
from . import get_ctdcal_config
from . import oxy_fitting as oxy_fitting

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

BOTTLE_FIRE_COL = "btl_fire"
BOTTLE_FIRE_NUM_COL = "btl_fire_num"


# Retrieve the bottle data from a converted file.
def retrieveBottleDataFromFile(converted_file):

    converted_df = pd.read_pickle(converted_file)

    return retrieveBottleData(converted_df)


# Retrieve the bottle data from a dataframe created from a converted file.
def retrieveBottleData(converted_df):
    if BOTTLE_FIRE_COL in converted_df.columns:
        converted_df[BOTTLE_FIRE_NUM_COL] = (
            (
                (converted_df[BOTTLE_FIRE_COL])
                & (
                    converted_df[BOTTLE_FIRE_COL]
                    != converted_df[BOTTLE_FIRE_COL].shift(1)
                )
            )
            .astype(int)
            .cumsum()
        )
        # converted_df['bottle_fire_num'] = ((converted_df[BOTTLE_FIRE_COL] == False)).astype(int).cumsum()
        return converted_df.loc[converted_df[BOTTLE_FIRE_COL]]
        # return converted_df
    else:
        log.error("Bottle fire column:", BOTTLE_FIRE_COL, "not found")

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


def _load_btl_data(btl_file, cols=None):
    """
    Loads "bottle mean" CTD data from .pkl file. Function will return all data unless
    cols is specified (as a list of column names)
    """

    btl_data = pd.read_pickle(btl_file)
    if cols is not None:
        btl_data = btl_data[cols]
    btl_data["SSSCC"] = Path(btl_file).stem.split("_")[0]

    return btl_data


def _load_reft_data(reft_file, index_name="btl_fire_num"):
    """
    Loads reft_file to dataframe and reindexes to match bottle data dataframe
    """
    reft_data = pd.read_csv(reft_file, usecols=["btl_fire_num", "T90", "REFTMP_FLAG_W"])
    reft_data.set_index(index_name)
    reft_data["SSSCC_TEMP"] = Path(reft_file).stem.split("_")[0]
    reft_data["REFTMP"] = reft_data["T90"]

    return reft_data


def _load_salt_data(salt_file, index_name="SAMPNO"):
    """
    Loads salt_file to dataframe and reindexes to match bottle data dataframe
    """
    salt_data = pd.read_csv(
        salt_file, usecols=["SAMPNO", "SALNTY", "BathTEMP", "CRavg"]
    )
    salt_data.set_index(index_name)
    salt_data["SSSCC_SALT"] = Path(salt_file).stem.split("_")[0]
    salt_data.rename(columns={"SAMPNO": "SAMPNO_SALT"}, inplace=True)

    return salt_data


def _add_btl_bottom_data(df, cast, lat_col="LATITUDE", lon_col="LONGITUDE", decimals=4):
    cast_details = pd.read_csv(
        # cfg.dirs["logs"] + "cast_details.csv", dtype={"SSSCC": str}
        cfg.dirs["logs"] + "bottom_bottle_details.csv",
        dtype={"SSSCC": str},
    )
    cast_details = cast_details[cast_details["SSSCC"] == cast]
    # df[lat_col] = np.round(cast_details["latitude"].iat[0], decimals)
    # df[lon_col] = np.round(cast_details["longitude"].iat[0], decimals)
    df[lat_col] = cast_details["latitude"].iat[0]
    df[lon_col] = cast_details["longitude"].iat[0]

    ts = pd.to_datetime(cast_details["bottom_time"].iat[0], unit="s")
    date = ts.strftime("%Y%m%d")
    hour = ts.strftime("%H%M")
    df["DATE"] = date
    df["TIME"] = hour
    return df


def load_all_btl_files(ssscc_list, cols=None):
    """
    Load bottle and secondary (e.g. reference temperature, bottle salts, bottle oxygen)
    files for station/cast list and merge into a dataframe.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to load
    cols : list of str, optional
        Subset of columns to load, defaults to loading all

    Returns
    -------
    df_data_all : DataFrame
        Merged dataframe containing all loaded data

    """
    df_data_all = pd.DataFrame()

    for ssscc in ssscc_list:
        log.info("Loading BTL data for station: " + ssscc + "...")
        btl_file = cfg.dirs["bottle"] + ssscc + "_btl_mean.pkl"
        btl_data = _load_btl_data(btl_file, cols)

        ### load REFT data
        reft_file = cfg.dirs["reft"] + ssscc + "_reft.csv"
        try:
            reft_data = _load_reft_data(reft_file)
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFT Data Station: "
                + ssscc
                + "...filling with NaNs"
            )
            reft_data = pd.DataFrame(index=btl_data.index, columns=["T90"], dtype=float)
            reft_data["btl_fire_num"] = btl_data["btl_fire_num"].astype(int)
            reft_data["SSSCC_TEMP"] = ssscc  # TODO: is this ever used?

        ### load REFC data
        refc_file = cfg.dirs["salt"] + ssscc + "_salts.csv"
        try:
            refc_data = _load_salt_data(refc_file, index_name="SAMPNO")
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFC Data Station: "
                + ssscc
                + "...filling with NaNs"
            )
            refc_data = pd.DataFrame(
                index=btl_data.index,
                columns=["CRavg", "BathTEMP", "BTLCOND"],
                dtype=float,
            )
            refc_data["SAMPNO_SALT"] = btl_data["btl_fire_num"].astype(int)

        ### load OXY data
        oxy_file = cfg.dirs["oxygen"] + ssscc
        try:
            oxy_data, params = oxy_fitting.load_winkler_oxy(oxy_file)
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFO Data Station: "
                + ssscc
                + "...filling with NaNs"
            )
            oxy_data = pd.DataFrame(
                index=btl_data.index,
                columns=[
                    "FLASKNO",
                    "TITR_VOL",
                    "TITR_TEMP",
                    "DRAW_TEMP",
                    "TITR_TIME",
                    "END_VOLTS",
                ],
                dtype=float,
            )
            oxy_data["STNNO_OXY"] = ssscc[:3]  # TODO: are these values
            oxy_data["CASTNO_OXY"] = ssscc[3:]  # ever used?
            oxy_data["BOTTLENO_OXY"] = btl_data["btl_fire_num"].astype(int)

        ### clean up dataframe
        # Horizontally concat DFs to have all data in one DF
        btl_data = pd.merge(btl_data, reft_data, on="btl_fire_num", how="outer")
        btl_data = pd.merge(
            btl_data,
            refc_data,
            left_on="btl_fire_num",
            right_on="SAMPNO_SALT",
            how="outer",
        )
        btl_data = pd.merge(
            btl_data,
            oxy_data,
            left_on="btl_fire_num",
            right_on="BOTTLENO_OXY",
            how="outer",
        )

        if len(btl_data) > 36:
            log.error(
                f"""Length of bottle data for {ssscc} is > 36, check for errors in reference parameter files"""
            )

        # Add bottom of cast information (date,time,lat,lon,etc.)
        btl_data = _add_btl_bottom_data(btl_data, ssscc)

        # Merge cast into df_data_all
        try:
            df_data_all = pd.concat([df_data_all, btl_data], sort=False)
        except AssertionError:
            raise AssertionError(
                "Columns of " + ssscc + " do not match those of previous columns"
            )
        # print("* Finished BTL data station: " + ssscc + " *")

    # Drop duplicated columns generated by concatenation
    df_data_all = df_data_all.loc[:, ~df_data_all.columns.duplicated()]

    df_data_all["master_index"] = range(len(df_data_all))

    return df_data_all


def _reft_loader(ssscc, reft_dir):
    # semi-flexible search for reft file (in the form of *ssscc.cap)
    try:
        reft_path = sorted(Path(reft_dir).glob(f"*{ssscc}.cap"))[0]
    except IndexError:
        raise FileNotFoundError

    # this works better than pd.read_csv as format is semi-inconsistent (cf .cap files)
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
    reftDF[1] = reftDF[[1, 2, 3, 4]].agg(" ".join, axis=1)  # dd/mm/yy/time cols are
    reftDF.drop(columns=[2, 3, 4], inplace=True)  # read separately; combine into one

    columns = OrderedDict(  # having this as a dict streamlines next steps
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

    # assign initial flags (large "diff" = unstable reading, flag questionable)
    reftDF["REFTMP_FLAG_W"] = 2
    reftDF.loc[reftDF["diff"].abs() >= 3000, "REFTMP_FLAG_W"] = 3

    # add in STNNBR, CASTNO columns
    # TODO: should these be objects or floats? be consistent!
    # string prob better for other sta/cast formats (names, letters, etc.)
    reftDF["STNNBR"] = ssscc[0:3]
    reftDF["CASTNO"] = ssscc[3:5]
    return reftDF


def process_reft(ssscc_list, reft_dir=cfg.dirs["reft"]):
    # TODO: import reft_dir from a config file
    """
    SBE35 reference thermometer processing function. Load in .cap files for given
    station/cast list, perform basic flagging, and export to .csv files.

    Inputs
    ------
    ssscc_list : list of str
        List of stations to process
    reft_dir : str, optional
        Path to folder containing raw salt files (defaults to data/reft/)

    """
    for ssscc in ssscc_list:
        if not Path(reft_dir + ssscc + "_reft.csv").exists():
            try:
                reftDF = _reft_loader(ssscc, reft_dir)
                reftDF.to_csv(reft_dir + ssscc + "_reft.csv", index=False)
            except FileNotFoundError:
                log.warning(
                    "refT file for cast " + ssscc + " does not exist... skipping"
                )
                return


def add_btlnbr_cols(df, btl_num_col):
    df["BTLNBR"] = df[btl_num_col].astype(int)
    # default to everything being good
    df["BTLNBR_FLAG_W"] = 2
    return df


def load_hy_file(path_to_hyfile):
    df = pd.read_csv(path_to_hyfile, comment="#", skiprows=[0])
    df = df[df["EXPOCODE"] != "END_DATA"]
    return df


def export_report_data(df):

    df["STNNBR"] = [int(x[0:3]) for x in df["SSSCC"]]
    df["CTDPRS"] = df["CTDPRS"].round(1)
    cruise_report_cols = [
        "STNNBR",
        "CTDPRS",
        "CTDTMP1",
        "CTDTMP1_FLAG_W",
        "CTDTMP2",
        "CTDTMP2_FLAG_W",
        "REFTMP",
        "CTDCOND1",
        "CTDCOND1_FLAG_W",
        "CTDCOND2",
        "CTDCOND2_FLAG_W",
        "BTLCOND",
        "CTDSAL",
        "CTDSAL_FLAG_W",
        "SALNTY",
        "CTDOXY",
        "CTDOXY_FLAG_W",
        "CTDRINKO",
        "CTDRINKO_FLAG_W",
        "OXYGEN",
    ]

    # add in missing flags
    df["CTDTMP1_FLAG_W"] = flagging.by_residual(
        df["CTDTMP1"], df["REFTMP"], df["CTDPRS"]
    )
    df["CTDTMP2_FLAG_W"] = flagging.by_residual(
        df["CTDTMP1"], df["REFTMP"], df["CTDPRS"]
    )
    df["CTDCOND1_FLAG_W"] = flagging.by_residual(
        df["CTDCOND1"], df["BTLCOND"], df["CTDPRS"]
    )
    df["CTDCOND2_FLAG_W"] = flagging.by_residual(
        df["CTDCOND2"], df["BTLCOND"], df["CTDPRS"]
    )
    df["CTDOXY_FLAG_W"] = flagging.by_percent_diff(df["CTDOXY"], df["OXYGEN"])
    df["CTDRINKO_FLAG_W"] = flagging.by_percent_diff(df["CTDRINKO"], df["OXYGEN"])

    df[cruise_report_cols].to_csv("data/scratch_folder/report_data.csv", index=False)

    return


def export_hy1(df, out_dir=cfg.dirs["pressure"], org="ODF"):
    log.info("Exporting bottle file")
    btl_data = df.copy()
    now = datetime.now()
    file_datetime = now.strftime("%Y%m%d")

    # TODO: move to config; integrate Barna's "params" package instead?
    btl_columns = {
        "EXPOCODE": "",
        "SECT_ID": "",
        "STNNBR": "",
        "CASTNO": "",
        "SAMPNO": "",
        "BTLNBR": "",
        "BTLNBR_FLAG_W": "",
        "DATE": "",
        "TIME": "",
        "LATITUDE": "",
        "LONGITUDE": "",
        "DEPTH": "METERS",
        "CTDPRS": "DBAR",
        "CTDTMP": "ITS-90",
        "CTDSAL": "PSS-78",
        "CTDSAL_FLAG_W": "",
        "SALNTY": "PSS-78",
        "SALNTY_FLAG_W": "",
        # "CTDOXY": "UMOL/KG",
        # "CTDOXY_FLAG_W": "",
        # "CTDRINKO": "UMOL/KG",
        # "CTDRINKO_FLAG_W": "",
        "CTDOXY": "UMOL/KG",
        "CTDOXY_FLAG_W": "",
        "OXYGEN": "UMOL/KG",
        "OXYGEN_FLAG_W": "",
        "REFTMP": "ITS-90",
        "REFTMP_FLAG_W": "",
    }

    # rename outputs as defined in user_settings.yaml
    for param, attrs in cfg.ctd_outputs.items():
        if param not in btl_data.columns:
            btl_data.rename(columns={attrs["sensor"]: param}, inplace=True)

    btl_data["EXPOCODE"] = cfg.expocode
    btl_data["SECT_ID"] = cfg.section_id
    btl_data["STNNBR"] = [int(x[0:3]) for x in btl_data["SSSCC"]]
    btl_data["CASTNO"] = [int(x[3:]) for x in btl_data["SSSCC"]]
    btl_data["SAMPNO"] = btl_data["btl_fire_num"].astype(int)
    btl_data = add_btlnbr_cols(btl_data, btl_num_col="btl_fire_num")

    # sort by decreasing sample number (increasing pressure) and reindex
    btl_data = btl_data.sort_values(
        by=["STNNBR", "SAMPNO"], ascending=[True, False], ignore_index=True
    )

    # switch oxygen primary sensor to rinko
    btl_data["CTDOXY"] = btl_data.loc[:, "CTDRINKO"]
    btl_data["CTDOXY_FLAG_W"] = btl_data.loc[:, "CTDRINKO_FLAG_W"]

    # round data
    # for col in ["CTDTMP", "CTDSAL", "SALNTY", "REFTMP"]:
    #     btl_data[col] = btl_data[col].round(4)
    # for col in ["CTDPRS", "CTDOXY", "OXYGEN"]:
    #     btl_data[col] = btl_data[col].round(1)

    # add depth
    depth_df = pd.read_csv(
        cfg.dirs["logs"] + "depth_log.csv", dtype={"SSSCC": str}, na_values=-999
    ).dropna()
    manual_depth_df = pd.read_csv(
        cfg.dirs["logs"] + "manual_depth_log.csv", dtype={"SSSCC": str}
    )
    full_depth_df = pd.concat([depth_df, manual_depth_df])
    full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)
    btl_data["DEPTH"] = -999
    for index, row in full_depth_df.iterrows():
        btl_data.loc[btl_data["SSSCC"] == row["SSSCC"], "DEPTH"] = int(row["DEPTH"])

    # deal with nans
    # TODO: missing REFTMP not obvious til loading data - where to put this?
    # _reft_loader() is not the right place
    # maybe during loading step flag missing OXYGEN, REFTMP, BTLCOND?
    btl_data["REFTMP_FLAG_W"] = flagging.nan_values(
        btl_data["REFTMP_FLAG_W"], old_flags=btl_data["REFTMP_FLAG_W"]
    )
    btl_data = btl_data.where(~btl_data.isnull(), -999)

    # check columns
    try:
        btl_data[btl_columns.keys()]
        # this is lazy, do better
    except KeyError as err:
        log.info("Column names not configured properly... attempting to correct")
        bad_cols = err.args[0].split("'")[1::2]  # every other str is a column name
        for col in bad_cols:
            if col.endswith("FLAG_W"):
                log.warning(col + " missing, flagging with 9s")
                btl_data[col] = 9
            else:
                log.warning(col + " missing, filling with -999s")
                btl_data[col] = -999

    btl_data = btl_data[btl_columns.keys()]
    time_stamp = file_datetime + org
    with open(out_dir + cfg.expocode + "_hy1.csv", mode="w+") as f:
        f.write("BOTTLE, %s\n" % (time_stamp))
        f.write(",".join(btl_columns.keys()) + "\n")
        f.write(",".join(btl_columns.values()) + "\n")
        btl_data.to_csv(f, header=False, index=False)
        f.write("\n" + "END_DATA")

    return

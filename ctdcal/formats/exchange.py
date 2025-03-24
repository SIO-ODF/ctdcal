"""
Functions to support importing from, converting to and exporting exchange files.
"""
import logging
from datetime import datetime as dt

import pandas as pd

from ctdcal import get_ctdcal_config
from ctdcal.flagging.flag_common import nan_values

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def load_hy_file(path_to_hyfile):
    """
    Read in an exchange-formatted bottle file as a Pandas DataFrame.

    Inputs
    ------
    path_to_hyfile : String or Path object
        The path to the bottle file.

    Returns
    -------
    df : Pandas DataFrame
        The bottle file without the lead/end rows, comments, or units
    """

    df = pd.read_csv(path_to_hyfile, comment="#", skiprows=[0])
    df = df.drop(df.index[0])  #   Drop the units
    df = df[df["EXPOCODE"] != "END_DATA"]  #   Drop the final row
    return df


def merge_hy1(
    df1,
    df2,
):
    """
    Merges two hy1 files, returning the combined Pandas DataFrame.
    If the hy1 file has not been loaded yet, use load_hy_file.

    Inputs
    -------
    df1 : Pandas DataFrame
        First hy1 file for concatination
    df2 : Pandas DataFrame
        Second hy1 file for concatination

    Returns
    df: Pandas DataFrame
        Merged bottle file as a DataFrame
    """

    if set(df1.columns) != set(df2.columns):
        print("Bottle file columns do not match. Concatenating with NaNs.")

    df = pd.concat([df1, df2], axis=0, ignore_index=True)  #   Staple df2 onto there

    sorting_cols = {"STNNBR", "CASTNO", "SAMPNO"}
    if sorting_cols.issubset(df):
        if df[list(sorting_cols)].isna().any().any():
            print("NaNs found in station/cast/sample number. Check source files.")

        else:
            df = df.sort_values(
                by=["STNNBR", "CASTNO", "SAMPNO"],
                ascending=[True, True, False],
                ignore_index=True,
            )
    return df


def add_btlnbr_cols(df, btl_num_col):
    """
    Initialize bottle number column and initialize WOCE bottle flags.

    Parameters
    ----------
    df : Pandas DataFrame
        Bottle DataFrame containing a defined rosette bottle number
    btl_num_col : String
        String of bottle column to be reassigned

    Returns
    -------
    df : Pandas DataFrame
        Bottle DataFrame with BTLNBR and flag columns as type int
    """
    df["BTLNBR"] = df[btl_num_col].astype(int)
    # default to everything being good
    df["BTLNBR_FLAG_W"] = 2
    return df


def export_hy1(df, out_dir=cfg.dirs["pressure"], org="ODF"):
    """
    Write out the exchange-lite formatted hy1 bottle file.

    Params
    ------
    df : Pandas DataFrame
        Fit bottle data
    out_dir = String or Path object, optional
        The path for where to write the hy1 file
    org : String, optional
        The organization or group used to determine subroutines

    """
    log.info("Exporting bottle file")
    btl_data = df.copy()
    now = dt.now()
    file_datetime = now.strftime("%Y%m%d")

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
    btl_data["REFTMP_FLAG_W"] = nan_values(
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

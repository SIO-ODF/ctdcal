"""
Functions to support importing from, converting to and exporting exchange files.
"""
import logging
from datetime import datetime as dt, timezone

import numpy as np
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

def roll_filter(df, p_col="CTDPRS", direction="down"):
    """
    Filter out heaving in CTD data due to ship rolls.

    Parameters
    ----------
    df : DataFrame
        CTD data
    p_col : str, optional
        Name of pressure column
    direction : str, optional
        Direction of cast (i.e. "down" or "up" cast)

    Returns
    -------
    var2 : dtype
        var description

    """
    if direction == "down":
        monotonic_sequence = df[p_col].expanding().max()
    elif direction == "up":
        monotonic_sequence = df[p_col].expanding().min()
    else:
        raise ValueError("direction must be one of (up, down)")

    return df[df[p_col] == monotonic_sequence]


def pressure_sequence(df, p_col="CTDPRS", direction="down"):
    """
    Convert CTD time series to a pressure series.

    Parameters
    ----------
    df : DataFrame
        CTD time series data
    p_col : str, optional
        Name of pressure column
    direction : str, optional
        Direction to sequence data

    Returns
    -------
    df_binned : DataFrame
        Pressure binned CTD data
    """
    # change to take dataframe with the following properties
    # * in water data only (no need to find cast start/end)
    # * The full down and up time series (not already split since this method will do it)
    # New "algorithm"
    # * if direction is "down", use the input as is
    # * if direction is "up", invert the row order of the input dataframe
    # Use the "roll filter" method to get only the rows to be binned
    # * the roll filter will treat the "up" part of the cast as a giant roll to be filtered out
    # * the reversed dataframe will ensure we get the "up" or "down" part of the cast
    # * there is no need to reverse the dataframe again as the pressure binning process will remove any "order" information (it doesn't care about the order)
    # That's basically all I (barna) have so far

    df_filtered = roll_filter(df, p_col, direction=direction)

    # 04/11/21 MK: this is not behaving properly or the order is wrong?
    # Currently this function fills the top-most good CTD value *before* bin avg,
    # so those bin averages are not the same as the first good binned value.
    # df_filled = _fill_surface_data(df_filtered, bin_size=2)

    df_binned = binning_df(df_filtered, bin_size=2)
    fill_rows = df_binned["CTDPRS"].isna()
    df_binned.loc[fill_rows, "CTDPRS"] = df_binned[fill_rows].index.to_numpy()
    df_binned.bfill(inplace=True)
    df_binned.loc[:, "interp_bool"] = False
    df_binned.loc[fill_rows, "interp_bool"] = True
    df_filled = _flag_backfill_data(df_binned).drop(columns="interp_bool")

    return df_filled.reset_index(drop=True)


def binning_df(df, p_column="CTDPRS", bin_size=2):
    """Calculate the bin-mean of each column in input dataframe

    Parameters
    ----------
    df : DataFrame
        Data to be bin-meaned
    p_column : str, optional
        Pressure column name to use for binning
    bin_size : int, optional
        Width of bins (in decibars)

    Returns
    -------
    df_out : DataFrame
        Bin-meaned data

    """
    if p_column not in df.columns:
        raise KeyError(f"{p_column} column missing from dataframe")

    p_max = np.ceil(df[p_column].max())
    labels = np.arange(0, p_max, bin_size)
    bin_edges = np.arange(0, p_max + bin_size, bin_size)
    df_out = df.copy()
    df_out.loc[:, "bins"] = pd.cut(
        df[p_column], bins=bin_edges, right=False, include_lowest=True, labels=labels
    )
    df_out.loc[:, p_column] = df_out["bins"].astype(float)

    try:
        #   Python 3.8
        return df_out.groupby("bins", observed=False).mean()
    except TypeError:
        return df_out.groupby("bins", observed=False).mean(numeric_only=True)


def _fill_surface_data(df, bin_size=2):
    """Copy first scan from top of cast and propagate up to surface at bin centers"""
    df = df.copy(deep=True)
    p_min = np.floor(df["CTDPRS"].iloc[0])
    df_surface = pd.DataFrame({"CTDPRS": np.arange(bin_size / 2, p_min, bin_size)})
    df_surface["interp_bool"] = True
    df = df_surface.merge(df, on="CTDPRS", how="outer")
    df["interp_bool"].fillna(False, inplace=True)
    df = _flag_backfill_data(df).drop(columns="interp_bool")

    return df.bfill()


def manual_backfill(df, p_cutoff, p_col="CTDPRS", flag_suffix="_FLAG_W"):
    """
    Overwrite values below cutoff pressure by backfilling the first data point past
    threshold upward to the surface. Backfilled data are flagged 6.

    Parameters
    ----------
    df : DataFrame
        Input data
    p_cutoff : float
        Cutoff pressure for backfilling
    p_col : str, optional
        Name of pressure column in df
    flag_suffix : str, optional
        Parameter suffix for data flags

    Returns
    -------
    df : DataFrame
        Input DataFrame with backfilled data
    """
    df = df.copy(deep=True)
    cols = df.columns.drop(p_col)
    df["interp_bool"] = df[p_col] < p_cutoff
    df.loc[df["interp_bool"], cols] = np.nan
    df = _flag_backfill_data(df).drop(columns="interp_bool")

    return df.bfill()


def _flag_backfill_data(
    df, p_col="CTDPRS", flag_bool_col="interp_bool", flag_suffix="_FLAG_W"
):
    """Flag data columns which have been interpolated with flag 6."""
    for col in df.columns:
        if flag_suffix in col:
            df.loc[df[flag_bool_col], col] = 6

    return df


def export_ct1(df, ssscc_list):
    """
    Export continuous CTD (i.e. time) data to data/pressure/ directory as well as
    adding quality flags and removing unneeded columns.

    Parameters
    ----------
    df : DataFrame
        Continuous CTD data
    ssscc_list : list of str
        List of stations to export

    Returns
    -------

    Notes
    -----
    Needs depth_log.csv and manual_depth_log.csv to run successfully

    """
    log.info("Exporting CTD files")

    # initial flagging (some of this should be moved)
    df["CTDFLUOR_FLAG_W"] = 1
    df["CTDXMISS_FLAG_W"] = 1
    # df["CTDBACKSCATTER_FLAG_W"] = 1

    # rename outputs as defined in user_settings.yaml
    for param, attrs in cfg.ctd_outputs.items():
        if param not in df.columns:
            df.rename(columns={attrs["sensor"]: param}, inplace=True)

    # check that all columns are there
    try:
        df[cfg.ctd_col_names]
        # this is lazy, do better
    except KeyError as err:
        log.info("Column names not configured properly... attempting to correct")
        bad_cols = err.args[0].split("'")[1::2]  # every other str is a column name
        for col in bad_cols:
            if col.endswith("FLAG_W"):
                log.warning(col + " missing, flagging with 9s")
                df[col] = 9
            else:
                log.warning(col + " missing, filling with -999s")
                df[col] = -999

    df["SSSCC"] = df["SSSCC"].astype(str).copy()
    cast_details = pd.read_csv(
        # cfg.dirs["logs"] + "cast_details.csv", dtype={"SSSCC": str}
        cfg.dirs["logs"] + "bottom_bottle_details.csv",
        dtype={"SSSCC": str},
    )
    depth_df = pd.read_csv(
        cfg.dirs["logs"] + "depth_log.csv", dtype={"SSSCC": str}, na_values=-999
    ).dropna()
    try:
        manual_depth_df = pd.read_csv(
            cfg.dirs["logs"] + "manual_depth_log.csv", dtype={"SSSCC": str}
        )
    except FileNotFoundError:
        log.warning("manual_depth_log.csv not found... duplicating depth_log.csv")
        manual_depth_df = depth_df.copy()  # write manual_depth_log as copy of depth_log
        manual_depth_df.to_csv(cfg.dirs["logs"] + "manual_depth_log.csv", index=False)
    full_depth_df = pd.concat([depth_df, manual_depth_df])
    full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)

    for ssscc in ssscc_list:

        time_data = df[df["SSSCC"] == ssscc].copy()
        time_data = pressure_sequence(time_data)
        # switch oxygen primary sensor to rinko
        # if int(ssscc[:3]) > 35:
        print(f"Using Rinko as CTDOXY for {ssscc}")
        time_data.loc[:, "CTDOXY"] = time_data["CTDRINKO"]
        time_data.loc[:, "CTDOXY_FLAG_W"] = time_data["CTDRINKO_FLAG_W"]
        time_data = time_data[cfg.ctd_col_names]
        # time_data = time_data.round(4)
        time_data = time_data.where(~time_data.isnull(), -999)  # replace NaNs with -999

        # force flags back to int
        for col in time_data.columns:
            if col.endswith("FLAG_W"):
                time_data[col] = time_data[col].astype(int)

        try:
            depth = full_depth_df.loc[full_depth_df["SSSCC"] == ssscc, "DEPTH"].iloc[0]
        except IndexError:
            log.warning(f"No depth logged for {ssscc}, setting to -999")
            depth = -999

        # get cast_details for current SSSCC
        cast_dict = cast_details[cast_details["SSSCC"] == ssscc].to_dict("records")[0]
        b_datetime = (
            dt.fromtimestamp(cast_dict["bottom_time"], tz=timezone.utc)
            .strftime("%Y%m%d %H%M")
            .split(" ")
        )
        btm_lat = cast_dict["latitude"]
        btm_lon = cast_dict["longitude"]

        now = dt.now(timezone.utc)
        file_datetime = now.strftime("%Y%m%d")  # %H:%M")
        file_datetime = file_datetime + "ODFSIO"
        with open(f"{cfg.dirs['pressure']}{ssscc}_ct1.csv", "w+") as f:
            # put in logic to check columns?
            # number_headers should be calculated, not defined
            ctd_header = (  # this is ugly but prevents tabs before label
                f"CTD,{file_datetime}\n"
                f"NUMBER_HEADERS = 11\n"
                f"EXPOCODE = {cfg.expocode}\n"
                f"SECT_ID = {cfg.section_id}\n"
                f"STNNBR = {ssscc[:3]}\n"  # STNNBR = SSS
                f"CASTNO = {ssscc[3:]}\n"  # CASTNO = CC
                f"DATE = {b_datetime[0]}\n"
                f"TIME = {b_datetime[1]}\n"
                f"LATITUDE = {btm_lat:.4f}\n"
                f"LONGITUDE = {btm_lon:.4f}\n"
                f"INSTRUMENT_ID = {cfg.ctd_serial}\n"
                f"DEPTH = {depth:.0f}\n"
            )
            f.write(ctd_header)
            np.asarray(cfg.ctd_col_names).tofile(f, sep=",", format="%s")
            f.write("\n")
            np.asarray(cfg.ctd_col_units).tofile(f, sep=",", format="%s")
            f.write("\n")
            time_data.to_csv(f, header=False, index=False)
            f.write("END_DATA")

#!/usr/bin/env python
import csv
import logging
import math
import warnings
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import scipy.signal as sig

from . import flagging as flagging
from . import get_ctdcal_config
from . import oxy_fitting as oxy_fitting
from . import report_ctd as report_ctd

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

warnings.filterwarnings("ignore", "Mean of empty slice.")


def cast_details(df, ssscc, log_file=None):
    """
    We determine the cast details using pandas magic.
    First find alternating periods of pumps on and pumps off, then select the
    pumps on period with the highest pressure. Get values from the row with the
    highest pressure, and return all values to be sent to log.

    Parameters
    ----------
    df : DataFrame
        Filtered CTD data
    ssscc : integer
        The station and cast, as SSSCC format
    log_file : file handle or string
        File destination for cast details

    Returns
    -------
    df_downcast : DataFrame
        CTD data with the soak period and upcast trimmed off

    Notes
    -----
    The following (float) variables are output to log_file:
    time_start : Time at start of cast (in unix epoch time)
    time_end : Time at end of cast (in unix epoch time)
    time_bottom : Time at bottom of cast (in unix epoch time)
    p_start : Pressure at which cast started
    p_max : Bottom of the cast pressure
    b_lat : Latitude at bottom of cast
    b_lon : Longitude at bottom of cast
    b_alt : Altimeter reading at bottom of cast
    """
    df_cast = _trim_soak_period(df)

    # TODO: call parameters from config file instead
    p_start = float(np.around(df_cast["CTDPRS"].head(1), 4))
    p_max_ind = df_cast["CTDPRS"].argmax()
    p_max = float(np.around(df_cast["CTDPRS"].max(), 4))
    time_start = float(df_cast["scan_datetime"].head(1))
    time_end = float(df_cast["scan_datetime"].tail(1))
    time_bottom = float(df_cast["scan_datetime"][p_max_ind])
    b_lat = float(np.around(df_cast["GPSLAT"][p_max_ind], 4))
    b_lon = float(np.around(df_cast["GPSLON"][p_max_ind], 4))
    b_alt = float(np.around(df_cast["ALT"][p_max_ind], 4))

    report_ctd.report_cast_details(
        ssscc,
        log_file,
        time_start,
        time_end,
        time_bottom,
        p_start,
        p_max,
        b_alt,
        b_lat,
        b_lon,
    )

    # remove upcast
    df_downcast = df_cast[:p_max_ind].copy()

    return df_downcast


def _trim_soak_period(df=None):
    """
    1) Find pump on/off patterns
    2) Select pump_on=True group with largest pressure recording
    3) Find soak period before start of downcast
    4) Trim cast, return everything after top of cast (i.e. minimum pressure)
    """
    df_list = [
        g for i, g in df.groupby(df["pump_on"].ne(df["pump_on"].shift()).cumsum())
    ]
    df_pump_on_list = [df for df in df_list if df["pump_on"].all()]
    df_cast = df_pump_on_list[np.argmax([df["CTDPRS"].max() for df in df_pump_on_list])]
    df_cast = df_cast.reset_index(drop=True)
    # next fn deals w/ edge cases, leave as is for now
    df_cast = _find_last_soak_period(df_cast)
    start_ind = df_cast.loc[: len(df) // 4, "CTDPRS"].argmin()
    df_trimmed = df_cast[start_ind:].reset_index(drop=True).copy()

    return df_trimmed


def _find_last_soak_period(df_cast, time_bin=8, P_surface=2, P_downcast=50):
    """
    Find the soak period before the downcast starts.

    The algorithm is tuned for repeat hydrography work, specifically US GO-SHIP
    parameters. This assumes the soak depth will be somewhere between 10 and 30
    meters, the package will sit at the soak depth for at least 20 to 30 seconds
    before starting ascent to the surface and descent to target depth.

    The algorithm is not guaranteed to catch the exact start of the soak period,
    but within a minimum period of time_bin seconds(?) from end of the soak if
    the soak period assumption is valid. This should be shorter than the total
    soak period time, and able to catch the following rise and descent of the
    package that signals the start of the cast.

    The algorithm has been designed to handle four general cases of casts:
        * A routine cast with pumps turning on in water and normal soak
        * A cast where the pumps turn on in air/on deck
        * A cast where the pumps turn on and off due to rosette coming out of water
        * A cast where there are multiple stops on the downcast to the target depth

    Parameters
    ----------
    df_cast : DataFrame
        DataFrame of the entire cast, from deckbox on to deckbox off
    time_bin : integer, optional
        Number of seconds to bin average for descent rate calculation
    P_surface : integer, optional
        Minimum surface pressure threshold required to look for soak depth
        (2 dbar was chosen as an average rosette is roughly 1.5 to 2 meters tall)
    P_downcast : integer, optional
        Minimum pressure threshold required to assume downcast has started
        (50 dbar has been chosen as double the deep soak depth of 20-30 dbar)

    Returns
    -------
    df_cast_trimmed : DataFrame
        DataFrame starting within time_bin seconds of the last soak period.
    """
    # Validate user input
    if time_bin <= 0:
        raise ValueError("Time bin value should be positive whole seconds.")
    if P_downcast <= 0:
        raise ValueError(
            "Starting downcast pressure threshold must be positive integers."
        )
    if P_downcast < P_surface:
        raise ValueError(
            "Starting downcast pressure threshold must be greater \
                        than surface pressure threshold."
        )

    # If pumps have not turned on until in water, return DataFrame
    if df_cast.iloc[0]["CTDPRS"] > P_surface:
        return df_cast

    # Bin the data by time, and compute the average rate of descent
    df_cast["index"] = df_cast.index  # needed at end to identify start_idx
    df_cast["bin"] = pd.cut(
        df_cast.index,
        np.arange(df_cast.index[0], df_cast.index[-1], time_bin * 24),
        labels=False,
        include_lowest=True,
    )
    df_binned = df_cast.groupby("bin").mean()

    # Compute difference of descent rates and label bins
    df_binned["dP"] = df_binned["CTDPRS"].diff().fillna(0).round(0)
    df_binned["movement"] = pd.cut(
        df_binned["dP"], [-1000, -0.5, 0.5, 1000], labels=["up", "stop", "down"]
    )

    # Find all periods where the rosette is not moving
    df_group = df_binned.groupby(
        df_binned["movement"].ne(df_binned["movement"].shift()).cumsum()
    )
    df_list = [g for i, g in df_group]

    # Find last soak period before starting descent to target depth
    def find_last(df_list, P_downcast):
        for idx, df in enumerate(df_list):
            if df["CTDPRS"].max() < P_downcast:
                # make sure it's soak, not a stop to switch to autocast (i.e. A20 2021)
                # TODO: try instead finding max depth then working backwards?
                if df.max()["movement"] == "stop" and len(df) > 1:
                    last_idx = idx
            else:
                return last_idx
        return last_idx

    # Trim off everything before last soak
    start_idx = int(df_list[find_last(df_list, P_downcast)].head(1)["index"])
    df_cast_trimmed = df_cast.loc[start_idx:].reset_index()

    return df_cast_trimmed


def ctd_align(inMat=None, col=None, time=0.0):
    """ctd_align function

    Function takes full NUMPY ndarray with predefined dtype array
    and adjusts time of sensor responce and water flow relative to
    the time frame of temperature sensor.

    Originally written by Courtney Schatzman, docstring by Joseph Gum.
    Need to generate alignment plots in order to properly use ctd_align.

    Args:
        param1 (ndarray): inMat, numpy ndarray with dtype array
        param2 (float): col, column to apply time advance to.
        param3 (float): time, advance in seconds to apply to raw data.

    Returns:
        Narray: The return value is ndarray with adjusted time of parameter
          specified.

    """
    # Num of frames per second.
    fl = 24

    if (inMat is not None) & (col is not None) & (time > 0.0):
        # Time to advance
        advnc = int(fl * time)
        tmp = np.arange(advnc, dtype=np.float)
        last = inMat[col][len(inMat) - 1]
        tmp.fill(float(last))
        inMat[col] = np.concatenate((inMat[col][advnc:], tmp))

    return inMat


def raw_ctd_filter(df=None, window="triangle", win_size=24, parameters=None):
    """
    Filter raw CTD data using one of three window types (boxcar, hanning, triangle).

    Parameters
    ----------
    df : DataFrame
        Raw CTD data
    window : str, optional
        Type of filter window
    win_size : int, optional
        Length of window in number of samples
    parameters : list of str, optional
        List of DataFrame columns to be filtered

    Returns
    -------
    filtered_df : DataFrame
        CTD data with filtered parameters
    """

    filter_df = df.copy()
    if parameters is not None:
        for p in parameters:
            if window == "boxcar":
                win = sig.boxcar(win_size)
            elif window == "hanning":
                win = sig.hann(win_size)
            elif window == "triangle":
                win = sig.triang(win_size)
            filter_df[p] = sig.convolve(filter_df[p], win, mode="same") / np.sum(win)

    return filter_df


def remove_on_deck(df, stacast, cond_startup=20.0, log_file=None):
    """
    Find and remove times when rosette is on deck.
    Optionally log average pressure at start and end of cast.

    Parameters
    ----------
    df : DataFrame
        Raw CTD data
    stacast : str
        Station/cast name
    cond_startup : float, optional
        Minimum conductivity (units?) threshold indicating rosette is in water
    log_file : str, optional
        Path and filename to save start/end deck pressure values

    Returns
    -------
    trimmed_df : DataFrame
        Raw CTD data trimmed to times when rosette is in water
    """
    # TODO: move these to config file
    # Frequency
    fl = 24
    fl2 = fl * 2
    # Half minute
    ms = 30
    time_delay = fl * ms  # time to let CTD pressure reading settle/sit on deck

    # split dataframe into upcast/downcast
    downcast = df.iloc[: (df["CTDPRS"].argmax() + 1)]
    upcast = df.iloc[(df["CTDPRS"].argmax() + 1) :]

    # Search each half of df for minimum conductivity
    # threshold to identify when rosette is out of water
    start_df = downcast.loc[
        (downcast[cfg.column["c1"]] < cond_startup)
        & (downcast[cfg.column["c2"]] < cond_startup),
        cfg.column["p"],
    ]
    end_df = upcast.loc[
        (upcast[cfg.column["c1"]] < cond_startup)
        & (upcast[cfg.column["c2"]] < cond_startup),
        cfg.column["p"],
    ]

    # Evaluate starting and ending pressures
    start_samples = len(start_df)
    if start_samples > time_delay:
        start_p = np.average(start_df.iloc[fl2 : (start_samples - time_delay)])
    else:
        start_seconds = start_samples / fl
        log.warning(
            f"{stacast}: Only {start_seconds:0.1f} seconds of start pressure averaged."
        )
        start_p = np.average(start_df.iloc[fl2:start_samples])

    end_samples = len(end_df)
    if end_samples > time_delay:
        end_p = np.average(end_df.iloc[(time_delay):])
    else:
        end_seconds = end_samples / fl
        log.warning(
            f"{stacast}: Only {end_seconds:0.1f} seconds of end pressure averaged."
        )
        end_p = np.average(end_df)  # just average whatever there is

    # Remove ondeck start and end pressures
    if len(start_df) == 0:
        log.warning("Failed to find starting deck pressure.")
        for n in [1, 2]:
            try:
                (downcast[cfg.column[f"c{n}"]] < cond_startup).value_counts()[True]
            except KeyError:
                log.warning(
                    f"No values below {cond_startup} found for {cfg.column[f'c{n}']}"
                )
        breakpoint()
    if len(end_df) == 0:
        log.warning("Failed to find ending deck pressure.")
        for n in [1, 2]:
            try:
                (upcast[cfg.column[f"c{n}"]] < cond_startup).value_counts()[True]
            except KeyError:
                log.warning(
                    f"No values below {cond_startup} found for {cfg.column[f'c{n}']}"
                )
        breakpoint()
    # MK (3/23/20, 11am):
    # auto end calculation failed bc cond2 is still >30
    # may have to do manually or just use cond1 for station 00901
    trimmed_df = df.iloc[start_df.index.max() : end_df.index.min()].copy()

    # Log ondeck pressures
    if log_file is not None:
        report_ctd.report_pressure_details(stacast, log_file, start_p, end_p)

    return trimmed_df


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
    # TODO: optional start/end pressure values?

    # change to take dataframe with the following properties
    # * in water data only (no need to find cast start/end)
    # * The full down and up time series (not already split since this method will do it)
    # New "algorithm" (TODO spell this right)
    # * if direction is "down", use the input as is
    # * if direction is "up", invert the row order of the input dataframe
    # Use the "roll filter" method to get only the rows to be binned
    # * the roll filter will treat the "up" part of the cast as a giant roll to be filtered out
    # * the reversed dataframe will ensure we get the "up" or "down" part of the cast
    # * there is no need to reverse the dataframe again as the pressure binning process will remove any "order" information (it doesn't care about the order)
    # That's basically all I (barna) have so far TODO Binning, etc...
    # pandas.cut() to do binning

    df_filtered = roll_filter(df, p_col, direction=direction)

    # 04/11/21 MK: this is not behaving properly or the order is wrong?
    # Currently this function fills the top-most good CTD value *before* bin avg,
    # so those bin averages are not the same as the first good binned value.
    # df_filled = _fill_surface_data(df_filtered, bin_size=2)

    df_binned = binning_df(df_filtered, bin_size=2)  # TODO: abstract bin_size in config
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

    return df_out.groupby("bins").mean()


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


def _get_pressure_offset(start_vals, end_vals):
    """
    Finds unique values and calculate mean for pressure offset

    Parameters
    ----------
    start_vals : array_like
        Array of initial ondeck pressure values

    end_vals : array_like
        Array of ending ondeck pressure values
    Returns
    -------
    p_off : float
         Average pressure offset

    """
    p_start = pd.Series(np.unique(start_vals))
    p_end = pd.Series(np.unique(end_vals))
    p_start = p_start[p_start.notnull()]
    p_end = p_end[p_end.notnull()]
    p_off = p_start.mean() - p_end.mean()

    # JACKSON THINKS THIS METHOD SHOULD BE USED TO KEEP START END PAIRS
    #    p_df = pd.DataFrame()
    #    p_df['p_start'] = p_start
    #    p_df['p_end'] = p_end
    #    p_df = p_df[p_df['p_end'].notnull()]
    #    p_df = p_df[p_df['p_start'].notnull()]
    #    p_off = p_df['p_start'].mean() - p_df['p_end'].mean()
    ##########################################################

    p_off = np.around(p_off, decimals=4)

    return p_off


def apply_pressure_offset(df, p_col="CTDPRS"):
    # TODO: import p_col from config file
    """
    Calculate pressure offset using deck pressure log and apply it to the data.
    Pressure flag column is added with value 2, indicating the data are calibrated.

    Parameters
    ----------
    df : DataFrame
        DataFrame containing column with pressure values
    p_col : str, optional
        Pressure column name in DataFrame (defaults to CTDPRS)

    Returns
    -------
    df : DataFrame
        DataFrame containing updated pressure values and a new flag column

    """
    p_log = pd.read_csv(
        cfg.dirs["logs"] + "ondeck_pressure.csv",
        dtype={"SSSCC": str},
        na_values="Started in Water",
    )
    p_offset = _get_pressure_offset(p_log.ondeck_start_p, p_log.ondeck_end_p)
    df[p_col] += p_offset
    df[p_col + "_FLAG_W"] = 2

    return df


def make_depth_log(time_df, threshold=80):
    # TODO: get column names from config file
    """
    Create depth log file from maximum depth of each station/cast in time DataFrame.
    If rosette does not get within the threshold distance of the bottom, returns NaN.

    Parameters
    ----------
    time_df : DataFrame
        DataFrame containing continuous CTD data
    threshold : int, optional
        Maximum altimeter reading to consider cast "at the bottom" (defaults to 80)

    """
    # TODO: make inputs be arraylike rather than dataframe
    df = time_df[["SSSCC", "CTDPRS", "GPSLAT", "ALT"]].copy().reset_index()
    df_group = df.groupby("SSSCC", sort=False)
    idx_p_max = df_group["CTDPRS"].idxmax()
    bottom_df = pd.DataFrame(
        data={
            "SSSCC": df["SSSCC"].unique(),
            "max_p": df.loc[idx_p_max, "CTDPRS"],
            "lat": df.loc[idx_p_max, "GPSLAT"],
            "alt": df.loc[idx_p_max, "ALT"],
        }
    )
    bottom_df.loc[bottom_df["alt"] > threshold, "alt"] = np.nan
    # pandas 1.2.1 ufunc issue workaround with pd.to_numpy()
    bottom_df["DEPTH"] = (
        (
            bottom_df["alt"]
            + np.abs(gsw.z_from_p(bottom_df["max_p"], bottom_df["lat"].to_numpy()))
        )
        .fillna(value=-999)
        .round()
        .astype(int)
    )
    bottom_df[["SSSCC", "DEPTH"]].to_csv(
        cfg.dirs["logs"] + "depth_log.csv", index=False
    )

    return True


def make_ssscc_list(fname="data/ssscc.csv"):
    """
    Attempt to automatically generate list of station/casts from raw files.
    """
    raw_files = Path(cfg.dirs["raw"]).glob("*.hex")
    ssscc_list = sorted([f.stem for f in raw_files])
    pd.Series(ssscc_list).to_csv(fname, header=None, index=False, mode="x")

    return ssscc_list


def get_ssscc_list(fname="data/ssscc.csv"):
    """
    Load in list of stations/casts to process.
    """
    ssscc_list = []
    with open(fname, "r") as filename:
        ssscc_list = [line.strip() for line in filename]

    return ssscc_list


def load_all_ctd_files(ssscc_list):
    """
    Load CTD files for station/cast list and merge into a dataframe.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to load

    Returns
    -------
    df_data_all : DataFrame
        Merged dataframe containing all loaded data

    """
    df_list = []
    for ssscc in ssscc_list:
        log.info("Loading TIME data for station: " + ssscc + "...")
        time_file = cfg.dirs["time"] + ssscc + "_time.pkl"
        time_data = pd.read_pickle(time_file)
        time_data["SSSCC"] = str(ssscc)
        time_data["dv_dt"] = oxy_fitting.calculate_dV_dt(
            time_data["CTDOXYVOLTS"], time_data["scan_datetime"]
        )
        df_list.append(time_data)
        # print("** Finished TIME data station: " + ssscc + " **")
    df_data_all = pd.concat(df_list, axis=0, sort=False)

    df_data_all["master_index"] = range(len(df_data_all))

    return df_data_all


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
    # TODO: lump all uncalibrated together; smart flagging like ["CTD*_FLAG_W"] = 1
    # TODO: may not always have these channels so don't hardcode them in!
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
        # TODO: add logging; look into inheriting/extending a class to add features
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

        # force flags back to int (TODO: make flags categorical)
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
            datetime.fromtimestamp(cast_dict["bottom_time"], tz=timezone.utc)
            .strftime("%Y%m%d %H%M")
            .split(" ")
        )
        # TODO: yo-yo casts are an edge case where this may be different
        btm_lat = cast_dict["latitude"]
        btm_lon = cast_dict["longitude"]

        now = datetime.now(timezone.utc)
        file_datetime = now.strftime("%Y%m%d")  # %H:%M")
        file_datetime = file_datetime + "ODFSIO"
        # TODO: only "cast" needs to be int; "station" is explicitly allowed to incl.
        # letters/etc. Moving from SSSCC to station & cast fields will be beneficial
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

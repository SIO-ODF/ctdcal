#!/usr/bin/env python
import csv
import math
import warnings
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path

import config as cfg
import gsw
import numpy as np
import pandas as pd
import scipy.signal as sig

import ctdcal.oxy_fitting as oxy_fitting
import ctdcal.report_ctd as report_ctd

warnings.filterwarnings("ignore", 'Mean of empty slice.')

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
    df_downcast = df_cast[: p_max_ind].copy()

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


def _find_last_soak_period(df_cast, surface_pressure=2, time_bin=8, downcast_pressure=50):
    """Find the soak period before the downcast starts.

    The algorithm is tuned for repeat hydrography work, specifically US GO-SHIP
    parameters. This assumes the soak depth will be somewhere between 10 and 30
    meters, the package will sit at the soak depth for at least 20 to 30 seconds
    before starting ascent to the surface and descent to target depth.

    Parameters
    ----------
    df_cast : DataFrame
        DataFrame of the entire cast
    surface_pressure : integer
        Minimum surface pressure threshold required to look for soak depth.
        2 dbar was chosen as an average rosette is roughly 1.5 to 2 meters tall.
    time_bin : integer
        Time, in whole seconds.
    downcast_pressure : integer
        Minimum pressure threshold required to assume downcast has started.
        50 dbar has been chosen as double the deep soak depth of 20-30 dbar.

    Returns
    -------
    df_cast_ret : DataFrame
        DataFrame starting within time_bin seconds of the last soak period.

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

    """
    #Validate user input
    if time_bin <= 0:
        raise ValueError('Time bin value should be positive whole seconds.')
    if downcast_pressure <=0:
        raise ValueError('Starting downcast pressure threshold must be positive integers.')
    if downcast_pressure < surface_pressure:
        raise ValueError(f'Starting downcast pressure threshold must be greater \
                        than surface pressure threshold.')

    # If pumps have not turned on until in water, return DataFrame
    if df_cast.iloc[0]['CTDPRS'] > surface_pressure:
        return df_cast

    '''code_pruning: variables should always have relevant names, even if it needs to be changed later with find/replace'''
    #Bin the data by time, and compute the average rate of descent
    df_blah = df_cast.loc[:,:]
    df_blah['index'] = df_blah.index
    df_blah['bin'] = pd.cut(df_blah.loc[:,'index'],
                            range(df_blah.iloc[0]['index'],df_blah.iloc[-1]['index'],time_bin*24),
                            labels=False, include_lowest=True)
    df_blah2 = df_blah.groupby('bin').mean()

    #Compute difference of descent rates and label bins
    df_blah2['prs_diff'] = df_blah2['CTDPRS'].diff().fillna(0).round(0)
    df_blah2['movement'] = pd.cut(df_blah2['prs_diff'], [-1000,-0.5,0.5,1000], labels=['up','stop','down'])

    #Find all periods where the rosette is not moving
    df_stop = df_blah2.groupby('movement').get_group('stop')
    groupby_test = df_blah2.groupby(df_blah2['movement'].ne(df_blah2['movement'].shift()).cumsum())
    list_test = [g for i,g in groupby_test]

    #Find a dataframe index of the last soak period before starting descent
    def poop(list_obj, downcast_pressure):
        """ Return dataframe index in the last soak period before starting
            descent to target depth.
        """
        for i, x in zip(range(len(list_test)),list_test):
            if x['CTDPRS'].max() < downcast_pressure:
                if x.max()['movement'] == 'stop':
                    index = i
            if x['CTDPRS'].max() > downcast_pressure:
                return index
        return index

    #Truncate dataframe to new starting index : end of dataframe
    start_index = np.around(list_test[poop(list_test, downcast_pressure)].head(1)['index'])
    df_cast = df_cast.set_index('index')
    df_cast = df_cast.loc[int(start_index):,:]
    df_cast_ret = df_cast.reset_index()
    return df_cast_ret

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

    if (inMat is not None) & (col is not None) & ( time > 0.0):
        # Time to advance
        advnc = int(fl * time)
        tmp = np.arange(advnc, dtype=np.float)
        last = inMat[col][len(inMat)-1]
        tmp.fill(float(last))
        inMat[col] = np.concatenate((inMat[col][advnc:],tmp))

    return inMat


'''code_pruning: do timing exercises to see if this is vectorized or not'''
def hysteresis_correction(H1=-0.033, H2=5000, H3=1450, inMat = None):
    """Hysteresis Correction function

    Function takes data ndarray and hysteresis coefficiants
    and returns hysteresis corrected oxygen data.

    Args:
        param1 (float): H1, hysteresis correction coefficiant 1
        param2 (float): H2, hysteresis correction coefficiant 2
        param3 (float): H3, hysteresis correction coefficiant 3
        param5 (array): inMat, raw ctd data.

    Returns:
        array: Return dissolved oxygen hysteresis corrected data.

    .. REF PAGE:
       http://http://www.seabird.com/document/an64-3-sbe-43-dissolved-oxygen-do-sensor-hysteresis-corrections
    """
    Oxnewconc = np.arange(0,len(inMat),1)

    Oxnewconc[0] = inMat['o1_mll'][1]

    if inMat is None:
       print("Hysteresis Correction function: No data")
       return
    else:
        for i in range(1,len(inMat)-1):
            D = 1 + H1 * (math.exp(inMat['p_dbar'][i] / H2) - 1)
            C = math.exp(-1 * 0.04167/ H3)
            Oxnewconc[i] = ((inMat['o1_mll'][i] + (Oxnewconc[i-1] * C * D)) - (inMat['o1_mll'][i-1] * C)) / D

        inMat['o1_mll'][:] = Oxnewconc[:]
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
    fl2 = fl*2
    # One minute
    mt = 60
    # Half minute
    ms = 30
    time_delay = fl*ms

    # split dataframe into upcast/downcast
    downcast = df.iloc[:(df["CTDPRS"].argmax() + 1)]
    upcast = df.iloc[(df["CTDPRS"].argmax() + 1):]

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
    if (start_samples > time_delay):
        start_p = np.average(start_df.iloc[fl2:(start_samples - time_delay)])
    else:
        start_p = np.average(start_df.iloc[fl2:start_samples])

    end_samples = len(end_df)
    if (end_samples > time_delay):
        end_p = np.average(end_df.iloc[(time_delay):])
    else:
        try:
            end_p = np.average(end_df.iloc[(end_samples):])
        except ZeroDivisionError:
            end_p = np.NaN

    # Remove ondeck start and end pressures
    trimmed_df = df.iloc[start_df.index.max():end_df.index.min()].copy()

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
    if direction == 'down':
        monotonic_sequence = df[p_col].expanding().max()
    elif direction == 'up':
        monotonic_sequence = df[p_col].expanding().min()
    else:
        raise ValueError("direction must be one of (up, down)")

    return df[df[p_col] == monotonic_sequence]


def pressure_sequence(df, p_col='CTDPRS', direction='down'):
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
    df_filled = _fill_surface_data(df_filtered, bin_size=2)
    df_binned = binning_df(df_filled, bin_size=2)  # TODO: abstract bin_size in config

    return df_binned.reset_index(drop=True)


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
    df["bins"] = pd.cut(
        df[p_column], bins=bin_edges, right=False, include_lowest=True, labels=labels
    )
    df[p_column] = df["bins"].astype(float)
    df_out = df.groupby("bins").mean()
    return df_out


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
    # assign initial qality flags
    reftDF.loc[:, "REFTMP_FLAG_W"] = 2
    reftDF.loc[abs(reftDF["diff"]) >= 3000, "REFTMP_FLAG_W"] = 3
    # add in STNNBR, CASTNO columns
    # TODO: should these be objects or floats? be consistent!
    # string prob better for other sta/cast formats (names, letters, etc.)
    reftDF["STNNBR"] = ssscc[0:3]
    reftDF["CASTNO"] = ssscc[3:5]
    return reftDF


def process_reft(ssscc_list, reft_dir=cfg.directory["reft"]):
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
                print("refT file for cast " + ssscc + " does not exist... skipping")
                return

def _load_reft_data(reft_file, index_name="btl_fire_num"):
    """
    Loads reft_file to dataframe and reindexes to match bottle data dataframe
    """
    reft_data = pd.read_csv(reft_file, usecols=["btl_fire_num", "T90", "REFTMP_FLAG_W"])
    reft_data.set_index(index_name)
    reft_data['SSSCC_TEMP'] = Path(reft_file).stem.split("_")[0]
    reft_data['REFTMP'] = reft_data['T90']

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



def _load_btl_data(btl_file, cols=None):
    """
    Loads "bottle mean" CTD data from .pkl file. Function will return all data unless
    cols is specified (as a list of column names)
    """

    btl_data = pd.read_pickle(btl_file)
    if cols != None:
        btl_data = btl_data[cols]
    btl_data["SSSCC"] = Path(btl_file).stem.split("_")[0]

    return btl_data


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

    p_off = np.around(p_off,decimals=4)

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
    p_log = pd.read_csv(cfg.directory["logs"] + "ondeck_pressure.csv", dtype={"SSSCC":str}, na_values="Started in Water")
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
    bottom_df = pd.DataFrame(data={
        "SSSCC": df["SSSCC"].unique(),
        "max_p": df.loc[idx_p_max, "CTDPRS"],
        "lat": df.loc[idx_p_max, "GPSLAT"],
        "alt": df.loc[idx_p_max, "ALT"],
    })
    bottom_df.loc[bottom_df["alt"] > threshold, "alt"] = np.nan
    bottom_df["DEPTH"] = (
        (bottom_df["alt"] + np.abs(gsw.z_from_p(bottom_df["max_p"], bottom_df["lat"])))
        .fillna(value=-999)
        .round()
        .astype(int)
    )
    bottom_df[["SSSCC", "DEPTH"]].to_csv(
        cfg.directory["logs"] + "depth_log.csv", index=False
    )

    return True

def get_ssscc_list():
    """
    Load in list of stations/casts to process.
    """
    ssscc_list = []
    with open(cfg.directory["ssscc_file"], "r") as filename:
        ssscc_list = [line.strip() for line in filename]

    return ssscc_list

def load_hy_file(path_to_hyfile):
    df = pd.read_csv(path_to_hyfile, comment='#', skiprows=[0])
    df = df[df['EXPOCODE'] != 'END_DATA']
    return df

def load_all_ctd_files(ssscc_list, series, cols=None):
    """
    Load CTD and secondary (e.g. reference temperature, bottle salts, bottle oxygen)
    files for station/cast list and merge into a dataframe.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to load
    series : str
        Data series to load ("bottle" or "time")
    cols : list of str, optional
        Subset of columns to load, defaults to loading all

    Returns
    -------
    df_data_all : DataFrame
        Merged dataframe containing all loaded data
    
    """
    df_data_all = pd.DataFrame()

    if series == 'bottle':
        for ssscc in ssscc_list:
            print('Loading BTL data for station: ' + ssscc + '...')
            btl_file = cfg.directory["bottle"] + ssscc + '_btl_mean.pkl'
            btl_data = _load_btl_data(btl_file,cols)

            ### load REFT data
            reft_file = cfg.directory["reft"] + ssscc + '_reft.csv'
            try:
                reft_data = _load_reft_data(reft_file)
            except FileNotFoundError:
                print('Missing (or misnamed) REFT Data Station: ' + ssscc + '...filling with NaNs')
                reft_data = pd.DataFrame(index=btl_data.index, columns=["T90"], dtype=float)
                reft_data["btl_fire_num"] = btl_data["btl_fire_num"].astype(int)
                reft_data["SSSCC_TEMP"] = ssscc

            ### load REFC data
            refc_file = cfg.directory["salt"] + ssscc + '_salts.csv'
            try:
                refc_data = _load_salt_data(refc_file, index_name='SAMPNO')
            except FileNotFoundError:
                print('Missing (or misnamed) REFC Data Station: ' + ssscc + '...filling with NaNs')
                refc_data = pd.DataFrame(
                    index=btl_data.index,
                    columns=["CRavg", "BathTEMP", "BTLCOND"],
                )
                refc_data['SAMPNO_SALT'] = btl_data['btl_fire_num'].astype(int)

            ### load OXY data
            oxy_file = cfg.directory["oxy"] + ssscc
            try:
                oxy_data,params = oxy_fitting.oxy_loader(oxy_file)
            except FileNotFoundError:
                print('Missing (or misnamed) REFO Data Station: ' + ssscc + '...filling with NaNs')
                oxy_data = pd.DataFrame(
                    index=btl_data.index,
                    columns=[
                        "STNNO_OXY",
                        "CASTNO_OXY",
                        "FLASKNO",
                        "TITR_VOL",
                        "TITR_TEMP",
                        "DRAW_TEMP",
                        "TITR_TIME",
                        "END_VOLTS",
                    ],
                )
                oxy_data['BOTTLENO_OXY'] = btl_data['btl_fire_num'].astype(int)

            ### clean up dataframe
            # Horizontally concat DFs to have all data in one DF
            btl_data = pd.merge(btl_data,reft_data,on='btl_fire_num',how='outer')
            btl_data = pd.merge(btl_data,refc_data,left_on='btl_fire_num',right_on='SAMPNO_SALT',how='outer')
            btl_data = pd.merge(btl_data,oxy_data,left_on='btl_fire_num',right_on='BOTTLENO_OXY',how='outer')

            if len(btl_data) > 36:
                print("***** Len of btl data for station: ",ssscc,' is > 36, check for multiple stations/casts in reference parameter files *****')

            # Add bottom of cast information (date,time,lat,lon,etc.)
            btl_data = _add_btl_bottom_data(btl_data, ssscc)

            # Merge cast into df_data_all
            try:
                df_data_all = pd.concat([df_data_all,btl_data],sort=False)
            except AssertionError:
                raise AssertionError('Columns of ' + ssscc + ' do not match those of previous columns')
            print('* Finished BTL data station: ' + ssscc + ' *')

        # Drop duplicated columns generated by concatenation
        df_data_all = df_data_all.loc[:,~df_data_all.columns.duplicated()]
        
    elif series == 'time':
        df_data_all = []
        for ssscc in ssscc_list:
            print('Loading TIME data for station: ' + ssscc + '...')
            time_file = cfg.directory["time"] + ssscc + '_time.pkl'
            time_data = pd.read_pickle(time_file)
            time_data['SSSCC'] = str(ssscc)
            time_data['dv_dt'] = oxy_fitting.calculate_dVdT(time_data['CTDOXYVOLTS'],time_data['scan_datetime'])
            df_data_all.append(time_data)
            print('** Finished TIME data station: ' + ssscc + ' **')
        df_data_all = pd.concat(df_data_all, axis=0, sort=False)

    df_data_all['master_index'] = range(len(df_data_all))

    return df_data_all


def add_btlnbr_cols(df,btl_num_col):
    df['BTLNBR'] = df[btl_num_col].astype(int)
    # default to everything being good
    df['BTLNBR_FLAG_W'] = 2
    return df

def _add_btl_bottom_data(df, cast, lat_col='LATITUDE', lon_col='LONGITUDE', decimals=4):
    cast_details = pd.read_csv(cfg.directory["logs"] + "cast_details.csv", dtype={"SSSCC": str})
    cast_details = cast_details[cast_details["SSSCC"] == cast]
    df[lat_col] = np.round(cast_details['latitude'].iat[0], decimals)
    df[lon_col] = np.round(cast_details['longitude'].iat[0], decimals)

    ts = pd.to_datetime(cast_details['bottom_time'].iat[0], unit="s")
    date = ts.strftime('%Y%m%d')
    hour= ts.strftime('%H%M')
    df['DATE'] = date
    df['TIME'] = hour
    return df

def flag_missing_values(df,flag_columns,flag_suffix='_FLAG_W'):
    for column in flag_columns:
        flag_name = column + flag_suffix
        df.loc[df[column].isna(),flag_name] = 9
        df[column].fillna(value=int(-999), inplace=True)
        df.loc[df[column].astype(int) == -999, flag_name] = 9
    return df


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


def _flag_backfill_data(df,p_col='CTDPRS',flag_bool_col='interp_bool',flag_suffix='_FLAG_W'):
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
    print("Exporting *_ct1.csv files")
    # clean up columns
    p_column_names = cfg.ctd_time_output["col_names"]
    p_column_units = cfg.ctd_time_output["col_units"]
    # initial flagging (some of this should be moved)
    # TODO: actual CTDOXY flagging (oxy_fitting.calibrate_oxy)
    df["CTDOXY_FLAG_W"] = 2
    # TODO: flag bad based on cond/temp and handcoded salt
    df["CTDSAL_FLAG_W"] = 2
    df["CTDTMP_FLAG_W"] = 2
    # TODO: lump all uncalibrated together; smart flagging like ["CTD*_FLAG_W"] = 1
    # TODO: may not always have these channels so don't hardcode them in!
    df["CTDFLUOR_FLAG_W"] = 1
    df["CTDXMISS_FLAG_W"] = 1
    df["CTDBACKSCATTER_FLAG_W"] = 1
    # renames
    df = df.rename(columns={"CTDTMP1": "CTDTMP", "FLUOR": "CTDFLUOR"})
    # check that all columns are there
    # TODO: make this better... 
    # #should it fail and return list of bad cols or just use fill values?
    try:
        df[p_column_names];  # this is lazy, do better
    except KeyError:
        print("Column names not configured properly... attempting to correct")
        for col in p_column_names:
            try:
                df[col];
            except KeyError:
                if col.endswith("FLAG_W"):
                    print(col + " missing, flagging with 9s")
                    df[col] = 9
                else:
                    print(col + " missing, filling with -999s")
                    df[col] = -999

    df["SSSCC"] = df["SSSCC"].astype(str).copy()
    cast_details = pd.read_csv(cfg.directory["logs"] + "cast_details.csv", dtype={"SSSCC": str})
    depth_df = pd.read_csv(cfg.directory["logs"] + 'depth_log.csv', dtype={"SSSCC": str}, na_values=-999).dropna()
    try:
        manual_depth_df = pd.read_csv(cfg.directory["logs"] + 'manual_depth_log.csv', dtype={"SSSCC": str})
    except FileNotFoundError:
        # TODO: add logging; look into inheriting/extending a class to add features
        print("manual_depth_log.csv not found... duplicating depth_log.csv")
        manual_depth_df = depth_df.copy()  # write manual_depth_log as copy of depth_log
        manual_depth_df.to_csv(cfg.directory["logs"] + 'manual_depth_log.csv', index=False)
    full_depth_df = pd.concat([depth_df,manual_depth_df])
    full_depth_df.drop_duplicates(subset='SSSCC', keep='first',inplace=True)

    for ssscc in ssscc_list:

        time_data = df[df['SSSCC'] == ssscc].copy()
        time_data = pressure_sequence(time_data)
        time_data = time_data[p_column_names]
        time_data = time_data.round(4)
        time_data = time_data.where(~time_data.isnull(), -999)  #replace NaNs with -999

        depth = full_depth_df.loc[full_depth_df['SSSCC'] == ssscc,'DEPTH'].iloc[0]
        # get cast_details for current SSSCC
        cast_dict = cast_details[cast_details["SSSCC"] == ssscc].to_dict("records")[0]
        b_datetime = datetime.fromtimestamp(cast_dict["bottom_time"], tz=timezone.utc).strftime('%Y%m%d %H%M').split(" ")
        # TODO: yo-yo casts are an edge case where this may be different
        btm_lat = cast_dict["latitude"]
        btm_lon = cast_dict["longitude"]
        btm_alt = cast_dict["altimeter_bottom"]

        now = datetime.now(timezone.utc)
        file_datetime = now.strftime("%Y%m%d") #%H:%M")
        file_datetime = file_datetime + 'ODFSIO'
        # TODO: only "cast" needs to be int; "station" is explicitly allowed to incl.
        # letters/etc. Moving from SSSCC to station & cast fields will be beneficial
        with open(f"{cfg.directory['pressure']}{ssscc}_ct1.csv", "w+") as f:
            # put in logic to check columns?
            # number_headers should be calculated, not defined
            ctd_header = (  # this is ugly but prevents tabs before label
                f"CTD,{file_datetime}\n"
                f"NUMBER_HEADERS = 11\n"
                f"EXPOCODE = {cfg.cruise['expocode']}\n"
                f"SECT_ID = {cfg.cruise['sectionid']}\n"
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
            np.asarray(p_column_names).tofile(f, sep=',', format='%s')
            f.write('\n')
            np.asarray(p_column_units).tofile(f, sep=',', format='%s')
            f.write('\n')
            time_data.to_csv(f, header=False, index=False)
            f.write("END_DATA")


def export_btl_data(df, out_dir=cfg.directory["pressure"], org='ODF'):

    btl_data = df.copy()
    now = datetime.now()
    file_datetime = now.strftime("%Y%m%d")
    btl_columns = [
        "EXPOCODE",
        "SECT_ID",
        "STNNBR",
        "CASTNO",
        "SAMPNO",
        "BTLNBR",
        "BTLNBR_FLAG_W",
        "DATE",
        "TIME",
        "LATITUDE",
        "LONGITUDE",
        "DEPTH",
        "CTDPRS",
        "CTDTMP",
        "CTDSAL",
        "CTDSAL_FLAG_W",
        "SALNTY",
        "SALNTY_FLAG_W",
        "CTDOXY",
        "CTDOXY_FLAG_W",
        "OXYGEN",
        "OXYGEN_FLAG_W",
        "REFTMP",
        "REFTMP_FLAG_W",
    ]
    btl_units = [
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "METERS",
        "DBAR",
        "ITS-90",
        "PSS-78",
        "",
        "PSS-78",
        "",
        "UMOL/KG",
        "",
        "UMOL/KG",
        "",
        "ITS-90",
        "",
    ]

    # rename
    btl_data = btl_data.rename(columns={"CTDTMP2": "CTDTMP"})
    btl_data["EXPOCODE"] = cfg.cruise["expocode"]
    btl_data["SECT_ID"] = cfg.cruise["sectionid"]
    btl_data["STNNBR"] = [int(x[0:3]) for x in btl_data["SSSCC"]]
    btl_data["CASTNO"] = [int(x[3:]) for x in btl_data["SSSCC"]]
    btl_data["SAMPNO"] = btl_data["btl_fire_num"].astype(int)
    btl_data = add_btlnbr_cols(btl_data, btl_num_col="btl_fire_num")

    # sort by decreasing sample number (increasing pressure)
    btl_data = btl_data.sort_values(by=["STNNBR", "SAMPNO"], ascending=[True, False])

    # round data
    for col in ["CTDTMP", "CTDSAL", "SALNTY", "REFTMP"]:
        btl_data[col] = btl_data[col].round(4)
    for col in ["CTDPRS", "CTDOXY", "OXYGEN"]:
        btl_data[col] = btl_data[col].round(1)

    # add depth
    depth_df = pd.read_csv(cfg.directory["logs"] + 'depth_log.csv', dtype={"SSSCC": str}, na_values=-999).dropna()
    manual_depth_df = pd.read_csv(cfg.directory["logs"] + 'manual_depth_log.csv', dtype={"SSSCC": str})
    full_depth_df = pd.concat([depth_df,manual_depth_df])
    full_depth_df.drop_duplicates(subset='SSSCC', keep='first',inplace=True)
    btl_data["DEPTH"] = -999
    for index, row in full_depth_df.iterrows():
        btl_data.loc[btl_data["SSSCC"] == row["SSSCC"], "DEPTH"] = int(row["DEPTH"])

    # deal with nans
    btl_data.loc[btl_data["REFTMP_FLAG_W"].isnull(), "REFTMP_FLAG_W"] = 9
    btl_data["REFTMP_FLAG_W"] = btl_data["REFTMP_FLAG_W"].astype(int)
    btl_data = btl_data.where(~btl_data.isnull(), -999)

    # flag CTDOXY with more than 1% difference
    btl_data["CTDOXY_FLAG_W"] = 2
    diff = btl_data["CTDOXY"] - btl_data["OXYGEN"]
    btl_data.loc[diff.abs()*100 > btl_data["CTDOXY"], "CTDOXY_FLAG_W"] = 3

    # flag CTDSAL using stepped filter
    thresh = [0.002, 0.005, 0.010, 0.020]
    btl_data["CTDSAL_FLAG_W"] = 2
    diff = btl_data["CTDSAL"] - btl_data["SALNTY"]
    btl_data.loc[(btl_data["CTDPRS"] > 2000) & (diff.abs() > thresh[0]), "CTDSAL_FLAG_W"] = 3
    btl_data.loc[(btl_data["CTDPRS"] <= 2000) & (btl_data["CTDPRS"] > 1000) & (diff.abs() > thresh[1]), "CTDSAL_FLAG_W"] = 3
    btl_data.loc[(btl_data["CTDPRS"] <= 1000) & (btl_data["CTDPRS"] > 500) & (diff.abs() > thresh[2]), "CTDSAL_FLAG_W"] = 3
    btl_data.loc[(btl_data["CTDPRS"] <= 500) & (diff.abs() > thresh[3]), "CTDSAL_FLAG_W"] = 3

    # check columns
    try:
        btl_data[btl_columns];  # this is lazy, do better
    except KeyError:
        print("Column names not configured properly... attempting to correct")
        for col in btl_columns:
            try:
                btl_data[col];
            except KeyError:
                if col.endswith("FLAG_W"):
                    print(col + " missing, flagging with 9s")
                    btl_data[col] = 9
                else:
                    print(col + " missing, filling with -999s")
                    btl_data[col] = -999

    btl_data = btl_data[btl_columns]
    time_stamp = file_datetime+org
    with open(out_dir + cfg.cruise["expocode"] + "_hy1.csv", mode="w+") as f:
        f.write("BOTTLE, %s\n" % (time_stamp))
        f.write(",".join(btl_columns) + "\n")
        f.write(",".join(btl_units) + "\n")
        btl_data.to_csv(f, header=False, index=False)
        f.write("\n" + "END_DATA")

    return

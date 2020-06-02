#!/usr/bin/env python
from scipy import interpolate
import numpy as np
from numpy.lib.recfunctions import append_fields
import scipy.signal as signal
import scipy.stats as st
import time, os
import pandas as pd
import xarray as xr
import math
#import report_ctd
import ctdcal.report_ctd as report_ctd
import warnings
import ctdcal.fit_ctd as fit_ctd
from datetime import datetime, timezone
from decimal import Decimal
import config as cfg

import sys
sys.path.append('ctdcal/')
import settings
import oxy_fitting
import gsw
import csv
from collections import OrderedDict
from pathlib import Path

warnings.filterwarnings("ignore", 'Mean of empty slice.')


def cast_details(stacast, log_file, p_col, time_col, b_lat_col, b_lon_col, alt_col, ds=None):
    """
    We determine the cast details using pandas magic.
    First find alternating periods of pumps on and pumps off, then select the
    pumps on period with the highest pressure. Get values from the row with the
    highest pressure, and return all values to be sent to log.

    Parameters
    ----------
    stacast : integer
        The station and cast, as SSSCC format
    log_file : file handle or string
        File destination for cast details
    p_col : str
        Name of the pressure column
    time_col : str
        Name of the time column
    b_lat_col : str
        Name of the latitude column
    b_lon_col : str
        Name of the longitude column
    alt_col : str
        Name of the altimeter column
    ds : xarray Dataset
        The filtered Dataset

    Returns
    -------
    ds_cast : xarray Dataset
        The Dataset that came in with the soak period trimmed off

    Notes
    -----
    The following variables are output to log_file:
    time_start : float
        Unix epoch seconds?, start of cast time, to be reported to log file
    time_end : float
        Unix epoch seconds?, end of cast time, to be reported to log file
    time_bottom : float
        Unix epoch seconds?, bottom of cast time, to be reported to log file
    p_start : float
        Pressure at which cast started, to be reported to log file
    p_max : float
        Bottom of the cast pressure, to be reported to log file
    b_lat : float
        Latitude at bottom of cast
    b_lon : float
        Longitude at bottom of cast
    b_alti : float
        Altimeter reading at bottom of cast - volts only!
    """

    ds_cast = _trim_soak_period(ds)
    # TODO: call parameters from config file instead
    p_start = float(np.around(ds_cast['CTDPRS'].head(1),4))
    p_max_ind = ds_cast['CTDPRS'].argmax()
    p_max = float(np.around(ds_cast['CTDPRS'].max(),4))
    time_start = float(ds_cast['scan_datetime'].head(1))
    time_end = float(ds_cast['scan_datetime'].tail(1))
    time_bottom = float(ds_cast['scan_datetime'][p_max_ind])
    b_lat = float(np.around(ds_cast['GPSLAT'][p_max_ind],4))
    b_lon = float(np.around(ds_cast['GPSLON'][p_max_ind],4))
    b_alt = float(np.around(ds_cast['ALT'][p_max_ind],4))

    # last two lines must be in to return the same as old - change to slices of df later
    report_ctd.report_cast_details(stacast, log_file, time_start, time_end,
                                   time_bottom, p_start, p_max, b_alt,
                                   b_lat, b_lon)
    # remove upcast
    p_max_ind = int(ds_cast['CTDPRS'].argmax())
    ds_cast = ds_cast.sel(index=np.arange(0, p_max_ind + 1))  # exclusive indexing

    return ds_cast


def _trim_soak_period(ds=None):
    """
    1) Find pump on/off patterns
    2) Select pump_on=True group with largest pressure recording
    3) Find soak period before start of downcast
    4) Trim cast, return everything after top of cast (i.e. minimum pressure)
    """
    ds_list = [
        g
        for i, g in ds.groupby((ds["pump_on"] != ds["pump_on"].shift(index=1)).cumsum())
    ]
    ds_pump_on_list = [ds for ds in ds_list if ds["pump_on"].all()]
    ds_cast = ds_pump_on_list[np.argmax([ds["CTDPRS"].max() for ds in ds_pump_on_list])]
    ds_cast["index"] = np.arange(0, len(ds_cast["index"]))
    ds_cast = _find_last_soak_period(ds_cast)  # deals w/ edge cases, leave as is for now
    start_ind = ds_cast["CTDPRS"][: len(ds["CTDPRS"]) // 4].argmin().values
    ds_trimmed = ds_cast.sel(index=np.arange(start_ind, len(ds_cast["index"])))
    ds_trimmed["index"] = np.arange(0, len(ds_trimmed["index"]))

    return ds_trimmed


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
    if df_cast['CTDPRS'][0] > surface_pressure:
        return df_cast
    breakpoint()
    # TODO: this still needs to be converted to xarray
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

#End move four functions
# def cast_details_old(stacast, log_file, p_col, time_col, b_lat_col, b_lon_col, alt_col, inMat=None):
#     """cast_details function
#
#     Function takes full NUMPY ndarray with predefined dtype array
#     and adjusts ndarray to remove all extraneous surface data.
#     Function returns cast start time, end time, bottom time and
#     cleaned up matrix.
#
#     Args:
#         param1 (str): stacast, station cast input
#         param2 (str): log_file, log file to write cast data.
#         param3 (str): p_col, pressure data column name
#         param4 (str): time_col, time data column name
#         param5 (ndarray): inMat, numpy ndarray with dtype array
#
#     Returns:
#         Narray: The return value is ndarray with adjusted time of parameter
#           specified.
#
#     """
#
#
#     if inMat is None:
#        print("In cast_details: No data")
#        return
#     else:
#         # Top of cast time, bottom of cast time, end of cast time,
#         start_cast_time = 0.0
#         bottom_cast_time = 0.0
#         end_cast_time = 0.0
#         # Test cycle time constant
#         fl = 24
#         # starting P
#         start_pressure = 2.0
#         # Max P
#         max_pressure = 10000.0
#         lm = len(inMat)-1
#         rev = np.arange(int(lm/4),0,-1)
#
#         # Find starting top of cast
#         # Smallest P from reverse array search
#         for i in rev:
#             if start_pressure < inMat[p_col][i]:
#                tmp = i
#             elif start_pressure > inMat[p_col][i]:
#                start_pressure = inMat[p_col][i]
#                tmp = abs(i - 24) #patched to not break through the c(sea)-floor, can be made cleaner
#                break
#         start_cast_time = inMat[time_col][tmp]
#
#         # Remove everything before cast start
#         inMat = inMat[tmp:]
#
#         # Max P and bottom time
#         max_pressure = max(inMat[p_col])
#         tmp = np.argmax((inMat[p_col]))
#         bottom_cast_time = inMat[time_col][tmp]
#         b_lat = inMat[b_lat_col][tmp]
#         b_lon = inMat[b_lon_col][tmp]
#         b_alti = inMat[alt_col][tmp]
#
#         tmp = len(inMat)
#         # Find ending top of cast time
#         for i in range(int(tmp/2),tmp):
#             if start_pressure > inMat[p_col][i]:
#                 end_cast_time = inMat[time_col][i]
#                 if i < tmp: tmp = i + 24
#                 break
#
#         # Remove everything after cast end
#         inMat = inMat[:tmp]
#
#     report_ctd.report_cast_details(stacast, log_file, start_cast_time, end_cast_time, bottom_cast_time, start_pressure, max_pressure, b_alti, b_lat, b_lon)
#
#     return start_cast_time, end_cast_time, bottom_cast_time, start_pressure, max_pressure, b_lat, b_lon, b_alti, inMat

def ctd_align(raw_ds=None, col=None, time=0.0):
    """
    Adjust time of sensor response and water flow relative to the time frame of
    temperature sensor. Last value is duplicated and appended to time series.

    Parameters
    ----------
    raw_ds : xarray Dataset
        Dataset containing raw CTD data
    col : str
        Name of variable being offset
    time : float
        Time offset (in seconds) to apply to variable

    Returns
    -------
    raw_ds : xarray Dataset
        Dataset with adjusted variable

    """
    fl = 24  # sampling frequency (TODO: put in sensor attrs?)

    if (raw_ds is not None) & (col is not None) & (time > 0.0):
        advnc = int(fl * time)  # number of samples to advance
        last = raw_ds[col][-1]
        tmp_col = raw_ds[col].pad(index=(0, advnc), constant_values=last)[advnc:]
        tmp_col["index"] = np.arange(0, len(tmp_col))  # trim and reindex to match len
        raw_ds[col].values = tmp_col

    return raw_ds

def ctd_quality_codes(column=None, p_range=None, qual_code=None, oxy_fit=False, p_qual_col=None, qual_one=None, inMat=None):
    """ctd_quality_codes function

    Function takes full NUMPY ndarray with predefined dtype array

    Args:
        param1 (ndarray):
        param2 (float):

    Returns:
        Narray: The return value is ndarray with adjusted time of parameter
          specified.

    """
    #If p_range set apply qual codes to part of array and return
    if p_range is not None:
        print("Some algoirythm for formatting qual codes per pressure range")
        return
    else:
        q_df = pd.DataFrame(index=np.arange(len(inMat)), columns=p_qual_col)
        for pq in p_qual_col:
            if pq in list(qual_one):
                q_df[pq] = q_df[pq].fillna(1)
            elif oxy_fit and pq is column:
                q_df[pq] = q_df[pq].fillna(2)
            else:
                q_df[pq] = q_df[pq].fillna(2)

    return q_df.values  # ndarray format

def formatTimeEpoc(time_zone='UTC', time_pattern='%Y-%m-%d %H:%M:%S', input_time = None):
    """formatTimeEpoc function

    Function takes pattern of time input, relative time zone, and
    date time data array and returns array of epoc time.

    title and the second row are the units for each column.
    Args:
        param1 (str): relative time zone for data.
        param2 (str): pattern of incoming data.
        param3 (ndarray): input_time, numpy 1d ndarray time array

    Returns:
        1D ndarray: The return array of epoch time
    """
    if input_time is None:
        print("In formatTimeEpoc: No data entered.")
        return
    else:
        os.environ['TZ'] = 'UTC'
        epoch_time = input_time
        for i in range(0,len(input_time)):
            epoch_time[i] = int(time.mktime(time.strptime(str(input_time[i], "utf-8"), time_pattern)))

    return epoch_time

def dataToDataFrame(inFile):
    """dataToDataFrame function

    Function takes full file path to csv type data file and returns a
    PANDAS dataframe for data treatment with a two row header.

    Data file should have a two row header. The first row being the column
    title and the second row are the units for each column.
    Args:
        param1 (str): Full path to data file.

    Returns:
        DataFrame: The return value is a full dataframe with header.

    .. REF PAGE:
       http://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html#pandas.read_csv
    """
    #df = pd.read_csv(inFile, header=[0,2])
    df = pd.read_csv(inFile)
    return df

def dataToNDarray(inFile, dtype=None, names=None, separator=',', skip=None):
    """dataToNDarray function

    Function takes full file path to csv type data file and returns NUMPY
    ndarray type ndarray for data manipulation with a two row header.

    Data file should have a two row header. The first row being the column
    title and the second row are the units for each column.
    Args:
        param1 (str): inFile, full path to csv file
        param2 (arr): dtype list
        param3 (str): separator, default comma ','

    Returns:
        Narray: The return value is a full data ndarray with two row header.

    Reference Page:
        https://scipy.github.io/old-wiki/pages/Cookbook/InputOutput.html
    """
    try:
        return pd.read_pickle(inFile).to_records()
    except:
        if skip is None:
            arr = np.genfromtxt(inFile, delimiter=separator, dtype=dtype, names=names)
        else:
            arr = np.genfromtxt(inFile, delimiter=separator, dtype=dtype, names=names, skip_header=skip)

    return arr

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


def data_interpolater(inArr):
    """data_interpolater to handle indices and logical indices of NaNs.

    Input:
        - inArr, 1d numpy array with return True np.isnans()
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
        - interpolated array
    Example:
        >>> # linear interpolation of NaNs
        >>> outArray = data_interpolater(inArr)
    """
    nans, tmp= np.isnan(inArr), lambda z: z.nonzero()[0]
    inArr[nans] = np.interp(tmp(nans), tmp(~nans), inArr[~nans])
    return inArr

def o2pl2pkg(p_col, t_col, sal_col, dopl_col, dopkg_col, lat_col, lon_col, inMat):
    """o2pl2pkg convert ml/l dissolved oxygen to umol/kg

    Input:
        - t_col, temperature column header deg c.
        - sal_col, salinity column header psu.
        - dopl_col, dissolved column header ml/l.
        - dopkg_col, dissolved column header umol/kg
        - lat_col, latitude for entire cast deg.
        - lon_col, longitude for entire cast deg.
        - inMat, dtype ndarray processed ctd time data.
    Output:
        - Converted Oxygen column umol/kg
    Example:
        >>> # linear interpolation of NaNs
        >>> outArray = o2pl2kg(inArr)
    """
    pkg = np.ndarray(shape=len(inMat), dtype=[(dopkg_col, np.float)])
    # Absolute sailinity from Practical salinity.
    SA = gsw.SA_from_SP(inMat[sal_col], inMat[p_col], inMat[lat_col], inMat[lon_col])

    # Conservative temperature from insitu temperature.
    CT = gsw.CT_from_t(SA, inMat[t_col], inMat[p_col])
    s0 = gsw.sigma0(SA, CT) # Potential density from Absolute Salinity g/Kg Conservative temperature deg C.

    # Convert DO ml/l to umol/kg
    for i in range(0,len(inMat[dopl_col])):
        pkg[i] = inMat[dopl_col][i] * 44660 / (s0[i] + 1000)
    return pkg

def oxy_to_umolkg(df_sal, df_pressure, df_lat, df_lon, df_temp, df_oxy):
    '''Rewritten from Courtney's method to use array-likes (aka use dataframes and ndarrays).
    '''
    # Absolute salinity from Practical salinity.
    SA = gsw.SA_from_SP(df_sal, df_pressure, df_lat, df_lon)

    # Conservative temperature from insitu temperature.
    CT = gsw.CT_from_t(SA, df_temp, df_pressure)
    s0 = gsw.sigma0(SA, CT) # Potential density from Absolute Salinity g/Kg Conservative temperature deg C.
    series = df_oxy * 44660 / (s0 + 1000)
    return series

def raw_ctd_filter(raw_ds=None, filter_type='triangle', win_size=24, parameters=None):
    """raw_ctd_filter function

    Function takes NUMPY array
    of raw ctd data and returns filtered data. This function also needs
    one of three filter types (boxcar, gaussian, triangle) as well as
    window size.

    Args:
        param1 (ndarray): Numpy ndarray with predefined header with at
        param2 (str): One of three tested filter types
          boxcar, gaussian_std, triangle.
          default is triangle
        param3 (int): A window size for the filter. Default is 24, which
          is the number of frames per second from a SBE9+/11 CTD/Dech unit.
        param4 (ndarray): parameters the dtype names used in filtering the
          analytical inputs.

    Returns:
        Narray: The return value is a matrix of filtered ctd data with
          the above listed header values.

    """

    if raw_ds is None:
        print("In raw_ctd_filter: No data array.")
        return
    elif parameters is None:
        print("In raw_ctd_filter: Empty parameter list.")
        return
    else:
        filt_ds = raw_ds.copy()
        for p in parameters:
            if p not in raw_ds.variables:
                print(f"'{p}' not found in dataset, filter not applied")
                continue

            filt_ds[p].attrs["filter_type"] = filter_type
            filt_ds[p].attrs["filter_window_size"] = win_size

            if filter_type == "boxcar":
                win = signal.boxcar(win_size)
                filt_ds[p].values = signal.convolve(raw_ds[p], win, mode="same")/len(win)
            elif filter_type == "gaussian":
                sigma = np.std(arr)
                win = signal.general_gaussian(win_size, 1.0, sigma)
                filt_ds[p].values = signal.convolve(raw_ds[p], win, mode="same")/len(win)
            elif filter_type == "triangle":
                win = signal.triang(win_size)
                filt_ds[p].values = 2*signal.convolve(raw_ds[p], win, mode="same")/len(win)

    return filt_ds


def ondeck_pressure(stacast, p_col, c1_col, c2_col, time_col, raw_ds=None, cond_startup=20.0, log_file=None):
    """ondeck_pressure function
    Function takes full NUMPY ndarray with predefined dtype array
    of filtered ctd raw data the stores, analizes and removes ondeck
    values from data.
    Args:
        param1 (str): stacast, station cast info
        param1 (str): p_col, pressure data column name
        param2 (str): c1_col, cond1 data column name
        param3 (str): c2_col, cond2 data column name
        param4 (str): time_col, time data column name
        param5 (ndarray): numpy ndarray with dtype array
        param6 (float): conductivity_startup, threshold value
        param7 (str): log_file, log file name
    Returns:
        Narray: The return ndarray with ondeck data removed.
        Also output start/end ondeck pressure.
    """
    start_pressure = []

    # Frequency
    fl = 24
    fl2 = fl*2
    # One minute
    mt = 60
    # Half minute
    ms = 30
    time_delay = fl*ms

    if raw_ds is None:
        print("Ondeck_pressure function: No data.")
        return
    else:
        # Search first quarter of matrix, using conductivity
        # threshold min to capture startup pressure
        c1 = raw_ds[c1_col].to_series()
        c2 = raw_ds[c2_col].to_series()
        p = raw_ds[p_col].to_series()
        n = len(p)//4
        start_pressure = p[
            np.flatnonzero((c1[:n] < cond_startup) & (c2[:n] < cond_startup))
        ]

        # Evaluate starting pressures
        if start_pressure is None:
            start_p = "Started in Water"
        else:
            n_start = len(start_pressure)
            if (n > time_delay):
                start_p = np.average(start_pressure[fl2:(n_start-time_delay)]).round(4)
            else:
                start_p = np.average(start_pressure[fl2:n]).round(4)

        # Remove on-deck startup
        raw_ds = raw_ds.sel(index=np.arange(n_start, len(p)))  # reset_index differs
        raw_ds["index"] = np.arange(0, len(raw_ds.index))  # from how it works in pandas

        # Searches last half of NDarray for conductivity threshold
        c1 = raw_ds[c1_col].to_series()
        c2 = raw_ds[c2_col].to_series()
        p = raw_ds[p_col].to_series()
        n = len(p)//2
        end_pressure = p[  # flatnonzero effectively resets index, so re-add n
            np.flatnonzero((c1[n:] < cond_startup) & (c2[n:] < cond_startup)) + n
        ]

        # Evaluate ending pressures
        n_end = len(end_pressure)
        if (n_end > time_delay):
            end_p = np.average(end_pressure[time_delay:]).round(4)
        else:
            end_p = np.average(end_pressure[n:]).round(4)

        # Remove on-deck ending
        raw_ds = raw_ds.sel(index=np.arange(0, len(p)-n_end))

        # Store ending on-deck pressure
        if log_file is not None:
            report_ctd.report_pressure_details(stacast, log_file, start_p, end_p)

    return raw_ds

def ondeck_pressure_2(df, stacast, p_col, c1_col, c2_col, conductivity_startup=20.0, log_file=None):
    """ondeck_pressure function
    Function takes pandas Dataframe of filtered ctd raw data the stores, analizes and removes ondeck
    values from data.
    Args:
        param1 (str): stacast, station cast info
        param1 (str): p_col, pressure data column name
        param2 (str): c1_col, cond1 data column name
        param3 (str): c2_col, cond2 data column name
        param4 (str): time_col, time data column name
        param5 (ndarray): numpy ndarray with dtype array
        param6 (float): conductivity_startup, threshold value
        param7 (str): log_file, log file name
    Returns:
        Narray: The return ndarray with ondeck data removed.
        Also output start/end ondeck pressure.
    """
    # Frequency
    fl = 24
    fl2 = fl*2
    # One minute
    mt = 60
    # Half minute
    ms = 30
    time_delay = fl*ms

    #split dataframe into upcast/downcast

    down, up = df.split(p_col=p_col)


    # Searches each half of df, uses conductivity
    # threshold min to capture startup pressure

    start_df = down.loc[(down[c1_col] < 20) & (down[c2_col] < 20)]
    end_df = up.loc[(up[c1_col] < 20) & (up[c2_col] < 20)]

    # Evaluate starting and ending pressures
    sp = len(start_df)

    if (sp > time_delay):
        start_p = np.average(start_df.iloc[fl2:sp-(time_delay)][p_col])
    else:
        start_p = np.average(start_df[fl2:sp])

    ep = len(end_df)

    ep = len(end_df)

    if (ep > time_delay):
        end_p = np.average(end_df.iloc[(time_delay):])
    else:
        try:
            end_p = np.average(end_df.iloc[(ep):])
        except ZeroDivisionError:
            end_p = np.NaN

    # Remove ondeck start and end pressures

    df.iloc[start_df.index.max():end_df.index.min()][p_col].max().copy()
        # Store ending on-deck pressure
    if log_file != None:
        report_ctd.report_pressure_details(stacast, log_file, start_p, end_p)

    return outMat



def _roll_filter(df, pressure_column="CTDPRS", direction="down"):
    #fix/remove try/except once serialization is fixed
    try:
        if direction == 'down':
            monotonic_sequence = df[pressure_column].expanding().max()
        elif direction == 'up':
            monotonic_sequence = df[pressure_column].expanding().min()
        else:
            raise ValueError("direction must be one of (up, down)")
    except KeyError:
        pressure_column = 'CTDPRS'
        if direction == 'down':
            monotonic_sequence = df[pressure_column].expanding().max()
        elif direction == 'up':
            monotonic_sequence = df[pressure_column].expanding().min()
        else:
            raise ValueError("direction must be one of (up, down)")

    return df[df[pressure_column] == monotonic_sequence]


def roll_filter(inMat, p_col, up='down', frames_per_sec=24, search_time=15, **kwargs):
    """roll_filter function
    Function takes full NUMPY ndarray with predefined dtype array
    and subsample arguments to return a roll filtered ndarray.
    Args:
        param1 (str): stacast, station cast info
        param2 (ndarray): inMat, numpy ndarray with dtype array
        param3 (str): up, direction to filter cast (up vs down)
        param4 (int): frames_per_sec, subsample selection rate
        param5 (int): seach_time, search time past pressure inversion
    Returns:
        Narray: The return value ndarray of data with ship roll removed
    """
    #When the "pressure sequence" code is fixed, uncomment and use this instead
    start = kwargs.get("start", 0)
    end = kwargs.get("end", -1)
    full_matrix = kwargs.get("full_matrix", inMat)
    tmp_df = pd.DataFrame.from_records(full_matrix[start:end])
    tmp_df = _roll_filter(tmp_df)
    #return tmp_df.to_records(index=False)
    return tmp_df

    remove = []
    frequency = 24 # Hz of package

    if (frames_per_sec > 0) & (frames_per_sec <= 24):
        sample = int(frequency/frames_per_sec) # establish subsample rate to time ratio
    else: sample = frequency

    # Adjusted search time with subsample rate
    search_time = int(sample*frequency*int(search_time))

    if inMat is None:
        print("Roll filter function: No input data.")
        return
    else:
        P = inMat[p_col]
        dP = np.diff(P,1)

        if up == 'down':
            index_to_remove = np.where(dP < 0)[0] # Differential filter
            subMat = np.delete(inMat, index_to_remove, axis=0)

            P = subMat[p_col]
            tmp = np.array([])
            for i in range(0,len(P)-1):
               if P[i] > P[i+1]:
                   deltaP = P[i+1] + abs(P[i] - P[i+1])
                   # Remove aliasing
                   k = np.where(P == min(P[i+1:i+search_time], key=lambda x:abs(x-deltaP)))[0]
                   tmp = np.arange(i+1,k[0]+1,1)
               remove = np.append(remove,tmp)
               deltaP = 0
        elif up == 'up':
            index_to_remove = np.where(dP > 0)[0] # Differential filter
            subMat = np.delete(inMat, index_to_remove, axis=0)

            P = subMat[p_col]
            tmp = np.array([])
            for i in range(0,len(P)-1):
               if P[i] < P[i+1]:
                   deltaP = P[i+1] - abs(P[i] - P[i+1])
                   # Remove aliasing
                   k = np.where(P == min(P[i+1:i+search_time], key=lambda x:abs(x-deltaP)))[0]
                   tmp = np.arange(i+1,k[0]+1,1)
               remove = np.append(remove,tmp)
               deltaP = 0

        subMat = np.delete(subMat,remove,axis=0)

    return subMat

def pressure_sequence(df, p_col='CTDPRS', intP=2.0, startT=-1.0, startP=0.0, up='down', sample_rate=12, search_time=15):
    """pressure_sequence function

    Function takes a dataframe and several arguments to return a pressure
    sequenced data ndarray.

    Pressure sequencing includes rollfilter.

    Necessary inputs are input Matrix (inMat) and pressure interval (intP).
    The other inputs have default settings. The program will figure out
    specifics for those settings if left blank.
    Start time (startT), start pressure (startP) and up are mutually exclusive.
    If sensors are not not fully functional when ctd starts down cast
    analyst can select a later start time or start pressure but not both.
    There is no interpolation to the surface for other sensor values.
    'up' indicates direction for pressure sequence. If up is set startT and startP
    are void.

    Args:
        param1 (Dataframe: Dataframe containing measurement data
        param2 (str): p_col, pressure column name
        param3 (float): starting pressure interval
        param5 (float): start time (startT) for pressure sequence
        param6 (float): start pressure (startP) for pressure sequence
        param7 (str): pressure sequence direction (down/up)
        param8 (int): sample_rate, sub sample rate for roll_filter. Cleans & speeds processing.
        param9 (int): search_time, truncate search index for the aliasing part of ship roll.
        param10 (ndarray): inMat, input data ndarray

    Returns:
        Narray: The return value is a matrix of pressure sequenced data

    todo: deep data bin interpolation to manage empty slices
    """
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

    #lenP, prvPrs not used
    # Passed Time-Series, Create Pressure Series

    start = 0

    # Roll Filter
    roll_filter_matrix = roll_filter(df, p_col, up, sample_rate, search_time, start=start)

    df_roll_surface = fill_surface_data(roll_filter_matrix, bin_size=2)
    #bin_size should be moved into config
    binned_df = binning_df(df_roll_surface, bin_size=2)
    binned_df = binned_df.reset_index(drop=True)
    return binned_df

def binning_df(df, **kwargs):
    '''Bins records according to bin_size, then finds the mean of each bin and returns a df.
    '''
    bin_size = kwargs.get("bin_size", 2)
    try:
        labels_in = [x for x in range(0,int(np.ceil(df['CTDPRS_DBAR'].max())),2)]
        df['bins'] = pd.cut(df['CTDPRS_DBAR'], range(0,int(np.ceil(df['CTDPRS_DBAR'].max()))+bin_size,bin_size), right=False, include_lowest=True, labels=labels_in)
        df['CTDPRS_DBAR'] = df['bins'].astype('float64')
        df_out = df.groupby('bins').mean()
        return df_out
    except KeyError:
        labels_in = [x for x in range(0,int(np.ceil(df['CTDPRS'].max())),2)]
        df['bins'] = pd.cut(df['CTDPRS'], range(0,int(np.ceil(df['CTDPRS'].max()))+bin_size,bin_size), right=False, include_lowest=True, labels=labels_in)
        df['CTDPRS'] = df['bins'].astype('float64')
        df_out = df.groupby('bins').mean()
        return df_out

def fill_surface_data(df, **kwargs):
    '''Copy first scan from top of cast, and propgate up to surface
    '''
    surface_values = []
    bin_size = kwargs.get("bin_size", 2)
    try:
        for x in range(1, int(np.floor(df.iloc[0]['CTDPRS_DBAR'])), bin_size):
            surface_values.append(x)
        df_surface = pd.DataFrame({'CTDPRS_DBAR': surface_values})
        df_surface['interp_bol'] = 1
        df_merged = pd.merge(df_surface, df, on='CTDPRS_DBAR', how='outer')
    except KeyError:
        for x in range(1, int(np.floor(df.iloc[0]['CTDPRS'])), bin_size):
            surface_values.append(x)
        df_surface = pd.DataFrame({'CTDPRS': surface_values})
        # Added by KJ to keep track of backfilled values
        df_surface['interp_bol'] = 1
        if len(df_surface['interp_bol']) == 1:
            df_surface['interp_bol'] = 0
        df_merged = pd.merge(df_surface.astype('float64'), df, on='CTDPRS', how='outer')
    if 'interp_bol' not in df_merged.columns:
        df_merged['interp_bol'] = np.NaN
    df_merged['interp_bol'].fillna(0,inplace=True)

    return df_merged.fillna(method='bfill')


def _reft_loader(ssscc, reft_dir):
    # semi-flexible search for reft file (in the form of *ssscc.cap)
    reft_path = sorted(Path(reft_dir).glob(f"*{ssscc}.cap"))[0]
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

    reft_df = pd.DataFrame.from_records(reftArray)
    reft_df = reft_df.replace(  # remove text columns, only need numbers and dates
        to_replace=["bn", "diff", "val", "t90", "="], value=np.nan
    )
    reft_df = reft_df.dropna(axis=1)
    reft_df[1] = reft_df[[1, 2, 3, 4]].agg(" ".join, axis=1)  # dd/mm/yy/time cols are
    reft_df.drop(columns=[2, 3, 4], inplace=True)  # read separately; combine into one
    cols = OrderedDict(  # having this as a dict streamlines next steps
        [
            ("index_memory", int),
            ("datetime", object),
            ("btl_fire_num", int),
            ("diff", int),
            ("raw_value", float),
            ("T90", float),
        ]
    )
    reft_df.columns = list(cols.keys())  # name columns
    reft_df = reft_df.astype(cols)  # force dtypes
    # assign initial qality flags
    reft_df.loc[:, "REFTMP_FLAG_W"] = 2
    reft_df.loc[abs(reft_df["diff"]) >= 3000, "REFTMP_FLAG_W"] = 3
    # add in STNNBR, CASTNO columns
    # TODO: should these be objects or floats? be consistent!
    # string prob better for other sta/cast formats (names, letters, etc.)
    reft_df["STNNBR"] = int(ssscc[0:3])  # int for now while spinning up xarray
    reft_df["CASTNO"] = int(ssscc[3:5])  # will need to sort out salt naming (only int?)

    # convert to dataset and add attrs
    reft_ds = xr.Dataset.from_dataframe(reft_df.set_index("btl_fire_num"))
    reft_ds = reft_ds.set_coords(["STNNBR", "CASTNO"])
    reft_ds["T90"].attrs = {
        "sensor_type": "sbe_35",  # TODO: double check this info
        "standard_name": "Sea-Bird SBE 35 thermometer",
        "averaging_period_seconds": 15,
        "ancillary_variables": "REFTMP_FLAG_W",  # do this here or at end?
    }

    return reft_ds


def process_reft(ssscc_list, reft_dir=cfg.directory["reft"]):
    # TODO: import reft_dir from a config file
    """
    SBE35 reference thermometer processing function. Load in .cap files for given
    station/cast list, perform basic flagging, and export to .nc files.

    Inputs
    ------
    ssscc_list : list of str
        List of stations to process
    reft_dir : str, optional
        Path to folder containing raw salt files (defaults to data/reft/)

    """
    print("Processing reft files")
    for ssscc in ssscc_list:
        if not Path(reft_dir + ssscc + "_reft.nc").exists():
            try:
                reft_ds = _reft_loader(ssscc, reft_dir)
                reft_ds.to_netcdf(reft_dir + ssscc + "_reft.nc")
            except FileNotFoundError:
                print("refT file for cast " + ssscc + " does not exist... skipping")
                return

def _load_reft_data(reft_file, index_name="btl_fire_num"):
    """
    Loads reft_file to dataframe and reindexes to match bottle data dataframe

    Note: loading in REFTMP_FLAG_W here will conflict with the REFTMP_FLAG_W
    determined during temperature calibration.
    """
    reft_data = pd.read_csv(reft_file, usecols=["btl_fire_num","T90"])
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
    btl_data = dataToNDarray(btl_file,float,True,',',0)
    btl_data = pd.DataFrame.from_records(btl_data)
    if cols != None:
        btl_data = btl_data[cols]
    btl_data["SSSCC"] = Path(btl_file).stem.split("_")[0]

    return btl_data


# MK: deprecated 04/28/20, use load_all_ctd_files
# def load_time_data(time_file):

#     time_data = dataToNDarray(time_file,float,True,',',1)
#     time_data = pd.DataFrame.from_records(time_data)

#     return time_data


# def calibrate_param(param,ref_param,press,calib,order,ssscc,btl_num,xRange=None,):
def calibrate_param(param,ref_param,press,ssscc,btl_num,xRange=None,):
### NOTE: REF VALUES DEEMED QUESTIONABLE ARE STILL BEING USED FOR CALIBRATION


    df_good = quality_check(param,ref_param,press,ssscc,btl_num,find='good')
    df_ques = quality_check(param,ref_param,press,ssscc,btl_num,find='quest')

    df_ques['Parameter'] = param.name


    #report questionable data to a csv file

    #constrain pressure to within limits of xRange

    if xRange != None:
        x0 = int(xRange.split(":")[0])
        x1 = int(xRange.split(":")[1])

        df_good_cons = df_good[(df_good[press.name] >= x0) & (df_good[press.name] <= x1)]


    else:
        #Take full range of temperature values
        x0 = df_good[param.name].min()
        x1 = df_good[param.name].max()

        df_good_cons = df_good[(df_good[param.name] >= x0) & (df_good[param.name] <= x1)]

    return df_good_cons,df_ques
    # if 'P' in calib:
    #     coef = get_param_coef(df_good_cons[press.name],df_good_cons['Diff'],order,calib)
    # elif 'T' or 'C' in calib:
    #     coef = get_param_coef(df_good_cons[param.name],df_good_cons['Diff'],order,calib)
    # else:
    #     print('calib argument not valid, use CP TP T or C')

    # return coef,df_ques

def quality_check(param,param_2,press,ssscc,btl_num,find,thresh=[0.002, 0.005, 0.010, 0.020]):


    param = fit_ctd.array_like_to_series(param)
    param_2 = fit_ctd.array_like_to_series(param_2)
    press = fit_ctd.array_like_to_series(press)
    ssscc = fit_ctd.array_like_to_series(ssscc)
    btl_num = fit_ctd.array_like_to_series(btl_num)

    diff = param_2 - param

    df = pd.concat([ssscc,btl_num.rename('Bottle'),param.rename('Param_1'),param_2.rename('Param_2'),press.rename('CTDPRS'),diff.rename('Diff')],axis=1)

    if find == 'good':
    # Find data values for each sensor that are below the threshold (good)
        df['Flag'] = 1
        #df_range_comp = df_range[(df_range[diff].abs() < threshold)]# & (df_range[d_2].abs() < threshold) & (df_range[d_12].abs() < threshold)]
        df.loc[(df.CTDPRS > 2000) & (df.Diff.abs() < thresh[0]), 'Flag'] = 2
        df.loc[(df.CTDPRS <= 2000) & (df.CTDPRS >1000) & (df.Diff.abs() < thresh[1]), 'Flag'] = 2
        df.loc[(df.CTDPRS <= 1000) & (df.CTDPRS >500) & (df.Diff.abs() < thresh[2]), 'Flag'] = 2
        df.loc[(df.CTDPRS <= 500)  & (df.Diff.abs() < thresh[3]), 'Flag'] = 2
#
        # Filter out bad values

        df = df[df['Flag'] == 2]

        # Rename Columns back to what they were

        if param.name != None:
            df.rename(columns = {'Param_1' : param.name}, inplace=True)

        if param_2.name != None:
            df.rename(columns = {'Param_2' : param_2.name},inplace=True)

        if press.name != None:
            df.rename(columns = {'CTDPRS' : press.name}, inplace=True )

    elif find == 'quest':
    # Find data values for each sensor that are above the threshold (questionable)

        df['Flag'] = 1
        df.loc[(df.CTDPRS > 2000) & (df.Diff.abs() > thresh[0]), 'Flag'] = 3
        df.loc[(df.CTDPRS <= 2000) & (df.CTDPRS >1000) & (df.Diff.abs() > thresh[1]), 'Flag'] = 3
        df.loc[(df.CTDPRS <= 1000) & (df.CTDPRS >500) & (df.Diff.abs() > thresh[2]), 'Flag'] = 3
        df.loc[(df.CTDPRS <= 500)  & (df.Diff.abs() > thresh[3]), 'Flag'] = 3

        # Filter out good values

        df = df[df['Flag'] == 3]

        # Remove unneeded columns

        df = df.drop(['Param_1','Param_2'],axis=1)

        # Re-Order Columns for better readability

        df = df[[ssscc.name,'Bottle',press.name,'Flag','Diff']]


    else:
        print('Find argument not valid, please enter "good" or "quest" to find good or questionable values')

    return df

def get_param_coef(calib_param,diff,order,calib):


    cf1 = np.polyfit(calib_param, diff, order)


    if 'T' in calib:
        coef = np.zeros(shape=5)

        if order == 0:
            coef[4] = cf1[0]

        elif (order == 1) and (calib == 'TP'):
            coef[1] = cf1[0]
            coef[4] = cf1[1]

        elif (order == 2) and (calib == 'TP'):
            coef[0] = cf1[0]
            coef[1] = cf1[1]
            coef[4] = cf1[2]

        elif (order == 1) and (calib == 'T'):
            coef[3] = cf1[0]
            coef[4] = cf1[1]

        elif (order == 2) and (calib == 'T'):
            coef[2] = cf1[0]
            coef[3] = cf1[1]
            coef[4] = cf1[2]

    if 'C' in calib:
        coef = np.zeros(shape=7)
        if order == 0:
            coef[6] = cf1[0]
        elif (order == 1) and (calib == 'CP'):
            coef[1] = cf1[0]
            coef[6] = cf1[1]
        elif (order == 2) and (calib == 'CP'):
            coef[0] = cf1[0]
            coef[1] = cf1[1]
            coef[6] = cf1[2]
        elif (order == 1) and (calib == 'C'):
            coef[5] = cf1[0]
            coef[6] = cf1[1]
        elif (order == 2) and (calib == 'C'):
            coef[4] = cf1[0]
            coef[5] = cf1[1]
            coef[6] = cf1[2]

    return coef

def combine_quality_flags(df_list):

    combined_df = pd.concat(df_list)
    combined_df = combined_df.sort_values(['SSSCC','Bottle'])

    combined_df = combined_df.round(4)

    return combined_df

    #Combine these three into a dataframe and write out to a csv
    #Sort by sta/cast, bottle number, rev. press


def calibrate_conductivity(df,order,calib_param,sensor,xRange=None,
                           refc_col='BTLCOND',cond_col_1='CTDCOND1',cond_col_2='CTDCOND2',
                           p_col='CTDPRS'):#refc_data
### NOTE: REF VALUES DEEMED QUESTIONABLE ARE STILL BEING USED FOR CALIBRATION
    if sensor == 1:
        postfix = 'c1'
        cond_col = 'CTDCOND1'
        t_col = 'CTDTMP1'
    elif sensor ==2:
        postfix = 'c2'
        cond_col = 'CTDCOND2'
        t_col = 'CTDTMP2'
    else:
        print('No sensor name supplied, difference column name will be: diff')

    if calib_param == 'P':
        calib_col = p_col
    elif calib_param == 'T':
        calib_col = t_col
    elif calib_param == 'C':
        calib_col = cond_col
    else:
        print('No calib_param supplied')

    diff = 'd_'+postfix #Difference between ref and prim sensor

    # Calculate absolute differences between sensors and salt sample data

    #df[diff] = refc_data[refc_col] - df[cond_col]
    df[diff] = df[refc_col] - df[cond_col]

    #df['primary_diff'] = refc_data[refc_col] - df[cond_col_1]
    df['primary_diff'] = df[refc_col] - df[cond_col_1]

    #df['secondary_diff'] = refc_data[refc_col] - df[cond_col_2]
    df['secondary_diff'] = df[refc_col] - df[cond_col_2]

    df['P-S'] = df[cond_col_1] - df[cond_col_2]



    #Greater than 2000 dBar
    lower_lim = 2000
    upper_lim = df[p_col].max()
    threshold = 0.002


    df_deep_good = quality_check(df,diff,lower_lim,upper_lim,threshold)
    df_deep_ques = quality_check(df,diff,lower_lim,upper_lim,threshold,find='quest')
    df_deep_ref = quality_check(df,diff,lower_lim,upper_lim,threshold,find='ref')


    #Between 2000 and 1000
    lower_lim = 1000
    upper_lim = 2000
    threshold = 0.005


    df_lmid_good = quality_check(df,diff,lower_lim,upper_lim,threshold)
    df_lmid_ques = quality_check(df,diff,lower_lim,upper_lim,threshold,find='quest')
    df_lmid_ref = quality_check(df,diff,lower_lim,upper_lim,threshold,find='ref')


    #Between 1000 and 500
    lower_lim = 500
    upper_lim = 1000
    threshold = 0.010



    df_umid_good = quality_check(df,diff,lower_lim,upper_lim,threshold)
    df_umid_ques = quality_check(df,diff,lower_lim,upper_lim,threshold,find='quest')
    df_umid_ref = quality_check(df,diff,lower_lim,upper_lim,threshold,find='ref')


    #Less than 500
    lower_lim = df[p_col].min() - 1
    upper_lim = 500
    threshold = 0.020

    df_shal_good = quality_check(df,diff,lower_lim,upper_lim,threshold)
    df_shal_ques = quality_check(df,diff,lower_lim,upper_lim,threshold,find='quest')
    df_shal_ref = quality_check(df,diff,lower_lim,upper_lim,threshold,find='ref')

    #concat dataframes into two main dfs
    df_good = pd.concat([df_deep_good,df_lmid_good,df_umid_good,df_shal_good])
    df_ques = pd.concat([df_deep_ques,df_lmid_ques,df_umid_ques,df_shal_ques])
    df_ref = pd.concat([df_deep_ref,df_lmid_ref,df_umid_ref,df_shal_ref])

    if sensor == 1:
        df_ques['Parameter'] = 'C1'
        df_ques['Flag'] = 3

        df_ref['Parameter'] = 'C'
        df_ref['Flag'] = 3

    elif sensor == 2:
        df_ques['Parameter'] = 'C2'
        df_ques['Flag'] = 3

        df_ref['Flag'] = 3

    if xRange != None:
        x0 = int(xRange.split(":")[0])
        x1 = int(xRange.split(":")[1])

        df_good_cons = df_good[(df_good[calib_col] >= x0) & (df_good[calib_col] <= x1)]


    else:
        #Take full range of temperature values
#        x0 = df_good[t_col].min()
#        x1 = df_good[t_col].max()

        df_good_cons = df_good#[(df_good[calib_col] >= x0) & (df_good[calib_col] <= x1)]

    cf = np.polyfit(df_good_cons[calib_col], df_good_cons[diff], order)

    sensor = '_c'+str(sensor)
    coef = np.zeros(shape=7)

    if order == 0:
        coef[6] = cf[0]
    elif (order == 1) and (calib_param == 'P'):
        coef[1] = cf[0]
        coef[6] = cf[1]
    elif (order == 2) and (calib_param == 'P'):
        coef[0] = cf[0]
        coef[1] = cf[1]
        coef[6] = cf[2]
    elif (order == 1) and (calib_param == 'T'):
        coef[3] = cf[0]
        coef[6] = cf[1]
    elif (order == 2) and (calib_param == 'T'):
        coef[2] = cf[0]
        coef[3] = cf[1]
        coef[6] = cf[2]
    elif (order == 1) and (calib_param == 'C'):
        coef[5] = cf[0]
        coef[6] = cf[1]
    elif (order == 2) and (calib_param == 'C'):
        coef[4] = cf[0]
        coef[5] = cf[1]
        coef[6] = cf[2]
    return coef,df_ques,df_ref


def prepare_fit_data(df,ref_col):

    good_data = df.copy()
    good_data = good_data[np.isfinite(good_data[ref_col])]

    return good_data



def prepare_conductivity_data(ssscc,df,refc,ssscc_col = 'SSSCC',index_col = 'btl_fire_num'):

    btl_concat = pd.DataFrame()
    for x in ssscc:
        btl_data = df[df[ssscc_col] == x]
        refc_data = refc[refc[ssscc_col] == x]
        btl_data_clean = prepare_fit_data(btl_data,refc_data,'C')
        btl_concat = pd.concat([btl_concat,btl_data_clean])
    refc = refc[refc[index_col] != 0]
    refc = refc.reset_index(drop=True)
    btl_concat = btl_concat.reset_index(drop=True)

    return btl_concat, refc

def prepare_all_fit_data(ssscc,df,ref_data,param):

    data_concat = pd.DataFrame()


    for x in ssscc:
        btl_data = df[df['SSSCC']==x]
        ref_data_stn= ref_data[ref_data['SSSCC']==x]
        btl_data_good = prepare_fit_data(btl_data,ref_data_stn,param)
        data_concat = pd.concat([data_concat,btl_data_good])

    return data_concat


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
    p_log = pd.read_csv(cfg.directory["logs"] + "ondeck_pressure.csv", dtype={"SSSCC":str})
    p_offset = _get_pressure_offset(p_log.ondeck_start_p, p_log.ondeck_end_p)
    df[p_col] += p_offset
    df[p_col + "_FLAG_W"] = 2

    return df


def make_depth_log(time_df):
    # TODO: get column names from config file
    """
    Create depth log file from maximum depth of each station/cast in time DataFrame.

    Parameters
    ----------
    time_df : DataFrame
        DataFrame containing continuous CTD data

    """
    ssscc_list = time_df["SSSCC"].unique()
    depth_dict = {}
    for ssscc in ssscc_list:
        time_rows = time_df["SSSCC"] == ssscc
        # TODO: improve error handling s.t. SSSCC is reported with error message for
        # stations with altimeter readings below threshold (in_find_cast_depth)
        max_depth = _find_cast_depth(
            ssscc,
            time_df.loc[time_rows, "CTDPRS"],
            time_df.loc[time_rows, "GPSLAT"],
            time_df.loc[time_rows, "ALT"],
        )
        depth_dict[ssscc] = max_depth

    depth_df = pd.DataFrame.from_dict(depth_dict, orient="index")
    depth_df.reset_index(inplace=True)
    depth_df.rename(columns={0: "DEPTH", "index": "SSSCC"}, inplace=True)
    depth_df.to_csv(cfg.directory["logs"] + "depth_log.csv", index=False)
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
                reft_data = pd.DataFrame(index=btl_data.index, columns=["T90"])
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

            # Calculate dv/dt for oxygen fitting
            btl_data['dv_dt'] = oxy_fitting.calculate_dVdT(btl_data['CTDOXYVOLTS'],btl_data['scan_datetime'])

            # Add bottom of cast information (date,time,lat,lon,etc.)
            btl_data = _add_btl_bottom_data(btl_data, ssscc, cfg.directory["logs"] + "cast_details.csv")

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

def merge_refcond_flags(btl_data, qual_flag_cond):
    # Merge df
    mask = qual_flag_cond[qual_flag_cond['Parameter'] == 'REF_COND'].copy()
    mask['SSSCC'] = mask['SSSCC'].astype(str)
    btl_data = btl_data.merge(mask,left_on=['SSSCC','btl_fire_num'], right_on=['SSSCC','Bottle'],how='left')
    # Rename Columns
    btl_data.rename(columns={'CTDPRS_x':'CTDPRS','SSSCC_x':'SSSCC','Flag':'SALNTY_FLAG_W'},inplace=True)
    btl_data.drop(columns=['Parameter','CTDPRS_y','Bottle','Diff'],inplace=True)
    btl_data['SALNTY_FLAG_W'].fillna(value=2,inplace=True)
    try:
        btl_data.loc[btl_data['BTLCOND'].isna(),'SALNTY_FLAG_W'] = 9
    except:
        btl_data[btl_data['SALNTY'].isna(),'SALNTY_FLAG_W'] = 9
    btl_data['SALNTY_FLAG_W'] = btl_data['SALNTY_FLAG_W'].astype(int)

    return btl_data

def merge_cond_flags(btl_data, qual_flag_cond,parameter):
    # Merge df
    #if sensor == 1:
    #    parameter = 'CTDCOND1'
    #elif sensor == 2:
    #    parameter = 'CTDCOND2'
    mask = qual_flag_cond[qual_flag_cond['Parameter'] == parameter].copy()
    mask['SSSCC'] = mask['SSSCC'].astype(str)
    btl_data = btl_data.merge(mask,left_on=['SSSCC','btl_fire_num'], right_on=['SSSCC','Bottle'],how='left')
    # Rename Columns
    btl_data.rename(columns={'CTDPRS_x':'CTDPRS','SSSCC_x':'SSSCC','Flag':'CTDSAL_FLAG_W'},inplace=True)
    btl_data.drop(columns=['Parameter','CTDPRS_y','Bottle','Diff'],inplace=True)
    btl_data['CTDSAL_FLAG_W'].fillna(value=2,inplace=True)
    btl_data.loc[btl_data[parameter].isna(),'CTDSAL_FLAG_W'] = 9
    btl_data['CTDSAL_FLAG_W'] = btl_data['CTDSAL_FLAG_W'].astype(int)

    return btl_data

def merged_reftemp_flags(btl_data, qual_flag_temp):

    mask = qual_flag_temp[qual_flag_temp['Parameter'] == 'REF_TEMP'].copy()
    mask['SSSCC'] = mask['SSSCC'].astype(str)
    btl_data = btl_data.merge(mask,left_on=['SSSCC','btl_fire_num'], right_on=['SSSCC','Bottle'],how='left')
    # Rename Columns
    btl_data.rename(columns={'CTDPRS_x':'CTDPRS','SSSCC_x':'SSSCC','Flag':'REFTMP_FLAG_W'},inplace=True)
    btl_data.drop(columns=['Parameter','CTDPRS_y','Bottle','Diff'],inplace=True)
    btl_data['REFTMP_FLAG_W'].fillna(value=2,inplace=True)
    try:
        btl_data.loc[btl_data['T90'].isna(),'REFTMP_FLAG_W'] = 9
    except:
        btl_data.loc[btl_data['REFTMP'].isna(),'REFTMP_FLAG_W'] = 9
        btl_data.loc[btl_data['REFTMP'].isna(),'REFTMP'] = -999
    #btl_data['REFTMP_FLAG_W'] = btl_data['REFTMP_FLAG_W'].astype(int)

    return btl_data

def merge_temp_flags(btl_data, qual_flag_temp, parameter):

    mask = qual_flag_temp[qual_flag_temp['Parameter'] == parameter].copy()
    mask['SSSCC'] = mask['SSSCC'].astype(str)
    btl_data = btl_data.merge(mask,left_on=['SSSCC','btl_fire_num'], right_on=['SSSCC','Bottle'],how='left')
    # Rename Columns
    btl_data.rename(columns={'CTDPRS_x':'CTDPRS','SSSCC_x':'SSSCC','Flag':'CTDTMP_FLAG_W'},inplace=True)
    btl_data.drop(columns=['Parameter','CTDPRS_y','Bottle','Diff'],inplace=True)
    btl_data['CTDTMP_FLAG_W'] = btl_data['CTDTMP_FLAG_W'].fillna(value=2)
    btl_data['CTDTMP_FLAG_W'] = btl_data['CTDTMP_FLAG_W'].astype(int)

    return btl_data

def merge_oxy_flags(btl_data):

    mask = (btl_data['OXYGEN'].isna())
    btl_data.loc[mask,'OXYGEN_FLAG_W'] = 9

def _find_cast_depth(ssscc,press,lat,alt,threshold=80):
    """
    Calculate the depth of a given cast. If rosette does not get within the threshold
    distance of the bottom, returns NaN.

    Parameters
    -----------
    ssscc : str
        Current station/cast name
    press : array-like
        CTD pressure
    lat : array-like
        Ship latitude
    alt : array-like
        CTD altimeter reading
    threshold : int, optional
        Maximum altimeter reading to consider cast "at the bottom" (defaults to 80)

    Returns
    --------
    max_depth : int
        Maximum depth from the cast

    """
    # Create Dataframe containing args
    df = pd.DataFrame({"CTDPRS": press, "LAT": lat, "ALT": alt})
    # Calculate DEPTH using gsw
    df['DEPTH'] = np.abs(gsw.z_from_p(df['CTDPRS'],df['LAT']))

    # Find max depth and see if ALT has locked in
    bottom_alt = df.loc[df['CTDPRS'] == df['CTDPRS'].max(),'ALT']
    if bottom_alt.values[0] <= threshold:
        max_depth = bottom_alt + df['CTDPRS'].max()
        max_depth = int(max_depth.values[0])
    else:
        print(
            ssscc,
            "- minimum altimeter reading above",
            threshold,
            "m threshold, setting max depth to -999",
        )
        max_depth = -999

    return max_depth

def format_time_data(df):

    format_columns = settings.pressure_series_output['column_names'].copy()
    if 'SSSCC' not in format_columns:
        print('Adding SSSCC')
        format_columns.append('SSSCC')
    if 'DEPTH' not in format_columns:
        print('Adding DEPTH')
        format_columns.append('DEPTH')
    try:
        df = df[format_columns]
    except:
        print('missing required pressure series output columns!')
        df_columns = list(df.keys())
        missing_list = list(np.setdiff1d(format_columns,df_columns))
        raise KeyError('missing columns: ',missing_list)

    return df

def add_btlnbr_cols(df,btl_num_col):
    df['BTLNBR'] = df[btl_num_col].astype(int)
    # default to everything being good
    df['BTLNBR_FLAG_W'] = 2
    return df

def castno_from_ssscc(ssscc):
    # ssscc: column (pandas series) containing station and cast numbers in the SSSCC format
    ssscc = pd.Series(ssscc)
    castno = ssscc.str[3:].astype(int)
    return castno

def stnnbr_from_ssscc(ssscc):
    # ssscc: column (pandas series) containing station and cast numbers in the SSSCC format
    ssscc = pd.Series(ssscc)
    stnno = ssscc.str[0:3].astype(int)
    return stnno

def add_sampno_col(df,btl_num_col):
    df['SAMPNO'] = df[btl_num_col].astype(int)
    return df

def get_btl_time(df,btl_num_col,time_col):
    # Get time for first btl fire
    time = df[df[btl_num_col]== df[btl_num_col].min()][time_col].values
    ts = pd.to_datetime(time,unit='s')
    date = ts.strftime('%Y%m%d')
    hour= ts.strftime('%H%M')
    df['DATE'] = date[0]
    df['TIME'] = hour[0]

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

def flag_missing_btl_values(df,flag_columns,flag_suffix='_FLAG_W'):
    for column in flag_columns:
        flag_name = column + flag_suffix
        df.loc[df[column].isna(),flag_name] = 9
        df[column].fillna(value=int(-999), inplace=True)
        df.loc[df[column].astype(int) == -999, flag_name] = 9
    return df

def format_btl_data(df,data_cols, prcn=4):
    """
    data_cols :(list) list containing the "numerical" columns to be rounded to prcn decimal places
    """
    # Choose correct sensors
    #df['CTDTMP'] = df[settings.bottle_inputs['t']]
    format_columns = settings.btl_series_output['btl_column_names'].copy()
#    if 'SSSCC' not in format_columns:
#        print('Adding SSSCC')
#        format_columns.append('SSSCC')
    if 'DEPTH' not in format_columns:
        print('Adding DEPTH')
        format_columns.append('DEPTH')
    try:
        df = df[format_columns].copy()
    except:
        print('missing required pressure series output columns!')
        df_columns = list(df.keys())
        missing_list = list(np.setdiff1d(format_columns,df_columns))
        raise KeyError('missing columns: ',missing_list)
        #print('missing columns: ',missing_list)

    for i in df.columns:
        if '_FLAG_W' in i:
            df[i] = df[i].astype(int)
        elif i in data_cols:
            df[i] = df[i].round(prcn)
    return df

def manual_backfill_values(df,bfill_prs,p_col='CTDPRS',flag_suffix='_FLAG_W'):
    col_list = df.columns.tolist()
    col_list.remove(p_col)
    df.loc[df[p_col] < bfill_prs,col_list] = np.nan
    press_flag_column = p_col + flag_suffix
    for column in col_list:
        if flag_suffix in column:
            if press_flag_column in column:
                pass
            else:
                df.loc[df[p_col] < bfill_prs,column] = 6
    df.bfill(inplace=True)
    return df

def flag_backfill_data(df,p_col='CTDPRS',flag_bol_col='interp_bol',flag_suffix='_FLAG_W'):
    col_list = df.columns.tolist()
    for column in col_list:

        if flag_suffix in column:
            df.loc[df[flag_bol_col] == 1, column] = 6
    return df

def flag_missing_values(df,flag_suffix='_FLAG_W'):
    full_col_list = df.columns.tolist()
    col_list = full_col_list

    for column in full_col_list:
        if flag_suffix in column:
            col_list.remove(column)
    for column in col_list:
        flag_name = column + flag_suffix

        df.loc[df[column].isna(),flag_name] = 9
        df.loc[df[column].astype(int) == -999, flag_name] = 9
    return df

def export_bin_data(df, ssscc, sample_rate, search_time, p_column_names, p_col='CTDPRS', ssscc_col='SSSCC', bin_size=2, direction='down'):
    # remove
    df_binned = pd.DataFrame()
    for cast in ssscc:
        time_data = df.loc[df[ssscc_col] == cast].copy()
        time_orig = time_data.copy()
        time_data = pressure_sequence(time_data,p_col,2.0,-1.0,0.0,'down',sample_rate,search_time)
        if time_data[p_col].hasnans:
            time_orig['CTDOXY'] = pd.to_numeric(time_orig['CTDOXY'])
            time_data = binning_df(time_orig, bin_size=2)
            time_data['interp_bol'] = 0
            time_data.loc[time_data['CTDPRS'].isnull(),'interp_bol'] = 1
            time_data['CTDPRS'] = time_data.index.astype('float')
        time_data = flag_backfill_data(time_data)
        time_data = fill_surface_data(time_data)
        time_data = time_data[p_column_names]
        time_data = time_data.round(4)
        try:
            time_data = flag_missing_values(time_data)
        except KeyError:
            raise KeyError('missing columns: ',missing_list)
        for i in time_data.columns:
            if '_FLAG_W' in i:
                time_data[i] = time_data[i].astype(int)
            else:
                time_data[i] = time_data[i].round(4)
        time_data[ssscc_col] = cast

        df_binned = pd.concat([df_binned,time_data])
    return df_binned

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
    # TODO: make this better... (see process_ctd.format_time_data())
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
    depth_df = pd.read_csv(cfg.directory["logs"] + 'depth_log.csv', dtype={"SSSCC": str}).dropna()
    try:
        manual_depth_df = pd.read_csv(cfg.directory["logs"] + 'manual_depth_log.csv', dtype={"SSSCC": str})
    except FileNotFoundError:
        # TODO: add logging; look into inheriting/extending a class to add features
        print("manual_depth_log.csv not found... duplicating depth_log.csv")
        manual_df = depth_df.copy()  # write manual_depth_log as copy of depth_log
        manual_df.to_csv(cfg.directory["logs"] + 'manual_depth_log.csv', index=False)
    full_depth_df = pd.concat([depth_df,manual_depth_df])
    full_depth_df.drop_duplicates(subset='SSSCC', keep='first',inplace=True)

    for ssscc in ssscc_list:

        time_data = df[df['SSSCC'] == ssscc].copy()
        time_data = pressure_sequence(time_data)
        time_data = flag_backfill_data(time_data)
        time_data = fill_surface_data(time_data)
        time_data = time_data[p_column_names]
        time_data = time_data.round(4)

        depth = full_depth_df.loc[full_depth_df['SSSCC'] == ssscc,'DEPTH'].iloc[0]
        # get cast_details for current SSSCC
        cast_dict = cast_details[cast_details["SSSCC"] == ssscc].to_dict("r")[0]
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
                f"DEPTH = {depth}\n"
            )
            f.write(ctd_header)
            np.asarray(p_column_names).tofile(f, sep=',', format='%s')
            f.write('\n')
            np.asarray(p_column_units).tofile(f, sep=',', format='%s')
            f.write('\n')
            time_data.to_csv(f, header=False, index=False)
            f.write("END_DATA")


def export_btl_data(df,expocode,btl_columns, btl_units, sectionID,out_dir=cfg.directory["pressure"],org='ODF'):

    btl_data = df.copy()
    now = datetime.now()
    file_datetime = now.strftime("%Y%m%d")

    time_stamp = file_datetime+org

    outfile = open(out_dir + expocode + '_hy1.csv', mode='w+')
    outfile.write("BOTTLE, %s\n" % (time_stamp))
    cn = np.asarray(btl_columns)
    cn.tofile(outfile,sep=',', format='%s')
    outfile.write('\n')
    cu = np.asarray(btl_units)
    cu.tofile(outfile,sep=',', format='%s')
    outfile.write('\n')
    outfile.close()

    file = out_dir + expocode + '_hy1.csv'
    with open(file,'a') as f:
        btl_data.to_csv(f, header=False,index=False)
    f.close()

    outfile = open(out_dir + expocode + '_hy1.csv', "a")
#    outfile.write('\n')
    outfile.write('END_DATA')
    outfile.close()

    return

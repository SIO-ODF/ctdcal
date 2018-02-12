#!/usr/bin/env python
from scipy import interpolate
import numpy as np
from numpy.lib.recfunctions import append_fields
import scipy.signal as sig
import scipy.stats as st
import time, os
import pandas as pd
import math
#import report_ctd
import ctdcal.report_ctd as report_ctd
import warnings

import gsw

warnings.filterwarnings("ignore", 'Mean of empty slice.')

def cast_details(stacast, log_file, p_col, time_col, b_lat_col, b_lon_col, alt_col, inMat=None):
    '''
    We determine the cast details using pandas magic.
    First find alternating periods of pumps on and pumps off, then select the
    pumps on period with the highest pressure. Get values from the row with the
    highest pressure, and return all values to be sent to log.

    Input:
    stacast - integer, the station and cast, as SSSCC format
    log_file - file handle or string, log_file
    p_col - string, name of the pressure column
    time_col - string, name of the time column
    b_lat_col - string, name of the latitude column
    b_lon_col - string, name of the longitude column
    alt_col - string, name of the altimeter column
    inMat - pandas dataframe, the dataframe to come in

    Output:
    start_cast_time - float, unix epoch seconds?, start of cast time, to be reported to log file
    end_cast_time - float, unix epoch seconds?, end of cast time, to be reported to log file
    bottom_cast_time - float, unix epoch seconds?, bottom of cast time, to be reported to log file
    start_pressure - float, pressure at which cast started, to be reported to log file
    max_pressure - float, bottom of the cast pressure, to be reported to log file
    b_lat - float, latitude at bottom of cast
    b_lon - float, longitude at bottom of cast
    b_alti - float, altimeter reading at bottom of cast - volts only!
    inMat - the dataframe that came in, with soak period trimmed off

    don't need end_cast_time, max_pressure
    inMat is trimmed to start and end of cast
    '''

    df_test = pd.DataFrame.from_records(inMat)

    dfs = find_pump_on_off_dfs(df_test)
    dfs_1 = find_pumps_on_dfs(dfs)
    df_cast = find_max_pressure_df(dfs_1)
    df_cast2 = trim_soak_period_from_df(df_cast)

    start_cast_time = float(df_cast2['scan_datetime'].head(1))
    start_pressure = float(df_cast2['CTDPRS'].head(1))
    end_cast_time = float(df_cast2['scan_datetime'].tail(1))
    max_pressure = float(df_cast2['CTDPRS'].max())
    bottom_cast_time = float(df_cast2.loc[df_cast2['CTDPRS'].idxmax()]['scan_datetime'])
    b_lat = float(df_cast2.loc[df_cast2['CTDPRS'].idxmax()]['GPSLAT'])
    b_lon = float(df_cast2.loc[df_cast2['CTDPRS'].idxmax()]['GPSLON'])
    b_alti = float(df_cast2.loc[df_cast2['CTDPRS'].idxmax()]['ALT'])

    #last two lines must be in to return the same as old - change to slices of df later
    report_ctd.report_cast_details(stacast, log_file, start_cast_time, end_cast_time,
                                   bottom_cast_time, start_pressure, max_pressure, b_alti,
                                   b_lat, b_lon)
    #reconvert to ndarray - might need to be altered to remove second index
    # inMat = df_cast2.loc[:df_cast2['CTDPRS'].idxmax()].to_records(index=False)
    inMat = df_cast2.loc[:df_cast2['CTDPRS'].idxmax()]

    return start_cast_time, end_cast_time, bottom_cast_time, start_pressure, max_pressure, b_lat, b_lon, b_alti, inMat
#Move next four functions to a library or class(?) Clean up module
def find_pump_on_off_dfs(df):
    '''Find pump_on patterns of dataframes, and return a list(?) of dataframes to iterate over.
    '''
    return [g for i,g in df.groupby(df['pump_on'].ne(df['pump_on'].shift()).cumsum())]

def find_max_pressure_df(dfs):
    '''Giving a list of data frames, return a reference to the frame with which contians the highest pressure value
    '''
    max_pressure_df = dfs[0]
    max_pressure = max_pressure_df['CTDPRS'].max() #TODO make into config var
    for df in dfs:
        if df['CTDPRS'].max() > max_pressure:
            max_pressure_df = df
    return max_pressure_df

def find_pumps_on_dfs(dfs):
    '''given a list of dataframes, remove all the frames with one or more rows containing a "false" pump on flag
    '''
    return list(filter(lambda df: df['pump_on'].all(), dfs))

def trim_soak_period_from_df(df):
    '''Look for minimum pressure in dataframe, then return everything after minimum pressure/top of cast.
    '''
    test = int(df.loc[1:(len(df)/4),['CTDPRS']].idxmin())
    return df.loc[test:]
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

def ctd_align(inMat=None, col=None, time=0.0):
    """ctd_align function

    Function takes full NUMPY ndarray with predefined dtype array
    and adjusts time of sensor responce and water flow relative to
    the time frame of temperature sensor.

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

        q_nd = q_df.as_matrix(columns=q_df.columns)
    return q_nd

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

def raw_ctd_filter(input_array=None, filter_type='triangle', win_size=24, parameters=None):
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

    if input_array is None:
        print("In raw_ctd_filter: No data array.")
        return
    else:
        return_array = input_array
        if parameters is None:
            print("In raw_ctd_filter: Empty parameter list.")
        else:
            for p in parameters:
                if filter_type is 'boxcar':
                    win = sig.boxcar(win_size)
                    return_array[str(p)] = sig.convolve(input_array[str(p)], win, mode='same')/len(win)
                elif filter_type is 'gaussian':
                    sigma = np.std(arr)
                    win = sig.general_gaussian(win_size, 1.0, sigma)
                    return_array[str(p)] = sig.convolve(input_array[str(p)], win, mode='same')/(len(win))
                elif filter_type is 'triangle':
                    win = sig.triang(win_size)
                    return_array[p] = 2*sig.convolve(input_array[p], win, mode='same')/len(win)
    return return_array


def ondeck_pressure(stacast, p_col, c1_col, c2_col, time_col, inMat=None, conductivity_startup=20.0, log_file=None):
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
    tmpMat = []
    outMat = []
    tmp = 0
    start_p = 0.0
    n = 0
    ep = []
    end_p = 0.0

    # Frequency
    fl = 24
    fl2 = fl*2
    # One minute
    mt = 60
    # Half minute
    ms = 30
    time_delay = fl*ms

    if inMat is None:
        print("Ondeck_pressure function: No data.")
        return
    else:
        # Searches first quarter of matrix, uses conductivity
        # threshold min to capture startup pressure
        for j in range(0,int(len(inMat)/4)):
            if ((inMat[c1_col][j] < conductivity_startup) and (inMat[c2_col][j] < conductivity_startup)):
                tmp = j
                start_pressure.append(inMat[p_col][j])

        # Evaluate starting pressures
        if not start_pressure: start_p = "Started in Water"
        else:
            n = len(start_pressure)
            if (n > time_delay): start_p = np.average(start_pressure[fl2:n-(time_delay)])
            else: start_p = np.average(start_pressure[fl2:n])

        # Remove on-deck startup
        inMat = inMat[tmp:]

        tmp = len(inMat);
        # Searches last half of NDarray for conductivity threshold
        for j in range(int(len(inMat)*0.5), len(inMat)):
            if ((inMat[c1_col][j] < conductivity_startup) and (inMat[c2_col][j] < conductivity_startup)):
                ep.append(inMat[p_col][j])
                if (tmp > j): tmp = j

        # Evaluate ending pressures
        if (len(ep) > (time_delay)): end_p = np.average(ep[(time_delay):])
        else: end_p = np.average(ep[(len(ep)):])

        # Remove on-deck ending
        outMat = inMat[:tmp]

        # Store ending on-deck pressure
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


def roll_filter(p_col, inMat=None, up='down', frames_per_sec=24, search_time=15, **kwargs):
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

        if up is 'down':
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
        elif up is 'up':
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

def pressure_sequence(stacast, p_col, time_col, intP=2.0, startT=-1.0, startP=0.0, up='down', sample_rate=12, search_time=15, inMat=None,):
    """pressure_sequence function
    Function takes full NUMPY ndarray with predefined dtype array
    and several arguments to return a pressure sequenced data ndarray.
    Pressure sequencing includes rollfilter.
    Necissary inputs are input Matrix (inMat) and pressure interval (intP).
    The other inputs have default settings. The program will figure out
    specifics for those settings if left blank.
    Start time (startT), start pressure (startP) and up are mutually exclusive.
    If sensors are not not fully functional when ctd starts down cast
    analyst can select a later start time or start pressure but not both.
    There is no interpolation to the surface for other sensor values.
    'up' indicates direction for pressure sequence. If up is set startT and startP
    are void.
    Args:
        param1 (str): stacast, station cast input
        param2 (str): p_col, pressure column name
        param3 (str): time_col, time column name
        param4 (float): starting pressure interval
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

    ###forget the check, let it crash
    # if inMat is None:
    #     print("Pressure sequence function: No input data.")
    #     return
    # else:
    #     pF = inMat[p_col]
    #     full_length = len(pF)-1

        #btm = max(pF) # bottom max P ###not needed once we change bin code
        ### Following two lines are only for btmTime, which is used in edge case check
        #indBtm = np.argmax(pF) # bottom index
        #btmTime = inMat[time_col][indBtm] # bottom time

        ##  Initialise input parameters OLD
        # if (startT > 0.0).bool() and (startT > inMat[time_col][0]).bool(): #patched startT with int() cast, complaining about it being a series
        #     start = (np.abs(inMat[time_col] - startT)).argmin()
        #     lenP = np.arange(start,indBtm,1)
        #     end = indBtm
        #     prvPrs = inMat[p_col][start]
        #     if btmTime <= startT:
        #         print("-startT start time is greater than down cast time. Cast issue.")
        #         return
        # elif ((startP > 0.0) and (startP > pF[0])):
        #     start = (np.abs(inMat[p_col] - startP)).argmin()
        #     lenP = np.arange(start,indBtm,1)
        #     end = indBtm
        #     prvPrs = inMat[p_col][start]
        #     if btm <= startP:
        #         print("-startP start pressure is greater than bottom pressure. Cast issue.")
        #         import pdb; pdb.set_trace() ####
        #         return
        # elif up is 'up':
        #     start = indBtm
        #     end = full_length
        #     lenP = np.arange(start,end,1)
        #     prvPrs = btm
        # else:
        #     lenP = np.arange(0,indBtm,1)
        #     start = 0
        #     end = indBtm
        #     prvPrs = 0.0

    df = pd.DataFrame.from_records(inMat)
    start = 0
#try except to be patched out when serialization changes
    try:
        end = df['CTDPRS_DBAR'].idxmax()
        btm = df['CTDPRS_DBAR'].max()
    except KeyError:
        end = df['CTDPRS'].idxmax()
        btm = df['CTDPRS'].max()
    # Roll Filter

    roll_filter_matrix = roll_filter(p_col, inMat[:][start:end:sample_rate], up, sample_rate, search_time, start=start, end=end, full_matrix=inMat)
    ### Needs to be removed and replaced with pd.fillna()
    # Treat surface data.
    #roll_filter_matrix = treat_surface_data(p_col, sample_rate, roll_filter_matrix)
    df_roll_surface = fill_surface_data(roll_filter_matrix, bin_size=2)
    #bin_size should be moved into config
    binned_df = binning_df(df_roll_surface, bin_size=2)
    ### OLD code - bin data and find mean of each bin
    # # Frame Pressure Bins
    # pressure_bins = np.arange(0,int(btm),intP)
    # p_bin_index = np.digitize(roll_filter_matrix[p_col],pressure_bins)
    #
    # # Define binned output array
    # binned_matrix = np.empty(shape=(len(pressure_bins),), dtype=inMat.dtype)
    #
    # # Iterate over input data by column, sort bins and find mean of binned data.
    # for col in binned_matrix.dtype.names:
    #     if col == p_col:
    #         binned_matrix[col] = pressure_bins
    #     elif binned_matrix[col].dtype is np.dtype(np.float64):
    #         binned_matrix[col] = [roll_filter_matrix[col][p_bin_index == i].mean() for i in range(1,len(pressure_bins)+1)]
    #         # Interpolate over NaN or missing data
    #         if np.isnan(binned_matrix[col]).any():
    #             binned_matrix[col] = data_interpolater(binned_matrix[col])

    binned_matrix = binned_df.to_records(index=False)
    return binned_matrix
### Once serialization has been fixed, fix try/except to compact code

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
        df_merged = pd.merge(df_surface, df, on='CTDPRS_DBAR', how='outer')
    except KeyError:
        for x in range(1, int(np.floor(df.iloc[0]['CTDPRS'])), bin_size):
            surface_values.append(x)
        df_surface = pd.DataFrame({'CTDPRS': surface_values})
        df_merged = pd.merge(df_surface, df, on='CTDPRS', how='outer')

    return df_merged.fillna(method='bfill')

def load_reft_data(reft_file,index_name = 'index_memory'):
    """ Loads reft_file to dataframe and reindexes to match bottle data dataframe"""
    
    reft_data = pd.read_csv(reft_file)
    reft_data.set_index(index_name)
    
    return reft_data

def load_btl_data(btl_file):
    
    """ex. '/Users/k3jackson/p06e/data/bottle/00201_btl_mean.csv'"""
    
    btl_data = dataToNDarray(btl_file,float,True,',',0) 

    btl_data = pd.DataFrame.from_records(btl_data)
    
    return btl_data


def calibrate_temperature(df,order,reft_data,calib_param,sensor,xRange=None,
                          t_col_1 = 'CTDTMP1', t_col_2 = 'CTDTMP2', reft_col = 'T90',
                          p_col = 'CTDPRS'):
    
    d_1 = 'd_t1' #Difference between ref and prim sensor
    d_2 = 'd_t2' #Difference between ref and second sensor
    d_12 = 'd_t1_t2' #Difference between prim and sec sensor
    
    # Calculate absolute differences between sensors and reference thermom
    
    df['d_t1'] = reft_data[reft_col] - df[t_col_1]
    #df['d_t1'] = df['d_t1'].abs()
    df['d_t2'] = reft_data[reft_col] - df[t_col_2]
    #df['d_t2'] = df['d_t2'].abs()
    df['d_t1_t2'] = df[t_col_1] - df[t_col_2]
    #df['d_t1_t2'] = df['d_t1_t2'].abs()
    
    #split dataframes by pressure ranges
    
    #Greater than 2000 dBar
    lower_lim = 2000
    upper_lim = df[p_col].max()
    threshold = 0.002
    
    df_deep_good = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold)
    df_deep_ques = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold,find='quest')
    
    #Between 2000 and 1000
    lower_lim = 1000
    upper_lim = 2000
    threshold = 0.005
    
    df_lmid_good = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold)
    df_lmid_ques = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold,find='quest')
    
    #Between 1000 and 500
    lower_lim = 500
    upper_lim = 1000
    threshold = 0.010
    
    df_umid_good = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold)
    df_umid_ques = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold,find='quest')
    
    #Less than 500
    lower_lim = df[p_col].min()
    upper_lim = 500
    threshold = 0.020

    df_shal_good = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold)
    df_shal_ques = quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold,find='quest')
    
    #concat dataframes into two main dfs
    df_good = pd.concat([df_deep_good,df_lmid_good,df_umid_good,df_shal_good])
    df_ques = pd.concat([df_deep_ques,df_lmid_ques,df_umid_ques,df_shal_ques])
    
    x0 = int(xRange.split(":")[0])
    x1 = int(xRange.split(":")[1])
    
    #report questionable data to a csv file
    
    #constrain dataframes to within limits of xRange
    
    if xRange != None:
        x0 = int(xRange.split(":")[0])
        x1 = int(xRange.split(":")[1])
        
        df_good_cons = df_good[(df_good[p_col] >= x0) & (df_good[p_col] <= x1)]
     
        #Add here is planning on using for other calibrate code
#    else: 
#        if order == 1:
#            
#            x0 = 
#            x1 = 
#        
#        elif:
#            
#            x0 = 
#            x1 = 
#            
    else:
        print('Invalid xRange')
    
    # Determine fitting ranges
    
    fit = np.arange(x0,x1,(x1-x0)/50)
    
    cf1 = np.polyfit(df_good_cons[p_col], df_good_cons[d_1], order)
    cf2 = np.polyfit(df_good_cons[p_col], df_good_cons[d_2], order)
    
   
    sensor = '_t'+str(sensor)
    coef1 = np.zeros(shape=5)
    coef2 = np.zeros(shape=5)
    
    if order is 0:
        coef1[4] = cf1[0]
        
        coef2[4] = cf2[0]
        
    elif (order is 1) and (calib_param == 'P'):
        coef1[1] = cf1[0]
        coef1[4] = cf1[1]
        
        coef2[1] = cf2[0]
        coef2[4] = cf2[1]
        
    elif (order is 2) and (calib_param == 'P'):
        coef1[0] = cf1[0]
        coef1[1] = cf1[1]
        coef1[4] = cf1[2]
        
        coef2[0] = cf2[0]
        coef2[1] = cf2[1]
        coef2[4] = cf2[2]
    elif (order is 1) and (calib_param == 'T'):
        coef1[3] = cf1[0]
        coef1[4] = cf1[1]
        
        coef2[3] = cf2[0]
        coef2[4] = cf2[1]
    elif (order is 2) and (calib_param == 'T'):
        coef1[2] = cf1[0]
        coef1[3] = cf1[1]
        coef1[4] = cf1[2]
    
        coef2[2] = cf2[0]
        coef2[3] = cf2[1]
        coef2[4] = cf2[2]
        
#    Y = fit_ctd.conductivity_polyfit(coef, fit, fit, np.full(len(fit), 0.0))
#
#
#    fitfile = str('fitting'+sensor+'.' + FILE_EXT)
#    fitfilePath = os.path.join(log_directory, fitfile)
#    report_ctd.report_polyfit(coef, file_base_arr, fitfilePath)
        
    return df
    
def quality_check(df,d_1,d_2,d_12,lower_lim,upper_lim,threshold,find='good',col_name = 'CTDPRS'):
    
    #Choose Data range to compare with
    df_range = df[(df[col_name] > lower_lim) & (df[col_name] <= upper_lim)]
    
    
    if find == 'good':
    # Find data values for each sensor that are below the threshold (good)
        df_range_comp_1 = df_range[df_range[d_1].abs() < threshold]
        df_range_comp_2 = df_range[df_range[d_2].abs() < threshold]
        df_range_comp_3 = df_range[df_range[d_12].abs() < threshold]
    
    elif find == 'quest':
    # Find data values for each sensor that are above the threshold (questionable)
        df_range_comp_1 = df_range[df_range[d_1].abs() > threshold]
        df_range_comp_2 = df_range[df_range[d_2].abs() > threshold]
        df_range_comp_3 = df_range[df_range[d_12].abs() > threshold]
   
    else:
        print('Find argument not valid, please enter "good" or "quest" to find good or questionable values')
    
    #concatenate dataframe to merge all values together
    df_concat = pd.concat([df_range_comp_1,df_range_comp_2,df_range_comp_3])
        
    # Remove duplicate values
    df_concat = df_concat.drop_duplicates(subset=[col_name],keep='first')
    
    return df_concat
    

    #Combine these three into a dataframe and write out to a csv 
    #Sort by sta/cast, bottle number, rev. press
    
    
###End try/except fix

### OLD UNUSED
# def treat_surface_data(p_col,sample_rate,inMat):
#     """treat surface data function
#
#     Function takes full NUMPY ndarray with predefined dtype array
#     and several arguments to treat missing sirface bin data. It
#     basically takes the first valid data row and repeats that row
#     with the sample rate defined used in roll filter and constructs
#     an interval based on a normal surface decent rate back to the surface.
#
#     implement the following as fail safe surface treatment,
#     Args:
#         param1 (str): p_col,
#         param1 (int): sample_rate,
#         param2 (int):
#         param3 (array):
#         param4 (ndarray):
#
#     Returns:
#         Narray: The return value is a matrix of pressure sequenced data
#     """
#
#     fl = 24
#     fps = fl / sample_rate # Number of frames per second
#     dr = 2 # dbar/sec
#     fpdb = fps / dr # Frames per dbar
#     start_pressure = inMat[p_col][0] # Start p
#     fn = math.ceil(start_pressure * fpdb)
#     dp = start_pressure / fn  # Delta pressure
#
#     surface_vals = np.empty(shape=(fn,), dtype=inMat.dtype)
#
#     for i in range(0,fn):
#         surface_vals[i] = inMat[0]
#         surface_vals[p_col][i] = surface_vals[p_col][i] - (start_pressure - i*dp)
#
#     inMat = np.concatenate((surface_vals, inMat), axis=0)
#     return inMat

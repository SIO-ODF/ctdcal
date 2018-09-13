#!/usr/bin/env python
from scipy import interpolate
import numpy as np
from numpy.lib.recfunctions import append_fields
import scipy.signal as sig
import scipy.stats as st
import time, os
import pandas as pd
import math
import ctdcal.report_ctd as report_ctd
import warnings
import ctdcal.fit_ctd as fit_ctd
import datetime

import sys
sys.path.append('ctdcal/')
import oxy_fitting
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
        
        if len(inMat) % 2 == 0:
            inMat_2 = inMat.copy()
        else:
            inMat_2 = inMat.iloc[1:].copy()
            
        inMat_half1, inMat_half2 = np.split(inMat_2,2)
        ep = inMat_half2[(inMat_half2[c1_col] < conductivity_startup) & (inMat_half2[c2_col] < conductivity_startup)][p_col]
            
#        for j in range(int(len(inMat)*0.5), len(inMat)):
#            if ((inMat[c1_col][j] < conductivity_startup) and (inMat[c2_col][j] < conductivity_startup)):
#                ep.append(inMat[p_col][j])
#                if (tmp > j): tmp = j

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


def roll_filter(df, p_col='CTDPRS', up='down', frames_per_sec=24, search_time=15, start=0):
    """roll_filter function

    Function takes full NUMPY ndarray with predefined dtype array
    and subsample arguments to return a roll filtered ndarray.

    
    Jackson Notes: calculate sampling frequency 
        Two types of filters: differential filter, alias filt
        Get index of second diff function and do discrete analysis

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
    #start = kwargs.get("start", 0)
    if df is None:
        print("Roll filter function: No input data.")
        return
    
    end = df[p_col].idxmax()
    full_matrix = df
    tmp_df = df[start:end]
    tmp_df = _roll_filter(tmp_df)

    return tmp_df

#    remove = []
#    frequency = 24 # Hz of package
#
#    if (frames_per_sec > 0) & (frames_per_sec <= 24):
#        sample = int(frequency/frames_per_sec) # establish subsample rate to time ratio
#    else: sample = frequency
#
#    # Adjusted search time with subsample rate
#    search_time = int(sample*frequency*int(search_time))

#    else:
#        P = df[p_col]
#        dP = P.diff()
#        
#        if up is 'down':
#            #index_to_remove = np.where(dP < 0)[0] # Differential filter Use DIff command
#            #subMat = np.delete(df, index_to_remove, axis=0)# use dataframe boolean
#            
#            P = P[(dP>0) | (dP.isna()==True)]#Remove pressure value increases and recover first element with or
#            P = P.reset_index(drop=True)
#            dP2 = P.diff()
#            P2 = P[(dP2<0)]# Questionable data points
#            indicies = P2.index
#
#            tmp = np.array([])
#            for i in range(0,len(P)-1):#Lose If Statement
#               if P[i] > P[i+1]:# Use another diff command to find indicies
#                   deltaP = P[i+1] + abs(P[i] - P[i+1])
#                   # Remove aliasing
#                   k = np.where(P == min(P[i+1:i+search_time], key=lambda x:abs(x-deltaP)))[0]
#                   tmp = np.arange(i+1,k[0]+1,1)
#               remove = np.append(remove,tmp)
#               deltaP = 0
#        elif up is 'up':
#            index_to_remove = np.where(dP > 0)[0] # Differential filter
#            subMat = np.delete(inMat, index_to_remove, axis=0)
#
#            P = subMat[p_col]
#            tmp = np.array([])
#            for i in range(0,len(P)-1):
#               if P[i] < P[i+1]:
#                   deltaP = P[i+1] - abs(P[i] - P[i+1])
#                   # Remove aliasing
#                   k = np.where(P == min(P[i+1:i+search_time], key=lambda x:abs(x-deltaP)))[0]
#                   tmp = np.arange(i+1,k[0]+1,1)
#               remove = np.append(remove,tmp)
#               deltaP = 0
#
#        subMat = np.delete(subMat,remove,axis=0)
#
#    return subMat

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
        df_merged = pd.merge(df_surface.astype('float64'), df, on='CTDPRS', how='outer')

    return df_merged.fillna(method='bfill')

def load_reft_data(reft_file,index_name = 'btl_fire_num'):
    """ Loads reft_file to dataframe and reindexes to match bottle data dataframe"""
    
    reft_data = pd.read_csv(reft_file,usecols=['btl_fire_num','T90'])
    reft_data.set_index(index_name)
    
    reft_data['SSSCC_TEMP'] = reft_file[-14:-9]
    
    return reft_data

    
    
def load_btl_data(btl_file,cols=None):
    
    """ex. '/Users/k3jackson/p06e/data/bottle/00201_btl_mean.pkl'"""
    
    btl_data = dataToNDarray(btl_file,float,True,',',0) 
    
    btl_data = pd.DataFrame.from_records(btl_data)
    if cols != None:
        btl_data = btl_data[cols]
    
    ssscc = btl_file[-18:-13]
    
    btl_data['SSSCC'] = ssscc
    
    return btl_data


def load_time_data(time_file):
    
    time_data = dataToNDarray(time_file,float,True,',',1)
    time_data = pd.DataFrame.from_records(time_data)
    
    return time_data


def calibrate_param(param,ref_param,press,calib,order,ssscc,btl_num,xRange=None,):
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
        
    if 'P' in calib:    
        coef = get_param_coef(df_good_cons[press.name],df_good_cons['Diff'],order,calib)
    elif 'T' or 'C' in calib:
        coef = get_param_coef(df_good_cons[param.name],df_good_cons['Diff'],order,calib)
    else:
        print('calib argument not valid, use CP TP T or C')
                
    return coef,df_ques
    
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
    
        if order is 0:
            coef[4] = cf1[0]
        
        elif (order is 1) and (calib == 'TP'):
            coef[1] = cf1[0]
            coef[4] = cf1[1]
        
        elif (order is 2) and (calib == 'TP'):
            coef[0] = cf1[0]
            coef[1] = cf1[1]
            coef[4] = cf1[2]
        
        elif (order is 1) and (calib == 'T'):
            coef[3] = cf1[0]
            coef[4] = cf1[1]
        
        elif (order is 2) and (calib == 'T'):
            coef[2] = cf1[0]
            coef[3] = cf1[1]
            coef[4] = cf1[2]
            
    if 'C' in calib:
        coef = np.zeros(shape=7)
        if order is 0:
            coef[6] = cf1[0]
        elif (order is 1) and (calib == 'CP'):
            coef[1] = cf1[0]
            coef[6] = cf1[1]
        elif (order is 2) and (calib == 'CP'):
            coef[0] = cf1[0]
            coef[1] = cf1[1]
            coef[6] = cf1[2]
        elif (order is 1) and (calib == 'C'):
            coef[5] = cf1[0]
            coef[6] = cf1[1]
        elif (order is 2) and (calib == 'C'):
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
    
    if order is 0:
        coef[6] = cf[0]
    elif (order is 1) and (calib_param == 'P'):
        coef[1] = cf[0]
        coef[6] = cf[1]
    elif (order is 2) and (calib_param == 'P'):
        coef[0] = cf[0]
        coef[1] = cf[1]
        coef[6] = cf[2]
    elif (order is 1) and (calib_param == 'T'):
        coef[3] = cf[0]
        coef[6] = cf[1]
    elif (order is 2) and (calib_param == 'T'):
        coef[2] = cf[0]
        coef[3] = cf[1]
        coef[6] = cf[2]
    elif (order is 1) and (calib_param == 'C'):
        coef[5] = cf[0]
        coef[6] = cf[1]
    elif (order is 2) and (calib_param == 'C'):
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
              

def get_pressure_offset(start_vals,end_vals):
    """
    Finds unique values and calclates mean for pressure offset
            
    Parameters
    ----------
    
    start_vals :array_like
                Array containing initial ondeck pressure values
           
    end_vals :array_like
              Array containing ending ondeck pressure values
    Returns
    -------
    
    p_off :float
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

def load_pressure_logs(file):
    """
        Loads pressure offset file from logs.
        
    Parameters
    ----------
    
    file : string
           Path to ondeck_pressure log
         
    Returns
    -------
    
    df : DataFrame
         Pandas DataFrame containing ondeck start and end pressure values
         
    """
    
    df = pd.read_csv(file,names=['SSSCC','ondeck_start_p','ondeck_end_p'])
    
    # Change vaules in each row by removing non-number parts
    for i in range(len(df['SSSCC'])):
        df['SSSCC'][i] = int(df['SSSCC'].loc[i][-5:])
        df['ondeck_start_p'][i] = float(df['ondeck_start_p'].loc[i][16:])
        df['ondeck_end_p'][i] = float(df['ondeck_end_p'].loc[i][14:])
        
    return df

def write_offset_file(df,p_off,write_file='data/logs/poffset_test.csv'):
    """
    
    """
    df_out = pd.DataFrame()
    df_out['SSSCC'] = df['SSSCC']
    df_out['offset'] = p_off
    
    df_out.to_csv(write_file,index=False)
    
    return
    
def pressure_calibrate(file):
    
    pressure_log = load_pressure_logs(file)
    p_off = get_pressure_offset(pressure_log)
    
    return p_off

    
def load_all_ctd_files(ssscc,prefix,postfix,series,cols,reft_prefix='data/reft/',reft_postfix='_reft.csv',
                       refc_prefix='data/salt/',refc_postfix='',press_file='data/logs/ondeck_pressure.csv',
                       oxy_prefix='data/oxygen/', oxy_postfix='',index_col='btl_fire_num',t_col='CTDTMP1',
                       p_col='CTDPRS',ssscc_col='SSSCC'):
    """
    LOAD ALL CTD FILES was changed (commented out)
    Lines 1324-1328,1335,1337, 1338,345
    """
    df_data_all = pd.DataFrame()
   
    if series == 'bottle':
        for x in ssscc:
            print('Loading BTL data for station: ' + x + '...')
            btl_file = prefix + x + postfix
            btl_data = load_btl_data(btl_file,cols)
            
            
            reft_file = reft_prefix + x + reft_postfix
            try:
                reft_data = load_reft_data(reft_file)
            except FileNotFoundError:
                print('Missing (or misnamed) REFT Data Station: ' + x + '...filling with NaNs')
                reft_data = pd.DataFrame()
                reft_data[index_col] = pd.Series(btl_data[index_col].values.astype(int))
                reft_data['T90'] = pd.Series([np.nan]*len(btl_data))
                ref_ssscc = ssscc_col + '_TEMP'
                reft_data[ref_ssscc] = x
                reft_data.index = btl_data.index
            
            refc_file = refc_prefix + x + refc_postfix
            try:
                refc_data = fit_ctd.salt_calc(refc_file,index_col,t_col,p_col,btl_data)

            except FileNotFoundError:
                print('Missing (or misnamed) REFC Data Station: ' + x + '...filling with NaNs')
                refc_data = pd.DataFrame()
                refc_data['SAMPNO'] = pd.Series(btl_data[index_col].values.astype(int))
                refc_data['CRavg'] = pd.Series([np.nan]*len(btl_data))
                refc_data['BathTEMP'] = pd.Series([np.nan]*len(btl_data))
                refc_data['BTLCOND'] = pd.Series([np.nan]*len(btl_data))
                refc_data.index = btl_data.index
                
            
            #Fix Index for each parameter to bottle number
            
#            btl_data[index_col] = btl_data[index_col].astype(int)
#            btl_data=btl_data.set_index(btl_data[index_col].values)
#            
#            reft_data = reft_data.set_index(reft_data[index_col].values)
            
            oxy_file = oxy_prefix + x + oxy_postfix
            oxy_data,params = oxy_fitting.oxy_loader(oxy_file)
            
#            #Horizontally concat DFs to have all data in one DF
#            btl_data_full = pd.concat([btl_data,reft_data,refc_data,oxy_data],axis=1)
            
            btl_data = pd.merge(btl_data,reft_data,on='btl_fire_num',how='outer')
            btl_data = pd.merge(btl_data,refc_data,left_on='btl_fire_num',right_on='SAMPNO',how='outer')
            btl_data = pd.merge(btl_data,oxy_data,left_on='btl_fire_num',right_on='BOTTLENO_OXY',how='outer')
            
            
#            #Drop columns that have no CTD data
#            btl_data_full = btl_data_full.dropna(subset=cols)
            #btl_data = btl_data.set_index(['SSSCC','GPSLAT','GPSLON','CTDPRS'],drop=True)
            try:
                df_data_all = pd.concat([df_data_all,btl_data],sort=False)
            except AssertionError:
                raise AssertionError('Colums of ' + x + ' do not match those of previous columns')
            print('* Finished BTL data station: ' + x + ' *')
            
        #Drops duplicated columns generated by concatenation
        df_data_all = df_data_all.loc[:,~df_data_all.columns.duplicated()]
    elif series == 'time':
        for x in ssscc:
            
            print('Loading TIME data for station: ' + x + '...')
            file = prefix + x + postfix
            time_data = load_time_data(file)
            time_data['SSSCC'] = str(x)
            df_data_all = pd.concat([df_data_all,time_data], sort=False)
            print('** Finished TIME data station: ' + x + ' **')
   
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

def merge_cond_flags(btl_data, qual_flag_cond,sensor=1):
    # Merge df
    if sensor == 1:
        parameter = 'CTDCOND1'
    elif sensor == 2:
        parameter = 'CTDCOND2'
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
        btl_data[btl_data['REFTMP'].isna()]['REFTMP_FLAG_W'] = 9
    btl_data['REFTMP_FLAG_W'] = btl_data['REFTMP_FLAG_W'].astype(int)
    
    return btl_data

def get_btl_time(df,btl_num_col,time_col):
    # Get time for first btl fire
    time = df[df[btl_num_col]== df[btl_num_col].min()][time_col].values
    ts = pd.to_datetime(time,unit='s')
    date = ts.strftime('%Y%m%d')
    hour= ts.strftime('%H%M')
    df['DATE'] = date[0]
    df['TIME'] = hour[0]
    
    return df
    
    

def export_time_data(df,ssscc,sample_rate,search_time,expocode,section_id,ctd,p_column_names,p_column_units,
                     t_sensor=1,out_dir='data/pressure/',p_col='CTDPRS',stacst_col='SSSCC',
                     logFile='data/logs/cast_details.csv'):
    """ Export Time data to pressure directory as well as adding qual_flags and 
    removing unneeded columns"""
     
#    pressure_seq_data = pressure_sequence(df,p_col,2.0,-1.0,0.0,'down',sample_rate,search_time)
    
    df[stacst_col] = df[stacst_col].astype(int)
    df[stacst_col] = df[stacst_col].astype(str)
    
    # Choose Sensors
    if t_sensor == 1:
        df['CTDTMP'] = df['CTDTMP1']
    elif t_sensor ==2:
       df['CTDTMP'] = df['CTDTMP1']
    
    # Change column names
    
    df['CTDFLUOR'] = df['FLUOR']
    #df['CTDOXY'] = np.NaN
    df['CTDRINKO'] = df['FREE1']
    
    # Add Flagged colummns
    
    df['CTDPRS_FLAG_W'] = 2
    df['CTDTMP_FLAG_W'] = 2
    df['CTDSAL_FLAG_W'] = 2
    df['CTDOXY_FLAG_W'] = 2
    df['CTDXMISS_FLAG_W'] = 1
    df['CTDFLUOR_FLAG_W'] = 1
    df['CTDBACKSCATTER_FLAG_W'] = 1
    df['CTDRINKO_FLAG_W'] = 1
    
    # Round to 4 decimal places
    
    df = df.round(4)
    
    #Remove unwanted columns
    
#    pressure_seq_data = pressure_seq_data[['CTDPRS','CTDPRS_FLAG_W','CTDTMP','CTDTMP_FLAG_W','CTDSAL','CTDSAL_FLAG_W','CTDOXY','CTDOXY_FLAG_W',
#                                           'CTDXMISS','CTDXMISS_FLAG_W','CTDFLUOR','CTDFLUOR_FLAG_W','CTDBACKSCATTER','CTDBACKSCATTER_FLAG_W',
#                                           'CTDRINKO','CTDRINKO_FLAG_W']]

    cast_details = dataToNDarray(logFile,str,None,',',0)
    
    for cast in ssscc:
        #time_data = pressure_seq_data.copy()
        #time_data = time_data[pressure_seq_data['SSSCC'] == cast]
        time_data = df[df['SSSCC'] == cast]
        time_data = pressure_sequence(time_data,p_col,2.0,-1.0,0.0,'down',sample_rate,search_time)
        time_data = time_data[['CTDPRS','CTDPRS_FLAG_W','CTDTMP','CTDTMP_FLAG_W','CTDSAL','CTDSAL_FLAG_W','CTDOXY','CTDOXY_FLAG_W',
                               'CTDXMISS','CTDXMISS_FLAG_W','CTDFLUOR','CTDFLUOR_FLAG_W','CTDBACKSCATTER','CTDBACKSCATTER_FLAG_W',
                               'CTDRINKO','CTDRINKO_FLAG_W']]
        
        time_data=time_data.round(4)
        
        
        s_num = cast[-5:-2]
        c_num = cast[-2:]
#        line_ID = 'stacast:'+sss
#        print(line_ID)
        for line in cast_details:
            if cast in line[0]:
                for val in line:
                    if 'at_depth' in val: btime = float(str.split(val, ':')[1])
                    if 'latitude' in val: btm_lat = float(str.split(val, ':')[1])
                    if 'longitude' in val: btm_lon = float(str.split(val, ':')[1])
                    if 'altimeter_bottom' in val: btm_alt = float(str.split(val, ':')[1])
                break
        
        bdt = datetime.datetime.fromtimestamp(btime).strftime('%Y%m%d %H%M').split(" ")
        b_date = bdt[0]
        b_time = bdt[1]
        depth = -99
        now = datetime.datetime.now()
        file_datetime = now.strftime("%Y%m%d") #%H:%M")
        file_datetime = file_datetime + 'ODFSIO'
        outfile = open(out_dir+cast+'_ct1.csv', "w+")
        outfile.write("CTD, %s\nNUMBER_HEADERS = %s \nEXPOCODE = %s \nSECT_ID = %s\nSTNNBR = %s\nCASTNO = %s\n DATE = %s\nTIME = %s\nLATITUDE = %f\nLONGITUDE = %f\nINSTRUMENT_ID = %s\n" % (file_datetime, 10, expocode, section_id, s_num, c_num[1], b_date, b_time, btm_lat, btm_lon, ctd))
        cn = np.asarray(p_column_names)
        cn.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')
        cu = np.asarray(p_column_units)
        cu.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')
        outfile.close()

        file = out_dir+cast+'_ct1.csv'
        with open(file,'a') as f:
            time_data.to_csv(f, header=False,index=False)
        f.close()
        
        outfile = open(out_dir+cast+'_ct1.csv', "a")
        outfile.write('END_DATA')
        outfile.close()
        
def export_btl_data(df,expocode,sectionID,cruise_line,out_dir='data/pressure/',t_sensor=1,org='ODF'):
    
    btl_columns = ['EXPOCODE','SECT_ID','STNNBR','CASTNO','SAMPNO','BTLNBR','BTLNBR_FLAG_W','DATE','TIME','LATITUDE','LONGITUDE','DEPTH',
                   'CTDPRS','CTDTMP','REFTMP','REFTMP_FLAG_W','CTDSAL','CTDSAL_FLAG_W','SALNTY','SALNTY_FLAG_W',
                   'CTDOXY','CTDOXY_FLAG_W','OXYGEN','OXYGEN_FLAG_W']
    btl_units = ['','','','','','','','','','','','METERS','DBARS','ITS-90','C','','PSS-78','','PSS-78','','UMOL/KG','','UMOL/KG','']
    
    btl_data = df.copy()
    
    if t_sensor ==1:
        btl_data['CTDTMP'] = btl_data['CTDTMP1']
    elif t_sensor ==2:
        btl_data['CTDTMP'] = btl_data['CTDTMP2']
        
    
    btl_data['EXPOCODE'] = expocode
    btl_data['SECT_ID'] = sectionID
    btl_data['STNNBR'] = btl_data['SSSCC'].str[0:3]
    btl_data['CASTNO'] = btl_data['SSSCC'].str[3:5]
    btl_data['SAMPNO'] = btl_data['btl_fire_num']
    btl_data['BTLNBR'] = btl_data['btl_fire_num']
    btl_data['BTLNBR_FLAG_W'] = '2'
#    btl_data['DATE'] = np.NaN
#    btl_data['TIME'] = np.NaN
    btl_data['LATITUDE'] = btl_data['GPSLAT']
    btl_data['LONGITUDE'] = btl_data['GPSLON']
    btl_data['DEPTH'] = np.NaN
    btl_data['REFTMP'] = btl_data['T90']
#    btl_data['REFTMP_FLAG'] = '2'
#    btl_data['CTDSAL_FLAG'] = '2_HARDCODE'
    btl_data['SALNTY'] = gsw.SP_from_C(btl_data['BTLCOND'],btl_data['CTDTMP'],btl_data['CTDPRS'])
#    btl_data['SALNTY_FLAG'] = '2_HARDCODE'
#    btl_data['CTDOXY'] = np.NaN
#    btl_data['CTDOXY_FLAG'] = '2'
#    btl_data['OXYGEN'] = np.NaN
    btl_data['OXYGEN_FLAG_W'] = '2'
    
    
    
    btl_data = btl_data[btl_columns]
    btl_data = btl_data.round(4)
    
    btl_data = btl_data.fillna(value=-999)
    btl_data = btl_data.round(4)    
    
    now = datetime.datetime.now()
    file_datetime = now.strftime("%Y%m%d")
    
    time_stamp = file_datetime+org
    
    outfile = open(out_dir+cruise_line+'_hy1.csv', "w+")
    outfile.write("BOTTLE, %s\n" % (time_stamp))
    cn = np.asarray(btl_columns)
    cn.tofile(outfile,sep=',', format='%s')
    outfile.write('\n')
    cu = np.asarray(btl_units)
    cu.tofile(outfile,sep=',', format='%s')
    outfile.write('\n')
    outfile.close()
    
    file = out_dir+cruise_line+'_hy1.csv'
    with open(file,'a') as f:
        btl_data.to_csv(f, header=False,index=False)
    f.close()
    
    outfile = open(out_dir+cruise_line+'_hy1.csv', "a")    
#    outfile.write('\n')
    outfile.write('END_DATA')
    outfile.close()
    
    return
        
        
   
    

    
    
    
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

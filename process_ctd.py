#!/usr/bin/env python
# Python script with tested CTD signal filters
# Ideal default CTD filter is the Triangle filter
# Other filters may be added if they are theoretically sound and tested first
# 
# EXAMPLE:
# import filters  
# data_array = []
# filters.traingle(data_array, 24)
# ** This means triangle filter with a 12 frame window size, filtered twice with over passed array
#  
import numpy as np
import scipy.signal as sig
import scipy.stats as st
import matplotlib.pyplot as plt
import pandas as pd


# This function manages the start time of cast post standard 10m start-up soak and return to surface.
def cast_details(inMat):
    """cast_details function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and adjusts ndarray to remove all extraneous surface data. 
    Function returns cast start time, end time, bottom time and 
    cleaned up matrix.

    Args:
        param1 (ndarray): inMat, numpy ndarray with dtype array 

    Returns:
        Narray: The return value is ndarray with adjusted time of parameter 
          specified. 

    """
    # Top of cast time, bottom of cast time, end of cast time, 
    s = 0.0
    b = 0.0
    e = 0.0
    # Test cycle time constant
    fl = 24
    # starting P 
    sp = 2.0
    # Max P 
    mp = 10000.0
    lm = len(inMat)
    rev = np.arange(int(lm/4),0,-1)
    
    # Find starting top of cast
    # Smallest P from reverse array search 
    for i in rev: 
        if sp > inMat['Pdbar'][i]:
           sp = inMat['Pdbar'][i]
           tmp = i
           break

    s = inMat['TIMEs'][tmp]

    # Remove everything before cast start
    inMat = inMat[tmp:]
  
    # Max P and bottom time
    mp = max(inMat['Pdbar'])    
    tmp = np.argmax((inMat['Pdbar'])) 
    b = inMat['TIMEs'][tmp]

    tmp = len(inMat)
    # Find ending top of cast time
    for i in range(int(lm/2),lm):
        if sp > inMat['Pdbar'][i]:
            e = inMat['TIMEs'][i]
            tmp = i
            break
     
    # Remove everything after cast end
    inMat = inMat[:tmp]
  
    print('start time '+str(s))
    print('start pres '+str(sp))
    print('bottom time '+str(b))
    print('max press '+str(mp))
    print('end time' +str(e))
    print(inMat)

    return s, e, b, sp, mp, inMat

def ctd_align(inMat, col, time=0.0):
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
    
    # Time to advance 
    advnc = int(fl * time)
    tmp = np.arange(advnc, dtype=np.float)
    last = inMat[col][len(inMat)-1]
    tmp.fill(float(last))
    inMat[col] = np.concatenate((inMat[col][advnc:],tmp))

    return inMat

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
    df = pd.read_csv(inFile, header=[0,1])
    return df

# Read in csv file to data object
# dtype is defined below
# https://scipy.github.io/old-wiki/pages/Cookbook/InputOutput.html
def dataToMatrix(inFile, dtype=None, separator=','):
    """dataToMatrix function 

    Function takes full file path to csv type data file and returns NUMPY
    ndarray type matrix for data manipulation with a two row header. 

    Data file should have a two row header. The first row being the column 
    title and the second row are the units for each column.
    Args:
        param1 (str): inFile, full path to csv file 
        param2 (arr): dtype list 
        param3 (str): separator, default comma ','  

    Returns:
        Narray: The return value is a full data matrix with two row header.

    Reference Page:
        https://scipy.github.io/old-wiki/pages/Cookbook/InputOutput.html
    print(dtype)
    if dtype is None: 
        return 
    """

    arr = np.genfromtxt(inFile, delimiter=separator, dtype=dtype, skip_header=2)
    
    return arr 

def raw_ctd_filter(arr = None, filter_type='triangle',win_size=24):
    """raw_ctd_filter function 

    Function takes NUMPY array 
    of raw ctd data and returns filtered data. This function also needs 
    one of three filter types (boxcar, gaussian, triangle) as well as 
    window size. 

    Args:
        param1 (ndarray): Numpy ndarray with predefined header with at 
          "TIMEs, Pdbar, T1C, T2C, C1mScm, C2mScm, SBE43FV, V1V ..., V7V" 
        param2 (str): One of three tested filter types
          boxcar, gaussian_std, triangle.
          default is triangle
        param3 (int): A window size for the filter. Default is 24, which 
          is the number of frames per second from a SBE9+/11 CTD/Dech unit.

    Returns:
        Narray: The return value is a matrix of filtered ctd data with 
          the above listed header values.

    """
    
    if arr is None:
        return
    else:  
        if filter_type is 'boxcar':
            win = sig.boxcar(win_size)
            rtn = sig.convolve(arr, win, mode='same')/len(win)
        elif filter_type is 'gaussian':
            sigma = np.std(arr)
            win = sig.general_gaussian(win_size, 1.0, sigma)
            rtn = sig.convolve(arr, win, mode='same')/(len(win))
        elif filter_type is 'triangle':
            win = sig.triang(win_size)
            rtn = 2*sig.convolve(arr, win, mode='same')/len(win)
    return rtn 

def ondeck_pressure(inMat, scond):
    """ondeck_pressure function 

    Function takes full NUMPY ndarray with predefined dtype array 
    of filtered ctd raw data the stores, analizes and removes ondeck 
    values from data. 

    Args:
        param1 (ndarray): numpy ndarray with dtype array 
        param2 (float): in water startup conductivity value 

    Returns:
        Narray: The return value is a matrix of filtered ctd data with 
          the above listed header values.

    """
    sp = []
    tmpMat = []
    outMat = []
    tmp = 0
    start_p = 0.0
    n = 0 
    ep = []
    end_p = 0.0 

    # Frame Length
    fl = 24 
    fl2 = fl*2 
    # One minute
    mt = 60
    # Half minute 
    ms = 30
    sdelay = fl*ms

    # Searches first quarter of matrix, uses start conductivity  
    # condition min to capture startup Press
    for j in range(0,int(len(inMat)/4)):
        if ((inMat['TIMEs'][j] > 0.0) and (inMat['C1mScm'][j] < scond) and (inMat['C2mScm'][j] < scond)):
            tmp = j
            sp.append(inMat['Pdbar'][j])
    
    # Evaluate starting pressures
    if not sp: start_p = "Started in Water"
    else: 
        n = len(sp)
        if (n > sdelay): start_p = np.average(sp[fl2:n-(sdelay)])
        else: start_p = np.average(sp[fl2:n])

    # Remove on-deck startup
    inMat = inMat[tmp:]

    tmp = len(inMat); 
    # Searches last half of Matrix for conductivity threshold 
    for j in range(int(len(inMat)*0.5), len(inMat)):
        if ((inMat['C1mScm'][j] < scond) and (inMat['C2mScm'][j] < scond)):
            ep.append(inMat['Pdbar'][j])
            if (tmp > j): tmp = j 

    # Evaluate ending pressures
    if (len(ep) > (sdelay)): end_p = np.average(ep[(sdelay):])
    else: end_p = np.average(ep[(len(ep)):])

    # Remove on-deck ending 
    inMat = inMat[:tmp]

    # Store ending on-deck pressure
    print("Sta/Cast ondeck start "+str(start_p)+" "+str(end_p))

    return inMat

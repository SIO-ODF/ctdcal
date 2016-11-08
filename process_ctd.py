#!/usr/bin/env python
import numpy as np
import scipy.signal as sig
import scipy.stats as st
import matplotlib.pyplot as plt
import pandas as pd

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
    # condition min to capture startup pressure
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
    outMat = inMat[:tmp]

    # Store ending on-deck pressure
    print("Sta/Cast ondeck start "+str(start_p)+" "+str(end_p))

    return outMat

def roll_filter(inMat=None, up='down', frames_per_sec=24, search_time=15):
    """roll_filter function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and subsample arguments to return a roll filtered ndarray. 

    Args:
        param1 (ndarray): inMat, numpy ndarray with dtype array 
        param2 (str): up, direction to filter cast (up vs down)
        param3 (int): frames_per_sec, subsample selection rate
        param4 (int): seach_time, search time past pressure inversion  

    Returns:
        Narray: The return value ndarray of data with ship roll removed
    """
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
        P = inMat['Pdbar']
        dP = np.diff(P,1) 

        if up is 'down':
            index_to_remove = np.where(dP < 0)[0] # Differential filter
            subMat = np.delete(inMat, index_to_remove, axis=0)

            P = subMat['Pdbar']
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

            P = subMat['Pdbar']
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

def pressure_sequence(inMat=None, intP=2.0, startT=-1.0, startP=0.0, up='down', sample_rate=1, search_time=15):
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
        param1 (ndarray): input matrix (inMat), numpy ndarray with dtype array 
        param2 (float): pressure interval (intP), default 2.0 dbar 
        param3 (float): starting pressure interval 
        param4 (float): start time (startT) for pressure sequence 
        param5 (float): start pressure (startP) for pressure sequence
        param5 (str): pressure sequence direction (down/up)
        param6 (int): sample_rate, sub sample rate for roll_filter. Cleans & speeds processing.
        param7 (int): search_time, truncate search index for the aliasing part of ship roll. 

    Returns:
        Narray: The return value is a matrix of pressure sequenced data 

    todo: implement the following as fail safe surface treatment, 
          repeatPt: the point in the matrix that is repeated back to the surface.
          repeatVal: the values in the matrix that is repeated back to the surface.
          
    todo: deep data bin interpolation to manage empty slices
    """

    # Passed Time-Series, Create Pressure Series
    if inMat is None:
        print("Pressure sequence function: No input data.")
        return
    else:
        pF = inMat['Pdbar']
        full_length = len(pF)-1

        btm = max(pF) # bottom max P
        indBtm = np.argmax(pF) # bottom index
        btmTime = inMat['TIMEs'][indBtm] # bottom time 

        # Initialise input parameters
        if ((startT > 0.0) and (startT > inMat['TIMEs'][0])):
            repeatPt = (np.abs(inMat['TIMEs'] - startT)).argmin()
            repeatVal = inMat[:][repeatPt]
            lenP = np.arange(repeatPt,indBtm,1)
            start = repeatPt
            end = len(lenP)
            prvPrs = inMat['Pdbar'][repeatPt]
            if btmTime <= startT:
                print("-startT start time is greater than down cast time. Cast issue.")
                return
        elif ((startP > 0.0) and (startP > pF[0])):
            repeatPt = (np.abs(inMat['Pdbar'] - startP)).argmin()
            repeatVal = inMat[:][repeatPt]
            lenP = np.arange(repeatPt,indBtm,1)
            start = repeatPt
            end = len(lenP) 
            prvPrs = inMat['Pdbar'][repeatPt]
            if btm <= startP:
                print("-startP start pressure is greater than bottom pressure. Cast issue.")
                return
        elif up is 'up':
            repeatPt = full_length
            repeatVal = inMat[:][repeatPt]
            lenP = np.arange(indBtm,repeatPt,1) 
            start = indBtm
            end = full_length 
            prvPrs = btm
        else:
            lenP = np.arange(0,indBtm,1)
            start = 0
            end = len(lenP)
            prvPrs = 0.0
        
        # Roll Filter 
        roll_filter_matrix = roll_filter(inMat[:][start:end:sample_rate], up, sample_rate, search_time)

        # Frame Pressure Bins
        pressure_bins = np.arange(0,int(btm),2)
        p_bin_index = np.digitize(roll_filter_matrix['Pdbar'],pressure_bins)

        # Define output array
        inMat_dtype = inMat.dtype
        binned_matrix = np.empty(shape=(len(pressure_bins),), dtype=inMat_dtype)

        # todo: remove explicit column 'Pdbar' call here
        for col in binned_matrix.dtype.names:
            if col == 'Pdbar':
                binned_matrix[col] = pressure_bins
            elif binned_matrix[col].dtype is np.dtype(np.float64):
                binned_matrix[col] = [roll_filter_matrix[col][p_bin_index == i].mean() for i in range(0,len(pressure_bins))]

    return binned_matrix

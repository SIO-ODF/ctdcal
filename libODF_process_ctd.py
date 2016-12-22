#!/usr/bin/env python
from scipy import interpolate
import numpy as np
from numpy.lib.recfunctions import append_fields
import scipy.signal as sig
import scipy.stats as st
import time, os
import pandas as pd
import math
import libODF_report_ctd as report_ctd
import warnings

warnings.filterwarnings("ignore", 'Mean of empty slice.')

def cast_details(stacast, log_file, p_col, time_col, blat_col, blon_col, alt_col, inMat=None):
    """cast_details function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and adjusts ndarray to remove all extraneous surface data. 
    Function returns cast start time, end time, bottom time and 
    cleaned up matrix.

    Args:
        param1 (str): stacast, station cast input 
        param2 (str): log_file, log file to write cast data.
        param3 (str): p_col, pressure data column name 
        param4 (str): time_col, time data column name 
        param5 (ndarray): inMat, numpy ndarray with dtype array 

    Returns:
        Narray: The return value is ndarray with adjusted time of parameter 
          specified. 

    """

    if inMat is None:
       print("In cast_details: No data")
       return
    else:
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
        lm = len(inMat)-1
        rev = np.arange(int(lm/4),0,-1)
    
        # Find starting top of cast
        # Smallest P from reverse array search 
        for i in rev:
            if sp < inMat[p_col][i]:
               tmp = i
            elif sp > inMat[p_col][i]:
               sp = inMat[p_col][i]
               tmp = i - 24
               break

        s = inMat[time_col][tmp]

        # Remove everything before cast start
        inMat = inMat[tmp:]
  
        # Max P and bottom time
        mp = max(inMat[p_col])    
        tmp = np.argmax((inMat[p_col])) 
        b = inMat[time_col][tmp]
        b_lat = inMat[blat_col][tmp]
        b_lon = inMat[blon_col][tmp]
        b_a = inMat[alt_col][tmp]

        tmp = len(inMat)
        # Find ending top of cast time
        for i in range(int(lm/2),lm):
            if sp > inMat[p_col][i]:
                e = inMat[time_col][i]
                tmp = i + 24
                break
     
        # Remove everything after cast end
        inMat = inMat[:tmp]
  
    report_ctd.report_cast_details(stacast, log_file, s, e, b, sp, mp, b_a, b_lat, b_lon)
    
    return s, e, b, sp, mp, b_lat, b_lon, b_a, inMat

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
    #print(p_qual_col)
    #If p_range set apply qual codes to part of array and return
    if p_range is not None:
        print("Some algoirythm for formatting qual codes per pressure range")
        return
    else: 
        q_df = pd.DataFrame(index=np.arange(len(inMat)), columns=p_qual_col)
        for pq in p_qual_col:
            if pq in list(qual_one):
                q_df[pq] = q_df[pq].fillna(1)
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

def dataToNDarray(inFile, dtype=None, names=None, separator=',', ):
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

    arr = np.genfromtxt(inFile, delimiter=separator, dtype=dtype, names=names, skip_header=2)
    
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
                #indices = [i for i, s in input_array.dtype.names if p in s]
                #print(indices)
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
    sp = []
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
                sp.append(inMat[p_col][j])
    
        # Evaluate starting pressures
        if not sp: start_p = "Started in Water"
        else: 
            n = len(sp)
            if (n > time_delay): start_p = np.average(sp[fl2:n-(time_delay)])
            else: start_p = np.average(sp[fl2:n])

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

def roll_filter(p_col, inMat=None, up='down', frames_per_sec=24, search_time=15):
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

def pressure_sequence(stacast, p_col, time_col, intP=2.0, startT=-1.0, startP=0.0, up='down', sample_rate=12, search_time=15, inMat=None):
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

    # Passed Time-Series, Create Pressure Series
    if inMat is None:
        print("Pressure sequence function: No input data.")
        return
    else:
        pF = inMat[p_col]
        full_length = len(pF)-1

        btm = max(pF) # bottom max P
        indBtm = np.argmax(pF) # bottom index
        btmTime = inMat[time_col][indBtm] # bottom time 

        # Initialise input parameters
        if ((startT > 0.0) and (startT > inMat[time_col][0])):
            start = (np.abs(inMat[time_col] - startT)).argmin()
            lenP = np.arange(start,indBtm,1)
            end = len(lenP)
            prvPrs = inMat[p_col][start]
            if btmTime <= startT:
                print("-startT start time is greater than down cast time. Cast issue.")
                return
        elif ((startP > 0.0) and (startP > pF[0])):
            start = (np.abs(inMat[p_col] - startP)).argmin()
            lenP = np.arange(start,indBtm,1)
            end = len(lenP) 
            prvPrs = inMat[p_col][start]
            if btm <= startP:
                print("-startP start pressure is greater than bottom pressure. Cast issue.")
                return
        elif up is 'up':
            start = indBtm
            end = full_length 
            lenP = np.arange(start,end,1) 
            prvPrs = btm
        else:
            lenP = np.arange(0,indBtm,1)
            start = 0
            end = len(lenP)
            prvPrs = 0.0
        
        # Roll Filter 
        roll_filter_matrix = roll_filter(p_col, inMat[:][start:end:sample_rate], up, sample_rate, search_time)

        # Treat surface data.
        roll_filter_matrix = treat_surface_data(p_col, sample_rate,roll_filter_matrix)

        # Frame Pressure Bins
        pressure_bins = np.arange(0,int(btm),intP)
        p_bin_index = np.digitize(roll_filter_matrix[p_col],pressure_bins)

        # Define binned output array
        binned_matrix = np.empty(shape=(len(pressure_bins),), dtype=inMat.dtype)

        # Iterate over input data by column, sort bins and find mean of binned data. 
        for col in binned_matrix.dtype.names:
            if col == p_col:
                binned_matrix[col] = pressure_bins
            elif binned_matrix[col].dtype is np.dtype(np.float64):
                binned_matrix[col] = [roll_filter_matrix[col][p_bin_index == i].mean() for i in range(1,len(pressure_bins)+1)]
                # Interpolate over NaN or missing data
                if np.isnan(binned_matrix[col]).any(): 
                    binned_matrix[col] = data_interpolater(binned_matrix[col])

    return binned_matrix


def treat_surface_data(p_col,sample_rate,inMat):
    """pressure_sequence function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and several arguments to treat missing sirface bin data. It
    basically takes the first valid data row and repeats that row 
    with the sample rate defined used in roll filter and constructs 
    an interval based on a normal surface decent rate back to the surface. 

    implement the following as fail safe surface treatment, 
    Args:
        param1 (str): p_col,
        param1 (int): sample_rate,
        param2 (int): 
        param3 (array): 
        param4 (ndarray): 

    Returns:
        Narray: The return value is a matrix of pressure sequenced data 
    """

    fl = 24
    fps = fl / sample_rate # Number of frames per second
    dr = 2 # dbar/sec 
    fpdb = fps / dr # Frames per dbar
    sp = inMat[p_col][0] # Start p
    fn = math.ceil(sp * fpdb)
    dp = sp / fn  # Delta pressure
    
    surface_vals = np.empty(shape=(fn,), dtype=inMat.dtype)
   
    for i in range(0,fn):
        surface_vals[i] = inMat[0] 
        surface_vals[p_col][i] = surface_vals[p_col][i] - (sp - i*dp)

    inMat = np.concatenate((surface_vals, inMat), axis=0)
    return inMat

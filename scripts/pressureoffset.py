#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 15:10:39 2018

@author: k3jackson
"""
import pandas as pd
import numpy as np

def pressure_offset(ondeck_file):
    """ Calculates average pressure offset
    Needs ondeck pressure csv file
    
    Variables - ondeck file: './ondeck_pressure.csv'
    """
    ondeck_press = pd.read_csv(ondeck_file,header=None,names=['STN','Start','End'])
    ondeck_press['p_start'] = ondeck_press['Start'].str.split(':').str.get(1)
    ondeck_press['p_start'] = pd.to_numeric(ondeck_press['p_start'],errors='coerce') 
    ondeck_press['p_end'] = ondeck_press['End'].str.split(':').str.get(1)
    ondeck_press['p_end'] = pd.to_numeric(ondeck_press['p_end'],errors='coerce') 
    
    p_off_start = ondeck_press['p_start'].unique()
    p_off_end = ondeck_press['p_end'].unique()
    
    p_offset = np.nanmean(p_off_start) - np.nanmean(p_off_end)
    #
    #Add reporting here
    return p_offset

def pressure_fit(offset,df,column='CTDPRS'):
    """
    Apply Pressure offset to dataframe
    
    """
    
    df[column] = df[column] + offset
    
    return df
    
pressure_seq_data = process_ctd.pressure_sequence(filename_base, p_col, timedate, 2.0, -1.0, 0.0, 'down', int(sample_rate), int(search_time), time_data) 
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
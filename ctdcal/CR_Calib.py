#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 06:04:34 2017

A set of functions to analyze autosal conductivity files/data

@author: Kenneth Jackson
"""
# break into two
#docstrings
# keyword argument in calibration default = worm

import csv
import numpy as np
import pandas as pd
import sys
import os

def SaltLoad(saltFile):
    
    """ Converts a autosal salinometer output file to a Pandas Dataframe.

    Input:
        - saltFile (file), an unextended file containing the output file
        from the autosal salinometer. Contains columns/values such
        as STATION NUMBER, CAST NUMBER,SAMPLENUMBER, CONDUCTIVITY RATIO, 
        etc. 
        Ex. saltFile = '/data/salt/ssscc'
        
    Output:
        - saltDF (Pandas Dataframe),Dataframe with 15 Columns containing the input data with 
        appropriate column names for the data.

    Usage:
        >>> saltDF = SaltLoad(saltFile)
    """
    
    f = open(saltFile, newline='')
    saltF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
    
    saltArray = []
    for row in saltF:
        saltArray.append(row)
    del saltArray[0]
         
    header = ['STNNBR','CASTNO','SAMPNO','BathTEMP','CRavg','autosalSAMPNO',\
              'Unknown','StartTime','EndTime','Attempts','Reading1','Reading2',\
              'Reading3', 'Reading4', 'Reading5']
    f.close()
    # make all rows of Salt files the same length as header   
    for row in saltArray:
        if len(row) < len(header):
            row.extend([np.NaN]*(len(header)-len(row)))
            
    saltArray = np.array(saltArray) # change to np array
    
    saltDF = pd.DataFrame(saltArray,columns=header) # change to DataFrame
    saltDF = saltDF.apply(pd.to_numeric, errors='ignore')
    
    return saltDF
        
def Cr_Calibration(saltDF,StandName='worm'):
    
    """
        Cleans up a salinity dataframe as well as applied a time-dependent correction
        for the offset associated with prolonged use of the autosal
        
        Input:
            -saltDF (pandas dataframe), a Dataframe containing the data from
            the output of the autosal salinometer (usually the output of the 
            SaltLoad function)
            -StandName, the value given to the standard seawater calibration
            standard used to calibrate the Conductivity Ratio values. 
            Default Value: worm
            Ex.shown in 5th column:
                Column:     0   1  2  3     4     5                    
                Standard: 0001 01 00 24 1.99967 worm 3807 04:40:13  04:40:13  01 1.99967
                Sample 1: 0001 01 01 24 2.03956    1 3808 04:44:55  04:45:38  03 2.03938
            
        Output:
            -outputDF (pandas dataframe), a corrected Dataframe containing 4 columns: 
            Station Number, Cast Number, Sample Number, and a new Conductivity
            Ratio which has the time-dependent offset from the salinometer
            removed from it.
            
        Usage:
            outputDF = Cr_Calibration(saltDF)
    
    """
    
    CrStrt = saltDF['CRavg'][saltDF.autosalSAMPNO==StandName][0] #First Standard Cr Value
    CrEnd = saltDF['CRavg'][saltDF.autosalSAMPNO==StandName][len(saltDF['CRavg'])-1] #Second Standard Cr Value
   
    #calculate start and end times (endtimes are when measurements are taken)    
    saltDF['EndTime'] = saltDF['EndTime'].apply(pd.Timedelta)

    startTime = saltDF['EndTime'][0]
    saltDF['ElapsedTime_(s)'] = (saltDF['EndTime']-startTime) / np.timedelta64(1,'s')

    duration = saltDF['ElapsedTime_(s)'][len(saltDF['ElapsedTime_(s)'])-1]
    Offset = CrEnd-CrStrt
    Offset = Offset/duration # Offset per second
       
    #Apply Offsets to Measured Data
    saltDF['CRavg_Corr'] = saltDF['CRavg']-(Offset*saltDF['ElapsedTime_(s)'])       
            
    saltDF = saltDF[(saltDF['autosalSAMPNO']!=StandName)]  #remove calibration entries
    #Create Export Dataframe
    outputDF = pd.DataFrame()
    outputDF = saltDF.loc[:,['STNNBR','CASTNO','SAMPNO','CRavg_Corr']] #copy wanted columns to new DF     
    return outputDF

def saltCat(saltDir):
    
    """
        Concatenates all corrected salt dataframes in the user-specified directory
        and writes dataframe to a master salt .csv/.pkl file in that directory
        
        Inputs:
            - saltDir (Directory), the directory containing the corrected dataframes
            for each file that is to be concatenated.
            
        Outputs:
            - a master salt .csv/.pkl file containing all a master dataframe saved
            to the input directory
            
        Usage:
            saltCat('/data/salt')
        
    """
    
    fileName = 'master_salt_DF.csv'  
    fileList = os.listdir(path=saltDir) #Creates list of files in the salt Directory
    exten = '_corr.csv' # type of file to be parsed out into dataframe
    extFiles = []
    for i in range(len(fileList)-1): # Parse out files that have the wanted extension
        if fileList[i][-9:] == exten:
            extFiles.append(fileList[i])
    masterDF = pd.DataFrame()
    for i in extFiles: #concatenate all Dataframes together
        catFrame = pd.read_csv(i)
        masterDF = pd.concat([masterDF,catFrame])
    #Clean up data before saving
    if 'Count' in masterDF.columns: #Remove extra 'Count' column in DF
        del masterDF['Count']
    
    masterDF.index = range(len(masterDF))
            
    f = open(fileName,'w')
    masterDF.to_csv(fileName,index_label='Count')
    f.close()

#if __name__ == '__main__':
##   Calibration Analysis:    
#    outputDF=Cr_Calibration(sys.argv[1:][0])
##    #save data into pickle format
#    exten = '_corr.csv'
#    fpath = os.path.split(sys.argv[1:][0])
#    file = fpath[1]+exten    
#    f = open(file,'w')
#    outputDF.to_csv(file,index_label='Count')
#    f.close()

#   Concatenation Analysis:
#    saltCat(sys.argv[1:][0])
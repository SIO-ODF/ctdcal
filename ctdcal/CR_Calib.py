#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 06:04:34 2017

This program loads in a salt file (no extension) from the autosal salinometer,
converts it into a dataframe and appropriately attaches headers based on the 
parameters of the measurement. From there, it determines the amount of offset
from the salinometer as a function of time and applies a correction to each value
in the file.

INPUT:
    unextended salt file/Pandas Dataframe.
    
OUTPUT: 4-column dataframe with Station Number (STNNBR), Cast Number (CASTNO),
        Sample Number (SAMPNO) and the corrected average Conductivity Ratio 
        (CRavg).


@author: Kenneth Jackson
"""

import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import os
import pdb

def Cr_Calibration(saltFile):
  
#  IF isfile statement    
#   if os.path.isfile('./00101')==True:
    #Load in salt file
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
        
    CrStrt = saltDF['CRavg'][saltDF.autosalSAMPNO=='worm'][0] #First Standard Cr Value
    CrEnd = saltDF['CRavg'][saltDF.autosalSAMPNO=='worm'][len(saltDF['CRavg'])-1] #Second Standard Cr Value
   
    #calculate start and end times (endtimes are when measurements are taken)    
    saltDF['EndTime'] = saltDF['EndTime'].apply(pd.Timedelta)

    #start = pd.Timedelta('00:00:00')    
    saltDF['ElapsedTime_(s)'] = 0 
    t = pd.Timedelta('00:00:00')
    for i in range(1,len(saltDF['EndTime'])):        
        t = t+saltDF['EndTime'][i] - saltDF['EndTime'][i-1]
        saltDF['ElapsedTime_(s)'].iat[i] = t.seconds

    duration = saltDF['ElapsedTime_(s)'][len(saltDF['ElapsedTime_(s)'])-1]
    Offset = CrEnd-CrStrt
    Offset = Offset/duration # Offset per second
       
    #Apply Offsets to Measured Data
    saltDF['CRavg_Corr']=saltDF['CRavg'][0]
      
    for i in range(1,len(saltDF['CRavg'])-1):
        if Offset>=0:   #Positive or 0 Drift (Samples Hi or no drift)
            saltDF['CRavg_Corr'].iat[i] = saltDF['CRavg'][i]-(Offset*saltDF['ElapsedTime_(s)'][i])
        else:     #Negative Drift (Samples Low)
            saltDF['CRavg_Corr'].iat[i] = saltDF['CRavg'][i]+(Offset*saltDF['ElapsedTime_(s)'][i])    
            
#    #Plots Original Data vs. Corrected
#    plt.figure(999)
#    plt.plot(saltDF['ElapsedTime_(s)'][1:len(saltDF['ElapsedTime_(s)'])-1],\
#                    saltDF['CRavg'][1:len(saltDF['CRavg'])-1],"bo",label='Original Data')
#    plt.plot(saltDF['ElapsedTime_(s)'][1:len(saltDF['ElapsedTime_(s)'])-1],\
#                    saltDF['CRavg_Corr'][1:len(saltDF['CRavg_Corr'])-1],"rx",label='Corrected Data')
#    plt.xlabel('Elapsed Time (s)')
#    plt.ylabel('CR Value')
#    plt.title('Conductivity Ratio Comparison')
#    plt.legend()
#    plt.show()
    
    #Create Export Dataframe
    outputDF = pd.DataFrame()
    outputDF = saltDF.loc[:,['STNNBR','CASTNO','SAMPNO','CRavg_Corr']] #copy wanted columns to new DF
       
    #DATAframe Concatenate Dataframe to a master DF in a file
    
    return outputDF

if __name__ == '__main__':
    pdb.set_trace()
    outputDF=Cr_Calibration(sys.argv[1:][0])
    pdb.set_trace()
    #save data into pickle format
    exten = '.pkl'
    fpath = os.path.split(sys.argv[1:][0])
    file = fpath[1]+exten
    
    f = open(file,'w')
    outputDF.to_pickle(file)
    f.close()

    #f=open(file,'r')
    #data = pd.read_pickle(file)
    
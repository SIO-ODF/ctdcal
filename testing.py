#! /usr/bin/env python
# Python script fit RINKOIII Data to DO bottle data 
# extracts dissolved o2 data from sensors and bottles
# and fits dat to botlle data using least squares fit
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as sig
import configparser 
import process_ctd 

def CalcPressureAverages():
  return

def roll_filter(up, inMat):
    """roll_filter function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and serveral arguments and returns inMat will ship roll removed. 

    Args:
        param1 (ndarray): numpy ndarray with dtype array 
        param2 (float): intP, pressure interval, default 2.0 dbar 

    Returns:
        Narray: The return value ndarray of data with ship roll removed

    """
    # p0, is the option record
    # rollfilter is an option, but should not be - should be manditory
    # up, is an option to start from bottom and pressure sequence up. Used in cases where the down cast record was damaged. 
    # pF is the pressure filter array pF[0] record starts 
    # prvPrs previous pressure ???
    # startP starting P
    # pPs pressure sequence
    # nInt number of intervals
    # pD is the an array of pressure
    # pN an integers of pF contents
    # m is an int = to pN[7] units and mutiplied by 0xffff
    
    #if up is True or not()

    return inMat

def pressure_sequence(inMat, intP=2.0, startP=0.0, endP=10000.0, Pdir='down', fromP=-1.0 , toP=-10, delay=-1.0):
    """pressure_sequence function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and several arguments to return a pressure sequenced data ndarray. 

    Args:
        param1 (ndarray): numpy ndarray with dtype array 
        param2 (float): intP, pressure interval, default 2.0 dbar 
        param3 (float): starting pressure interval 
        param4 (float): endng pressure for interval
        param5 (str): pressure sequence direction (down/up)
        param6 (float): from time for pressure sequence 
        param7 (float): to time for pressure sequence 
        param8 (float): delay time for start up pressure sequence. 

    Returns:
        Narray: The return value is a matrix of pressure sequenced data 

    """

# Initialise Options
    if (intP <= 0.0):
        print("Pressure interval must be larger that 0.0")
        return
    if (startP >= endP):
        v = startP
        startP = endP
        endP = v
        if (startP is endP):
            endP += intP
        if ((Pdir is 'up') or (fromP >= 0.0)):
            delay = 0.0
        elif (delay < 0.0):
            delay = 2.0

  # Calculate half pressure sequence and number of intervals
    pSeq_half = intP*0.5
    nInt = int((endP - startP + pSeq_half)/intP) + 1
    pSeq_nInt = nInt

  # Allocate Buffers
  # Not sure we need to do this in python


  # psarr pressure sequenced array
    return psarr

# Import configuration file
config = configparser.ConfigParser()
config.read('configuration.ini')
rdir = config['ctd_processing']['raw_data_directory']
infile = rdir+config['ctd_processing']['raw_data_file']
tc1_align = config['ctd_processing']['TC_primary_align']
tc2_align = config['ctd_processing']['TC_secondary_align']
do_align = config['ctd_processing']['DO_align']

# Directory input
dframe = process_ctd.dataToDataFrame(infile)
dtlist = np.dtype(list(dframe.columns.values))
dmatrix = process_ctd.dataToMatrix(infile,dtlist,',') 

if tc1_align: dmatrix = process_ctd.ctd_align(dmatrix,'SBE43FV', float(tc1_align))
if tc2_align: dmatrix = process_ctd.ctd_align(dmatrix,'SBE43FV', float(tc2_align))
if do_align: dmatrix = process_ctd.ctd_align(dmatrix,'SBE43FV', float(do_align))
dmatrix = process_ctd.ondeck_pressure(dmatrix, float(config['ctd_processing']['conductivity_start'])) 
print(dmatrix.dtype.names)

# Filter data
dmatrix['Pdbar'] = process_ctd.raw_ctd_filter(dmatrix['Pdbar'], 'triangle',24)
dmatrix['T1C'] = process_ctd.raw_ctd_filter(dmatrix['T1C'], 'triangle',24)
dmatrix['T2C'] = process_ctd.raw_ctd_filter(dmatrix['T2C'], 'triangle',24)
dmatrix['C1mScm'] = process_ctd.raw_ctd_filter(dmatrix['C1mScm'], 'triangle',24)
dmatrix['C2mScm'] = process_ctd.raw_ctd_filter(dmatrix['C2mScm'], 'triangle',24)
dmatrix['SBE43FV'] = process_ctd.raw_ctd_filter(dmatrix['SBE43FV'], 'triangle',24)
dmatrix['V0V'] = process_ctd.raw_ctd_filter(dmatrix['V0V'], 'triangle',24)
dmatrix['V1V'] = process_ctd.raw_ctd_filter(dmatrix['V1V'], 'triangle',24)
dmatrix['V2V'] = process_ctd.raw_ctd_filter(dmatrix['V2V'], 'triangle',24)
dmatrix['V3V'] = process_ctd.raw_ctd_filter(dmatrix['V3V'], 'triangle',24)
dmatrix['V4V'] = process_ctd.raw_ctd_filter(dmatrix['V4V'], 'triangle',24)
dmatrix['V5V'] = process_ctd.raw_ctd_filter(dmatrix['V5V'], 'triangle',24)
dmatrix['V6V'] = process_ctd.raw_ctd_filter(dmatrix['V6V'], 'triangle',24)
dmatrix['V7V'] = process_ctd.raw_ctd_filter(dmatrix['V7V'], 'triangle',24)


# Leave only the actual cast data 
# Find start, end and bottom time
stime, etime, btime, startP, maxP, dmatrix = process_ctd.cast_details(dmatrix)
# Testing cast start time.

# Testing cast start time.
#start_time, smatrix = cast_start(dmatrix)

#plottitle = cruise+': Filter' 
#plotfile = maindir+'BoxCar/bc.00101.png' 
plt.plot(dmatrix['TIMEs'], dmatrix['Pdbar'], color='b', label='Temp1')
#plt.plot(dmatrix['T1C'], dmatrix['Pdbar'], color='r', label='Temp1')
#plt.plot(dmatrix['T2C'], dmatrix['Pdbar'], color='b', label='Temp2')
#plt.plot(dmatrix['C1mScm'], dmatrix['Pdbar'], color='r', label='Cond1')
#plt.plot(dmatrix['C2mScm'], dmatrix['Pdbar'], color='b', label='Cond2')
#plt.plot(dmatrix['SBE43FV'], dmatrix['Pdbar'], color='g', label='DO')
#plt.ylim([-10,max(dmatrix['Pdbar'])*1.1])
#plt.xlim([-5,40])
plt.gca().invert_yaxis()
##plt.title(plottitle)
#plt.xlabel('')
#plt.ylabel('P')
#plt.legend( loc='best')
plt.axis()
plt.show()
#plt.savefig(maindir+plotfile)
#plt.close()
#plt.show()

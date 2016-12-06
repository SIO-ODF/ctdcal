#!/usr/bin/env python
import numpy as np
import scipy.signal as sig
import scipy.stats as st
import datetime 
import pandas as pd
import math

def report_time_series_data(stacast, printdir, expocode, column_names, column_units, column_data, column_format, inMat=None):
    """report_time_series_data function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and writes time series data with quality codes to csv file.   

    Args:
        param1 (str): stacast, station cast information
        param2 (str): printdir, file output directory 
        param3 (str): expocode, exposition code for data repository 
        param4 (str): column_names, column header names 
        param5 (str): column_units, column header units 
        param6 (str): column_data, ndarray of qulity codes. 
        param7 (str): column_format, ndarray of qulity codes. 
        param8 (ndarray): inMat, input data ndarray  

    Prints formatted csv file. 

    Returns:
        No return
    """

    if inMat is None:
       print("In report_time_series_data: No data")
       return
    else:
        out_col = []
        now = datetime.datetime.now()
        file_datetime = now.strftime("%Y%m%d %H:%M") 

        outfile = open(printdir+stacast+'_time.csv', "w+")
        outfile.write('expocode: '+expocode+', station cast: '+stacast+', printed: '+file_datetime+'\n')
        cn = np.asarray(column_names)
        cn.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')
        cu = np.asarray(column_units)
        cu.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')

        # Need to rewrite to remove hardcoded column data. 
        # No ideal method to print formatted output to csv in python ATM 
        # Other attemps to import formats from config file 
        # or use savetxt and csv have been buggy so far
        for i in range(0,len(inMat)-1):
            outfile.write("%8.1f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%d,%10.5f,%10.5f\n" % (inMat['p_dbar'][i], inMat['t1_C'][i], inMat['t2_C'][i], inMat['c1_mScm'][i], inMat['c1_mScm'][i], inMat['sal_PSU'][i], inMat['o1_mll'][i], inMat['cstar_ugl'][i], inMat['fluoro_ugl'][i], inMat['scan_datetime'][i], inMat['lat_ddeg'][i], inMat['lon_ddeg'][i]))
            #outfile.write(column_format % (column_data))

        #print(out_col)
        #np.savetxt(outfile, out_col, fmt=column_format, delimiter=',')
        #np.savetxt(outfile, out_col, delimiter=',')
  
    return

def report_pressure_series_data(inMat=None):
    """report_time_series_data function 

    Function takes full NUMPY ndarray with predefined dtype array 
    and writes time series data with quality codes to csv file.   

    Args:
        param1 (ndarray): inMat, input data ndarray  

    Prints formatted csv file. 

    Returns:
        No return
    """

    if input_array is None:
        print("In report_pressure_series_data: Not data array.")
        return
    #else:  

    return

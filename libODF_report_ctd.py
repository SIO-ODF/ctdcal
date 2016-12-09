#!/usr/bin/env python
import numpy as np
import scipy.signal as sig
import scipy.stats as st
import datetime 
import pandas as pd
import math

def report_pressure_details(stacast, log_file, start, end):
    """report_cast_details function 

    Function takes deck pressure and writes them to a file in log 
    directory 

    Args:
        param1 (str): stacast, station cast data for file 
        param2 (str): c_file, file name location for cast details 
        param3 (str): start, cast start time from top of cast after 10m soak 
        param4 (str): end, cast end time when instrument leaves water 

    Prints formatted csv file. 

    Returns:
        No return
    """
    outfile = open(log_file, "w+")
    outfile.write("stacast:%s, ondeck_start_p:%s, ondeck_end_p:%s\n" % (stacast, start, end))
    return

def report_cast_details(stacast, c_file, start, end, bottom, start_p, max_p):
    """report_cast_details function 

    Function takes cast details and writes them to a file in log 
    directory 

    Args:
        param1 (str): stacast, station cast data for file 
        param2 (str): c_file, file name location for cast details 
        param3 (str): start, cast start time from top of cast after 10m soak 
        param4 (str): end, cast end time when instrument leaves water 
        param5 (str): bottom, bottom of cast time when instrument reaches max depth 
        param6 (str): start_p, starting pressure at the time the cast begins 
        param7 (str): max_p, maximum pressure for entire cast 

    Prints formatted csv file. 

    Returns:
        No return
    """
    
    outfile = open(c_file, "w+")
    outfile.write("stacast:%s, , begin:%s, , bottom:%s, end:%s, start_pressure:%s, max_pressure:%s\n" % (stacast, start, bottom, end, start_p, max_p))
    outfile.close()

    return

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

def report_pressure_series_data(stacast, expocode, section_id, btm_lat, btm_lon, depth, btm_alt, ctd, inMat=None):
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
        print("In report_pressure_series_data: No data array.")
        return
    #else:  

    return

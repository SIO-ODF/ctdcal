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
    outfile = open(log_file, 'a')
    outfile.write("stacast:%s, ondeck_start_p:%s, ondeck_end_p:%s\n" % (stacast, start, end))
    return


def report_polyfit(coef, stacast_list, fitfile):
    """report_polyfit

    Function takes deck pressure offset and stores them for later fitting

    Args:
        p_off_file (str): path to file

    Prints formatted csv file.

    Returns:
        No return
    """
    outfile = open(fitfile, 'w+')
    coefTostr = ', '.join(map(str, coef))

    try:
        for stacast in stacast_list:
            outfile.write("stacast:%s, coef:%s\n" % (stacast, coefTostr))
    except:
        outfile.write("stacast:%s, coef:%s\n" % (stacast_list, coefTostr))
    return


def report_pressure_offset(p_off_file, p_offset, stacast_list):
    """report_pressure_offset function

    Function takes deck pressure offset and stores them for later fitting

    Args:
        p_off_file (str): path to file
        p_offset (float): offset to apply to specific station cast
        stacast (str): station cast for later file reference

    Prints formatted csv file.

    Returns:
        No return
    """
    outfile = open(p_off_file, 'a')
    try:
        for stacast in stacast_list:
            outfile.write("stacast:%s, offset:%9.4f\n" % (stacast, p_offset))
    except:
        outfile.write("stacast:%s, offset:%9.4f\n" % (stacast_list, p_offset))
    return


def report_cast_details(stacast, c_file, start, end, bottom, start_p, max_p, b_alt, b_lat, b_lon):
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
        param7 (str): max_p, maximum pressure for entire cast.
        param8 (str): b_alt, altimeter value at bottom of cast.
        param9 (str): b_lat, latitude at bottom of cast.
        param10 (str): b_lon, longitude at bottom of cast.

    Prints formatted csv file.

    Returns:
        No return
    """

    outfile = open(c_file, "a")
    outfile.write("stacast:%s, begin:%s, at_depth:%s, end:%s, start_pressure:%s, max_pressure:%s, altimeter_bottom:%s, latitude:%s, longitude:%s\n" % (stacast, start, bottom, end, start_p, max_p, b_alt, b_lat, b_lon))
    outfile.close()

    return


def report_quality_flags(qual_file, stacast, btl, pres, param, d1, d2, d12):
    """report_quality flags

    Prints formatted csv file.

    Returns:
        No return
    """
    #print(stacast+' '+str(btl)+' '+param+' '+str(pres)+' '+str(d1)+' '+str(d2)+' '+str(d12))
    outfile = open(qual_file, "a")
    outfile.write("stacast: %s, bottle: %s, parameter: %s, pressure: %8.3f, flag: 3, Primary Diff: %f, Secondary Diff: %f, P-S: %f\n" % (stacast, str(btl), param, pres, d1, d2, d12))
    outfile.close()

    return


def report_btl_data(btl_file, btl_dtype, btl_data):
    """report_btl_data function

    Args:
        btl_file (str): file path
        btl_dtype (dtype): list
        btl_data (ndarray):input data ndarray.

    Prints formatted csv file.

    Returns:
        No return
    """
    try:
        btl_data = pd.DataFrame.from_records(btl_data)
        #cut index from btl_data, i hope. Not guaranteed to work correctly...
        btl_data = btl_data.iloc[:, 1:]
        btl_data.to_pickle(btl_file)
    except:
        dtype_list = btl_data.dtype.names

        if btl_data is None:
            print("In report_btl_data: No data array.")
            return
        else:
            outfile = open(btl_file, "w+")

            # Print out headers
            outfile.write("%s" % (dtype_list[0]))
            for i in range(1, len(dtype_list)):
                outfile.write(",%s" % (dtype_list[i]))
            outfile.write('\n'+btl_dtype+'\n')

            # Print out data
            for i in range(0,len(btl_data)):
                outfile.write("%s" % (btl_data[dtype_list[0]][i]))
                for j in range(1, len(dtype_list)):
                    outfile.write(",%s" % (btl_data[dtype_list[j]][i]))
                outfile.write('\n')

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
    try:
        inMat = pd.DataFrame.from_records(inMat)
        #inMat = inMat.iloc[:, 1:]
        inMat.to_pickle(printdir+stacast+'_time.pkl')
    except:
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
                outfile.write('%9.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%12d,%10.5f,%10.5f\n' % (inMat[column_data[0]][i], inMat[column_data[1]][i], inMat[column_data[2]][i], inMat[column_data[3]][i], inMat[column_data[4]][i], inMat[column_data[5]][i], inMat[column_data[6]][i], inMat[column_data[7]][i], inMat[column_data[8]][i], inMat[column_data[9]][i], inMat[column_data[10]][i], inMat[column_data[11]][i], inMat[column_data[12]][i]))

    return

def report_pressure_series_data(stacast, expocode, section_id, btime=-999, btm_lat=-999, btm_lon=-999, depth=-999, btm_alt=-999, ctd=-999, p_dir=None, p_column_names=None, p_column_units=None, p_column_data=None, qualMat=None, doArr=None, inMat=None):
    """report_pressure_series_data function

    Function takes full NUMPY ndarray with predefined dtype array
    and writes time series data with quality codes to csv file.

    Args:
        param1 (str): stacast, actually file base name from hex data.
        param2 (str): expocode, vessel id and startd date cruise identification code
        param3 (str): section_id, US hydrographic section id
        param4 (str): btime, UTC date time stamp at bottom of cast
        param5 (str): btm_lat, latitude value at bottom of cast.
        param6 (str): btm_lon, longitude value at bottom of cast.
        param7 (str): depth, ctd depth value + altimeter value.
        param8 (str): btm_alt, altimeter value at bottom of cast.
        param9 (str): ctd, serial number of ctd inst.
        param10 (str): p_dir, directory to write file.
        param11 (str): p_column_names, header column names.
        param12 (str): p_column_units, header column units.
        param13 (dataframe): qualMat, input dataframe.
        param13 (ndarray): doArr, input do data.
        param14 (ndarray): inMat, input data ndarray.

    Prints formatted csv file.

    Returns:
        No return
    """

    if inMat is None:
        print("In report_pressure_series_data: No data array.")
        return
    else:
        out_col = []
        h_num = 11
        now = datetime.datetime.now()
        file_datetime = now.strftime("%Y%m%d %H:%M")
        bdt = datetime.datetime.fromtimestamp(btime).strftime('%Y%m%d %H%M').split(" ")
        b_date = bdt[0]
        b_time = bdt[1]

        s_num = stacast[-5:-2]
        c_num = stacast[-2:]

        outfile = open(p_dir+stacast+'_ct1.csv', "w+")
        outfile.write("CTD, %s\nNUMBER_HEADERS = %s \nEXPOCODE = %s \nSECT_ID = %s\nSTNNBR = %s\nCASTNO = %s\n DATE = %s\nTIME = %s\nLATITUDE = %f\nLONGITUDE = %f\nDEPTH = %s\nINSTRUMENT_ID = %s\n" % (file_datetime, h_num, expocode, section_id, s_num, c_num, b_date, b_time, btm_lat, btm_lon, depth, ctd))
        cn = np.asarray(p_column_names)
        cn.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')
        cu = np.asarray(p_column_units)
        cu.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')

        # Need to rewrite to remove hardcoded column data.
        # No ideal method to print formatted output to csv in python ATM
        for i in range(0,len(inMat)):
            ### Alter configuration.ini in pressure_series in order to change output
            ### Will need to be altered in order to put out RINKO in UMOL/KG
            ### CTDOXY is slipped in and piggybacks on CTDXMISS qual code
            outfile.write("%8.1f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d\n" % (inMat[p_column_data[0]][i], qualMat[i][0], inMat[p_column_data[1]][i], qualMat[i][1], inMat[p_column_data[2]][i], qualMat[i][2], doArr['CTDOXY1'][i], qualMat[i][3], inMat[p_column_data[3]][i], qualMat[i][3], inMat[p_column_data[4]][i], qualMat[i][4],inMat[p_column_data[5]][i], qualMat[i][5], inMat[p_column_data[6]][i], qualMat[i][6]))
        outfile.write('END_DATA')
    return

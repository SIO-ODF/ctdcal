#!/usr/bin/env python
import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def report_pressure_details(stacast, log_file, start, end):
    """
    Write start/end deck pressure to ondeck_pressure.csv log file.

    Parameters
    ----------
    stacast : str
        station cast data for file
    c_file : str
        file name location for cast details
    start : str
        cast start time from top of cast after 10m soak
    end : str
        cast end time when instrument leaves water

    Returns
    -------
    None
    """
    df = pd.DataFrame(
        {"SSSCC": stacast, "ondeck_start_p": start, "ondeck_end_p": end}, index=[0]
    )
    add_header = not Path(log_file).exists()  # add header iff file doesn't exist
    with open(log_file, "a") as f:
        df.to_csv(f, mode='a', header=add_header, index=False)

    return True


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
    """
    Write cast details to cast_details.csv log file.

    Parameters
    ----------
    stacast : str
        station cast data for file
    c_file : str
        file name location for cast details
    start : str
        cast start time from top of cast after 10m soak
    end : str
        cast end time when instrument leaves water
    bottom : str
        bottom of cast time when instrument reaches max depth
    start_p : str
        starting pressure at the time the cast begins
    max_p : str
        maximum pressure for entire cast
    b_alt : str
        altimeter value at bottom of cast
    b_lat : str
        latitude at bottom of cast
    b_lon : str
        longitude at bottom of cast

    Returns
    -------
    None
    """
    df = pd.DataFrame(
        {
            "SSSCC": stacast,
            "start_time": start,
            "bottom_time": bottom,
            "end_time": end,
            "start_pressure": start_p,
            "max_pressure": max_p,
            "altimeter_bottom": b_alt,
            "latitude": b_lat,
            "longitude": b_lon,
        },
        index=[0],
    )
    add_header = not Path(c_file).exists()  # add header iff file doesn't exist
    with open(c_file, "a") as f:
        df.to_csv(f, mode="a", header=add_header, index=False)

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
        # inMat is already a DataFrame
        # inMat = pd.DataFrame.from_records(inMat)
        # remove duplicate index column and reset to count from 0
        inMat = inMat.drop(labels='index', axis=1).reset_index(drop=True)
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
        bdt = datetime.datetime.fromtimestamp(btime, tz=datetime.timezone.utc).strftime('%Y%m%d %H%M').split(" ")
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

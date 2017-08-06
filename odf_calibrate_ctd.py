#! /usr/bin/env python
import matplotlib
matplotlib.use('agg')
import sys
import os
import math
import argparse
import fnmatch
import pylab
import numpy as np
import pandas as pd
import scipy as sp
import json
import libODF_process_ctd as process_ctd
import libODF_report_ctd as report_ctd
import libODF_fit_ctd as fit_ctd
import configparser
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'
PKL_EXT = 'pkl'

#File extension to use for raw output
REFT_SUFFIX = '_reft'

#File extension to use for raw output
BTL_SUFFIX = '_btl'

#File extension to use for raw output
FIT_SUFFIX = '_fit'

#File extension to use for raw output
MEAN_SUFFIX = '_mean'

#File extension to use for raw output
TIME_SUFFIX = '_time'

#File extension to use for converted output
CONVERTED_SUFFIX = '_cnv'

def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# -------------------------------------------------------------------------------------
# Main function of the script should it be run as a stand-alone utility.
# -------------------------------------------------------------------------------------
def main(argv):

    parser = argparse.ArgumentParser(description='Convert SBE raw data to a converted, csv-formatted text file')
    #parser.add_argument('castRange', metavar='cast_range', help='Range of cast values for shipboard calibration.')

    parser.add_argument('cast_file', metavar='cast_file', help='File of cast values for shipboard calibration.')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='Display debug messages')

    # update bottle data with down trace profile data
    parser.add_argument('-btl', '--bottle_data', action='store_true', help='Update bottle data with downtrace profile.')

    # Calibrate pressure
    parser.add_argument('-press', '--pressure', action='store_true', help='Calibrate pressure.')

    # Calibrate temperature
    parser.add_argument('-temp', '--temperature', action='store_true', help='Calibrate temperature.')

    # Calibrate conductivity
    parser.add_argument('-cond', '--conductivity', action='store_true', help='Calibrate conductivity.')

    # Calibrate primary side
    parser.add_argument('-primary', '--primary', action='store_true', help='Calibrate primary.')

    # Calibrate sensondary side
    parser.add_argument('-secondary', '--secondary', action='store_true', help='Calibrate secondary.')

    # Independent parameter to calibrate with
    parser.add_argument('-calib', metavar='calibrate', dest='calib', help='Independent parameter to calabrate with (P, T or C).')

    # Independent parameter range
    parser.add_argument('-xRange', metavar='xRange', dest='xRange', help='Independent parameter range.')

    # Order of fit
    parser.add_argument('-order', metavar='order', dest='order', help='Calibration order.')

    # Process Command-line args
    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")


    #Import Cruise Configuration File
    iniFile = 'data/ini-files/configuration.ini'
    config = configparser.RawConfigParser()
    config.read(iniFile)

    #Initialise Configuration Parameters
    expocode = config['cruise']['expocode']
    sectionID = config['cruise']['sectionid']
    time_directory = config['ctd_processing']['time_data_directory']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    reft_directory = config['ctd_processing']['reft_directory']
    salt_directory = config['ctd_processing']['salt_directory']
    btl_directory = config['ctd_processing']['bottle_directory']
    log_directory = config['ctd_processing']['log_directory']
    ctd = config['ctd_processing']['ctd_serial']

    time_zone = config['inputs']['time_zone']
    p_col = config['analytical_inputs']['p']
    p_btl_col = config['inputs']['p']
    t_col = config['analytical_inputs']['t']
    t_btl_col = config['inputs']['t']
    t1_col = config['analytical_inputs']['t1']
    t1_btl_col = config['inputs']['t1']
    t2_col = config['analytical_inputs']['t2']
    t2_btl_col = config['inputs']['t2']
    c_col = config['analytical_inputs']['c']
    c_btl_col = config['inputs']['c']
    c1_col = config['analytical_inputs']['c1']
    c1_btl_col = config['inputs']['c1']
    c2_col = config['analytical_inputs']['c2']
    c2_btl_col = config['inputs']['c2']
    sal_col = config['analytical_inputs']['salt']
    sal_btl_col = config['inputs']['salt']
    dov_col = config['inputs']['dov']
    dov_btl_col = config['inputs']['dov']
    btl_sal_col = config['analytical_inputs']['btl_salt']
    timedate = config['analytical_inputs']['datetime']
    lat_col = config['analytical_inputs']['lat']
    lat_btl_col = config['inputs']['lat']
    lon_col = config['analytical_inputs']['lon']
    lon_btl_col = config['inputs']['lon']
    reft_col = config['inputs']['reft']
    btl_num_col = config['inputs']['btl_num']
    btl_dtype = config['bottle_series_output']['dtype']

    file_base_arr = []

    logfileName = str('ondeck_pressure' + '.' + FILE_EXT)
    logfilePath = os.path.join(log_directory, logfileName)
    log_data = process_ctd.dataToNDarray(logfilePath,str,None,',',None)
    p_start = []
    p_end = []

    if args.cast_file:
        with open(args.cast_file, 'r') as filename:
            files = ['{}{}.{}'.format(line.strip(), TIME_SUFFIX, FILE_EXT) for line in filename]
    #Might be worthwhile to fix in future, but use a hand edited station list for now.
    #if args.castRange:
    #    casts = str.split(args.castRange, '-')
    #    if int(casts[0][2]) is int(casts[1][2]):
    #        regEx_fileList = '['+casts[0][0]+'-'+casts[1][0]+']['+casts[0][1]+'-'+casts[1][1]+']?0?'+TIME_SUFFIX+'.'+FILE_EXT
    #    else:
    #        regEx_fileList = '['+casts[0][0]+'-'+casts[1][0]+']['+casts[0][1]+'-'+casts[1][1]+']['+casts[0][2]+'-'+casts[1][2]+']0?'+TIME_SUFFIX+'.'+FILE_EXT

    #    files = [f for f in os.listdir(time_directory) if fnmatch.fnmatch(f, regEx_fileList)]

        for f in files:
            # Used later for building output file names
            filename_base = str.split(f,'_')[0] # original filename w/o ext
            file_base_arr.append(filename_base)

        # Update bottle file with isopycnal equivalnt from downtrace
        if args.bottle_data:
            for sc in file_base_arr:
                # Collect bottle data
                btlfileName = str(sc + BTL_SUFFIX + MEAN_SUFFIX + '.' + PKL_EXT)
                btlfilePath = os.path.join(btl_directory, btlfileName)
                if os.path.isfile(btlfilePath):
                    btl_data = process_ctd.dataToNDarray(btlfilePath,float,True,',',None)
                    btl_data = btl_data[:][1:]
                else:
                    print("Missing file: "+btlfilePath)
                    break

                timefileName = str(sc + TIME_SUFFIX + '.' + PKL_EXT)
                timefilePath = os.path.join(time_directory, timefileName)
                if os.path.isfile(timefilePath):
                    time_data = process_ctd.dataToNDarray(timefilePath,float,True,',',1)
                    time_data = time_data[:][1:]
                else:
                    print("Missing file: "+timefilePath)
                    break

                # Find Isopycnal Down Trace Bottle Trip Equivalent
                end = np.argmax(time_data[p_col])
                down_trace_btl = fit_ctd.find_isopycnals(p_btl_col, t1_btl_col, sal_btl_col, dov_btl_col, lat_btl_col, lon_btl_col, btl_data, p_col, t1_col, sal_col, dov_col, lat_col, lon_col, time_data)

                # Write bottle data to file
                report_ctd.report_btl_data(btlfilePath, btl_dtype, down_trace_btl)

            # Exit Bottle Data Update
            sys.exit()


        # Get Pressure Calibration Offset for castrange.
        # And save in Log Files
        if args.pressure:
            for sc in file_base_arr:
                for line in log_data:
                    # Search log file for station cast input
                    if sc in line[0]:
                        for val in line:
                            if 'ondeck_start_p' in val:
                                s_p = str.split(val, ':')
                                if ('Start' not in s_p[1]) and (s_p[1] != 'nan'):
                                    p_start.append(float(s_p[1]))
                            if 'ondeck_end_p' in val:
                                e_p = str.split(val, ':')
                                if ('End' not in e_p[1]) and (e_p[1] != 'nan'):
                                    p_end.append(float(e_p[1]))

            p_start = np.unique(p_start)
            p_end = np.unique(p_end)
            p_off = np.mean(p_start) - np.mean(p_end)

            # Record pressure offset for cast set list in log file for later fit
            report_ctd.report_pressure_offset(log_directory + 'poffset.' + FILE_EXT, p_off, file_base_arr)

            # Exit Pressure Offset
            sys.exit()

        if args.order:
            order = int(args.order)
        else: print('Set polyfit calibration order. Ex: ./odf_calibrate_ctd.py 03001-03901 -cond -secondary -param P -order 2 -prange 0:6000')

        btl_num = []
        P1 = []
        P2 = []
        P12 = []
        T1 = []
        T2 = []
        T12 = []
        C1 = []
        C2 = []
        C12 = []
        ref_data1 = []
        ref_data2 = []
        ref_data12 = []
        ctd_1 = []
        ctd_2 = []
        ctd_12 = []
        ctd_d1 = []
        ctd_d2 = []
        ctd_d12 = []
        calib1 = []
        calib2 = []
        calib12 = []
        qual = []
        indx = []

        #Determine parameter condition
        if 'P' in args.calib:
            calib1_btl_col = p_btl_col
            calib2_btl_col = p_btl_col
            calib12_btl_col = p_btl_col
        elif 'T' in args.calib:
            calib1_btl_col = t1_btl_col
            calib2_btl_col = t2_btl_col
            calib12_btl_col = t_btl_col
        elif 'C' in args.calib:
            calib1_btl_col = c1_btl_col
            calib2_btl_col = c2_btl_col
            calib12_btl_col = c_btl_col

        if args.temperature:
            param = 'T'
            qualfileName = str('quality_flag_temp.' + FILE_EXT)
            qualfilePath = os.path.join(log_directory, qualfileName)
            coef = np.zeros(shape=5)
            ref_col = reft_col
            ctd1_btl_col = t1_btl_col
            ctd2_btl_col = t2_btl_col
        elif args.conductivity:
            param = 'C'
            qualfileName = str('quality_flag_cond.' + FILE_EXT)
            qualfilePath = os.path.join(log_directory, qualfileName)
            coef = np.zeros(shape=7)
            ref_col = 'BTLCOND'
            ctd1_btl_col = c1_btl_col
            ctd2_btl_col = c2_btl_col

        if os.path.exists(qualfilePath): os.remove(qualfilePath)

        for filename_base in file_base_arr:
            btlfileName = str(filename_base + BTL_SUFFIX + MEAN_SUFFIX + '.' + PKL_EXT)
            btlfilePath = os.path.join(btl_directory, btlfileName)
            if os.path.isfile(btlfilePath):
                btl_data = process_ctd.dataToNDarray(btlfilePath,float,True,',',None)
                btl_data = btl_data[:][1:]
            else:
                print("Missing file: "+btlfilePath)
                sys.exit()

            timefileName = str(filename_base + TIME_SUFFIX + '.' + PKL_EXT)
            timefilePath = os.path.join(time_directory, timefileName)
            if os.path.isfile(timefilePath):
                time_data = process_ctd.dataToNDarray(timefilePath,float,True,',',1)
                time_data = time_data[:][1:]
            else:
                print("Missing file: "+timefilePath)
                sys.exit()

            if args.temperature:
                reftfileName = (filename_base + REFT_SUFFIX + '.' + FILE_EXT)
                reftfilePath = os.path.join(reft_directory, reftfileName)
                if os.path.isfile(reftfilePath):
                    print("Processing "+ filename_base)
                    ref_btl = process_ctd.dataToNDarray(reftfilePath,float,True,',',None)
                else:
                    print("Missing Reference temperature data for "+ filename_base)
                    continue
            elif args.conductivity:
                saltfileName = filename_base
                saltfilePath = os.path.join(salt_directory, saltfileName)
                if os.path.isfile(saltfilePath):
                    print("Processing "+ filename_base)
                    ref_btl, salt_btl = fit_ctd.salt_calc(saltfilePath,btl_num_col,t1_btl_col,p_btl_col,btl_data)
                else:
                    print("Missing reference salinity data for "+ filename_base + " in " + salt_directory)
                    continue

            for i in range(0,len(ref_btl[btl_num_col])):
            # Search bottle number
                j = ref_btl[btl_num_col][i]
                if j != 0:
                    k = int(np.where(btl_data[btl_num_col] == j)[0][0])
                    d1 = ref_btl[ref_col][i] - btl_data[ctd1_btl_col][k]
                    d2 = ref_btl[ref_col][i] - btl_data[ctd2_btl_col][k]
                    d12 = btl_data[ctd1_btl_col][k] - btl_data[ctd2_btl_col][k]
                    if btl_data[p_btl_col][k] > 2000:
                        if abs(d1) < 0.002:
                            ref_data1.append(ref_btl[ref_col][i])
                            P1.append(btl_data[p_btl_col][k])
                            T1.append(btl_data[t1_btl_col][k])
                            C1.append(btl_data[c1_btl_col][k])
                            calib1.append(btl_data[calib1_btl_col][k])
                            ctd_1.append(btl_data[ctd1_btl_col][k])
                            ctd_d1.append(d1)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'1',
                                                        d1, d2, d12)
                        if abs(d2) < 0.002:
                            ref_data2.append(ref_btl[ref_col][i])
                            P2.append(btl_data[p_btl_col][k])
                            T2.append(btl_data[t2_btl_col][k])
                            C2.append(btl_data[c2_btl_col][k])
                            calib2.append(btl_data[calib2_btl_col][k])
                            ctd_2.append(btl_data[ctd2_btl_col][k])
                            ctd_d2.append(d2)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'2',
                                                        d1, d2, d12)
                        if abs(d12) < 0.002:
                            ref_data12.append(ref_btl[ref_col][i])
                            P12.append(btl_data[p_btl_col][k])
                            T12.append(btl_data[t_btl_col][k])
                            C12.append(btl_data[c_btl_col][k])
                            calib12.append(btl_data[calib12_btl_col][k])
                            ctd_d12.append(d12)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param,
                                                        d1, d2, d12)
                    elif (btl_data[p_btl_col][k] < 2000) and (btl_data[p_btl_col][k] > 1000):
                        if abs(d1) < 0.005:
                            ref_data1.append(ref_btl[ref_col][i])
                            P1.append(btl_data[p_btl_col][k])
                            T1.append(btl_data[t1_btl_col][k])
                            C1.append(btl_data[c1_btl_col][k])
                            calib1.append(btl_data[calib1_btl_col][k])
                            ctd_1.append(btl_data[ctd1_btl_col][k])
                            ctd_d1.append(d1)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'1',
                                                        d1, d2, d12)
                        if abs(d2) < 0.005:
                            ref_data2.append(ref_btl[ref_col][i])
                            P2.append(btl_data[p_btl_col][k])
                            T2.append(btl_data[t2_btl_col][k])
                            C2.append(btl_data[c2_btl_col][k])
                            calib2.append(btl_data[calib2_btl_col][k])
                            ctd_2.append(btl_data[ctd2_btl_col][k])
                            ctd_d2.append(d2)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'2',
                                                        d1, d2, d12)
                        if abs(d12) < 0.005:
                            ref_data12.append(ref_btl[ref_col][i])
                            P12.append(btl_data[p_btl_col][k])
                            T12.append(btl_data[t_btl_col][k])
                            C12.append(btl_data[c_btl_col][k])
                            calib12.append(btl_data[calib12_btl_col][k])
                            ctd_d12.append(d12)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param,
                                                        d1, d2, d12)
                    elif (btl_data[p_btl_col][k] < 1000) and (btl_data[p_btl_col][k] > 500):
                        if abs(d1) < 0.010:
                            ref_data1.append(ref_btl[ref_col][i])
                            P1.append(btl_data[p_btl_col][k])
                            T1.append(btl_data[t1_btl_col][k])
                            C1.append(btl_data[c1_btl_col][k])
                            calib1.append(btl_data[calib1_btl_col][k])
                            ctd_1.append(btl_data[ctd1_btl_col][k])
                            ctd_d1.append(d1)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'1',
                                                        d1, d2, d12)
                        if abs(d2) < 0.010:
                            ref_data2.append(ref_btl[ref_col][i])
                            P2.append(btl_data[p_btl_col][k])
                            T2.append(btl_data[t2_btl_col][k])
                            C2.append(btl_data[c2_btl_col][k])
                            calib2.append(btl_data[calib2_btl_col][k])
                            ctd_2.append(btl_data[ctd2_btl_col][k])
                            ctd_d2.append(d2)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'2',
                                                        d1, d2, d12)
                        if abs(d12) < 0.010:
                            ref_data12.append(ref_btl[ref_col][i])
                            P12.append(btl_data[p_btl_col][k])
                            T12.append(btl_data[t_btl_col][k])
                            C12.append(btl_data[c_btl_col][k])
                            calib12.append(btl_data[calib12_btl_col][k])
                            ctd_d12.append(d12)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param,
                                                        d1, d2, d12)
                    elif btl_data[p_btl_col][k] < 500:
                        if abs(d1) < 0.020:
                            ref_data1.append(ref_btl[ref_col][i])
                            P1.append(btl_data[p_btl_col][k])
                            T1.append(btl_data[t1_btl_col][k])
                            C1.append(btl_data[c1_btl_col][k])
                            calib1.append(btl_data[calib1_btl_col][k])
                            ctd_1.append(btl_data[ctd1_btl_col][k])
                            ctd_d1.append(d1)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'1',
                                                        d1, d2, d12)
                        if abs(d2) < 0.020:
                            ref_data2.append(ref_btl[ref_col][i])
                            P2.append(btl_data[p_btl_col][k])
                            T2.append(btl_data[t2_btl_col][k])
                            C2.append(btl_data[c2_btl_col][k])
                            calib2.append(btl_data[calib2_btl_col][k])
                            ctd_2.append(btl_data[ctd2_btl_col][k])
                            ctd_d2.append(d2)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param+'2',
                                                        d1, d2, d12)
                        if abs(d12) < 0.020:
                            ref_data12.append(ref_btl[ref_col][i])
                            P12.append(btl_data[p_btl_col][k])
                            T12.append(btl_data[t_btl_col][k])
                            C12.append(btl_data[c_btl_col][k])
                            calib12.append(btl_data[calib12_btl_col][k])
                            ctd_d12.append(d12)
                        else:
                            report_ctd.report_quality_flags(qualfilePath,
                                                        filename_base, ref_btl[btl_num_col][i],
                                                        btl_data[p_btl_col][k], param,
                                                        d1, d2, d12)

        # Get Range
        if args.xRange:
            x0 = int(args.xRange.split(":")[0])
            x1 = int(args.xRange.split(":")[1])
        else:
            if args.primary:
                x0 = min(calib1)
                x1 = max(calib1)
            elif args.secondary:
                x0 = min(calib2)
                x1 = max(calib2)

        if args.primary:
            sensor = 1
            ctd_d = ctd_d1
            calib = calib1
            P = P1
        elif args.secondary:
            sensor = 2
            ctd_d = ctd_d2
            calib = calib2
            P = P2

        calib_indx = np.argsort(calib)
        for i in calib_indx:
            if (calib[i] >= x0) and (calib[i] <= x1):
                indx.append(i)
        calib = [calib[i] for i in indx]
        fit = np.arange(x0,x1,(x1-x0)/50)
        P = [P[i] for i in indx]
        ctd_d = [ctd_d[i] for i in indx]

        mean = np.mean(ctd_d)
        std = np.std(ctd_d)
        lstd = mean - std
        hstd = mean + std
        cf = np.polyfit(calib, ctd_d, order)

        indx12 = [index for index, value in enumerate(calib12)
                  if (value >= x0) and (value <= x1)]
        calib12 = [calib12[i] for i in indx12]
        P12 = [P12[i] for i in indx12]
        T12 = [T12[i] for i in indx12]
        C12 = [C12[i] for i in indx12]
        ref_data12 = [ref_data12[i] for i in indx12]
        ctd_d12 = [ctd_d12[i] for i in indx12]

        mean12 = np.mean(ctd_d12)
        std12 = np.std(ctd_d12)
        lstd12 = mean12 - std12
        hstd12 = mean12 + std12

        if args.temperature:
            sensor = '_t'+str(sensor)
            if order is 0:
                coef[4] = cf[0]
            elif (order is 1) and ('P' in args.calib):
                coef[1] = cf[0]
                coef[4] = cf[1]
            elif (order is 2) and ('P' in args.calib):
                coef[0] = cf[0]
                coef[1] = cf[1]
                coef[4] = cf[2]
            elif (order is 1) and ('T' in args.calib):
                coef[3] = cf[0]
                coef[4] = cf[1]
            elif (order is 2) and ('T' in args.calib):
                coef[2] = cf[0]
                coef[3] = cf[1]
                coef[4] = cf[2]

            Y = fit_ctd.temperature_polyfit(coef, fit, np.full(len(fit), 0.0))

        elif args.conductivity:
            sensor = '_c'+str(sensor)
            if order is 0:
                coef[6] = cf[0]
            elif (order is 1) and ('P' in args.calib):
                coef[1] = cf[0]
                coef[6] = cf[1]
            elif (order is 2) and ('P' in args.calib):
                coef[0] = cf[0]
                coef[1] = cf[1]
                coef[6] = cf[2]
            elif (order is 1) and ('T' in args.calib):
                coef[3] = cf[0]
                coef[6] = cf[1]
            elif (order is 2) and ('T' in args.calib):
                coef[2] = cf[0]
                coef[3] = cf[1]
                coef[6] = cf[2]
            elif (order is 1) and ('C' in args.calib):
                coef[5] = cf[0]
                coef[6] = cf[1]
            elif (order is 2) and ('C' in args.calib):
                coef[4] = cf[0]
                coef[5] = cf[1]
                coef[6] = cf[2]

            Y = fit_ctd.conductivity_polyfit(coef, fit, fit, np.full(len(fit), 0.0))


        fitfile = str('fitting'+sensor+'.' + FILE_EXT)
        fitfilePath = os.path.join(log_directory, fitfile)
        report_ctd.report_polyfit(coef, file_base_arr, fitfilePath)

        # plt.scatter(calib, ctd_d, c=P, cmap=plt.cm.gist_rainbow)
        # plt.plot(fit, Y, color='red')
        # pylab.axhline( mean, color="black", label='mean')
        # pylab.axhline( lstd, color="black", linestyle="--", label='std')
        # pylab.axhline( hstd, color="black", linestyle="--")
        # plt.axis()
        # plt.show()
        #
        # plt.scatter(calib12, ctd_d12, c=P12, cmap=plt.cm.gist_rainbow)
        # pylab.axhline( mean12, color="black", label='mean')
        # pylab.axhline( lstd12, color="black", linestyle="--", label='std')
        # pylab.axhline( hstd12, color="black", linestyle="--")
        # plt.axis()
        # plt.show()

    else:
        errPrint('ERROR: Input cast range:', args.castRange, 'not found\n')
        sys.exit(1)


    debugPrint('Done!')

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

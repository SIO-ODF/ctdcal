#! /usr/bin/env python
import sys
import os
import math
import argparse
import fnmatch
import numpy as np
import pandas as pd
import json
import libODF_process_ctd as process_ctd
import libODF_report_ctd as report_ctd
import libODF_fit_ctd as fit_ctd
import configparser
#import matplotlib.pyplot as plt
from scipy.optimize import leastsq

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

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
    parser.add_argument('castRange', metavar='cast_range', help='range of cast values for shipboard calibration.')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # raw output
    parser.add_argument('-press', '--pressure', action='store_true', help='Calibration pressure for fitting.')

    # raw output
    parser.add_argument('-temp', '--temperature', action='store_true', help='')

    # raw output
    parser.add_argument('-cond', '--conductivity', action='store_true', help='')

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
    c1_col = config['analytical_inputs']['c1']
    c1_btl_col = config['inputs']['c1']
    c2_col = config['analytical_inputs']['c2']
    c2_btl_col = config['inputs']['c2']
    sal_col = config['analytical_inputs']['salt']
    sal_btl_col = config['inputs']['salt']
    dov_col = config['inputs']['dov']
    btl_sal_col = config['analytical_inputs']['btl_salt']
    timedate = config['analytical_inputs']['datetime']
    lat_col = config['analytical_inputs']['lat']
    lon_col = config['analytical_inputs']['lon']
    reft_col = config['inputs']['reft']
    btl_num_col = config['inputs']['btl_num']
    
    file_base_arr = []

    logfileName = str('ondeck_pressure' + '.' + FILE_EXT)
    logfilePath = os.path.join(log_directory, logfileName)
    log_data = process_ctd.dataToNDarray(logfilePath,str,None,',',None)
    p_start = []
    p_end = []

    if args.castRange:
        casts = str.split(args.castRange, '-')
        if int(casts[0][2]) is int(casts[1][2]):
            regEx_fileList = '['+casts[0][0]+'-'+casts[1][0]+']['+casts[0][1]+'-'+casts[1][1]+']?0?'+TIME_SUFFIX+'.'+FILE_EXT
        else:
            regEx_fileList = '['+casts[0][0]+'-'+casts[1][0]+']['+casts[0][1]+'-'+casts[1][1]+']['+casts[0][2]+'-'+casts[1][2]+']0?'+TIME_SUFFIX+'.'+FILE_EXT

        files = [f for f in os.listdir(time_directory) if fnmatch.fnmatch(f, regEx_fileList)]

        for f in files:
            # Used later for building output file names
            filename_base = str.split(f,'_')[0] # original filename w/o ext
            file_base_arr.append(filename_base)

        # If Pressure fit selected
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
 
        # If Temperature fit selected
        if args.temperature:

            btl_num = []
            P1 = []
            P2 = []
            P12 = []
            reft1 = []
            reft2 = []
            reft12 = []
            ctd_btl_t1 = [] 
            ctd_btl_t2 = [] 
            ctd_btl_dt1 = [] 
            ctd_btl_dt2 = [] 
            ctd_btl_dt1t2 = [] 
            for filename_base in file_base_arr:
                reftfileName = (filename_base + REFT_SUFFIX + '.' + FILE_EXT) 
                reftfilePath = os.path.join(reft_directory, reftfileName)
                if os.path.isfile(reftfilePath):
                    print("Processing "+ filename_base)
                    reft_btl = process_ctd.dataToNDarray(reftfilePath,float,True,',',None)
                else: 
                    continue

                btlfileName = str(filename_base + BTL_SUFFIX + MEAN_SUFFIX + '.' + FILE_EXT)
                btlfilePath = os.path.join(btl_directory, btlfileName)
                if os.path.isfile(btlfilePath):
                    btl_data = process_ctd.dataToNDarray(btlfilePath,float,True,',',None)
                    btl_data = btl_data[:][1:]
                else:
                    print("Missing file: "+btlfilePath)
                    break

                timefileName = str(filename_base + TIME_SUFFIX + '.' + FILE_EXT)
                timefilePath = os.path.join(time_directory, timefileName)
                if os.path.isfile(timefilePath):
                    time_data = process_ctd.dataToNDarray(timefilePath,float,True,',',1)
                    time_data = time_data[:][1:]
                else:
                    print("Missing file: "+timefilePath)
                    break

                # Find Isopycnal Down Trace Bottle Trip Equivalent
                # Need to add density (sigma_theta) driven component to this search
                end = np.argmax(time_data[p_col])
                down_trace_btl = fit_ctd.find_isopycnals(p_btl_col, t1_btl_col, sal_btl_col, dov_col, btl_data, time_data[p_col], time_data[t1_col], time_data[sal_col], time_data[dov_col])

                for i in range(0,len(reft_btl[btl_num_col])):
                    j = reft_btl[btl_num_col][i] 
                    k = int(np.where(btl_data[btl_num_col] == j)[0][0])
                    dt1 = reft_btl[reft_col][i] - down_trace_btl[t1_btl_col][k]
                    dt2 = reft_btl[reft_col][i] - down_trace_btl[t2_btl_col][k]
                    dt1t2 = down_trace_btl[t1_btl_col][k] - down_trace_btl[t2_btl_col][k]
                    if down_trace_btl[p_btl_col][k] > 2000:
                        if abs(dt1) < 0.002:
                            reft1.append(reft_btl[reft_col][i])
                            P1.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t1.append(down_trace_btl[t1_btl_col][k])
                            ctd_btl_dt1.append(dt1)
                        if abs(dt2) < 0.002:
                            reft2.append(reft_btl[reft_col][i])
                            P2.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t2.append(down_trace_btl[t2_btl_col][k])
                            ctd_btl_dt2.append(dt2)
                        if abs(dt1t2) < 0.002:
                            reft12.append(reft_btl[reft_col][i])
                            P12.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_dt1t2.append(dt1t2)
                    elif (down_trace_btl[p_btl_col][k] < 2000) and (down_trace_btl[p_btl_col][k] > 1000):
                        if abs(dt1) < 0.005:
                            reft1.append(reft_btl[reft_col][i])
                            P1.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t1.append(down_trace_btl[t1_btl_col][k])
                            ctd_btl_dt1.append(dt1)
                        if abs(dt2) < 0.005:
                            reft2.append(reft_btl[reft_col][i])
                            P2.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t2.append(down_trace_btl[t2_btl_col][k])
                            ctd_btl_dt2.append(dt2)
                        if abs(dt1t2) < 0.005:
                            reft12.append(reft_btl[reft_col][i])
                            P12.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_dt1t2.append(dt1t2)
                    elif (down_trace_btl[p_btl_col][k] < 1000) and (down_trace_btl[p_btl_col][k] > 500):
                        if abs(dt1) < 0.010:
                            reft1.append(reft_btl[reft_col][i])
                            P1.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t1.append(down_trace_btl[t1_btl_col][k])
                            ctd_btl_dt1.append(dt1)
                        if abs(dt2) < 0.010:
                            reft2.append(reft_btl[reft_col][i])
                            P2.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t2.append(down_trace_btl[t2_btl_col][k])
                            ctd_btl_dt2.append(dt2)
                        if abs(dt1t2) < 0.010:
                            reft12.append(reft_btl[reft_col][i])
                            P12.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_dt1t2.append(dt1t2)
                    elif btl_data[p_btl_col][k] < 500:
                        if abs(dt1) < 0.020:
                            reft1.append(reft_btl[reft_col][i])
                            P1.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t1.append(down_trace_btl[t1_btl_col][k])
                            ctd_btl_dt1.append(dt1)
                        if abs(dt2) < 0.020:
                            reft2.append(reft_btl[reft_col][i])
                            P2.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_t2.append(down_trace_btl[t2_btl_col][k])
                            ctd_btl_dt2.append(dt2)
                        if abs(dt1t2) < 0.020:
                            reft12.append(reft_btl[reft_col][i])
                            P12.append(down_trace_btl[p_btl_col][k])
                            ctd_btl_dt1t2.append(dt1t2)


            # Get descrete ref temp data 
            coef1 = fit_ctd.find_temp_coef(reft1, P1, ctd_btl_t1)
            coef2 = fit_ctd.find_temp_coef(reft2, P2, ctd_btl_t2)

            fitfileT1 = str('fitting_t1.' + FILE_EXT)
            fitfileT1Path = os.path.join(log_directory, fitfileT1)
            fitfileT2 = str('fitting_t2.' + FILE_EXT)
            fitfileT2Path = os.path.join(log_directory, fitfileT2)

            report_ctd.report_polyfit(coef1, file_base_arr, fitfileT1Path)
            report_ctd.report_polyfit(coef2, file_base_arr, fitfileT2Path)

 
        # If Temperature fit selected
        if args.conductivity:

            btl_num = []
            P1 = []
            P2 = []
            P12 = []
            T1 = []
            T2 = []
            T12 = []
            cond_btl_col = 'BTLCOND'
            cond1 = []
            cond2 = []
            cond12 = []
            salt1 = []
            salt2 = []
            salt12 = []
            ctd_btl_c1 = [] 
            ctd_btl_c2 = [] 
            ctd_btl_dc1 = [] 
            ctd_btl_dc2 = [] 
            ctd_btl_dc1c2 = [] 
            for filename_base in file_base_arr:
                btlfileName = str(filename_base + BTL_SUFFIX + MEAN_SUFFIX + '.' + FILE_EXT)
                btlfilePath = os.path.join(btl_directory, btlfileName)
                if os.path.isfile(btlfilePath):
                    btl_data = process_ctd.dataToNDarray(btlfilePath,float,True,',',None)
                    btl_data = btl_data[:][1:]
                else:
                    print("Missing file: "+btlfilePath)
                    continue

                saltfileName = filename_base
                saltfilePath = os.path.join(salt_directory, saltfileName)
                if os.path.isfile(saltfilePath):
                    print("Processing "+ filename_base)
                    cond_btl, salt_btl = fit_ctd.salt_calc(saltfilePath,btl_num_col,t1_btl_col,p_btl_col,btl_data)
                else: 
                    continue

                timefileName = str(filename_base + TIME_SUFFIX + '.' + FILE_EXT)
                timefilePath = os.path.join(time_directory, timefileName)
                if os.path.isfile(timefilePath):
                    time_data = process_ctd.dataToNDarray(timefilePath,float,True,',',1)
                    time_data = time_data[:][1:]
                else:
                    print("Missing file: "+timefilePath)
                    break

                # Find Isopycnal Down Trace Bottle Trip Equivalent
                # Need to add density (sigma_theta) driven component to this search
                end = np.argmax(time_data[p_col])
                down_trace_btl = fit_ctd.find_isopycnals(p_btl_col, t1_btl_col, sal_btl_col, dov_col, btl_data, time_data[p_col], time_data[t1_col], time_data[sal_col], time_data[dov_col])

                for i in range(0,len(cond_btl[btl_num_col])):
                    j = cond_btl[btl_num_col][i] 
                    if j != 0:
                        k = int(np.where(btl_data[btl_num_col] == j)[0][0])
                        dc1 = cond_btl[cond_btl_col][i] - down_trace_btl[c1_btl_col][k]
                        dc2 = cond_btl[cond_btl_col][i] - down_trace_btl[c2_btl_col][k]
                        dc1c2 = down_trace_btl[c1_btl_col][k] - down_trace_btl[c2_btl_col][k]
                        if down_trace_btl[p_btl_col][k] > 2000:
                            if abs(dc1) < 0.002:
                                cond1.append(cond_btl[cond_btl_col][i])
                                P1.append(down_trace_btl[p_btl_col][k])
                                T1.append(down_trace_btl[t1_btl_col][k])
                                ctd_btl_c1.append(down_trace_btl[c1_btl_col][k])
                                ctd_btl_dc1.append(dc1)
                            if abs(dc2) < 0.002:
                                cond2.append(cond_btl[cond_btl_col][i])
                                P2.append(down_trace_btl[p_btl_col][k])
                                T2.append(down_trace_btl[t2_btl_col][k])
                                ctd_btl_c2.append(down_trace_btl[c2_btl_col][k])
                                ctd_btl_dc2.append(dc2)
                            if abs(dc1c2) < 0.002:
                                cond12.append(cond_btl[cond_btl_col][i])
                                P12.append(down_trace_btl[p_btl_col][k])
                                T12.append(down_trace_btl[t_btl_col][k])
                                ctd_btl_dc1c2.append(dc1c2)
                        elif (down_trace_btl[p_btl_col][k] < 2000) and (down_trace_btl[p_btl_col][k] > 1000):
                            if abs(dc1) < 0.005:
                                cond1.append(cond_btl[cond_btl_col][i])
                                P1.append(down_trace_btl[p_btl_col][k])
                                T1.append(down_trace_btl[t1_btl_col][k])
                                ctd_btl_c1.append(down_trace_btl[c1_btl_col][k])
                                ctd_btl_dc1.append(dc1)
                            if abs(dc2) < 0.005:
                                cond2.append(cond_btl[cond_btl_col][i])
                                P2.append(down_trace_btl[p_btl_col][k])
                                T2.append(down_trace_btl[t2_btl_col][k])
                                ctd_btl_c2.append(down_trace_btl[c2_btl_col][k])
                                ctd_btl_dc2.append(dc2)
                            if abs(dc1c2) < 0.005:
                                cond12.append(cond_btl[cond_btl_col][i])
                                P12.append(down_trace_btl[p_btl_col][k])
                                T12.append(down_trace_btl[t_btl_col][k])
                                ctd_btl_dc1c2.append(dc1c2)
                        elif (down_trace_btl[p_btl_col][k] < 1000) and (down_trace_btl[p_btl_col][k] > 500):
                            if abs(dc1) < 0.010:
                                cond1.append(cond_btl[cond_btl_col][i])
                                P1.append(down_trace_btl[p_btl_col][k])
                                T1.append(down_trace_btl[t1_btl_col][k])
                                ctd_btl_c1.append(down_trace_btl[c1_btl_col][k])
                                ctd_btl_dc1.append(dc1)
                            if abs(dc2) < 0.010:
                                cond2.append(cond_btl[cond_btl_col][i])
                                P2.append(down_trace_btl[p_btl_col][k])
                                T2.append(down_trace_btl[t2_btl_col][k])
                                ctd_btl_c2.append(down_trace_btl[c2_btl_col][k])
                                ctd_btl_dc2.append(dc2)
                            if abs(dc1c2) < 0.010:
                                cond12.append(cond_btl[cond_btl_col][i])
                                P12.append(down_trace_btl[p_btl_col][k])
                                T12.append(down_trace_btl[t_btl_col][k])
                                ctd_btl_dc1c2.append(dc1c2)
                        elif btl_data[p_btl_col][k] < 500:
                            if abs(dc1) < 0.020:
                                cond1.append(cond_btl[cond_btl_col][i])
                                P1.append(down_trace_btl[p_btl_col][k])
                                T1.append(down_trace_btl[t1_btl_col][k])
                                ctd_btl_c1.append(down_trace_btl[c1_btl_col][k])
                                ctd_btl_dc1.append(dc1)
                            if abs(dc2) < 0.020:
                                cond2.append(cond_btl[cond_btl_col][i])
                                P2.append(down_trace_btl[p_btl_col][k])
                                T2.append(down_trace_btl[t2_btl_col][k])
                                ctd_btl_c2.append(down_trace_btl[c2_btl_col][k])
                                ctd_btl_dc2.append(dc2)
                            if abs(dc1c2) < 0.020:
                                cond12.append(cond_btl[cond_btl_col][i])
                                P12.append(down_trace_btl[p_btl_col][k])
                                T12.append(down_trace_btl[t_btl_col][k])
                                ctd_btl_dc1c2.append(dc1c2)


            # Get descrete ref temp data 
            coef1 = fit_ctd.find_cond_coef(cond1, P1, T1, ctd_btl_c1)
            coef2 = fit_ctd.find_cond_coef(cond2, P2, T2, ctd_btl_c2)

            fitfileC1 = str('fitting_c1.' + FILE_EXT)
            fitfileC1Path = os.path.join(log_directory, fitfileC1)
            fitfileC2 = str('fitting_c2.' + FILE_EXT)
            fitfileC2Path = os.path.join(log_directory, fitfileC2)

            report_ctd.report_polyfit(coef1, file_base_arr, fitfileC1Path)
            report_ctd.report_polyfit(coef2, file_base_arr, fitfileC2Path)
    else:
        errPrint('ERROR: Input cast range:', args.castRange, 'not found\n')
        sys.exit(1)

    #plt.scatter(P1, ctd_btl_dc1, c=P1, cmap=plt.cm.gist_rainbow, label='SBE35 - T1')
    #plt.plot(P1, time_data[p_col], color='b', label='raw')
    #plt.plot(pressure_seq_data[dopl_col], pressure_seq_data[p_col], color='r', label='raw')
    #plt.plot(o2pl_btl['OXYGEN'], btl_data[p_btl_col], color='g', marker='o', label='raw')
    #plt.plot(pressure_seq_data[do_col], pressure_seq_data[p_col], color='b', label='raw')
    #plt.gca().invert_yaxis()
    #plt.axis()
    #plt.show()

    debugPrint('Done!')

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])

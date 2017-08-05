#! /usr/bin/env python
import sys
import os
import argparse
import numpy as np
import pandas as pd
import json
import libODF_process_ctd as process_ctd
import libODF_report_ctd as report_ctd
import libODF_fit_ctd as fit_ctd
import configparser
#import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import gsw

DEBUG = False

#File extension to use for output files (csv-formatted)
FILE_EXT = 'csv'

#File extension to use for output files (csv-formatted)
XML_EXT = 'XMLCON'

#File extension to use for output files (csv-formatted)
HEX_EXT = 'hex'

#File extension to use for raw output
BTL_SUFFIX = '_btl'

#File extension to use for raw output
FIT_SUFFIX = '_fit'

#File extension to use for raw output
MEAN_SUFFIX = '_mean'

#File extension to use for raw output
TIME_SUFFIX = '_time'

#File extension to use for output files (csv-formatted)
UPDT_SUFFIX = '_updt'

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
    parser.add_argument('timeFile', metavar='time_file', help='the .csv data file to fit by bottle data')

    # debug messages
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    # raw output
    parser.add_argument('-pres', '--pressure', action='store_true', help='Fit pressure data')

    # raw output
    parser.add_argument('-temp', '--temperature', action='store_true', help='Fit temperture data')

    # raw output
    parser.add_argument('-cond', '--conductivity', action='store_true', help='Fit conductivity data')

    # raw output
    parser.add_argument('-oxy', '--oxygen', type=argparse.FileType('r'), help='return the oxygen data file')

    # Process Command-line args
    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    # Verify hex file exists
    if not os.path.isfile(args.timeFile):
        errPrint('ERROR: Input time dependent .csv file:', args.timeFile, 'not found\n')
        sys.exit(1)

    # Used later for building output file names
    filename_ext = os.path.basename(args.timeFile) # original filename with ext
    filename_base = os.path.splitext(filename_ext)[0] # original filename w/o ext

    if '_' in filename_base:
        filename_base = filename_base.split('_')[0]

    #Import Cruise Configuration File
    iniFile = 'data/ini-files/configuration.ini'
    config = configparser.RawConfigParser()
    config.read(iniFile)

    #Initialise Configuration Parameters
    expocode = config['cruise']['expocode']
    sectionID = config['cruise']['sectionid']
    raw_directory = config['ctd_processing']['raw_data_directory']
    time_directory = config['ctd_processing']['time_data_directory']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    oxygen_directory = config['ctd_processing']['oxygen_directory']
    btl_directory = config['ctd_processing']['bottle_directory']
    o2flask_file = config['ctd_processing']['o2flask_file']
    log_directory = config['ctd_processing']['log_directory']
    sample_rate = config['ctd_processing']['sample_rate']
    search_time = config['ctd_processing']['roll_filter_time']
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
    btl_sal_col = config['analytical_inputs']['btl_salt']
    dov_col = config['analytical_inputs']['dov']
    dov_btl_col = config['inputs']['dov']
    dopl_col = config['analytical_inputs']['dopl']
    dopl_btl_col = config['inputs']['dopl']
    dopkg_col = config['analytical_inputs']['dopkg']
    btl_oxy_col = config['analytical_inputs']['btl_oxy']
    xmis_col = config['analytical_inputs']['xmis']
    fluor_col = config['analytical_inputs']['fluor']
    timedate = config['analytical_inputs']['datetime']
    lat_col = config['analytical_inputs']['lat']
    lat_btl_col = config['inputs']['lat']
    lon_col = config['analytical_inputs']['lon']
    lon_btl_col = config['inputs']['lon']
    reft_col = config['inputs']['reft']
    btl_num_col = config['inputs']['btl_num']

    #time_column_data = config['time_series_output']['data_names'].split(',')
    time_column_data = config['time_series_output']['data_output']
    time_column_names = config['time_series_output']['column_name'].split(',')
    time_column_units = config['time_series_output']['column_units'].split(',')
    time_column_format = config['time_series_output']['format']

    #pressure_column_data = config['time_series_output']['data_names'].split(',')
    p_column_data = config['pressure_series_output']['data'].split(',')
    p_column_names = config['pressure_series_output']['column_name'].split(',')
    p_column_units = config['pressure_series_output']['column_units'].split(',')
    p_column_format = config['pressure_series_output']['format']
    p_column_qual = config['pressure_series_output']['qual_columns'].split(',')
    p_column_one = list(config['pressure_series_output']['q1_columns'].split(','))

    #bottle_data outputs
    btl_dtype = config['bottle_series_output']['dtype']

    hexfileName = str(filename_base + '.' + HEX_EXT)
    hexfilePath = os.path.join(raw_directory, hexfileName)

    xmlfileName = str(filename_base + '.' + XML_EXT)
    xmlfilePath = os.path.join(raw_directory, xmlfileName)

    outtimefileName = str(filename_base + TIME_SUFFIX + '.' + FILE_EXT)
    outtimefilePath = os.path.join(time_directory, outtimefileName)

    pressfileName = str(filename_base + FIT_SUFFIX + '.' + FILE_EXT)
    pressfilePath = os.path.join(pressure_directory, pressfileName)

    btlfileName = str(filename_base + BTL_SUFFIX + MEAN_SUFFIX + '.' + FILE_EXT)
    btlfilePath = os.path.join(btl_directory, btlfileName)

    # Get bottle data
    btl_data = process_ctd.dataToNDarray(btlfilePath,float,True,',',0)
    btl_data = btl_data[:][1:]

    # Get procesed time data
    time_data = process_ctd.dataToNDarray(args.timeFile,float,True,',',1)
    btm = np.argmax(time_data[p_col][1:])
    time_data = time_data[:][1:btm]

    if args.pressure:
        print('In -pres flag fit condition')
        print(filename_base)
        pfileName = str('poffset' + '.' + FILE_EXT)
        pfilePath = os.path.join(log_directory, pfileName)
        poff_data = process_ctd.dataToNDarray(pfilePath,str,None,',',None)

        for line in poff_data:
            if filename_base in line[0]:
                for val in line:
                    if 'offset' in val:
                        offset = float(str.split(val, ':')[1])
            continue

        # Pressure offset
        btl_data[p_btl_col] = fit_ctd.offset(offset, btl_data[p_btl_col])
        time_data[p_col] = fit_ctd.offset(offset, time_data[p_col])
        # End pressure if condition

    if args.temperature:
        print('In -temp flag fit condition')
        print(filename_base)
        coef1 = []
        coef2 = []
        # Get descrete ref temp data
        t1fileName = str('fitting_t1' + '.' + FILE_EXT)
        t1filePath = os.path.join(log_directory, t1fileName)
        t1_coef = process_ctd.dataToNDarray(t1filePath,str,None,',',None)

        for line in t1_coef:
            if filename_base in line[0]:
                val = line[1:]
                for i in range(0,len(val)):
                    if 'coef' in val[i]:
                        coef1.append(float(str.split(val[i], ':')[1]))
                    else:
                        coef1.append(float(val[i]))
            continue
        btl_data[t1_btl_col] = fit_ctd.temperature_polyfit(coef1, btl_data[p_btl_col], btl_data[t1_btl_col])
        time_data[t1_col] = fit_ctd.temperature_polyfit(coef1, time_data[p_col], time_data[t1_col])

        t2fileName = str('fitting_t2' + '.' + FILE_EXT)
        t2filePath = os.path.join(log_directory, t2fileName)
        t2_coef = process_ctd.dataToNDarray(t2filePath,str,None,',',None)

        for line in t2_coef:
            if filename_base in line[0]:
                val = line[1:]
                for i in range(0,len(val)):
                    if 'coef' in val[i]:
                        coef2.append(float(str.split(val[i], ':')[1]))
                    else:
                        coef2.append(float(val[i]))
            continue
        btl_data[t2_btl_col] = fit_ctd.temperature_polyfit(coef2, btl_data[p_btl_col], btl_data[t2_btl_col])
        time_data[t2_col] = fit_ctd.temperature_polyfit(coef2, time_data[p_col], time_data[t2_col])

    if args.conductivity:
        print('In -cond flag fit condition')
        print(filename_base)
        coef1 = []
        coef2 = []

        # Get descrete cond data

        c1fileName = str('fitting_c1' + '.' + FILE_EXT)
        c1filePath = os.path.join(log_directory, c1fileName)
        if os.path.exists(c1filePath):
            c1_coef = process_ctd.dataToNDarray(c1filePath,str,None,',',None)

            for line in c1_coef:
                if filename_base in line[0]:
                    val = line[1:]
                    for i in range(0,len(val)):
                        if 'coef' in val[i]:
                            coef1.append(float(str.split(val[i], ':')[1]))
                        else:
                            coef1.append(float(val[i]))
                continue
            btl_data[c1_btl_col] = fit_ctd.conductivity_polyfit(coef1, btl_data[p_btl_col], btl_data[t1_btl_col], btl_data[c1_btl_col])
            time_data[c1_col] = fit_ctd.conductivity_polyfit(coef1, time_data[p_col], time_data[c1_col], time_data[c1_col])

        c2fileName = str('fitting_c2' + '.' + FILE_EXT)
        c2filePath = os.path.join(log_directory, c2fileName)
        if os.path.exists(c2filePath):
            c2_coef = process_ctd.dataToNDarray(c2filePath,str,None,',',None)

            for line in c2_coef:
                if filename_base in line[0]:
                    val = line[1:]
                    for i in range(0,len(val)):
                        if 'coef' in val[i]:
                            coef2.append(float(str.split(val[i], ':')[1]))
                        else:
                            coef2.append(float(val[i]))
                continue
            btl_data[c2_btl_col] = fit_ctd.conductivity_polyfit(coef2, btl_data[p_btl_col], btl_data[t2_btl_col], btl_data[c2_btl_col])
            time_data[c2_col] = fit_ctd.conductivity_polyfit(coef2, time_data[p_col], time_data[c2_col], time_data[c2_col])

        time_data[sal_col] = gsw.SP_from_C(time_data[c_col],time_data[t_col],time_data[p_col])

    if args.oxygen:
        print('In -oxy flag fit condition')
        print(filename_base)
        # Get Analytical Oxygen data
        o2pkg_btl, o2pl_btl = fit_ctd.o2_calc(o2flask_file,args.oxygen.name,btl_data[btl_num_col],btl_data[sal_btl_col])

        kelvin = []
        for i in range(0,len(time_data[t_col])):
            kelvin.append(time_data[t_col][i] + 273.15)

        # Find New Oxygen Coef
        oxy_coef = fit_ctd.find_oxy_coef(o2pl_btl['OXYGEN'], btl_data[p_btl_col], btl_data[t_btl_col], btl_data[sal_btl_col], btl_data[dov_btl_col], hexfilePath, xmlfilePath)

        # Convert CTD Oxygen Voltage Data with New DO Coef
        time_data[dopl_col] = fit_ctd.oxy_dict(oxy_coef, time_data[p_col], kelvin, time_data[t_col], time_data[sal_col], time_data[dov_col])
        # End oxygen flag fitting if condition

    # Find Isopycnal Down Trace Bottle Trip Equivalent
    # Write bottle data to file
    report_ctd.report_btl_data(btlfilePath, btl_dtype, btl_data)

    # Write time data to file
    report_ctd.report_time_series_data(filename_base, time_directory, expocode, time_column_names, time_column_units, time_column_names, time_column_format, time_data)

    # Pressure Sequence
    pressure_seq_data = process_ctd.pressure_sequence(filename_base, p_col, timedate, 2.0, -1.0, 0.0, 'down', int(sample_rate), int(search_time), time_data)

    # Convert dissolved oxygen from ml/l to umol/kg
    dopkg = process_ctd.o2pl2pkg(p_col, t1_col, sal_col, dopl_col, dopkg_col, lat_col, lon_col, pressure_seq_data)

    # Add quality codes to data
    qual_pseq_data = process_ctd.ctd_quality_codes(dopkg_col, None, None, True, p_column_qual, p_column_one, pressure_seq_data)

    # Collect Cast Details from Log
    logfileName = str('cast_details' + '.' + FILE_EXT)
    logfilePath = os.path.join(log_directory, logfileName)

    cast_details = process_ctd.dataToNDarray(logfilePath,str,None,',',0)
    for line in cast_details:
        if filename_base in line[0]:
            for val in line:
                if 'at_depth' in val: btime = float(str.split(val, ':')[1])
                if 'latitude' in val: btm_lat = float(str.split(val, ':')[1])
                if 'longitude' in val: btm_lon = float(str.split(val, ':')[1])
                if 'altimeter_bottom' in val: btm_alt = float(str.split(val, ':')[1])
            break

    # Write time data to file
    depth = -999
    #import pdb; pdb.set_trace()
    report_ctd.report_pressure_series_data(filename_base, expocode, sectionID, btime, btm_lat, btm_lon, depth, btm_alt, ctd, pressure_directory, p_column_names, p_column_units, p_column_data, qual_pseq_data, dopkg, pressure_seq_data)

    #plt.plot(o2pl_btl['OXYGEN'], btl_data[p_btl_col], color='b', marker='o')
    #plt.plot(tmpo2, time_data[p_col], color='g', label='raw')
    #plt.plot(time_data[dopl_col], time_data[p_col], color='b', label='raw')
    #plt.plot(pressure_seq_data[dopl_col], pressure_seq_data[p_col], color='r', label='raw')
    #plt.plot(pressure_seq_data[t2_col], pressure_seq_data[p_col], color='r', label='raw')
    #plt.plot(btl_data[t_btl_col], time_data[p_col], color='r', label='raw')
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

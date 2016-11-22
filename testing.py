#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as sig
import configparser 
import libODF_process_ctd as process_ctd 

# todo: Migrate testing.py to run.py next 

# Import configuration file
config = configparser.ConfigParser()
config.read('configuration.ini')
rdir = config['ctd_processing']['raw_data_directory']
infile = rdir+config['ctd_processing']['raw_data_file']

# Initialize config parameters
tc1_align = config['ctd_processing']['TC_primary_align']
tc2_align = config['ctd_processing']['TC_secondary_align']
do_align = config['ctd_processing']['DO_align']
sample_rate = config['ctd_processing']['sample_rate']
search_time = config['ctd_processing']['roll_filter_time']

# Data input to format
dframe = process_ctd.dataToDataFrame(infile)
dtlist = np.dtype(list(dframe.columns.values))
dmatrix = process_ctd.dataToMatrix(infile,dtlist,',') 

# Align dependent sensor/state var data 
if tc1_align: dmatrix = process_ctd.ctd_align(dmatrix,'C1mScm', float(tc1_align))
if tc2_align: dmatrix = process_ctd.ctd_align(dmatrix,'C2mScm', float(tc2_align))
if do_align: dmatrix = process_ctd.ctd_align(dmatrix,'SBE43FV', float(do_align))
dmatrix = process_ctd.ondeck_pressure(dmatrix, float(config['ctd_processing']['conductivity_start'])) 

# Filter data
# todo: filtered_matrix = process_ctd.raw_ctd_filter(dmatrix, 'triangle', 24)
# Use dtype assertion within filtering function to remove explicit column call here
# Ref: press_sequence func
dmatrix['Pdbar'] = process_ctd.raw_ctd_filter(dmatrix['Pdbar'], 'triangle',24)
dmatrix['T1C'] = process_ctd.raw_ctd_filter(dmatrix['T1C'], 'triangle',24)
dmatrix['T2C'] = process_ctd.raw_ctd_filter(dmatrix['T2C'], 'triangle',24)
dmatrix['C1mScm'] = process_ctd.raw_ctd_filter(dmatrix['C1mScm'], 'triangle',24)
dmatrix['C2mScm'] = process_ctd.raw_ctd_filter(dmatrix['C2mScm'], 'triangle',24)
dmatrix['SBE43FV'] = process_ctd.raw_ctd_filter(dmatrix['SBE43FV'], 'triangle',24)
dmatrix['SBE43UuMolKg'] = process_ctd.raw_ctd_filter(dmatrix['SBE43UuMolKg'], 'triangle',24)
dmatrix['V0V'] = process_ctd.raw_ctd_filter(dmatrix['V0V'], 'triangle',24)
dmatrix['V1V'] = process_ctd.raw_ctd_filter(dmatrix['V1V'], 'triangle',24)
dmatrix['V2V'] = process_ctd.raw_ctd_filter(dmatrix['V2V'], 'triangle',24)
dmatrix['V3V'] = process_ctd.raw_ctd_filter(dmatrix['V3V'], 'triangle',24)
dmatrix['V4V'] = process_ctd.raw_ctd_filter(dmatrix['V4V'], 'triangle',24)
dmatrix['V5V'] = process_ctd.raw_ctd_filter(dmatrix['V5V'], 'triangle',24)
dmatrix['V6V'] = process_ctd.raw_ctd_filter(dmatrix['V6V'], 'triangle',24)
dmatrix['V7V'] = process_ctd.raw_ctd_filter(dmatrix['V7V'], 'triangle',24)

# Cast Details 
stime, etime, btime, startP, maxP, dmatrix = process_ctd.cast_details(dmatrix)
# Pressure Sequence
pressure_matrix = process_ctd.pressure_sequence(dmatrix, 2.0, stime, startP, 'down', int(sample_rate), int(search_time))

#plottitle = cruise+': Filter' 
#plotfile = maindir+'BoxCar/bc.00101.png' 
plt.plot(pressure_matrix['T1C'], pressure_matrix['Pdbar'], color='r', label='Temp')
plt.plot(pressure_matrix['C1mScm'], pressure_matrix['Pdbar'], color='y', label='Cond')
plt.plot(pressure_matrix['SBE43UuMolKg'], pressure_matrix['Pdbar'], color='g', label='DOxy')
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

import sys
import os
import json

import numpy as np
import pandas as pd
import scipy as sp
import math
import pylab

import libODF_sbe_reader as sbe_reader
import libODF_convert as cnv
import libODF_process_ctd as process_ctd
import libODF_report_ctd as report_ctd
import libODF_fit_ctd as fit_ctd

import gsw

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
log_directory = config['ctd_processing']['log_directory']
conductivity_startup = config['ctd_processing']['conductivity_start']
tc1_align = config['ctd_processing']['TC_primary_align']
tc2_align = config['ctd_processing']['TC_secondary_align']
do_align = config['ctd_processing']['DO_align']
sample_rate = config['ctd_processing']['sample_rate']
search_time = config['ctd_processing']['roll_filter_time']
H1 = config['ctd_processing']['hysteresis_1']
H2 = config['ctd_processing']['hysteresis_2']
H3 = config['ctd_processing']['hysteresis_3']
ctd = config['ctd_processing']['ctd_serial']

lat_col = config['inputs']['lat']
lon_col = config['inputs']['lon']
alt_col = config['inputs']['alt']
input_parameters = config['inputs']['input_array'].split("\n")
p_col = config['inputs']['p']
t1_col = config['inputs']['t1']
t2_col = config['inputs']['t2']
c1_col = config['inputs']['c1']
c2_col = config['inputs']['c2']
sal_col = config['inputs']['salt']
dopl_col = config['inputs']['dopl']
dopkg_col = config['analytical_inputs']['dopkg']
xmis_col = config['inputs']['xmis']
fluor_col = config['inputs']['fluor']
v2_col = config['inputs']['backscatter']
v3_col = config['inputs']['rinko_oxy']
v4_col = config['inputs']['rinko_tmp']
time_zone = config['inputs']['time_zone']
nmea_time_col = config['inputs']['nmea_datetime']
scan_time_col = config['inputs']['scan_datetime']

#time_column_data = config['time_series_output']['data_names'].split(',')
time_column_data = config['time_series_output']['data_output'].split(',')
time_column_names = config['time_series_output']['column_name'].split(',')
time_column_units = config['time_series_output']['column_units'].split(',')
time_column_format = config['time_series_output']['format']

#pressure_column_data = config['time_series_output']['data_names'].split(',')
p_column_data = config['pressure_series_output']['data_output'].split(',')
p_column_names = config['pressure_series_output']['column_name'].split(',')
p_column_units = config['pressure_series_output']['column_units'].split(',')
p_column_format = config['pressure_series_output']['format']
p_column_qual = config['pressure_series_output']['qual_columns'].split(',')
p_column_one = list(config['pressure_series_output']['q1_columns'].split(','))

ext_xml = config['file_extensions']['ext_xml']
ext_hex = config['file_extensions']['ext_hex']
ext_ct1 = config['file_extensions']['ext_ct1']

                        ##### START LOGIC #####

all_casts = []
with open(file_ssscc, 'r') as filename:
    all_casts = [line.strip() for line in filename]

### odf_convert_sbe section
for cast in all_casts:
    sbeReader = sbe_reader.SBEReader.from_paths(raw_directory + cast + ext_hex,
                                                raw_directory + cast + ext_xml)
    ###foo
    converted_df = cnv.convertFromSBEReader(sbeReader, False)

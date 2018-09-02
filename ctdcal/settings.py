#!/usr/bin/env python
import os
import numpy as np
import pathlib



cruise = dict(
    cruise_title = '',
    cruise_name = '',
    cruisedb = '',
    vessel = '',
    ship_code = '',
    expocode = '',
    chief_sci = '',
    sectionid = '',
    start_date = '',
    start_port = '',
    start_latlon = '',
    end_date = '',
    end_port = '',
    end_latlon = ''
)

ctd_processing_dir = dict(
    raw_data_directory = 'data/raw/',
    time_data_directory = 'data/time/',
    pressure_data_directory = 'data/pressure/',
    fit_data_directory = 'data/fit/',
    log_directory = 'data/logs/',
    salt_directory = 'data/salt/',
    oxygen_directory = 'data/oxygen/',
    bottle_directory = 'data/bottle/',
    reft_directory = 'data/reft/',
    o2flask_file = 'data/oxygen/o2flasks.vol',
    pressure_log = 'data/logs/ondeck_pressure.csv',
    hex_prefix = 'data/raw/',
    hex_postfix ='.hex',
    xml_prefix = 'data/raw/',
    xml_postfix = '.XMLCON'
)

ctd_processing_constants = dict(
    conductivity_start = 20.0,
    pressure_interval = 2.0,
    sample_rate = 4,
    roll_filter_time = 15,
    TC_side = 2,
    roll_filter = 1,
    hysteresis_1 = -0.033,
    hysteresis_2 = 5000,
    hysteresis_3 = 1450,
    TC_primary_align = 0.0,
    TC_secondary_align = 0.0,
    DO_align = 2.0,
    ctd_serial = 1281
)


bottle_inputs = dict(
    p = 'CTDPRS', # Pressure column
    t = 'CTDTMP1', # CTD Temperature 
    t1 = 'CTDTMP1', # CTD Temperature 
    t2 = 'CTDTMP2', # CTD Temperature 
    c = 'CTDCOND1', # CTD Conductivity 
    c1 = 'CTDCOND1', # CTD Conductivity 1 
    c2 = 'CTDCOND2', # CTD Conductivity 2
    salt = 'CTDSAL', # CTD Salinity
    dopl = 'CTDOXY1', # CTD Oxygen (ml/L)
    dov = 'CTDOXYVOLTS', # CTD Oxygen volts
    xmis = 'CTDXMISS', # Transmissometer column
    fluor = 'FLUOR', # Fluorometer column
    backscatter = 'CTDBACKSCATTER', # Backscatter column
    rinko_oxy = 'FREE1', # Rinko Oxygen column (additional column)
    rinko_temp = 'FREE2', # Rinko Temp column (additional column)
    alt = 'ALT', # Altimeter 
    ref_par = 'REF_PAR',
    lat = 'GPSLAT', # Latitutde
    lon = 'GPSLON', # Longitude
    scan_datetime = 'scan_datetime', # Scan time
    nmea_datetime = 'nmea_datetime',
    btl_num = 'btl_fire_num', # Bottle fire number/bottle index column
    btl_fire = 'btl_fire', # Bottle fire check (0/1)
    pump_on = 'pump_on', # Pump on check (0/1)
    time_zone = 'UTC', # Time Zone
    reft = 'T90', # Reference Temp
    bath_temp = 'BathTEMP', # Salinometer bath temp
    cond_ratio = 'CRAVG', # Salinometer conductivity ratio
    btl_cond = 'BTLCOND', # Salinometer bottle conductivity
    btl_oxy = 'OXYGEN', # Bottle Oxygen

)

btl_input_array = [bottle_inputs['p'], bottle_inputs['t1'], bottle_inputs['t2'], bottle_inputs['c1'], bottle_inputs['c2'],
                   bottle_inputs['dov'], bottle_inputs['dopl'], bottle_inputs['xmis'], bottle_inputs['alt'],
                   bottle_inputs['rinko_oxy'], bottle_inputs['rinko_temp'], bottle_inputs['lat'], bottle_inputs['lon'], bottle_inputs['pump_on'], bottle_inputs['btl_fire'], 
                   bottle_inputs['scan_datetime'], bottle_inputs['btl_num']]

#btl_input_array = [bottle_inputs['p'], bottle_inputs['t1'], bottle_inputs['t2'], bottle_inputs['c1'], bottle_inputs['c2'],
#                   bottle_inputs['dov'], bottle_inputs['dopl'], bottle_inputs['xmis'], bottle_inputs['alt'], bottle_inputs['fluor'], bottle_inputs['backscatter'],
#                   bottle_inputs['rinko_oxy'], bottle_inputs['rinko_temp'], bottle_inputs['lat'], bottle_inputs['lon'], bottle_inputs['pump_on'], bottle_inputs['btl_fire'], 
#                   bottle_inputs['scan_datetime'], bottle_inputs['btl_num']]

ctd_inputs = dict(
    p = 'CTDPRS', # Pressure column
    t = 'CTDTMP1', # CTD Temp
    t1 = 'CTDTMP1', # CTD Temp 1
    t2 = 'CTDTMP2', # CTD Temp 2
    reft = 'REFTMP', # Reference Temp
    c = 'CTDCOND2', # CTD Conductivity
    c1 = 'CTDCOND1', # CTD Conductivity 1
    c2 = 'CTDCOND2', # CTD Conductivity 2
    salt = 'CTDSAL', # CTD Salinity
    btl_salt = 'SALNTY', # Bottle Salinity
    dov = 'CTDOXYVOLTS', # CTD Oxygen Volts
    dopl = 'CTDOXY1', # CTD Oxygen (ml/L)
    dopkg = 'CTDOXY1', # CTD Oxygem (umol/KG)
    xmis = 'TRANS', # Transmissometer
    fluor = 'FLUOR', # Fluorometer
    backscatter = 'CTDBACKSCATTER', # Backscatter
    rinko_oxy = 'FREE1', # Rinko Oxygen column (additional column)
    rinko_tmp = 'FREE2', # Rinko Temperature column (additional column)
    pump_on = 'pump_on', # Pump on check (0/1)
    alt = 'ALT', # Altimeter
    scan_datetime = 'scan_datetime', # Scan time
    datetime = 'GPSEPOCHTIME', # Time
    lat = 'GPSLAT', # Latitude
    lon = 'GPSLON', # Longitude    
)

ctd_input_array = [ctd_inputs['p'], ctd_inputs['t1'], ctd_inputs['t2'], ctd_inputs['c1'], ctd_inputs['c2'], ctd_inputs['dov'], ctd_inputs['dopl'],
                   ctd_inputs['xmis'], ctd_inputs['alt'], ctd_inputs['rinko_oxy'],ctd_inputs['rinko_tmp'],
                   ctd_inputs['lat'], ctd_inputs['lon'], ctd_inputs['pump_on'], ctd_inputs['scan_datetime']]

#ctd_input_array = [ctd_inputs['p'], ctd_inputs['t1'], ctd_inputs['t2'], ctd_inputs['c1'], ctd_inputs['c2'], ctd_inputs['dov'], ctd_inputs['dopl'],
#                   ctd_inputs['xmis'], ctd_inputs['alt'], ctd_inputs['fluor'], ctd_inputs['backscatter'], ctd_inputs['rinko_oxy'],ctd_inputs['rinko_tmp'],
#                   ctd_inputs['lat'], ctd_inputs['lon'], ctd_inputs['pump_on'], ctd_inputs['scan_datetime']]

# Create necessary files:
for key, val in  ctd_processing_dir.items():
    if key != 'o2flask_file' and key != 'pressure_log':
        pathlib.Path(val).mkdir(exist_ok=True)
    else:
        continue

# Generate ssscc from directory
ssscc = os.listdir(ctd_processing_dir['raw_data_directory'])
ssscc = [s[0:5] for s in ssscc]
for s in ssscc:
    if not s.startswith('.'):
        s = s[0:5]
    else:
        ssscc.remove(s)
ssscc = np.unique(ssscc)
ssscc = list(ssscc)

#### Use this to remove individual stations
#remove_list = ['14601','14402']
#ssscc = [x for x in ssscc if x not in remove_list]


### Create your own SSSCC:
# ssscc = []

# Choose Sensors to calibrate (primary or secondary)

do_primary = 1 #0
do_secondary = 1 #0


#!/usr/bin/env python
import os
import pathlib



cruise = dict(
    cruise_title = 'SOCCOM NBP1707',
    cruise_name = 'nbp1707',
    cruisedb = 'nbp1707',
    vessel = 'R/V Palmer',
    ship_code = '3206',
    expocode = '320620170820',
    chief_sci = 'Kevin Speer',
    sectionid = 'nbp1707',
    start_date = '2017-08-20 11:00:00',
    start_port = 'Papeete, French Polynesia',
    start_latlon = '17.535 S 149.5696 W',
    end_date = '2017-09-30 10:00:00',
    end_port = 'Valparaiso, Chile',
    end_latlon = '33.05 S 71.616 W'
)

ctd_processing_dir = dict(
    raw_data_directory = 'data/raw/',
    time_data_directory = 'data/time/',
    pressure_data_directory = 'data/pressure/',
    fit_data_directory = 'data/fit/',
    converted_directory = 'data/converted/',
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
    fluor = 'CTDFLUOR', # Fluorometer column
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
    cond_ratio = 'CRavg', # Salinometer conductivity ratio
    btl_cond = 'BTLCOND', # Salinometer bottle conductivity
    btl_oxy = 'OXYGEN', # Bottle Oxygen

)

btl_input_array = [bottle_inputs['p'], bottle_inputs['t1'], bottle_inputs['t2'], bottle_inputs['c1'], bottle_inputs['c2'],
                   bottle_inputs['dov'], bottle_inputs['dopl'], bottle_inputs['xmis'], bottle_inputs['alt'],
                   bottle_inputs['rinko_oxy'], bottle_inputs['rinko_temp'], bottle_inputs['lat'], bottle_inputs['lon'], bottle_inputs['pump_on'], bottle_inputs['btl_fire'], 
                   bottle_inputs['scan_datetime'], bottle_inputs['btl_num']]#bottle_inputs['fluor']

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
    xmis = 'CTDXMISS', # Transmissometer
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
                   ctd_inputs['lat'], ctd_inputs['lon'], ctd_inputs['pump_on'], ctd_inputs['scan_datetime']]#,ctd_inputs['fluor']

#ctd_input_array = [ctd_inputs['p'], ctd_inputs['t1'], ctd_inputs['t2'], ctd_inputs['c1'], ctd_inputs['c2'], ctd_inputs['dov'], ctd_inputs['dopl'],
#                   ctd_inputs['xmis'], ctd_inputs['alt'], ctd_inputs['fluor'], ctd_inputs['backscatter'], ctd_inputs['rinko_oxy'],ctd_inputs['rinko_tmp'],
#                   ctd_inputs['lat'], ctd_inputs['lon'], ctd_inputs['pump_on'], ctd_inputs['scan_datetime']]

time_series_output = dict(
data_output = ['CTDPRS_DBAR','CTDTMP1_ITS90','CTDTMP2_ITS90','CTDCOND1_MSPCM','CTDCOND2_MSPCM','CTDSAL_PSU','CTDOXY1_MLPL','CTDOXYVOLTS','CTDXMISS_05VDC','FLUOR_05VDC','scan_datetime','LATITUDE','LONGITUDE'],
column_name = ['CTDPRS','CTDTMP1','CTDTMP2','CTDCOND1','CTDCOND2','CTDSAL','CTDOXY','CTDOXYVOLTS','TRANS','FLUOR','GPSEPOCHTIME','GPSLAT','GPSLON'],
column_units =  ['DBAR', 'ITS90 degC', 'ITS90 degC', 'MS/CM', 'MS/CM', 'PSS78 PSU', 'ML/L', '0-5VDC', '0-5VDC', '0-5VDC', 'SECONDS', '-S', '-W'],
format = ['%8.1f','%10.4f','%10.4f','%10.4f','%10.4f','%10.4f','%10.4f','%10.4f','%10.4f','%d','%10.5f','%10.5f\n']
)

pressure_series_output = dict(
data = ['CTDPRS','CTDTMP1','CTDSAL','TRANS','FLUOR'],
data_output = ['CTDPRS_DBAR','CTDTMP1_ITS90','CTDSAL_PSU','CTDXMISS_05VDC','FLUOR_05VDC'],
column_name = ['CTDPRS','CTDPRS_FLAG_W','CTDTMP','CTDTMP_FLAG_W','CTDSAL','CTDSAL_FLAG_W','CTDOXY','CTDOXY_FLAG_W','TRANS','TRANS_FLAG_W','FLUOR','FLUOR_FLAG_W','CTDBACKSCATTER','CTDBACKSCATTER_FLAG_W','CTDRINKO','CTDRINKO_FLAG_W'],
column_units = ['DBAR','','ITS90 degC','','PSS78 PSU','','UMOL/KG','','0-5VDC','','0-5VDC','','0-5VDC','','0-5VDC',''],
qual_columns = ['CTDPRS_FLAG_W','CTDTMP_FLAG_W', 'CTDSAL_FLAG_W', 'CTDOXY_FLAG_W', 'TRANS_FLAG_W', 'FLUOR_FLAG_W'],
q1_columns = ['', 'CTDOXY_FLAG_W', 'TRANS_FLAG_W', 'FLUOR_FLAG_W']
)

input_array = [ctd_inputs['p'], ctd_inputs['t'],ctd_inputs['salt'],ctd_inputs['dopl'],ctd_inputs['dov'],ctd_inputs['xmis'],ctd_inputs['fluor']]
# Create necessary files:
for key, val in  ctd_processing_dir.items():
    if key != 'o2flask_file' and key != 'pressure_log':
        pathlib.Path(val).mkdir(exist_ok=True)
    else:
        continue

# Generate ssscc from directory
ssscc = os.listdir(ctd_processing_dir['raw_data_directory'])
xml_list = []
hex_list = []
for s in ssscc:
    if '.XMLCON' in s:
        xml_list.append(s[:-7])
    elif '.hex' in s:
        hex_list.append(s[:-4])
    else:
        pass
ssscc = list(set(xml_list) & set(hex_list))
ssscc.sort()
    


do_primary = 1 #0
do_secondary = 1 #0


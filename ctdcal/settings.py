#!/usr/bin/env python
import os
import pathlib

cruise = dict(
    cruise_title="GO-SHIP NBP1802",
    cruise_name="nbp1802",
    cruisedb="nbp1802",
    vessel="R/V Palmer",
    ship_code="3206",
    expocode="320620180309",
    chief_sci="Alison Macdonald",
    sectionid="S04P",
    start_date="2018-03-09 14:00:00",
    start_port="Hobart, Tasmania, Australia",
    start_latlon="42.8821 S 147.3272 W",
    end_date="2018-05-14 10:00:00",
    end_port="Punta Arenas, Chile",
    end_latlon="53.1638 S 70.9171 W",
)

ctd_processing_dir = dict(
    raw_data_directory="data/raw/",
    time_data_directory="data/time/",
    pressure_data_directory="data/pressure/",
    fit_data_directory="data/fit/",
    converted_directory="data/converted/",
    log_directory="data/logs/",
    salt_directory="data/salt/",
    oxygen_directory="data/oxygen/",
    bottle_directory="data/bottle/",
    reft_directory="data/reft/",
    o2flask_file="data/oxygen/o2flasks.vol",
    pressure_log="data/logs/ondeck_pressure.csv",
    hex_prefix="data/raw/",
    hex_postfix=".hex",
    xml_prefix="data/raw/",
    xml_postfix=".XMLCON",
)

ctd_processing_constants = dict(
    conductivity_start=20.0,
    pressure_interval=2.0,
    sample_rate=4,
    roll_filter_time=15,
    TC_side=2,
    roll_filter=1,
    hysteresis_1=-0.033,
    hysteresis_2=5000,
    hysteresis_3=1450,
    TC_primary_align=0.0,
    TC_secondary_align=0.0,
    DO_align=2.0,
    ctd_serial=1281,
)


bottle_inputs = dict(
    p="CTDPRS",  # Pressure column
    t="CTDTMP2",  # CTD Temperature
    t1="CTDTMP1",  # CTD Temperature
    t2="CTDTMP2",  # CTD Temperature
    c="CTDCOND1",  # CTD Conductivity
    c1="CTDCOND1",  # CTD Conductivity 1
    c2="CTDCOND2",  # CTD Conductivity 2
    salt="CTDSAL",  # CTD Salinity
    dopl="CTDOXY1",  # CTD Oxygen (ml/L)
    dov="CTDOXYVOLTS",  # CTD Oxygen volts
    xmis="CTDXMISS",  # Transmissometer column
    fluor="FLUOR",  # Fluorometer column
    backscatter="CTDBACKSCATTER",  # Backscatter column
    rinko_oxy="FREE1",  # Rinko Oxygen column (additional column)
    rinko_temp="FREE3",  # Rinko Temp column (additional column)
    alt="ALT",  # Altimeter
    ref_par="REF_PAR",
    lat="GPSLAT",  # Latitutde
    lon="GPSLON",  # Longitude
    scan_datetime="scan_datetime",  # Scan time
    nmea_datetime="nmea_datetime",
    btl_num="btl_fire_num",  # Bottle fire number/bottle index column
    btl_fire="btl_fire",  # Bottle fire check (0/1)
    pump_on="pump_on",  # Pump on check (0/1)
    time_zone="UTC",  # Time Zone
    #    T90 = 'T90', # Reference Temp
    reft="REFTMP",
    bath_temp="BathTEMP",  # Salinometer bath temp
    cond_ratio="CRavg",  # Salinometer conductivity ratio
    btl_cond="BTLCOND",  # Salinometer bottle conductivity
    btl_oxy="OXYGEN",  # Bottle Oxygen
)

btl_input_array = [
    bottle_inputs["p"],
    bottle_inputs["t1"],
    bottle_inputs["t2"],
    bottle_inputs["c1"],
    bottle_inputs["c2"],
    bottle_inputs["dov"],
    bottle_inputs["dopl"],
    bottle_inputs["xmis"],
    bottle_inputs["alt"],
    bottle_inputs["rinko_oxy"],
    bottle_inputs["rinko_temp"],
    bottle_inputs["lat"],
    bottle_inputs["lon"],
    bottle_inputs["pump_on"],
    bottle_inputs["btl_fire"],
    bottle_inputs["scan_datetime"],
    bottle_inputs["btl_num"],
    bottle_inputs["fluor"],
]  # bottle_inputs['fluor']

# Dimensioned/flagged bottle outputs:
# Form = ['PARAMETER', 'PARAMETER_FLAG_NAME', 'PARAMETER_UNITS','PARAMETER_FLAG_UNITS']
btl_outputs = dict(
    btl_nbr=["BTLNBR", "BTLNBR_FLAG_W", "", ""],
    press=["CTDPRS", "CTDPRS_FLAG_W", "DBAR", ""],
    temp=["CTDTMP", "CTDTMP_FLAG_W", "ITS-90", ""],
    ref_tmp=["REFTMP", "REFTMP_FLAG_W", "ITS-90", ""],
    salt=["CTDSAL", "CTDSAL_FLAG_W", "PSS-78", ""],
    slnty=["SALNTY", "SALNTY_FLAG_W", "PSS-78", ""],
    ctd_oxy=["CTDOXY", "CTDOXY_FLAG_W", "UMOL/KG", ""],
    oxy=["OXYGEN", "OXYGEN_FLAG_W", "UMOL/KG", ""],
    rinko=["CTDRINKO", "CTDRINKO_FLAG_W", "UMOL/KG", ""],
    xmiss=["CTDXMISS", "CTDXMISS_FLAG_W", "0-5VDC", ""],
    fluor=["CTDFLUOR", "CTDFLUOR_FLAG_W", "0-5VDC", ""],
    # bbp = ['CTDBBP700RAW', 'CTDBBP700RAW_FLAG_W', '0-5VDC', ''],
)
# Undimensioned/unflagged bottle outputs:
# Form = ['PARAMETER','PARAMETER NOTE/UNIT'] parameter note/unit is what will be displayed underneath the parameter name in the exchange file
btl_und_outputs = dict(
    expocd=["EXPOCODE", ""],
    sectid=["SECT_ID", ""],
    stnnbr=["STNNBR", ""],
    castno=["CASTNO", ""],
    sampno=["SAMPNO", ""],
    date=["DATE", ""],
    time=["TIME", ""],
    lat=["LATITUDE", ""],
    lon=["LONGITUDE", ""],
    dep=["DEPTH", ""],
)
btl_column_names = []
btl_column_units = []
btl_flagged_params = []

for i in range(len(btl_und_outputs)):
    param = list(btl_und_outputs.keys())[i]
    param_list = btl_und_outputs[param]
    btl_column_names.append(param_list[0])
    btl_column_units.append(param_list[1])

for i in range(len(btl_outputs)):
    param = list(btl_outputs.keys())[i]
    param_list = btl_outputs[param]
    btl_column_names.append(param_list[0])
    btl_column_names.append(param_list[1])
    btl_column_units.append(param_list[2])
    btl_column_units.append(param_list[3])
    btl_flagged_params.append(param_list[0])

btl_series_output = dict(
    btl_column_names=btl_column_names, btl_column_units=btl_column_units
)

ctd_inputs = dict(
    p="CTDPRS",  # Pressure column
    t="CTDTMP2",  # CTD Temp
    t1="CTDTMP1",  # CTD Temp 1
    t2="CTDTMP2",  # CTD Temp 2
    reft="REFTMP",  # Reference Temp
    c="CTDCOND2",  # CTD Conductivity
    c1="CTDCOND1",  # CTD Conductivity 1
    c2="CTDCOND2",  # CTD Conductivity 2
    salt="CTDSAL",  # CTD Salinity
    btl_salt="SALNTY",  # Bottle Salinity
    dov="CTDOXYVOLTS",  # CTD Oxygen Volts
    dopl="CTDOXY1",  # CTD Oxygen (ml/L)
    dopkg="CTDOXY1",  # CTD Oxygem (umol/KG)
    xmis="CTDXMISS",  # Transmissometer
    fluor="CTDFLUOR",  # Fluorometer
    backscatter="CTDBACKSCATTER",  # Backscatter
    rinko_oxy="FREE1",  # Rinko Oxygen column (additional column)
    rinko_tmp="FREE2",  # Rinko Temperature column (additional column)
    pump_on="pump_on",  # Pump on check (0/1)
    alt="ALT",  # Altimeter
    scan_datetime="scan_datetime",  # Scan time
    datetime="GPSEPOCHTIME",  # Time
    lat="GPSLAT",  # Latitude
    lon="GPSLON",  # Longitude
)
ctd_outputs = dict(
    press=["CTDPRS", "CTDPRS_FLAG_W", "DBAR", ""],
    temp=["CTDTMP", "CTDTMP_FLAG_W", "ITS-90", ""],
    salt=["CTDSAL", "CTDSAL_FLAG_W", "PSS-78", ""],
    doxy=["CTDOXY", "CTDOXY_FLAG_W", "UMOL/KG", ""],
    rinko=["CTDRINKO", "CTDRINKO_FLAG_W", "0-5VDC", ""],
    xmiss=["CTDXMISS", "CTDXMISS_FLAG_W", "0-5VDC", ""],
    fluor=["CTDFLUOR", "CTDFLUOR_FLAG_W", "0-5VDC", ""],
    # bbp = ['CTDBBP700RAW', 'CTDBBP700RAW_FLAG_W', '0-5VDC', ''],
)

ctd_input_array = [
    ctd_inputs["p"],
    ctd_inputs["t1"],
    ctd_inputs["t2"],
    ctd_inputs["c1"],
    ctd_inputs["c2"],
    ctd_inputs["dov"],
    ctd_inputs["dopl"],
    ctd_inputs["xmis"],
    ctd_inputs["alt"],
    ctd_inputs["rinko_oxy"],
    ctd_inputs["rinko_tmp"],
    ctd_inputs["lat"],
    ctd_inputs["lon"],
    ctd_inputs["pump_on"],
    ctd_inputs["scan_datetime"],
]  # ,ctd_inputs['fluor']

# ctd_input_array = [ctd_inputs['p'], ctd_inputs['t1'], ctd_inputs['t2'], ctd_inputs['c1'], ctd_inputs['c2'], ctd_inputs['dov'], ctd_inputs['dopl'],
#                   ctd_inputs['xmis'], ctd_inputs['alt'], ctd_inputs['fluor'], ctd_inputs['backscatter'], ctd_inputs['rinko_oxy'],ctd_inputs['rinko_tmp'],
#                   ctd_inputs['lat'], ctd_inputs['lon'], ctd_inputs['pump_on'], ctd_inputs['scan_datetime']]

time_series_output = dict(
    data_output=[
        "CTDPRS_DBAR",
        "CTDTMP1_ITS90",
        "CTDTMP2_ITS90",
        "CTDCOND1_MSPCM",
        "CTDCOND2_MSPCM",
        "CTDSAL_PSU",
        "CTDOXY1_MLPL",
        "CTDOXYVOLTS",
        "CTDXMISS_05VDC",
        "FLUOR_05VDC",
        "scan_datetime",
        "LATITUDE",
        "LONGITUDE",
    ],
    column_name=[
        ctd_outputs["press"][0],
        "CTDTMP1",
        "CTDTMP2",
        "CTDCOND1",
        "CTDCOND2",
        "CTDSAL",
        "CTDOXY",
        "CTDOXYVOLTS",
        "TRANS",
        "FLUOR",
        "GPSEPOCHTIME",
        "GPSLAT",
        "GPSLON",
    ],
    column_units=[
        "DBAR",
        "ITS90 degC",
        "ITS90 degC",
        "MS/CM",
        "MS/CM",
        "PSS78 PSU",
        "ML/L",
        "0-5VDC",
        "0-5VDC",
        "0-5VDC",
        "SECONDS",
        "-S",
        "-W",
    ],
    format=[
        "%8.1f",
        "%10.4f",
        "%10.4f",
        "%10.4f",
        "%10.4f",
        "%10.4f",
        "%10.4f",
        "%10.4f",
        "%10.4f",
        "%d",
        "%10.5f",
        "%10.5f\n",
    ],
)

ctd_column_names = []
ctd_column_units = []
for i in range(len(ctd_outputs)):
    param = list(ctd_outputs.keys())[i]
    param_list = ctd_outputs[param]
    ctd_column_names.append(param_list[0])
    ctd_column_names.append(param_list[1])
    ctd_column_units.append(param_list[2])
    ctd_column_units.append(param_list[3])

pressure_series_output = dict(
    column_names=ctd_column_names, column_units=ctd_column_units
)


input_array = [
    ctd_inputs["p"],
    ctd_inputs["t"],
    ctd_inputs["salt"],
    ctd_inputs["dopl"],
    ctd_inputs["dov"],
    ctd_inputs["xmis"],
    ctd_inputs["fluor"],
]
# Create necessary files:
for key, val in ctd_processing_dir.items():
    if key != "o2flask_file" and key != "pressure_log":
        pathlib.Path(val).mkdir(exist_ok=True)
    else:
        continue

# Generate ssscc from directory
ssscc = os.listdir(ctd_processing_dir["raw_data_directory"])
xml_list = []
hex_list = []
for s in ssscc:
    if ".XMLCON" in s:
        xml_list.append(s[:-7])
    elif ".hex" in s:
        hex_list.append(s[:-4])
    else:
        pass
ssscc = list(set(xml_list) & set(hex_list))
ssscc.sort()


do_primary = 1  # 0
do_secondary = 1  # 0

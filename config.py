### configuration file for odf_process_all_MK.py
# TODO: organize these by editable/fixed variables

# Cruise specifics
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

ctd_serial = 1281

# CTD variables/flags/units
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

# List of directories for I/O purposes
directory = {
    "ssscc_file": "data/ssscc.csv",
    "logs": "data/logs/",
    "qual_temp_primary": "quality_code/temp_primary/",
    "qual_temp_secondary": "quality_code/temp_secondary/",
    "qual_cond_primary": "quality_code/cond_primary/",
    "qual_cond_secondary": "quality_code/cond_secondary/",
}

# Labels for CTD columns
column = {
    "p": "CTDPRS",
    "p_btl": "CTDPRS",
    "t1": "CTDTMP1",
    "t2": "CTDTMP2",
    "t1_btl": "CTDTMP1",
    "t2_btl": "CTDTMP2",
    "c1": "CTDCOND1",
    "c2": "CTDCOND2",
    "c1_btl": "CTDCOND1",
    "c2_btl": "CTDCOND2",
    "reft": "T90",
    "refc": "BTLCOND",
    "sal": "CTDSAL",
    "sal_btl": "SALNTY",
    "rinko_oxy": "FREE1",
    "oxy_btl": "OXYGEN",
    "oxyvolts": "CTDOXYVOLTS",
    "lat": "GPSLAT",
    "lon": "GPSLON",
    "lat_btl": "GPSLAT",
    "lon_btl": "GPSLON",
}

# List of bottle columns to be read during calibration
btl_cols = [
    "index",
    "CTDTMP1",
    "CTDTMP2",
    "CTDPRS",
    "CTDCOND1",
    "CTDCOND2",
    "CTDSAL",
    "CTDOXY1",
    "CTDOXYVOLTS",
    "CTDXMISS",
    "ALT",
    "REF_PAR",
    "GPSLAT",
    "GPSLON",
    "new_fix",
    "pressure_temp_int",
    "pump_on",
    "btl_fire",
    "scan_datetime",
    "btl_fire_num",
]
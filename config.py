### configuration file for odf_process_all_MK.py

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

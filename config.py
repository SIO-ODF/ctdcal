### configuration file for odf_process_all_MK.py
# this is an example to fill with your own values
#
# TODO: organize these by editable/fixed variables

# Cruise specifics
expocode = "325020210316"
section_id = "A20"
ctd_serial = 914  # TODO: how to handle multiple CTDs on one cruise?

# CTD file (.ct1) variable outputs
# move elsewhere when xarray is implemented
# TODO: import dict from user_settings.yaml (working title), check against cchdo.params?
ctd_outputs = dict(
    press=["CTDPRS", "DBAR"],
    temp=["CTDTMP", "ITS-90"],
    salt=["CTDSAL", "PSS-78"],
    doxy=["CTDOXY", "UMOL/KG"],
    # rinko=["CTDRINKO", "UMOL/KG"],
    xmiss=["CTDXMISS", "0-5VDC"],
    fluor=["CTDFLUOR", "0-5VDC"],
    # backscatter=["CTDBACKSCATTER", "0-5VDC"],
)
ctd_col_names, ctd_col_units = [], []
for (param, unit) in ctd_outputs.values():
    if param == "CTDPRS":
        ctd_col_names += [param]
        ctd_col_units += [unit]
    else:
        ctd_col_names += [param, f"{param}_FLAG_W"]
        ctd_col_units += [unit, ""]

# List of directories for I/O purposes
directory = {
    "ssscc_file": "data/ssscc.csv",
    "ssscc": "data/ssscc/",
    "raw": "data/raw/",
    "converted": "data/converted/",
    "time": "data/time/",
    "pressure": "data/pressure/",
    "bottle": "data/bottle/",
    "salt": "data/salt/",
    "reft": "data/reft/",
    "oxy": "data/oxygen/",
    "logs": "data/logs/",
    "qual_temp_primary": "quality_code/temp_primary/",
    "qual_temp_secondary": "quality_code/temp_secondary/",
    "qual_cond_primary": "quality_code/cond_primary/",
    "qual_cond_secondary": "quality_code/cond_secondary/",
    "t1_fit_figs": "data/logs/fitting_figs/temp_primary/",
    "t2_fit_figs": "data/logs/fitting_figs/temp_secondary/",
    "c1_fit_figs": "data/logs/fitting_figs/cond_primary/",
    "c2_fit_figs": "data/logs/fitting_figs/cond_secondary/",
    "ox_fit_figs": "data/logs/fitting_figs/oxy_primary/",
    "rinko_fit_figs": "data/logs/fitting_figs/oxy_rinko/",
}

# remnant of old system, will be pushed into xarray metadata/attrs
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
    "rinko_oxy": "U_DEF_poly1",  # CHECK THIS!
    "oxy_btl": "OXYGEN",
    "oxyvolts": "CTDOXYVOLTS",
    "lat": "GPSLAT",
    "lon": "GPSLON",
    "lat_btl": "GPSLAT",
    "lon_btl": "GPSLON",
}

# List of columns to filter
filter_cols = []
for x in ["p", "t1", "t2", "c1", "c2", "sal", "rinko_oxy", "oxyvolts", "lat", "lon"]:
    filter_cols.append(column[x])

# List of bottle columns to be read during calibration
btl_cols = [
    "CTDTMP1",
    "CTDTMP2",
    "CTDPRS",
    "CTDCOND1",
    "CTDCOND2",
    "CTDSAL",
    "CTDOXY1",
    "CTDOXYVOLTS",
    "U_DEF_poly1",
    "CTDXMISS",
    "CTDFLUOR",
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

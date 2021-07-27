### configuration file for odf_process_all.py
#
# TODO: organize these by editable/fixed variables
import yaml

with open("user_settings.yaml", "r") as f:
    settings = yaml.safe_load(f)

# Unpack user settings (any sanitizing/checks needed? probably)
expocode = settings["expocode"]
section_id = settings["section_id"]
ctd_serial = settings["ctd_serial"]
ctd_outputs = settings["ctd_outputs"]

# CTD file (.ct1) variable outputs
# move elsewhere when xarray is implemented
# TODO: check against cchdo.params?
ctd_col_names, ctd_col_units = [], []
for (param, attrs) in ctd_outputs.items():
    if param == "CTDPRS":
        ctd_col_names += [param]
        ctd_col_units += [attrs["units"]]
    else:
        ctd_col_names += [param, f"{param}_FLAG_W"]
        ctd_col_units += [attrs["units"], ""]

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
    "rinko_oxy": "FREE1",  # CHECK THIS!
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

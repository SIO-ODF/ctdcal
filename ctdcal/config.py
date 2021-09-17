### configuration file for odf_process_all.py
#
# TODO: organize these by editable/fixed variables
from importlib import resources
import yaml

with resources.open_text("ctdcal", "user_settings.yaml") as f:
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
dirs = {
    "ssscc": "data/ssscc/",
    "raw": "data/raw/",
    "converted": "data/converted/",
    "time": "data/time/",
    "pressure": "data/pressure/",
    "bottle": "data/bottle/",
    "salt": "data/salt/",
    "reft": "data/reft/",
    "oxygen": "data/oxygen/",
    "logs": "data/logs/",
}
fig_dirs = {
    "t1": "data/logs/fitting_figs/temp_primary/",
    "t2": "data/logs/fitting_figs/temp_secondary/",
    "c1": "data/logs/fitting_figs/cond_primary/",
    "c2": "data/logs/fitting_figs/cond_secondary/",
    "ox": "data/logs/fitting_figs/oxy_primary/",
    "rinko": "data/logs/fitting_figs/oxy_rinko/",
}

# remnant of old system, will be pushed into xarray metadata/attrs
# Labels for CTD columns
column = {
    "p": "CTDPRS",
    "t1": "CTDTMP1",
    "t2": "CTDTMP2",
    "c1": "CTDCOND1",
    "c2": "CTDCOND2",
    "sal": "CTDSAL",
    # "s1": "CTDSAL1",  # TODO: calc salinity from primary and secondary sensors
    # "s2": "CTDSAL2",
    "rinko_oxy": "FREE1",  # CHECK THIS!
    "oxyvolts": "CTDOXYVOLTS",
    "refT": "T90",
    "refC": "BTLCOND",
    "refS": "SALNTY",
    "refO": "OXYGEN",
    "lat": "GPSLAT",
    "lon": "GPSLON",
}

# List of columns to filter
filter_cols = []
for x in ["p", "t1", "t2", "c1", "c2", "sal", "rinko_oxy", "oxyvolts", "lat", "lon"]:
    filter_cols.append(column[x])

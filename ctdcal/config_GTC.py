### configuration file for gp17_process_all.py

# Unpack user settings (any sanitizing/checks needed? probably)
platform = "GTC"
expocode = "33RR20221201"
section_id = "GP17-OCE"
ctd_serial = 888
ctd_outputs = {
    "CTDPRS": {"sensor": "CTDPRS", "units": "DBAR"},
    "CTDTMP": {"sensor": "CTDTMP1", "units": "ITS-90"},
    "CTDSAL": {"sensor": "CTDSAL1", "units": "PSS-78"},
    "CTDOXY": {"sensor": "CTDOXY1", "units": "UMOL/KG"},
    "CTDXMISS": {"sensor": "CTDXMISS1", "units": "0-5VDC"},
    "CTDFLUOR": {"sensor": "FLUOR", "units": "0-5VDC"},
}

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
    "ssscc": "data/ssscc_GTC/",
    "raw": "data/raw/",
    "converted": "data/converted/",
    "time": "data/time/",
    "pressure": "data/pressure/",
    "bottle": "data/bottle/",
    "salt": "data/salt/",
    "reft": "data/reft/",
    "oxygen": "data/oxygen/",
    "logs": "data/logs/",
    "scratch": "data/scratch_folder",
}
fig_dirs = {
    "t1": "data/logs/fitting_figs_GTC/temp_primary/",
    "t2": "data/logs/fitting_figs_GTC/temp_secondary/",
    "c1": "data/logs/fitting_figs_GTC/cond_primary/",
    "c2": "data/logs/fitting_figs_GTC/cond_secondary/",
    "ox": "data/logs/fitting_figs_GTC/oxy_primary/",
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
    "oxyvolts": "CTDOXYVOLTS",
    "refT": "REFTMP",
    "refC": "BTLCOND",
    "refS": "SALNTY",
    "refO": "OXYGEN",
    "lat": "GPSLAT",
    "lon": "GPSLON",
}

# List of columns to filter
filter_cols = []
for x in ["p", "t1", "t2", "c1", "c2", "sal", "oxyvolts", "lat", "lon"]:
    filter_cols.append(column[x])

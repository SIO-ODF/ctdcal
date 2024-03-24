# user configuration file
# TODO: consolidate and organize all the various config files
#
# TODO: organize these by editable/fixed variables
from importlib import resources
from munch import Munch

import yaml


class Parameters(object):
    """
    Base class to hold high level parameters for parsing and processing
    """
    def __init__(self, parameters=None):
        parameters = [] if parameters is None else parameters
        self.parameters = parameters

    def create_dict(self):
        """
        Create a Bunch class object to store the parameter names for the data
        files.
        """
        bunch = Munch()

        for name in self.parameters:
            bunch[name] = []

        return bunch

    def from_csv(self, fname):
        """
        Creates a list from a csv.
        """
        lines = []
        with open(fname, 'r', encoding='utf8') as fid:
            for line in fid:
                if not line.startswith('#'):
                    lines.append(line.strip())
        return lines


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
    "rinko_oxy": "U_DEF_poly1",  # CHECK THIS!
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
for x in ["p", "t1", "t2", "c1", "c2", "sal", "rinko_oxy", "oxyvolts", "lat", "lon"]:
    filter_cols.append(column[x])

# Fitting coefficients file
coeffs_file = '_fit_coefs.yaml'

# Cast Tools settings
win_size = 2     # default filter window size
max_soak = 20    # maximum soak pressure threshold
despike_deltas = {'CTDSAL': 0.001,
                  'CTDTMP1': 0.75,
                  'CTDTMP2': 0.75}   # despike filter thresholds

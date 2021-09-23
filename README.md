[![Documentation Status](https://readthedocs.org/projects/ctdcal/badge/?version=latest)](https://ctdcal.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

# ctdcal project

The ctdcal project is a library primarily designed to process data from CTD casts and calibrate
them against Niskin bottle samples.

In the future, parts of the ctdcal library will be split off into additional packages,
such as an "ocean sensors" package with Python implementations of conversion routines
for in-situ sensors used for ocean measurement.

---

## Installation
ctdcal can be installed using pip:

```
$ pip install ctdcal
```

---

## CLI usage
### Initialize data folders
Initialize default `/data/` folders by running:

```
$ ctdcal init
```

(Future versions of ctdcal are planned have more robust init options/flags/etc.)

### Import and process data
To process data, copy over raw `.hex` and `.xmlcon` files into `/data/raw/` and reference data to their appropriate folder (`oxygen`, `reft`, `salt`).

Users can process their data with individual ctdcal functions or try:

```
$ ctdcal process [--group ODF]
```

to process using ODF procedures.

---

## Package usage
### Explore user settings
Most ctdcal functions get settings from `user_settings.yaml` and subsequently `config.py`. Call the configuration loader to explore default settings:

```py
from ctdcal import get_ctdcal_config
cfg = get_ctdcal_config()

# directories for I/O purposes
print(cfg.dirs)
print(cfg.fig_dirs)

# experiment-specific settings (e.g., expocode, CTD serial number) from user_settings.yaml
print(cfg.settings)

# dictionary mapping of short/long column names
print(cfg.columns)
```

As ctdcal continues to be developed, more robust [tutorials](https://ctdcal.readthedocs.io/en/latest/tutorials.html) will be added to [our documentation](https://ctdcal.readthedocs.io/en/latest/).

---

## LICENSING
BSD 3-clause
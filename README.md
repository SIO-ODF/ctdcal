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
### Clone repository
Pull down the latest version of ctdcal:

```
git clone https://github.com/cchdo/ctdcal.git
```

### Install ctdcal and dependencies
Change directories to the top-level ctdcal:

```
cd ctdcal
```

Create a new virtual environment with your preferred environment manager and install with pip:

```
pip install .
```

Note: there is an occasional (conda?) bug where CLI tools are not immediately accessible after install – this can usually be remedied by deactivating and reactiving the virtual environment.

Initialize default `/data/` folders by running:

```
ctdcal init
```

(Future versions of ctdcal are planned have more robust init options/flags/etc.)

---

## Usage
### Import and process data
To process data, copy over raw `.hex` and `.xmlcon` files into `/data/raw/` and reference data to their appropriate folder (`oxygen`, `reft`, `salt`).

Users can process their data with individual ctdcal functions or try:

```
ctdcal process [--group ODF]
```

to process using ODF procedures.

---

## LICENSING
BSD 3-clause
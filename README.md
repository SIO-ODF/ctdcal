[![Documentation Status](https://readthedocs.org/projects/ctdcal/badge/?version=latest)](https://ctdcal.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

# ctdcal project

The ctdcal project is a library primarily designed to process data from CTD casts and calibrate
them against Niskin bottle samples.

In the future, parts of the ctdcal library will be split off into additional packages,
such as an "ocean sensors" package with Python implementations of conversion routines
for in-situ sensors used for ocean measurement.

***
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

***
## Usage

### Import and process data
To process data, copy over raw `.hex` and `.xmlcon` files into `/data/raw/` and reference data to their appropriate folder (`oxygen`, `reft`, `salt`).

Users can process their data with individual ctdcal functions or try:
```
ctdcal process [--group ODF]
```
to process using ODF procedures.

<!-- 
## Overview of available commands
this is outdated but a good template possibly

#### Utility for converting Seabird Electronics .hex/.XMLCON raw data into csv-formatted text files
```
usage: odf_convert_sbe.py [-h] [-d] [-r] [-o destDir] hex_file XMLCON_file

positional arguments:
  hex_file      the .hex data file to process
  XMLCON_file   the .XMLCON data file to process

optional arguments:
  -h, --help    show this help message and exit
  -d, --debug   display debug messages
  -r, --raw     return the raw data values
  -o dest_dir   location to save output files
```

#### Utility for extracting Niskin bottle firing related data from converted csv-formatted text files

```
usage: bottle.py [-h] [-o output file] [-d] cnv_file

positional arguments:
  cnv_file  the converted csv-formatted file to process

optional arguments:
  -h, --help      show this help message and exit
  -o output_file  name and location of output file
  -d, --debug     display debug messages
```

#### Utility for processing ctd sensor related data from converted csv-formatted text files

```
usage: odf_process_ctd.py [-h] [-d] [-i cnv_file] [-o dest_dir] ini_file

positional arguments:
  ini_file     the .ini file to use for processing

optional arguments:
  -h, --help   show this help message and exit
  -d, --debug  display debug messages
  -i cnv_file  the converted, csv-formatted ctd data to process, this
               argument overrides any data file specified in the .ini file
  -o dest_dir  location to save output files
```

#### Sample Script for importing converted SBE Data
```
usage: sampleImport.py [-h] [-d] converted_File

positional arguments:
  cnv_file  the converted, csv-formatted data file to process

optional arguments:
  -h, --help      show this help message and exit
  -d, --debug     display debug messages
```  

#### Sample Script for exporting Pandas dataframes to csv-formatted text files with the CLIVAR 2-row header record.

```
usage: sampleExport.py [-h] [-d] [-w] converted_file output_file

positional arguments:
  cnv_file     the converted, csv-formatted data file to import to a dataframe
  output_file  the filename to export the dataframe to

optional arguments:
  -h, --help       show this help message and exit
  -d, --debug      display debug messages
  -w, --overwrite  overwrite the pre-existing output file if found
```

#### Sample Script for converting raw SBE Data
```
usage: sampleConvert.py [-h] [-d] hex_file XMLCON_file

positional arguments:
  hex_file     the .hex data file to process
  XMLCON_file  the .XMLCON data file to process

optional arguments:
  -h, --help   show this help message and exit
  -d, --debug  display debug messages
``` 
-->

***
## LICENSING
BSD 3-clause
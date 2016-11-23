# General Utilities for processing SBE Data

### Dependencies:
 - python3 (need minimum version info)
 - numpy (need minimum version info)
 - scipy (need minimum version info)
 - gsw (need minimum version info)
 - matplotlib (need minimum version info)
 - pandas (>= v0.18.1)

###Overview of available scripts
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

### Additional information
 - Please take a look at [INSTALL.md](./INSTALL.md) for instructions on how to install the ORD-CTD-PROC Utilites.
 - Please take a look at [SETUP.md](./SETUP.md) for recommendation on how to setup the SBE CTD Hardware and SBE Data Acquisition software.
 - Please take a look at [COOKBOOK.md](./INSTALL.md) for suggestions on how to leverage the ORD-CTD-PROC Utilites in a real-world, at-sea scenario.

 

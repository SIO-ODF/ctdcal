# ODF-CTD-PROC Cookbook

This cookbook attempts to provide instrutions on how to use the ODF-CTD-PROC Utlities for interpreting and processing raw data files from Seabird Electronic SBE9, SBE9+ and similar CTD sensor packages.  This cookbook also provides some best practices for using ODF-CTD-PROC at-sea during CTD-heavy cruises where adopting a more assembly-line type operating model is required to maintain sanity.

For full usage statments for all the scripts please refer to [README.md](./README.md)

## Before you begin:
 - Please take a look at [INSTALL.md](./INSTALL.md) for instructions on how to install the ORD-CTD-PROC Utilites.
 - Please take a look at [SETUP.md](./SETUP.md) for recommendation on how to setup the SBE CTD Hardware and SBE Data Acquisition software.

## Defining the Challanges and End Goals:
For the purposes of this cookbook the challenge trying to be overcome is that we are on a 3-week cruise and the operational goals are to conduct a CTD cast every 4 hours, 24-hours a day.  For some of the casts there are specific locations where the cast is to be done.  For a percentage of the casts the location will be determined based on the data from previous casts.

Some of the data used for selecting a location of the opportunistic casts will be from in situ sensors on the CTD package.  Additional data will come from analysis of water samples (Surprise! CTD package includes a 24-bottle rosette!) and the water chemists will require the position, depth and environmental conditions at the location of collection.

The CTD operations team wants the resulting data to be organized in directories by cast.  The name of the directory is the same as the cast.  The naming convention for the cast is \<CRUISE_ID\>\_\<CAST###\>.  The cruise ID for this cruise is RR2001.  The cast number is a 3-digit incrementing number starting at 001.

## A sample post-cast processing workflow:
1. Convert the raw .hex data to scientific units:

   `python3 ./odf_convert_sbe.py -d -o ./procdata/RR2001_CAST001 ./rawdata/RR2001_CAST001.hex ./rawdata/RR2001.XMLCON`

   Breakdown:
  - odf_convert_sbe.py - The raw data conversion script.
  - -d - Show verbose messages (optional but possibly useful if there's a problem)
  - -o ./procdata/RR2001_CAST001 - where to save the converted data file.
  - ./rawdata/RR2001_CAST001.hex - the raw hex file
  - ./rawdata/RR2001.XMLCON - the cast configuration package.
  
  The output of this script is the RR2001_CAST001_cnv.csv file. This file is a csv-formatted text file with a 2-row header.  The first row is the data type (t = temp, p = pressure, latitude, longitude, etc) and the units of measure (C, mBar, ddeg, etc).  The second row is the datatype of the data from the standpoint of a computer (int, float, string, datetime, etc).
  
  The RR2001_CAST001_cnv.csv file is saved to ./procdata/RR2001_CAST001 per the specification of the -o command-line argument.

2. Process the cnv file to extract Niskin bottle firing data:

   `python3 ./odf_process_bottle.py -d -o ./procdata/RR2001_CAST001 ./procdata/RR2001_CAST001/RR2001_CAST001_cnv.csv`
   
   Breakdown:
  - odf_process_bottle.py - The data processing script.
  - -d - Show verbose messages (optional but possibly useful if there's a problem)
  - -o ./procdata/RR2001_CAST001 - where to save the converted data file.
  - ./procdata/RR2001_CAST001/RR2001_CAST001_cnv.csv - the cnv file to process

  The output of this script is comprised of 3 new files.  RR2001_CAST001_btl.csv, RR2001_CAST001_btl_mean.csv, RR2001_CAST001_btl_median.csv.  RR2001_CAST001_btl.csv includes all rows from the RR2001_CAST001_cnv.csv input file where the bottle firing command was detected and the sequential number of the firing (which hopefully matches the bottle num).  Given the high data collection rate (upto 24Hz) the data acquistion system will record multiple records for a single bottle firing event.  RR2001_CAST001_btl_mean.csv takes the average of those multiple recorded firing event and return a single record per bottle firing.  RR2001_CAST001_btl_median.csv takes the median of those multiple recorded firing event and return a single record per bottle firing.

  The RR2001_CAST001_btl.csv, RR2001_CAST001_btl_mean.csv and RR2001_CAST001_btl_median.csv files are saved to ./procdata/RR2001_CAST001 per the specification of the -o command-line argument.  In this particular case the -o command-line argument was unecessary as the default behavior of this script is to save the output in the same directory as the cnv file.

3. Process the sensor data in the cnv file:

   `python3 ./odf_process_ctd.py -d -i ./procdata/RR2001_CAST001/RR2001_CAST001_cnv.csv ./procdata/ini-files/configuration.ini`

   Breakdown:
  - odf_process_ctd.py - The data processing script.
  - -d - Show verbose messages (optional but possibly useful if there's a problem)
  - -i ./procdata/RR2001_CAST001/RR2001_CAST001_cnv.csv - the input cnv file to process
  - ./procdata/ini-files/configuration.ini - the ini file to use for processing
  
  Currently there is no output from this script as it is still under heavy development... stay tuned.
   
### Results of this workflow:
After running these three scripts there will be a new directory: ./procdata/RR2001_CAST001.  Within this directory will be the following files:
 - RR2001_CAST001_cnv.csv
 - RR2001_CAST001_btl.csv
 - RR2001_CAST001_btl_mean.csv
 - RR2001_CAST001_btl_median.csv
 
 ### Additional steps:
 
1. Copy the raw data files for the cast to the output directory
 Per the desired goals of the CTD team (all data organized by cast) the raw data needs to be copied to the resulting output directory.  This is a 2-command task:

 ```
 cp ./rawdata/RR2001_CAST001.hex ./procdata/RR2001_CAST001/ 
 cp ./rawdata/RR2001.XMLCON ./procdata/RR2001_CAST001/RR2001_CAST001.XMLCON
 ```

 Take note that the second copy command also renames the file.  The updated name reflects both the cruise and cast.  Taking the time to create per-cast XMLCON files can prevent some data management headaches down the line.  Should the CTD sensor loadout change over time, having per-cast versions of the sensor package loadout provides some great documentation.  It also allows the raw data to be reprocessed again if there were any problems the first time that weren't noticed until later in the cruise.

2. Copy the ini configuration file for the cast processing to the output directory
 As with the previous step, making a per-cast copy of the ini file helps keep all cast-related files in a single directory and also allows a cast to be reprocessed, tweaked and reprocessed again without disturbing the original ini file.

 ```
 cp ./ini-files/configuration.ini ./procdata/RR2001_CAST001/RR2001_CAST001_config.ini
 ```

 To reprocess the cnv file for this cast simply run the following command:

 ```
 python3 ./odf_process_ctd.py -d -i ./procdata/RR2001_CAST001/RR2001_CAST001_cnv.csv ./procdata/RR2001_CAST001/RR2001_CAST001_config.ini
 ```
 
## Final Thoughts
We hope this cookbook is helpful.  Please report any errors using GitHub's issue tracker mechanism. 

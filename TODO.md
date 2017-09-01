Setup script:
data directory outside of


generate (based on configuration.ini?):

- raw_data_directory : data/raw/
- time_data_directory: data/time/
- pressure_data_directory: data/pressure/
- fit_data_directory: data/fit/
- log_directory: data/logs/
- salt_directory: data/salt/
- oxygen_directory: data/oxygen/
- bottle_directory: data/bottle/
- reft_directory: data/reft/
- o2flask_file: data/oxygen/o2flasks.vol (needs to be moved into)

standardize flag names with directory names (oxy vs oxygen, cond vs conductivity, etc)
(standardize on WOCE parameter names/relationships where possible, but in a long form)

config and/or constants file

integrate following libraries:
- Software:
  - logging
  - pathlib
  - unittest
  - click* (third party)
- Science:
  - gsw (full integration across all modules)

Rereckon all module names to something better

Setup sphinx documentation

Separate all conversions (salts, oxygen) from merging/calibration routines

Binary search algorithm when fitting T/C - start looking for 2nd degree fit at 1000 db,
then try upwards with min=500 to start. minimum 50 db bins. Optimize based on residuals on all points
and on residuals against y=0 (differences are 0) for pressure range (1000:6000)
possibly fold in some of the data above the pressure range, just to make sure it's not overly throwing out the top

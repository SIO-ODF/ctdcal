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

integrate following libraries:
- Software:
  - logging
  - pathlib
  - unittest
  - click* (third party)
- Science:
  - gsw (full integration across all modules)

Rereckon all module names to something better

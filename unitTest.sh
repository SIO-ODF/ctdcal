#!/bin/bash

#General Testing
python3 ./odf_convert_sbe.py -h
python3 ./odf_process_bottle.py -h
python3 ./odf_process_ctd.py -h
python3 ./sampleConvert.py -h
python3 ./sampleImport.py -h
python3 ./sampleExport.py -h

python3 ./odf_convert_sbe.py -d -o ./unitTesting ./sample_data_dir/raw/GS3601101.hex ./sample_data_dir/raw/GS3601101.XMLCON

python3 ./odf_process_bottle.py -d -o ./unitTesting ./unitTesting/GS3601101_cnv.csv

python3 ./odf_process_ctd.py -d -i ./unitTesting/GS3601101_cnv.csv ./sample_data_dir/ini-files/configuration.ini

python3 ./sampleConvert.py -d ./sample_data_dir/raw/GS3601101.hex ./sample_data_dir/raw/GS3601101.XMLCON

python3 ./sampleImport.py -d -w ./unitTesting/GS3601101_cnv.csv

python3 ./sampleExport.py -d -w ./unitTesting/GS3601101_cnv.csv ./unitTesting/sampleExport_output.csv

ls ./unitTesting

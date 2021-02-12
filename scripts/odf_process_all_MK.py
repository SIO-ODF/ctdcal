"""
Attempt to write a cleaner processing script from scratch.
"""

# import necessary packages
import sys

from ctdcal import get_ctdcal_config
import ctdcal.convert as convert
import ctdcal.fit_ctd as fit_ctd
import ctdcal.odf_io as odf_io
import ctdcal.oxy_fitting as oxy_fitting
import ctdcal.process_bottle as process_bottle
import ctdcal.process_ctd as process_ctd


def process_all():

    # TODO: document which functions/scripts produce which files

    #####
    # Step 0: Load and define necessary variables
    #####

    cfg = get_ctdcal_config()

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####
    # load station/cast list from file
    try:
        ssscc_list = process_ctd.get_ssscc_list()
    except FileNotFoundError:
        print("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = process_ctd.make_ssscc_list()

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list)

    # process time files
    convert.make_time_files(ssscc_list)

    # process bottle file
    convert.make_btl_mean(ssscc_list)

    # generate salt .csv files
    odf_io.process_salts(ssscc_list)

    # generate reftemp .csv files
    process_bottle.process_reft(ssscc_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list, cfg.btl_cols)

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)

    # calibrate temperature against reference
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)

    # calibrate temperature against reference
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)

    # calculate params needs for oxy/rinko calibration
    oxy_fitting.prepare_oxy(btl_data_all, time_data_all, ssscc_list)

    # calibrate oxygen against reference
    oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)

    # TODO: calibrate rinko against reference, similar to oxy_fitting.calibrate_oxy()
    # rinko.calibrate_oxy()  # or something

    #####
    # Step 3: export data
    #####

    # export time data to _ct1.csv format
    # TODO: clean this up more
    process_ctd.export_ct1(time_data_all, ssscc_list)
    process_bottle.export_hy1(btl_data_all)

    # run: ctd_to_bottle.py


def main(argv):
    """Run everything."""
    process_all()


if __name__ == "__main__":
    main(sys.argv[1:])

"""
Process all ICES data in accordance with data requirements and different discrete/continuous data formats.

Call with cli option "ctdcal process -g ices"
"""

import logging

# import needed ctdcal modules
from ctdcal import (
    convert,
    fit_ctd,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
)
from ctdcal.common import load_user_config, validate_file

log = logging.getLogger(__name__)

USERCONFIG = 'ctdcal/cfg.yaml'
user_cfg = load_user_config(validate_file(USERCONFIG))


def ices_process_basic():
    """
    Basic processing using the beta methods.
    * No discrete fitting - ignore bullets 3, 4
    * Alternative file writeouts -> must submit as .CNV files
        * Final data (before QA and calibration?)
        * 1 db binned .CNV
    * Include all metadata and intermediate files (ew pickles, report data, ssscc lists)
    """

    print("If you can read this, then it means that the ices_process() script is running.")

    #####
    # Step 0: Load and define necessary variables
    #####

    # cfg = get_ctdcal_config()

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    prefix = "CE17007_" #   Could add this to user config file?
    try:
        ssscc_list = process_ctd.get_ssscc_list()
    except FileNotFoundError:
        log.info("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = process_ctd.make_ssscc_list(prefix=prefix)

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list, prefix=prefix)
    #   Should no longer need the prefix - files are saved as SSSCC pickles

    # process time files
    convert.make_time_files(ssscc_list, user_cfg.datadir, user_cfg)

    # process bottle file
    convert.make_btl_mean(ssscc_list)

    #   ICES not done until instructed
    # generate salt .csv files
    # odf_io.process_salts(ssscc_list, user_cfg)

    #   ICES not done until instructed
    # generate reftemp .csv files
    # process_bottle.process_reft(ssscc_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    #   ICES bottle data not done until instructed
    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    # btl_data_all = process_bottle.load_all_btl_files(ssscc_list)

    #   ICES CTD data begins/terminates in the water, do not do this
    # process pressure offset
    # TODO: these functions return an updated dataframe, which we aren't
    #   assigning or reassigning to anything. Instead we trust that the
    #   updates which happen in the other module are visible by this one
    #   too (they  indeed seem to be). Is this a safe assumption?
    # process_ctd.apply_pressure_offset(btl_data_all)
    # process_ctd.apply_pressure_offset(time_data_all)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)

    # calibrate temperature against reference
    # fit_ctd.calibrate_temp(btl_data_all, time_data_all)

    # calibrate conductivity against reference
    # btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all, user_cfg, 'salt')

    # calculate params needs for oxy/rinko calibration
    # oxy_fitting.prepare_oxy(btl_data_all, time_data_all, ssscc_list, user_cfg, 'oxygen')

    # calibrate oxygen against reference
    # oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    # rinko.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)

    #####
    # Step 3: export data
    #####

    #   ICES requires submissions as .CNV file, which has not been historically supported by ODF
    process_ctd.export_ct_as_cnv(time_data_all)

    # export files for making cruise report figs
    # process_bottle.export_report_data(btl_data_all)

    # export to Exchange format
    # process_ctd.export_ct1(time_data_all, ssscc_list)
    # process_bottle.export_hy1(btl_data_all)

    # run: ctd_to_bottle.py


if __name__ == "__main__":
    ices_process_basic()

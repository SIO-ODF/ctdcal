"""
CTD and bottle processing routine for use on WHOI OSNAP cruise 32 (2022).
Engineering-focused cruise, highlighting cal dips, mooring casts, and CTD 

This script is a modified version of odf_process_all

Aaron Mau, 2022
"""

# import needed ctdcal modules
from .. import (
    convert,
    fit_ctd,
    get_ctdcal_config,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
)

import logging

log = logging.getLogger(__name__)


def whoi_process_all():

    #   Considering changing the order of some of these steps

    #####
    # Step 0: Load and define necessary variables
    #####

    # cfg = get_ctdcal_config()

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    try:
        ssscc_list = process_ctd.get_ssscc_list()
    except FileNotFoundError:
        log.info("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = process_ctd.make_ssscc_list()

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list)

    # process time files
    convert.make_time_files(ssscc_list)

    # process bottle file
    convert.make_btl_mean(
        ssscc_list
    )  #   Note that bottles may not be fired when you expect them

    # generate salt .csv files
    #   Portasal files will require a new input method
    #   Unlikely that standards will be able to be run before and after a run

    odf_io.process_salts(ssscc_list)

    # generate reftemp .csv files
    #   NO SBE35: Perhaps flag by sensor differences from each other (flagging primary only?)
    # process_bottle.process_reft(ssscc_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list)

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)

    # calibrate temperature against reference
    #   No reference temperature: Generate additional figures and assign coeffs manually?
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)

    # calibrate temperature against reference
    #   Run through using the flags only
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)

    # calculate params needs for oxy/rinko calibration
    #   OXY titration data will be rare (2-4 points per cast), so fits may be looser
    oxy_fitting.prepare_oxy(btl_data_all, time_data_all, ssscc_list)

    # calibrate oxygen against reference
    oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    # rinko.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)  #   No rinko

    #####
    # Step 3: export data
    #####

    # export files for making cruise report figs
    process_bottle.export_report_data(btl_data_all)  #   Does this need to happen?

    #   Generate generic .csv/.netcdf of the data for others to use
    #       (bottle and continuous)
    #   Generate plots similar to what Leah has produced


if __name__ == "__main__":
    whoi_process_all()

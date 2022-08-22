"""
CTD and bottle processing routine for use on WHOI OSNAP cruise 32 (2022).
Engineering-focused cruise, highlighting cal dips, mooring casts, and CTD 

This script is a modified version of odf_process_all deveoped by SIO-ODF.

Aaron Mau, 2022
"""

# import needed ctdcal modules
from ctdcal import (
    convert,
    fit_ctd,
    get_ctdcal_config,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    salts_io,
    flagging,
)

import logging

log = logging.getLogger(__name__)


def whoi_process_all(group="WHOI"):
    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv)
    #####

    # load station/cast list from file
    try:
        ssscc_list = process_ctd.get_ssscc_list()
    except FileNotFoundError:
        log.info("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = process_ctd.make_ssscc_list()
        ssscc_list = [
            n.replace("ar69-03_", "") for n in ssscc_list
        ]  #   OSNAP default is just cast number CCC

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list)

    # process time files
    convert.make_time_files(ssscc_list, group)

    # process bottle file
    convert.make_btl_mean(ssscc_list)

    # generate salt .csv files
    salts_io.osnap_salts(ssscc_list)

    #####
    # Step 2: calibrate conductivity and oxygen
    #####

    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)

    #   Flag the temp data relative to itself, don't fit it w/o reference. Using ODF threshold defaults.
    time_data_all["CTDTMP_FLAG_W"] = flagging.by_residual(
        time_data_all.CTDTMP1, time_data_all.CTDTMP2, time_data_all.CTDPRS
    )
    btl_data_all["CTDTMP_FLAG_W"] = flagging.by_residual(
        btl_data_all.CTDTMP1, btl_data_all.CTDTMP2, btl_data_all.CTDPRS
    )  #   These columns should go in the output file.

    #   Calibrate conductivity
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)

    # calculate params needs for oxy calibration
    #   OXY titration data will be rare (2-4 points per cast), so fits may be looser
    oxy_fitting.prepare_oxy(btl_data_all, time_data_all, ssscc_list)

    # calibrate oxygen against reference
    oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)

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

"""
Process all CTD and bottle data using ODF routines.
"""

# import needed ctdcal modules
from ctdcal import (
    convert,
    fit_ctd,
    get_ctdcal_config,
    salts_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
)

import logging
cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def odf_process_all():

    #####
    # Step 0: Load and define necessary variables
    #####

    # cfg = get_ctdcal_config()

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####
    print("Beginning CTDCAL run...")
    # load station/cast list from file
    ssscc_odf = process_ctd.get_ssscc_list(fname="data/ssscc/ssscc_odf.csv")    #   For the sake of running different routines
    ssscc_gtc = process_ctd.get_ssscc_list(fname="data/ssscc/ssscc_gtc.csv")
    ssscc_all = sorted(ssscc_odf + ssscc_gtc)

    print("Converting new .hex files...")
    # convert raw .hex files
    convert.hex_to_ctd(ssscc_odf)
    convert.hex_to_ctd(ssscc_gtc)

    # process time files
    convert.make_time_files(ssscc_all)

    # process bottle file
    convert.make_btl_mean(ssscc_all)

    # generate reftemp .csv files
    process_bottle.process_reft(ssscc_all)

    # generate salt .csv files
    # odf_io.process_salts(ssscc_list)
    salts_io.portasal_salts(ssscc_all)

    # Generate oxygen .csv files (splitting up sources)
    oxy_fitting.winkler_to_csv(ssscc_all)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    print("Loading ODF ct and bottle data...")
    time_data_odf = process_ctd.load_all_ctd_files(ssscc_odf)
    btl_data_odf  = process_bottle.load_all_btl_files(ssscc_odf)
    print("Loading GTC ct and bottle data...")
    time_data_gtc = process_ctd.load_all_ctd_files(ssscc_gtc)
    btl_data_gtc  = process_bottle.load_all_btl_files(ssscc_gtc)

    import xarray as xr
    time_pre_xr = time_data_all.to_xarray()
    time_pre_xr.to_netcdf(path=cfg.dirs["converted"]+"all_ct1.nc")

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_odf, mode="by_ssscc")
    process_ctd.apply_pressure_offset(btl_data_gtc, mode="by_ssscc")
    process_ctd.apply_pressure_offset(time_data_odf, mode="by_ssscc")
    process_ctd.apply_pressure_offset(time_data_gtc, mode="by_ssscc")

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)
    print("Calibrating temperature...")
    # calibrate temperature against reference
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)
    print("Calibrating conductivity...")
    # calibrate temperature against reference
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)

    print("Preparing and calibrating oxygen...")
    # calculate params needs for oxy/rinko calibration
    # TODO: move density matching to prepare_oxy
    # oxy_fitting.prepare_oxy(btl_data_all, time_data_all, ssscc_list)

    # calibrate oxygen against reference
    # oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    print("And now the ODF RINKO...")
    # rinko.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)

    #####
    # Step 3: export data
    #####
    print("Exporting data products...")
    # export files for making cruise report figs
    process_bottle.export_report_data(btl_data_all)

    # export to Exchange format
    # TODO: clean this up more
    # process_ctd.export_ct1(time_data_all, ssscc_list)
    process_bottle.export_hy1(btl_data_all)

    # run: ctd_to_bottle.py


if __name__ == "__main__":
    odf_process_all()

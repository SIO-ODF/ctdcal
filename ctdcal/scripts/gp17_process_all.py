"""
Process all CTD and bottle data from ODF/GTC/etc. platforms.
"""

# import needed ctdcal modules
from .. import (
    get_ctdcal_config,
    convert,
    fit_ctd,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
)

import logging
from pathlib import Path

log = logging.getLogger(__name__)
cfg_GTC = get_ctdcal_config("config_GTC.py")


def gp17_process_all():

    # ctdcal --debug process -g GTC

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    ssscc_list = process_ctd.get_ssscc_list()
    GTC_ssscc_list = process_ctd.get_ssscc_list("data/ssscc_GTC.csv")
    # TODO: GTC 00502 was aborted cast, how to deal with this?

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list)
    convert.hex_to_ctd([f"GTC_{ssscc}" for ssscc in GTC_ssscc_list])

    # process time files
    convert.make_time_files(ssscc_list)
    convert.make_time_files(GTC_ssscc_list, cfg=cfg_GTC)

    # process bottle file
    convert.make_btl_mean(ssscc_list)
    convert.make_btl_mean(GTC_ssscc_list, cfg=cfg_GTC)

    # generate salt .csv files
    # salt file naming is inconsistent due to only ~12 sample per rosette
    # just convert whatever is in the salt folder
    salt_list = sorted([f.name for f in Path("data/salt/").glob("?????")])
    odf_io.process_salts(salt_list)

    # generate reftemp .csv files
    process_bottle.process_reft(ssscc_list)

    # generate oxygen .csv files
    oxy_list = sorted([f.name for f in Path("data/oxygen/").glob("?????")])
    oxy_fitting.convert_winkler_oxy(oxy_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list)
    time_data_GTC = process_ctd.load_all_ctd_files(GTC_ssscc_list, cfg=cfg_GTC)
    btl_data_GTC = process_bottle.load_all_btl_files(GTC_ssscc_list, cfg=cfg_GTC)

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_all, ssscc_list)
    process_ctd.apply_pressure_offset(time_data_all, ssscc_list)
    process_ctd.apply_pressure_offset(btl_data_GTC, GTC_ssscc_list)
    process_ctd.apply_pressure_offset(time_data_GTC, GTC_ssscc_list)

    # create cast depth log file
    import pandas as pd
    process_ctd.make_depth_log(pd.concat([time_data_all, time_data_GTC]))

    # calibrate temperature against reference
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)
    time_data_GTC["CTDTMP_FLAG_W"] = 2

    # calibrate temperature against reference
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)
    btl_data_GTC, time_data_GTC = fit_ctd.calibrate_cond(
        btl_data_GTC, time_data_GTC, cfg=cfg_GTC
    )

    # calculate params needs for oxy/rinko calibration
    # TODO: move density matching to prepare_oxy
    oxy_fitting.prepare_oxy(btl_data_all, time_data_all)
    oxy_fitting.prepare_oxy(btl_data_GTC, time_data_GTC, cfg=cfg_GTC)

    # calibrate oxygen against reference
    oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    rinko.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)

    #####
    # Step 3: export data
    #####

    # export files for making cruise report figs
    process_bottle.export_report_data(btl_data_all)

    # export to Exchange format
    # TODO: clean this up more
    process_ctd.export_ct1(time_data_all, ssscc_list)
    process_bottle.export_hy1(btl_data_all)
    process_ctd.export_ct1(time_data_GTC, GTC_ssscc_list, cfg=cfg_GTC)
    process_bottle.export_hy1(btl_data_GTC, cfg=cfg_GTC)

    # run: ctd_to_bottle.py


if __name__ == "__main__":
    gp17_process_all()

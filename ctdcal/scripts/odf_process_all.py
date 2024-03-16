"""
Process all CTD and bottle data using ODF routines.
"""

# import needed ctdcal modules
from .. import (
    convert,
    fit_ctd,
    # get_ctdcal_config,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
)

from ctdcal.config import Parameters

import logging

log = logging.getLogger(__name__)

SSSCC = 'data/ssscc.csv'
SSSCC_HAS_OXY = 'data/ssscc_oxy.csv'

parameter_names = [
    'ssscc',
    'ssscc_oxy'
]

def odf_process_all():

    #####
    # Step 0: Load and define necessary variables
    #####

    cfg = Parameters(parameter_names)
    cfg.create_dict()
    cfg.ssscc = cfg.from_csv(SSSCC)
    cfg.ssscc_oxy = cfg.from_csv(SSSCC_HAS_OXY)

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    # try:
    #     ssscc_list = process_ctd.get_ssscc_list()
    # except FileNotFoundError:
    #     log.info("No ssscc.csv file found, generating from .hex file list")
    #     ssscc_list = process_ctd.make_ssscc_list()

    # convert raw .hex files
    convert.hex_to_ctd(cfg.ssscc)

    # process time files
    convert.make_time_files(cfg.ssscc)

    # process bottle file
    convert.make_btl_mean(cfg.ssscc)

    # generate salt .csv files
    odf_io.process_salts(cfg.ssscc)

    # generate reftemp .csv files
    process_bottle.process_reft(cfg.ssscc)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(cfg.ssscc)
    btl_data_all = process_bottle.load_all_btl_files(cfg.ssscc)
    # spam = time_data_all.rolling(24, min_periods=24, step=24).mean(numeric_only=False)

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)

    # calibrate temperature against reference
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)

    # calibrate conductivity against reference
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)

    # calculate params needs for oxy/rinko calibration
    # TODO: move density matching to prepare_oxy
    oxy_fitting.prepare_oxy(btl_data_all, time_data_all, cfg.ssscc)

    # calibrate oxygen against reference
    oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, cfg)

    rinko.calibrate_oxy(btl_data_all, time_data_all, cfg)
    # spam = time_data_all[time_data_all['SSSCC'] == '00601']

    #####
    # Step 3: export data
    #####

    # export files for making cruise report figs
    # process_bottle.export_report_data(btl_data_all)

    # export to Exchange format
    # TODO: clean this up more
    process_ctd.export_ct1(time_data_all, cfg)
    process_bottle.export_hy1(btl_data_all)

    # run: ctd_to_bottle.py


if __name__ == "__main__":
    odf_process_all()

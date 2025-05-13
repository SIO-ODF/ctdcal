"""
Process all CTD and bottle data using ODF routines.
"""
import logging

from ctdcal.common import load_user_config, validate_file, get_ssscc_list, make_ssscc_list
from ctdcal.fitting.fit_ctd import calibrate_temp, calibrate_cond, load_all_ctd_files, apply_pressure_offset
from ctdcal.fitting.fit_oxy import calibrate_oxy
from ctdcal.fitting import fit_oxy_rinko as rinko
from ctdcal.formats.exchange import export_hy1, export_ct1
from ctdcal.processors.cast_tools import make_time_files
from ctdcal.processors.convert_legacy import hex_to_ctd
from ctdcal.processors.proc_bottle import make_btl_mean, load_all_btl_files
from ctdcal.processors.proc_oxy_ctd import prepare_oxy
from ctdcal.processors.proc_reft import process_reft
from ctdcal.processors.proc_salt_odf import process_salts
from ctdcal.reporting.report_odf import make_depth_log

log = logging.getLogger(__name__)

USERCONFIG = 'ctdcal/cfg.yaml'
user_cfg = load_user_config(validate_file(USERCONFIG))


def odf_process_all():

    #####
    # Step 0: Load and define necessary variables
    #####

    # cfg = get_ctdcal_config()

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    try:
        ssscc_list = get_ssscc_list()
    except FileNotFoundError:
        log.info("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = make_ssscc_list()

    # convert raw .hex files
    hex_to_ctd(ssscc_list)

    # process time files
    make_time_files(ssscc_list, user_cfg.datadir, user_cfg)

    # process bottle file
    make_btl_mean(ssscc_list)

    # generate salt .csv files
    process_salts(ssscc_list, user_cfg)

    # generate reftemp .csv files
    process_reft(ssscc_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    time_data_all = load_all_ctd_files(ssscc_list)
    btl_data_all = load_all_btl_files(ssscc_list)

    # process pressure offset
    # TODO: these functions return an updated dataframe, which we aren't
    #   assigning or reassigning to anything. Instead we trust that the
    #   updates which happen in the other module are visible by this one
    #   too (they  indeed seem to be). Is this a safe assumption?
    apply_pressure_offset(btl_data_all)
    apply_pressure_offset(time_data_all)

    # create cast depth log file
    make_depth_log(time_data_all)

    # calibrate temperature against reference
    calibrate_temp(btl_data_all, time_data_all)

    # calibrate conductivity against reference
    btl_data_all, time_data_all = calibrate_cond(btl_data_all, time_data_all, user_cfg, 'salt')

    # calculate params needs for oxy/rinko calibration
    prepare_oxy(btl_data_all, time_data_all, ssscc_list, user_cfg, 'oxygen')

    # calibrate oxygen against reference
    calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    rinko.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)

    #####
    # Step 3: export data
    #####

    # export files for making cruise report figs
    # process_bottle.export_report_data(btl_data_all)

    # export to Exchange format
    export_ct1(time_data_all, ssscc_list)
    export_hy1(btl_data_all)

    # run: ctd_to_bottle.py


if __name__ == "__main__":
    odf_process_all()

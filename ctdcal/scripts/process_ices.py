#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
from pathlib import Path

from ctdcal import (
    convert,
    process_ctd, process_bottle,
)
from ctdcal.common import load_user_config, validate_file
from ctdcal.fit_ctd import calibrate_temp, calibrate_cond
from ctdcal.oxy_fitting import prepare_oxy_ices, calibrate_oxy

from ctdcal.parsers.all_bottle_xlsx import parse_discrete
from ctdcal.process_bottle import export_hy1, export_report_data
from ctdcal.process_ctd import export_ct1

log = logging.getLogger(__name__)

USERCONFIG = '/Users/als026/data/ices/ices.yaml'
cfg = load_user_config(validate_file(USERCONFIG))
INST = 'ctd'
BOTTLEDATAFILE = '/Users/als026/data/ices/discrete/CE17007_Bottle_Data_with_Metadata.xlsx'


def main():
    # load station/cast list from file
    try:
        ssscc_list = process_ctd.get_ssscc_list(cfg.datadir)
    except FileNotFoundError:
        log.info("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = process_ctd.make_ssscc_list(cfg.datadir, INST)
    print('Station list loaded.')

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list, INST, cfg)
    print('All HEX files converted.')

    # process time files
    convert.make_time_files(ssscc_list, cfg.datadir, INST, cfg)
    print('All time files created.')

    # export "Stage 1" data
    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list, cfg.datadir, INST)
    process_ctd.export_ct_as_cnv(time_data_all, cfg.datadir, INST)
    print('All stage 1 data exported.')

    # make bottle files
    convert.make_btl_mean(ssscc_list, INST, cfg)
    print('All bottle files created.')

    # generate reftemp .csv files
    process_bottle.process_reft(ssscc_list, cfg.datadir, 'reft')
    print('All refT data parsed.')

    # parse bottle data excel file into salt csv files
    parse_discrete(BOTTLEDATAFILE, Path(cfg.datadir, 'converted', 'salt'), 'salt', ssscc_list,
                   'Bench Salinity', 'SALNTY', cast_id_col='Cast', btlnum_col='Niskin')
    print('All bottle salinity data parsed.')

    # parse bottle data excel file into oxygen csv files
    parse_discrete(BOTTLEDATAFILE, Path(cfg.datadir, 'converted', 'oxygen'), 'oxygen', ssscc_list,
                   'Winkler DO umol/kg', 'OXYGEN', cast_id_col='Cast', btlnum_col='Niskin')
    print('All bottle oxygen data parsed.')

    # load all bottle data
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list, cfg.datadir, INST, 'reft', 'salt', 'oxygen', )
    print('All bottle data loaded.')
    print('Found %s samples.' % len(btl_data_all))

    # calculate and apply average pressure offset
    process_ctd.apply_pressure_offset(btl_data_all, cfg.datadir)
    process_ctd.apply_pressure_offset(time_data_all, cfg.datadir)
    print('Pressure offsets applied')

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all, cfg.datadir)
    print('Depth log saved.')

    # calibrate temperature against reference
    calibrate_temp(btl_data_all, time_data_all, cfg.datadir, INST, ssscc_list)
    print('Temperature fitting complete.')

    # calibrate conductivity against reference
    btl_data_all, time_data_all = calibrate_cond(btl_data_all, time_data_all, cfg.datadir, INST, 'salt', ssscc_list, cfg.bottleflags_man)
    print('Conductivity fitting complete.')

    # OXY FIT HERE
    prepare_oxy_ices(btl_data_all, time_data_all, cfg.datadir, INST, 'oxygen', cfg.bottleflags_man)
    calibrate_oxy(btl_data_all, time_data_all, cfg.datadir, INST, 'oxygen', ssscc_list)
    print('Oxygen fitting complete.')

    # export to Exchange format
    print('Exporting to exchange...')
    export_ct1(time_data_all, cfg.datadir, INST, ssscc_list)
    export_hy1(btl_data_all, cfg.datadir, INST)

    # export files for making cruise report figs
    export_report_data(btl_data_all, cfg.datadir, INST)

if __name__ == "__main__":
    main()

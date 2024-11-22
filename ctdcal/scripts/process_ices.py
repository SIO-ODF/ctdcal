#!/usr/bin/env python
# -*- coding: utf-8 -*-
from ctdcal import (
    convert,
    process_ctd, process_bottle,
)
from ctdcal.common import load_user_config, validate_file

import logging

from ctdcal.parsers.all_bottle_xlsx import parse_salinity

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

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list, INST, cfg)

    # process time files
    convert.make_time_files(ssscc_list, cfg.datadir, INST, cfg)

    # export "Stage 1" data
    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list, cfg.datadir, INST)
    process_ctd.export_ct_as_cnv(time_data_all, cfg.datadir, INST)

    # make bottle files
    convert.make_btl_mean(ssscc_list, INST, cfg)

    # generate reftemp .csv files
    process_bottle.process_reft(ssscc_list, cfg.datadir, 'reft')

    # parse bottle data excel file into salt csv files
    parse_salinity(BOTTLEDATAFILE, cfg.datadir, 'salt', ssscc_list, cast_id_col='Cast', btlnum_col='Niskin', sal_col='Bench Salinity')


if __name__ == "__main__":
    main()

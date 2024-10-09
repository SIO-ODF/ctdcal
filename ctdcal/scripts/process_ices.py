#!/usr/bin/env python
# -*- coding: utf-8 -*-
from ctdcal import (
    convert,
    process_ctd,
)
from ctdcal.common import load_user_config, validate_file

import logging


log = logging.getLogger(__name__)

USERCONFIG = '/Users/als026/data/ices/ices.yaml'
cfg = load_user_config(validate_file(USERCONFIG))
INST = 'ctd'

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

    #####
    # Step 3: export data
    #####



if __name__ == "__main__":
    main()

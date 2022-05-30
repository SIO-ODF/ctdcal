"""
P02 Process the BIO casts, which run before the full cast, on
full cast coeffs.

Bio casts:
* Typically run about 1000 m
    * The altimeter is almost never used
* Have no Winkler/Autosal procedures done
* Tend to have most bottles fired at the surface in quick bursts
    * Therefore the refT can be sporadic in later bottles (<15 sec)

Order of operations:
* ctdcal process -> processes the full cast and generates fit coeffs
* ctdcal bio     -> processes bio casts using the previous ctdcal run
"""

from ctdcal import (
    convert,
    fit_ctd,
    get_ctdcal_config,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
    sbe_reader as sbe_rd,
)

from pathlib import Path
import logging
import pandas as pd

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

def odf_process_bio():
    """
    Rebuilding based on procedure outlined in odf_process_all.py
    """
    #   Step 0: Initialize
    
    #   Step 1: Import the fit coeffs for T, C, oxy, and Rinko
    #   Note: The refT for the bio casts exist and can be refit.
    fname_t1    = "data/logs/fit_coef_t1.csv"
    t1_coefs    = pd.read_csv(fname_t1)
    fname_t2    = "data/logs/fit_coef_t2.csv"
    t2_coefs    = pd.read_csv(fname_t2)
    fname_c1    = "data/logs/fit_coef_c1.csv"
    c1_coefs    = pd.read_csv(fname_c1)
    fname_c2    = "data/logs/fit_coef_c2.csv"
    c2_coefs    = pd.read_csv(fname_c2)
    fname_sbe4  = "data/logs/sbe43_coefs.csv"
    sbe43_coefs = pd.read_csv(fname_sbe4)
    fname_rinko = "data/logs/rinko_coefs.csv"
    rinko_coefs = pd.read_csv(fname_rinko)
    log.info("All coefs loaded.")

    #   Step 2: Make the bio_ssscc
    #   Should keep this separate from the full ssscc_list
    #   to prevent anything from trying to be fit without having
    #   data that can be used to constrain the fit
    fname_ssscc = "data/ssscc_bio.csv"
    ssscc_list = []
    with open(fname_ssscc, "r") as filename:
        ssscc_list = [line.strip() for line in filename]

    #   Step 3: Convert the hex into bottle averages and time .pkl
    log.info("Converting .hex files")
    for ssscc in ssscc_list:
        #   Stash .pkl files in the normal folders
        if not Path(cfg.dirs["converted"] + ssscc + ".pkl").exists():
            hexFile = "data/raw_bio/" + ssscc + ".hex"  #   Keep the raw files separated for safety
            xmlconFile = "data/raw_bio/" + ssscc + ".XMLCON"
            sbeReader = sbe_rd.SBEReader.from_paths(hexFile, xmlconFile)
            converted_df = convert.convertFromSBEReader(sbeReader, ssscc)
            converted_df.to_pickle(cfg.dirs["converted"] + ssscc + ".pkl")

    convert.make_time_files(ssscc_list)
    convert.make_btl_mean(ssscc_list)
    log.info('All .pkl files generated.')

    #   Step 4: Get the time and bottle data into the workspace
    time_data_all   = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all    = process_bottle.load_all_btl_files(ssscc_list) #   Extra cols should be full of NaNs
    log.info('Time and Bottle data loaded.')
    print(time_data_all.head())
    print(btl_data_all.head())

    #   Step 5: Apply the pressure offset

    #   Step 6: Make the depth log

    #   Step 7: Calibrate the temperature
    #   Note: The refT for the bio casts exist and can be refit.

    #   Step 8: Calibrate the cond w/ pre-existing coeffs

    #   Step 9: Calibrate the SBE43 w/ pre-existing coeffs

    #   Step 10: Calibrate the Rinko w/ pre-existing coeffs

    #   Step 11: Export the ct1 files

    #   Step 12: Export/merge the hy1 file

def load_all_bottle_files_bio():
    """
    A modified version of process_bottle.load_all_btl_files
    This skips REFT
    """
    
    pass

if __name__ == "__main__":
    odf_process_bio()
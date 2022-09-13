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
    oxy_fitting,
    process_bottle,
    process_ctd,
    salts_io,
    osnap_oxy,
    flagging,
    ctd_plots,
)

import pandas as pd

import logging

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def whoi_process_all(group="WHOI"):
    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv)
    #####
    print("Beginning processing. Checking files for conversion...")
    # load station/cast list from file
    try:
        #   Station 086 = 086B, do not use 086 hex/xmlcon
        ssscc_list = process_ctd.get_ssscc_list()
        microcat_list = process_ctd.get_ssscc_list(fname="data/ssscc_microcat.csv")
    except FileNotFoundError:
        log.info("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = process_ctd.make_ssscc_list()
        ssscc_list = [
            n.replace("ar69-03_", "") for n in ssscc_list
        ]  #   OSNAP default is just cast number CCC

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list, group)

    # process time files
    convert.make_time_files(ssscc_list, group, microcat_list)

    # process bottle file
    convert.make_btl_mean(ssscc_list)

    # generate salt .csv files
    salts_io.osnap_salts(ssscc_list)

    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list)
    print("Files loaded in to station", ssscc_list[-1])

    # time_microcat = process_ctd.load_all_ctd_files(microcat_list) #   Deal with microcats seperately

    #####
    # Step 2: calibrate conductivity and oxygen
    #####

    print("Beginning data processing...")

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

    btl_data_all, time_data_all = osnap_oxy.ctd_oxy_converter(
        btl_data_all, time_data_all
    )

    btl_data_prefit = (
        btl_data_all.copy()
    )  #   Stashing for prefit vs postfit (could reduce size)
    time_data_prefit = time_data_all.copy()

    btl_data_fit, time_data_fit = osnap_oxy.osnap_oxy_main(
        btl_data_all, time_data_all, ssscc_list
    )

    print("Data fitting and flagging complete.")

    #####
    # Step 3: export data
    #####

    print("Writing out data products...")

    try:
        import xarray as xr
        import datetime

        depth_df = pd.read_csv(
            cfg.dirs["logs"] + "depth_log.csv", dtype={"SSSCC": str}, na_values=-999
        ).dropna()
        manual_depth_df = pd.read_csv(
            cfg.dirs["logs"] + "manual_depth_log.csv", dtype={"SSSCC": str}
        )
        full_depth_df = pd.concat([depth_df, manual_depth_df])
        full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)
        btl_data_fit["DEPTH"] = -999
        for index, row in full_depth_df.iterrows():
            btl_data_fit.loc[btl_data_fit["SSSCC"] == row["SSSCC"], "DEPTH"] = int(
                row["DEPTH"]
            )

        btl_data_fit["DateTime"] = btl_data_fit.nmea_datetime  #   For MatLAB < 2020
        time_data_fit["DateTime"] = time_data_fit.nmea_datetime
        btl_data_fit = process_bottle.add_btlnbr_cols(
            btl_data_fit, btl_num_col="btl_fire_num"
        )  #   BTLNBR int
        btl_data_fit = btl_data_fit.rename(
            columns={"SALNTY": "BTL_SAL", "OXYGEN": "BTL_OXY"}
        )
        outfile = cfg.dirs["pressure"] + "bottle_data"
        btl_cols = {
            "SSSCC": "Station",
            "DateTime": "",
            "GPSLAT": "Dec Degrees",
            "GPSLON": "Dec Degrees",
            "BTLNBR": "",
            "CTDPRS": "DBAR",
            "CTDTMP1": "ITS-90",
            "CTDTMP2": "ITS-90",
            "CTDTMP_FLAG_W": "",
            "CTDCOND1": "mS/cm",
            "CTDCOND2": "mS/cm",
            "CTDSAL": "PSS-78",
            "CTDSAL_FLAG_W": "",
            "BTL_SAL": "PSS-78",
            "CTDOXY": "UMOL/KG",
            "CTDOXY_FLAG_W": "",
            "BTL_OXY": "UMOL/KG",
            "CTDOXY1": "ML/L",
            "CTDOXYVOLTS": "0-5VDC",
            "ALT": "M",
            "CTDFLUOR": "mg/m^3",
            "TURBIDITY": "0-5VDC",
            "CTDXMISS": "0-5VDC",
            "FLUOR_CDOM": "0-5VDC",
        }
        with open(outfile + ".csv", mode="w+") as f:
            f.write(",".join(btl_cols.keys()) + "\n")
            f.write(",".join(btl_cols.values()) + "\n")
            btl_data_fit[btl_cols.keys()].to_csv(f, header=False, index=False)
        save_btl = btl_data_fit[btl_cols.keys()].to_xarray()
        save_btl.to_netcdf(path=outfile + ".nc")
        # print("Exporting continuous .csv files...")
        # process_ctd.export_ct1(time_data_all, ssscc_list)
        time_cols = [
            "SSSCC",
            "DateTime",
            "GPSLAT",
            "GPSLON",
            "CTDPRS",
            "CTDTMP1",
            "CTDTMP_FLAG_W",
            "CTDCOND1",
            "CTDCOND2",
            "CTDSAL",
            "CTDSAL_FLAG_W",
            "CTDOXY",
            "CTDOXY_FLAG_W",
            "CTDOXY1",
            "CTDOXYVOLTS",
            "ALT",
            "CTDFLUOR",
            "TURBIDITY",
            "CTDXMISS",
            "FLUOR_CDOM",
        ]
        for ssscc in ssscc_list:
            print("Writing time file", ssscc)
            time_out = cfg.dirs["pressure"] + ssscc + "_profile.nc"
            time_ssscc = time_data_fit.loc[time_data_fit.SSSCC == ssscc]
            time_ssscc = time_ssscc[time_cols].to_xarray()
            time_ssscc.to_netcdf(path=time_out)

        print("Exporting OSNAP data suite figures...")
        ctd_plots.osnap_suite(
            btl_data_prefit, btl_data_fit, time_data_prefit, time_data_fit
        )
    except:
        print("Could not export final data.")


if __name__ == "__main__":
    whoi_process_all()

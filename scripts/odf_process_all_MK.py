"""
Attempt to write a cleaner processing script from scratch.
"""

# import necessary packages
import os
import glob  # replace with pathlib everywhere
import sys
import pathlib
import subprocess
import pandas as pd
import numpy as np
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import ctdcal.oxy_fitting as oxy_fitting
import gsw
import ctdcal.rinko as rinko
import ctdcal.odf_io as odf_io
import scripts.odf_sbe_metadata as odf_sbe_metadata


def process_all():

    # TODO: document which functions/scripts produce which files

    #####
    # Step 0: Load and define necessary variables
    #####

    import config as cfg

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    ssscc_list = []
    with open(cfg.directory["ssscc_file"], "r") as filename:
        ssscc_list = [line.strip() for line in filename]

    # convert hex to ctd (TODO: convert this to function form)
    cnv_dir_list = os.listdir("data/converted/")
    for ssscc in ssscc_list:
        if "{}.pkl".format(ssscc) not in cnv_dir_list:
            subprocess.run(
                [
                    "odf_convert_sbe.py",
                    "data/raw/" + ssscc + ".hex",
                    "data/raw/" + ssscc + ".XMLCON",
                    "-o",
                    "data/converted",
                ],
                stdout=subprocess.PIPE,
            )
            print("odf_convert_sbe.py SSSCC: " + ssscc + " done")

    # first half of CTD data processing
    # TODO: clean up and document better
    # TODO: export to odf_io instead of standalone script?
    # this produces #_ct1.csv (preliminary), #_time.pkl, and "ondeck_pressure.csv"
    time_dir_list = os.listdir("data/time/")
    for ssscc in ssscc_list:
        if "{}_time.pkl".format(ssscc) not in time_dir_list:
            odf_sbe_metadata.main("data/converted/" + ssscc + ".pkl")
            print("odf_sbe_metadata.py SSSCC: " + ssscc + " done")

    # process bottle file (TODO: convert this to function form)
    btl_dir_list = os.listdir("data/bottle/")
    for ssscc in ssscc_list:
        if "{}_btl_mean.pkl".format(ssscc) not in btl_dir_list:
            subprocess.run(
                [
                    "odf_process_bottle.py",
                    "data/converted/" + ssscc + ".pkl",
                    "-o",
                    "data/bottle/",
                ],
                stdout=subprocess.PIPE,
            )
            print("odf_process_bottle.py SSSCC: " + ssscc + " done")

    # generate salt .csv files
    odf_io.process_salts(ssscc_list)

    # generate reftemp .csv files
    process_ctd.process_reft(ssscc_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    # TODO: clean up process_ctd.load_all_ctd_files
    btl_data_all = process_ctd.load_all_ctd_files(ssscc_list, "bottle", cfg.btl_cols)
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list, "time", None)

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)

    # calibrate temperature against reference
    fit_ctd.fit_temp(btl_data_all, time_data_all)

    # calibrate temperature against reference
    fit_ctd.fit_cond(btl_data_all, time_data_all)

    ####################
    # oxygen calibration
    ####################
    # TODO: export to single line function
    # calculate sigma
    btl_data_all["sigma_btl"] = oxy_fitting.sigma_from_CTD(
        btl_data_all[cfg.column["sal_btl"]],
        btl_data_all[cfg.column["t1_btl"]],  # oxygen sensor is on primary line (ie t1)
        btl_data_all[cfg.column["p_btl"]],
        btl_data_all[cfg.column["lon_btl"]],
        btl_data_all[cfg.column["lat_btl"]],
    )
    time_data_all["sigma_ctd"] = oxy_fitting.sigma_from_CTD(
        time_data_all[cfg.column["sal"]],
        time_data_all[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_data_all[cfg.column["p"]],
        time_data_all[cfg.column["lon_btl"]],
        time_data_all[cfg.column["lat_btl"]],
    )

    # Calculate SA and CT
    btl_data_all["SA"] = gsw.SA_from_SP(
        btl_data_all[cfg.column["sal_btl"]],
        btl_data_all[cfg.column["p_btl"]],
        btl_data_all[cfg.column["lon_btl"]],
        btl_data_all[cfg.column["lat_btl"]],
    )
    btl_data_all["CT"] = gsw.CT_from_t(
        btl_data_all["SA"],
        btl_data_all[cfg.column["t1_btl"]],  # oxygen sensor is on primary line (ie t1)
        btl_data_all[cfg.column["p_btl"]],
    )
    time_data_all["SA"] = gsw.SA_from_SP(
        time_data_all[cfg.column["sal"]],
        time_data_all[cfg.column["p"]],
        time_data_all[cfg.column["lon_btl"]],
        time_data_all[cfg.column["lat_btl"]],
    )
    time_data_all["CT"] = gsw.CT_from_t(
        time_data_all["SA"],
        time_data_all[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_data_all[cfg.column["p"]],
    )

    # Calculate oxygen solubility in µmol/kg
    btl_data_all["OS_btl"] = gsw.O2sol(  # any reason to label as OS_btl? not really..
        btl_data_all["SA"],
        btl_data_all["CT"],
        btl_data_all[cfg.column["p_btl"]],
        btl_data_all[cfg.column["lon_btl"]],
        btl_data_all[cfg.column["lat_btl"]],
    )
    time_data_all["OS_ctd"] = gsw.O2sol(  # any reason to label as OS_ctd?
        time_data_all["SA"],
        time_data_all["CT"],
        time_data_all[cfg.column["p"]],
        time_data_all[cfg.column["lon"]],
        time_data_all[cfg.column["lat"]],
    )

    # Calculate bottle oxygen
    btl_data_all[cfg.column["oxy_btl"]] = oxy_fitting.calculate_bottle_oxygen(
        ssscc_list,
        btl_data_all["SSSCC"],
        btl_data_all["TITR_VOL"],
        btl_data_all["TITR_TEMP"],
        btl_data_all["FLASKNO"],
    )
    btl_data_all[cfg.column["oxy_btl"]] = oxy_fitting.oxy_ml_to_umolkg(
        btl_data_all[cfg.column["oxy_btl"]], btl_data_all["sigma_btl"]
    )
    btl_data_all["OXYGEN_FLAG_W"] = oxy_fitting.flag_winkler_oxygen(
        btl_data_all[cfg.column["oxy_btl"]]
    )

    # Prep vars, dfs, etc.
    all_sbe43_merged = pd.DataFrame()
    # all_rinko_merged = pd.DataFrame()
    # rinko_dict = {}
    sbe43_dict = {}
    all_sbe43_fit = pd.DataFrame()
    # rinko_flag = pd.DataFrame()

    # Density match time/btl oxy dataframes
    for ssscc in ssscc_list:

        time_data = time_data_all[time_data_all["SSSCC"] == ssscc].copy()
        btl_data = btl_data_all[btl_data_all["SSSCC"] == ssscc].copy()

        if (btl_data["OXYGEN_FLAG_W"] == 9).all():
            # rinko_dict[ssscc] = np.full(8, np.nan)
            sbe43_dict[ssscc] = np.full(5, np.nan)
            print(ssscc + " skipped, all oxy data is NaN")
            continue

        sbe43_merged = oxy_fitting.match_sigmas(
            btl_data[cfg.column["p_btl"]],
            btl_data[cfg.column["oxy_btl"]],
            btl_data["sigma_btl"],
            btl_data["btl_fire_num"],  # used for sorting later
            time_data["sigma_ctd"],
            time_data["OS_ctd"],
            time_data[cfg.column["p"]],
            time_data[cfg.column["t1"]],
            time_data[cfg.column["oxyvolts"]],
            time_data["scan_datetime"],
            time_data["SSSCC"],
        )

        # rinko_merged = rinko.match_sigmas(
        #     btl_data[cfg.column["p_btl"]],
        #     btl_data[cfg.column["oxy_btl"]],
        #     btl_data["sigma_btl"],
        #     time_data["sigma_ctd"],
        #     time_data["OS_ctd"],
        #     time_data[cfg.column["p"]],
        #     time_data[cfg.column["t1"]],
        #     time_data[cfg.column["rinko_oxy"]],
        #     btl_data["SSSCC"],
        # )

        all_sbe43_merged = pd.concat([all_sbe43_merged, sbe43_merged])
        # all_rinko_merged = pd.concat([all_rinko_merged, rinko_merged])
        print(ssscc + " density matching done")

    # Fit ALL oxygen stations together to get initial coefficient guess
    (sbe_coef0, _) = oxy_fitting.sbe43_oxy_fit(all_sbe43_merged)
    # (rinko_coef0, _) = rinko.rinko_oxygen_fit(all_rinko_merged)

    # Fit oxygen stations using SSSCC chunks to refine coefficients
    ssscc_files_ox = sorted(glob.glob("data/ssscc/ssscc_ox*.csv"))
    for f in ssscc_files_ox:
        ssscc_list_ox = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()

        sbe_coef, sbe_df = oxy_fitting.sbe43_oxy_fit(
            all_sbe43_merged.loc[all_sbe43_merged["SSSCC_sbe43"].isin(ssscc_list_ox)],
            sbe_coef0=sbe_coef0,
            f_out=f,
        )

        # TODO: calculate RINKO coefs

        # build coef dictionary
        for ssscc in ssscc_list_ox:
            if ssscc not in sbe43_dict.keys():  # don't overwrite NaN'd stations
                sbe43_dict[ssscc] = sbe_coef

        # all non-NaN oxygen data with flags
        all_sbe43_fit = pd.concat([all_sbe43_fit, sbe_df])

    # TODO: secondary oxygen flagging step (instead of just taking outliers from fit routine)
    # TODO: save outlier data from fits?
    # # TODO: abstract to oxy_fitting.py
    # TODO: figure out what this code is/was
    # breakpoint()
    # all_sbe43_fit["SSSCC_int"] = all_sbe43_fit["SSSCC_sbe43"].astype(int)
    # all_sbe43_fit = all_sbe43_fit.sort_values(
    #     by=["SSSCC_int", "btl_fire_num"], ascending=[True, True]
    # )
    # all_sbe43_fit["STNNBR"] = all_sbe43_fit["SSSCC_sbe43"].str[0:3]  # SSS from SSSCC
    # all_sbe43_fit["CASTNO"] = all_sbe43_fit["SSSCC_sbe43"].str[3:]  # CC from SSSCC
    # all_sbe43_fit = all_sbe43_fit.rename(columns={"btl_fire_num": "SAMPNO"})
    # all_sbe43_fit = all_sbe43_fit[
    #     ["STNNBR", "CASTNO", "SAMPNO", "CTDOXY", "CTDOXY_FLAG_W"]
    # ]
    # all_sbe43_fit.to_csv(cfg.directory["logs"] + "quality_flag_sbe43.csv", index=False)

    # apply coefs
    time_data_all["CTDOXY"] = -999
    time_data_all["RINKO"] = -999
    for ssscc in ssscc_list:

        if np.isnan(sbe43_dict[ssscc]).all():
            print(ssscc + " missing oxy data, leaving -999 values and flagging as 9")
            time_data_all.loc[time_data_all["SSSCC"] == ssscc, "CTDOXY_FLAG_W"] = 9
            time_data_all.loc[time_data_all["SSSCC"] == ssscc, "RINKO_FLAG_W"] = 9
            continue

        btl_rows = (btl_data_all["SSSCC"] == ssscc).values
        time_rows = (time_data_all["SSSCC"] == ssscc).values

        btl_data_all.loc[btl_rows, "CTDOXY"] = oxy_fitting.PMEL_oxy_eq(
            sbe43_dict[ssscc],
            (
                btl_data_all.loc[btl_rows, cfg.column["oxyvolts"]],
                btl_data_all.loc[btl_rows, cfg.column["p_btl"]],
                btl_data_all.loc[btl_rows, cfg.column["t1_btl"]],
                btl_data_all.loc[btl_rows, "dv_dt"],
                btl_data_all.loc[btl_rows, "OS_btl"],
            ),
        )
        print(ssscc + " btl data fitting done")
        time_data_all.loc[time_rows, "CTDOXY"] = oxy_fitting.PMEL_oxy_eq(
            sbe43_dict[ssscc],
            (
                time_data_all.loc[time_rows, cfg.column["oxyvolts"]],
                time_data_all.loc[time_rows, cfg.column["p"]],
                time_data_all.loc[time_rows, cfg.column["t1"]],
                time_data_all.loc[time_rows, "dv_dt"],
                time_data_all.loc[time_rows, "OS_ctd"],
            ),
        )
        print(ssscc + " time data fitting done")

    # TODO: flag oxy data here? compare w/ T/C routines

    # export fitting coefs
    sbe43_coefs = pd.DataFrame.from_dict(
        sbe43_dict, orient="index", columns=["Soc", "Voffset", "Tau20", "Tcorr", "E"]
    )
    sbe43_coefs.to_csv(cfg.directory["logs"] + "sbe43_coefs.csv")

    ###############################
    # final cleanup/prep for export
    ###############################

    # clean up columns
    p_column_names = cfg.ctd_time_output["col_names"]

    # initial flagging (some of this should be moved)
    # TODO: CTDOXY flags
    time_data_all["CTDOXY_FLAG_W"] = 2
    # TODO: flag bad based on cond/temp and handcoded salt
    time_data_all["CTDSAL_FLAG_W"] = 2
    time_data_all["CTDTMP_FLAG_W"] = 2
    # TODO: lump all uncalibrated together; smart flagging like ["CTD*_FLAG_W"] = 1
    # maybe not always have these channels so don't hardcode them in
    time_data_all["CTDFLUOR_FLAG_W"] = 1
    time_data_all["CTDXMISS_FLAG_W"] = 1
    time_data_all["CTDBACKSCATTER_FLAG_W"] = 1

    # renames
    time_data_all = time_data_all.rename(
        columns={"CTDTMP1": "CTDTMP", "FLUOR": "CTDFLUOR"}
    )

    # TODO: see process_ctd.format_time_data()
    try:
        df = time_data_all[p_column_names].copy()
    except KeyError:
        print("Column names not configured properly... attempting to correct")
        df = pd.DataFrame()
        df["SSSCC"] = time_data_all["SSSCC"].copy()
        for col in p_column_names:
            try:
                df[col] = time_data_all[col]
            except KeyError:
                if col.endswith("FLAG_W"):
                    print(col + " missing, flagging with 9s")
                    df[col] = 9
                else:
                    print(col + " missing, filling with -999s")
                    df[col] = -999

    # needs depth_log.csv, manual_depth_log.csv
    print("Exporting *_ct1.csv files")
    process_ctd.export_ct1(
        df,
        ssscc_list,
        cfg.cruise["expocode"],
        cfg.cruise["sectionid"],
        cfg.ctd_serial,
        cfg.ctd_time_output["col_names"],
        cfg.ctd_time_output["col_units"],
    )

    # run: ctd_to_bottle.py

    breakpoint()

    # TODO: write/abstract code to a function (in process_ctd?)
    # TODO: look at Argo .nc files for inspiration on structuring
    #
    # experiment with xarray
    # import xarray  # not needed apparently? at least as currently coded

    da_out = df.to_xarray()  # xarray calls them DataArrays instead of DataFrames

    # set attributes
    da_out["CTDTMP"].attrs["long_name"] = "Temperature"
    da_out["CTDTMP"].attrs["units"] = "ITS-90"
    da_out["CTDTMP"].attrs["description"] = "Continuous temperature from CTD downcast"

    # can set attrs from dict
    prs_attrs = dict(
        long_name="Pressure", units="dbar", description="Continuous pressure"
    )
    tmp_attrs = dict(
        long_name="Temperature", units="ITS-90", description="Continuous temperature"
    )
    sal_attrs = dict(
        long_name="Salinity", units="PSS-78", description="Continuous salinity",
    )

    # can do one at at time
    da_out["CTDSAL"].attrs = sal_attrs

    # maybe do nested dicts in config.py and loop? e.g.:
    ctd_attrs = dict(
        CTDPRS=prs_attrs,
        CTDTMP=tmp_attrs,
        CTDSAL=sal_attrs,
        # CTDOXY=oxy_attrs,
        # CTDRINKO=rinko_attrs,
        # CTDXMISS=xmiss_attrs,
        # CTDFLUOR=fluor_attrs,
    )

    for var in da_out.keys():
        if var == "SSSCC":
            continue
        if not var.endswith("_FLAG_W"):
            da_out[var].attrs = ctd_attrs[var]

    # output files
    # don't actually run bc this is 24Hz data...
    # da_out.to_netcdf('example.nc')


def main(argv):
    """Run everything.
    """
    process_all()


if __name__ == "__main__":
    main(sys.argv[1:])

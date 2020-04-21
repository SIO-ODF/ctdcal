"""
Attempt to write a cleaner processing script from scratch.
"""

# import necessary packages
import os
import glob
import sys
import subprocess
import pandas as pd
import numpy as np
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import ctdcal.oxy_fitting as oxy_fitting
import gsw
import ctdcal.rinko as rinko
import odf_salt_parser as salt_parser
import odf_reft_parser as reft_parser


def process_all():

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

    # ??? (TODO: convert this to function form)
    time_dir_list = os.listdir("data/time/")
    for ssscc in ssscc_list:
        if "{}_time.pkl".format(ssscc) not in time_dir_list:
            subprocess.run(
                ["odf_sbe_metadata.py", "data/converted/" + ssscc + ".pkl"],
                stdout=subprocess.PIPE,
            )
            print("odf_sbe_metadata.py SSSCC: " + ssscc + " done")

    # process bottle file (TODO: convert this to function form)
    # does this generate "ondeck_pressure.csv"?
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

    # generate salt files
    salt_dir_list = os.listdir("data/salt/")
    for ssscc in ssscc_list:
        if "{}_salts.csv".format(ssscc) not in salt_dir_list:
            salt_parser.process_salts(ssscc, "data/salt/")

    # generate ref temp files
    reft_dir_list = os.listdir("data/reft/")
    for ssscc in ssscc_list:
        if "{}_reft.csv".format(ssscc) not in reft_dir_list:
            reft_parser.process_reft(ssscc, "data/reft/")

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    btl_data_all = process_ctd.load_all_ctd_files(ssscc_list, "bottle", cfg.btl_cols)
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list, "time", None)

    # process pressure offset
    pressure_log = process_ctd.load_pressure_logs("data/logs/ondeck_pressure.csv")
    p_offset = process_ctd.get_pressure_offset(
        pressure_log.ondeck_start_p, pressure_log.ondeck_end_p
    )

    btl_data_all = fit_ctd.apply_pressure_offset(
        btl_data_all, cfg.column["p"], p_offset
    )
    time_data_all = fit_ctd.apply_pressure_offset(
        time_data_all, cfg.column["p"], p_offset
    )

    # create depth_log.csv
    depth_dict = {}
    for ssscc in ssscc_list:
        print(ssscc)
        time_rows = time_data_all["SSSCC"] == ssscc
        max_depth = process_ctd.find_cast_depth(
            time_data_all.loc[time_rows, "CTDPRS"],
            time_data_all.loc[time_rows, "GPSLAT"],
            time_data_all.loc[time_rows, "ALT"],
        )
        depth_dict[ssscc] = max_depth

    depth_df = pd.DataFrame.from_dict(depth_dict, orient="index")
    depth_df.reset_index(inplace=True)
    depth_df.rename(columns={0: "DEPTH", "index": "STNNBR"}, inplace=True)
    depth_df.to_csv("data/logs/depth_log.csv", index=False)

    #########################
    # temperature calibration
    #########################

    ssscc_files_t = sorted(glob.glob("data/ssscc/ssscc_t*.csv"))
    qual_flag_t1 = pd.DataFrame()
    qual_flag_t2 = pd.DataFrame()
    coef_t1_all = pd.DataFrame()
    coef_t2_all = pd.DataFrame()

    for f in ssscc_files_t:
        # 0) grab ssscc chunk to fit
        ssscc_list_t = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()
        btl_rows = btl_data_all["SSSCC"].isin(ssscc_list_t).values
        time_rows = time_data_all["SSSCC"].isin(ssscc_list_t).values

        # 1) remove non-finite data
        df_temp_good = process_ctd.prepare_fit_data(
            btl_data_all[btl_rows], cfg.column["reft"],
        )

        # TODO: allow for cast-by-cast T_order/P_order/xRange
        # 2 & 3) flag outliers and calculate fit params on flag 2s
        coef_t1, df_ques_t1 = fit_ctd.get_T_coefs(
            df_temp_good[cfg.column["t1_btl"]],
            df_temp_good[cfg.column["reft"]],
            df_temp_good[cfg.column["p_btl"]],
            df_temp_good["SSSCC"],
            df_temp_good["btl_fire_num"],
            T_order=2,
            P_order=2,
            xRange="1000:5000",
        )
        coef_t2, df_ques_t2 = fit_ctd.get_T_coefs(
            df_temp_good[cfg.column["t2_btl"]],
            df_temp_good[cfg.column["reft"]],
            df_temp_good[cfg.column["p_btl"]],
            df_temp_good["SSSCC"],
            df_temp_good["btl_fire_num"],
            T_order=2,
            P_order=2,
            xRange="1000:5000",
        )

        # 4) apply fit
        btl_data_all.loc[btl_rows, cfg.column["t1_btl"]] = fit_ctd.temperature_polyfit(
            btl_data_all.loc[btl_rows, cfg.column["t1_btl"]],
            btl_data_all.loc[btl_rows, cfg.column["p_btl"]],
            coef_t1,
        )
        btl_data_all.loc[btl_rows, cfg.column["t2_btl"]] = fit_ctd.temperature_polyfit(
            btl_data_all.loc[btl_rows, cfg.column["t2_btl"]],
            btl_data_all.loc[btl_rows, cfg.column["p_btl"]],
            coef_t2,
        )
        time_data_all.loc[time_rows, cfg.column["t1"]] = fit_ctd.temperature_polyfit(
            time_data_all.loc[time_rows, cfg.column["t1"]],
            time_data_all.loc[time_rows, cfg.column["p"]],
            coef_t1,
        )
        time_data_all.loc[time_rows, cfg.column["t2"]] = fit_ctd.temperature_polyfit(
            time_data_all.loc[time_rows, cfg.column["t2"]],
            time_data_all.loc[time_rows, cfg.column["p"]],
            coef_t2,
        )

        # 5) handle quality flags
        qual_flag_t1 = pd.concat([qual_flag_t1, df_ques_t1])
        qual_flag_t2 = pd.concat([qual_flag_t2, df_ques_t2])

        # 6) handle fit params
        coef_t1_df = pd.DataFrame()
        coef_t1_df["SSSCC"] = ssscc_list_t
        coef_t2_df = coef_t1_df.copy()

        for idx, val in enumerate(coef_t1):
            # build df w/ columns c0, c1, c2, etc.
            coef_t1_df["c" + str(idx)] = coef_t1[idx]
            coef_t2_df["c" + str(idx)] = coef_t2[idx]

        coef_t1_all = pd.concat([coef_t1_all, coef_t1_df])
        coef_t2_all = pd.concat([coef_t2_all, coef_t2_df])

    # export temp quality flags
    qual_flag_t1.to_csv(cfg.directory["logs"] + "qual_flag_t1.csv", index=False)
    qual_flag_t2.to_csv(cfg.directory["logs"] + "qual_flag_t2.csv", index=False)

    # export temp fit params
    coef_t1_all.to_csv(cfg.directory["logs"] + "fit_coef_t1.csv", index=False)
    coef_t2_all.to_csv(cfg.directory["logs"] + "fit_coef_t2.csv", index=False)

    ##########################
    # conductivity calibration
    ##########################

    ssscc_files_c = sorted(glob.glob("data/ssscc/ssscc_*c*.csv"))
    qual_flag_c1 = pd.DataFrame()
    qual_flag_c2 = pd.DataFrame()
    coef_c1_all = pd.DataFrame()
    coef_c2_all = pd.DataFrame()

    # calculate BTLCOND values from autosal data
    # TODO: what temp sensor to use? should cfg.py have a var for which sensor is used in final data?
    btl_data_all[cfg.column["refc"]] = fit_ctd.CR_to_cond(
        btl_data_all["CRavg"],
        btl_data_all["BathTEMP"],
        btl_data_all[cfg.column["t1_btl"]],
        btl_data_all[cfg.column["p_btl"]],
    )
    # could use REFTMP instead of T1; testing this is a good project

    # TODO: make salt flagger move .csv somewhere else? or just always have it
    # somewhere else and read it from that location (e.g. in data/scratch_folder/salts)
    salt_file = "tools/salt_flags_handcoded.csv"  # abstract to config.py
    if glob.glob(salt_file):
        handcoded_salts = pd.read_csv(
            salt_file, dtype={"SSSCC": str, "salinity_flag": int}
        )
        handcoded_salts = handcoded_salts.rename(
            columns={"SAMPNO": "btl_fire_num", "salinity_flag": "SALNTY_FLAG_W"}
        ).drop(columns=["diff", "Comments"])
        btl_data_all = btl_data_all.merge(
            handcoded_salts, on=["SSSCC", "btl_fire_num"], how="left"
        )
        btl_data_all.loc[btl_data_all["BTLCOND"].isnull(), "SALNTY_FLAG_W"] = 9
        btl_data_all["SALNTY_FLAG_W"] = btl_data_all["SALNTY_FLAG_W"].fillna(
            2, downcast="infer"  # fill remaining NaNs with 2s and cast to dtype int
        )
    else:
        btl_data_all["SALNTY_FLAG_W"] = 2

    for f in ssscc_files_c:
        # 0) grab ssscc chunk to fit
        ssscc_list_c = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()
        btl_rows = (btl_data_all["SSSCC"].isin(ssscc_list_c).values) & (
            btl_data_all["SALNTY_FLAG_W"] == 2  # only use salts flagged good (e.g. 2)
        )
        time_rows = time_data_all["SSSCC"].isin(ssscc_list_c).values

        # 1) remove non-finite data
        df_cond_good = process_ctd.prepare_fit_data(
            btl_data_all[btl_rows], cfg.column["refc"],
        )

        # 2 & 3) calculate fit params
        coef_c1, df_ques_c1 = fit_ctd.get_C_coefs(
            df_cond_good[cfg.column["c1_btl"]],
            df_cond_good[cfg.column["refc"]],
            df_cond_good[cfg.column["t1_btl"]],
            df_cond_good[cfg.column["p_btl"]],
            df_cond_good["SSSCC"],
            df_cond_good["btl_fire_num"],
            P_order=2,
            T_order=0,
            C_order=1,
            xRange="1000:5000",
        )
        coef_c2, df_ques_c2 = fit_ctd.get_C_coefs(
            df_cond_good[cfg.column["c2_btl"]],
            df_cond_good[cfg.column["refc"]],
            df_cond_good[cfg.column["t1_btl"]],
            df_cond_good[cfg.column["p_btl"]],
            df_cond_good["SSSCC"],
            df_cond_good["btl_fire_num"],
            P_order=2,
            T_order=0,
            C_order=1,
            xRange="1000:5000",
        )

        # 4) apply fit
        btl_data_all.loc[btl_rows, cfg.column["c1_btl"]] = fit_ctd.conductivity_polyfit(
            btl_data_all.loc[btl_rows, cfg.column["c1_btl"]],
            btl_data_all.loc[btl_rows, cfg.column["t1_btl"]],
            btl_data_all.loc[btl_rows, cfg.column["p_btl"]],
            coef_c1,
        )
        btl_data_all.loc[btl_rows, cfg.column["c2_btl"]] = fit_ctd.conductivity_polyfit(
            btl_data_all.loc[btl_rows, cfg.column["c2_btl"]],
            btl_data_all.loc[btl_rows, cfg.column["t2_btl"]],
            btl_data_all.loc[btl_rows, cfg.column["p_btl"]],
            coef_c2,
        )
        time_data_all.loc[time_rows, cfg.column["c1"]] = fit_ctd.conductivity_polyfit(
            time_data_all.loc[time_rows, cfg.column["c1"]],
            time_data_all.loc[time_rows, cfg.column["t1"]],
            time_data_all.loc[time_rows, cfg.column["p"]],
            coef_c1,
        )
        time_data_all.loc[time_rows, cfg.column["c2"]] = fit_ctd.conductivity_polyfit(
            time_data_all.loc[time_rows, cfg.column["c2"]],
            time_data_all.loc[time_rows, cfg.column["t2"]],
            time_data_all.loc[time_rows, cfg.column["p"]],
            coef_c2,
        )

        # 5) handle quality flags
        qual_flag_c1 = pd.concat([qual_flag_c1, df_ques_c1])
        qual_flag_c2 = pd.concat([qual_flag_c2, df_ques_c2])

        # 6) handle fit params
        coef_c1_df = pd.DataFrame()
        coef_c1_df["SSSCC"] = ssscc_list_c
        coef_c2_df = coef_c1_df.copy()

        for idx, val in enumerate(coef_c1):
            coef_c1_df["c" + str(idx)] = coef_c1[idx]
            coef_c2_df["c" + str(idx)] = coef_c2[idx]

        coef_c1_all = pd.concat([coef_c1_all, coef_c1_df])
        coef_c2_all = pd.concat([coef_c2_all, coef_c2_df])

    # export cond quality flags
    qual_flag_c1.to_csv(cfg.directory["logs"] + "qual_flag_c1.csv", index=False)
    qual_flag_c2.to_csv(cfg.directory["logs"] + "qual_flag_c2.csv", index=False)

    # export cond fit params
    coef_c1_all.to_csv(cfg.directory["logs"] + "fit_coef_c1.csv", index=False)
    coef_c2_all.to_csv(cfg.directory["logs"] + "fit_coef_c2.csv", index=False)

    # recalculate salinity with calibrated C/T
    time_data_all[cfg.column["sal"]] = gsw.SP_from_C(
        time_data_all[cfg.column["c1"]],
        time_data_all[cfg.column["t1"]],
        time_data_all[cfg.column["p"]],
    )

    ####################
    # oxygen calibration
    ####################

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
    btl_data_all["OS_btl"] = gsw.O2sol(  # any reason to label as OS_btl?
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
    ssscc_files_ox = sorted(glob.glob("data/ssscc/ssscc_*ox*.csv"))
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

    # TODO: save outlier data from fits?
    # # TODO: abstract to oxy_fitting.py
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
    time_data_all["CTDSAL_FLAG_W"] = 2
    time_data_all["CTDTMP_FLAG_W"] = 2
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

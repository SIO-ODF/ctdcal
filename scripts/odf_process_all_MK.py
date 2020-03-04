"""
Attempt to write a cleaner processing script from scratch.
"""

# import necessary packages
import os
import glob
import timeit
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

    # start_time = timeit.default_timer()

    #####
    # Step 0: Load and define necessary variables
    #####

    import config as cfg

    # define cruise and file information/extensions
    prefix = "nbp1802_"

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    ssscc_list = []
    with open(cfg.directory["ssscc_file"], "r") as filename:
        ssscc_list = [line.strip() for line in filename]

    # make list of already converted files to skip later
    cnv_dir_list = os.listdir("data/converted/")
    time_dir_list = os.listdir("data/time/")
    btl_dir_list = os.listdir("data/bottle/")
    salt_dir_list = os.listdir("data/salt/")
    reft_dir_list = os.listdir("data/reft/")

    # convert hex to ctd (TODO: convert this to function form)
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
    for ssscc in ssscc_list:
        if "{}_time.pkl".format(ssscc) not in time_dir_list:
            subprocess.run(
                ["odf_sbe_metadata.py", "data/converted/" + ssscc + ".pkl"],
                stdout=subprocess.PIPE,
            )
            print("odf_sbe_metadata.py SSSCC: " + ssscc + " done")

    # process bottle file (TODO: convert this to function form)
    # does this generate "ondeck_pressure.csv"?
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
    for ssscc in ssscc_list:
        if "{}_salts.csv".format(ssscc) not in salt_dir_list:
            salt_parser.process_salts(ssscc, "data/salt/")

    # generate ref temp files
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

    # temperature calibration
    # steps: 1) remove non-finite data
    #        2) flag points w/ large deviations
    #        3) calculate fit parameters (on data w/ flag 2) -> save them too!
    #        4) apply fit
    #        5) qualify flag file

    ssscc_files_t = sorted(glob.glob("data/ssscc/ssscc_*t*.csv"))
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

        # 2 & 3) calculate fit params
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

    # conductivity calibration

    ssscc_files_c = sorted(glob.glob("data/ssscc/ssscc_*c*.csv"))
    qual_flag_c1 = pd.DataFrame()
    qual_flag_c2 = pd.DataFrame()
    coef_c1_all = pd.DataFrame()
    coef_c2_all = pd.DataFrame()

    # calculate BTLCOND values from autosal data
    # TODO: what temp sensor to use? should cfg.py have a var for which sensor is used in final data?
    # TODO: clean up fit_ctd.CR_to_cond
    btl_data_all[cfg.column["refc"]] = fit_ctd.CR_to_cond(
        btl_data_all["CRavg"],
        btl_data_all["BathTEMP"],
        btl_data_all["CTDTMP1"],  # could use REFTMP; testing this is a good project
        btl_data_all["CTDPRS"],
    )

    for f in ssscc_files_c:
        # 0) grab ssscc chunk to fit
        ssscc_list_c = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()
        btl_rows = btl_data_all["SSSCC"].isin(ssscc_list_c).values
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

    # print(timeit.default_timer() - start_time)
    # breakpoint()

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

    # Calculate SA and PT
    btl_data_all["SA"] = gsw.SA_from_SP(
        btl_data_all[cfg.column["sal_btl"]],
        btl_data_all[cfg.column["p_btl"]],
        btl_data_all[cfg.column["lon_btl"]],
        btl_data_all[cfg.column["lat_btl"]],
    )
    btl_data_all["PT"] = gsw.pt0_from_t(
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
    time_data_all["PT"] = gsw.pt0_from_t(
        time_data_all["SA"],
        time_data_all[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_data_all[cfg.column["p"]],
    )

    # Calculate oxygen solubility in µmol/kg
    btl_data_all["OS_btl"] = oxy_fitting.os_umol_kg(  # any reason to label as OS_btl?
        btl_data_all["SA"], btl_data_all["PT"]
    )
    time_data_all["OS_ctd"] = oxy_fitting.os_umol_kg(  # any reason to label as OS_ctd?
        time_data_all["SA"], time_data_all["PT"]
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
    rinko_coef0 = rinko.rinko_o2_cal_parameters()
    all_rinko_df = pd.DataFrame()
    all_sbe43_df = pd.DataFrame()
    rinko_dict = {}
    sbe43_dict = {}

    for ssscc in ssscc_list:

        time_data = time_data_all[time_data_all["SSSCC"] == ssscc].copy()
        btl_data = btl_data_all[btl_data_all["SSSCC"] == ssscc].copy()

        if (btl_data["OXYGEN_FLAG_W"] == 9).all():
            rinko_dict[ssscc] = np.full(8, np.nan)
            sbe43_dict[ssscc] = np.full(7, np.nan)
            print(ssscc + " skipped, all oxy data is NaN")
            continue

        rinko_coef, rinko_oxy_df = rinko.rinko_oxygen_fit(
            btl_data[cfg.column["p_btl"]],
            btl_data[cfg.column["oxy_btl"]],
            btl_data["sigma_btl"],
            time_data["sigma_ctd"],
            time_data["OS_ctd"],
            time_data[cfg.column["p"]],
            time_data[cfg.column["t1"]],  # which T to use?
            time_data[cfg.column["rinko_oxy"]],
            rinko_coef0,
            btl_data["SSSCC"],
        )

        hex_file = "data/raw/" + ssscc + ".hex"  # move these into get_SB_coef?
        xml_file = "data/raw/" + ssscc + ".XMLCON"  # and just call with ssscc
        sbe_coef0 = oxy_fitting.get_SB_coef(hex_file, xml_file)
        sbe_coef, sbe_oxy_df = oxy_fitting.sbe43_oxy_fit(
            btl_data[cfg.column["p_btl"]],
            btl_data[cfg.column["oxy_btl"]],
            btl_data["sigma_btl"],
            time_data["sigma_ctd"],
            time_data["OS_ctd"],
            time_data[cfg.column["p"]],
            time_data[cfg.column["t1"]],
            time_data[cfg.column["oxyvolts"]],
            time_data["scan_datetime"],
            sbe_coef0,
            btl_data["SSSCC"],
        )

        rinko_dict[ssscc] = rinko_coef
        sbe43_dict[ssscc] = sbe_coef
        all_rinko_df = pd.concat([all_rinko_df, rinko_oxy_df])
        all_sbe43_df = pd.concat([all_sbe43_df, sbe_oxy_df])
        print(ssscc + " Done!")

    sbe43_coef_df = oxy_fitting.create_coef_df(sbe43_dict)
    rinko_coef_df = oxy_fitting.create_coef_df(rinko_dict)

    btl_data_all = btl_data_all.merge(
        all_rinko_df,
        left_on=["SSSCC", cfg.column["p_btl"]],
        right_on=["SSSCC_rinko", "CTDPRS_rinko_btl"],
        how="left",
    )
    btl_data_all = btl_data_all.merge(
        all_sbe43_df,
        left_on=["SSSCC", cfg.column["p_btl"]],
        right_on=["SSSCC_sbe43", "CTDPRS_sbe43_btl"],
        how="left",
    )

    btl_data_all.drop(list(btl_data_all.filter(regex="rinko")), axis=1, inplace=True)
    btl_data_all.drop(list(btl_data_all.filter(regex="sbe43")), axis=1, inplace=True)


def main(argv):
    """Run everything.
    """
    process_all()


if __name__ == "__main__":
    main(sys.argv[1:])

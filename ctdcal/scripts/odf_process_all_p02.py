"""
Process all CTD and bottle data using ODF routines.
"""

import logging

import gsw

# import needed ctdcal modules
import pandas as pd

from .. import (
    convert,
    fit_ctd,
    get_ctdcal_config,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
)

log = logging.getLogger(__name__)
cfg = get_ctdcal_config()


def odf_process_all_p02():

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    try:
        ssscc_list = process_ctd.get_ssscc_list()
        ssscc_bio_list = process_ctd.get_ssscc_list("data/ssscc_bio.csv")
    except FileNotFoundError:
        log.info("No ssscc.csv file found, generating from .hex file list")
        ssscc_list = process_ctd.make_ssscc_list()

    # convert raw .hex files
    convert.hex_to_ctd(ssscc_list)
    log.info("processing bio next:")
    convert.hex_to_ctd(ssscc_bio_list)

    # process time files
    convert.make_time_files(ssscc_list)
    log.info("processing bio next:")
    convert.make_time_files(ssscc_bio_list)

    # process bottle file
    convert.make_btl_mean(ssscc_list)
    log.info("processing bio next:")
    convert.make_btl_mean(ssscc_bio_list)

    # generate salt .csv files
    odf_io.process_salts(ssscc_list)

    # generate reftemp .csv files
    process_bottle.process_reft(ssscc_list)
    log.info("processing bio next:")
    process_bottle.process_reft(ssscc_bio_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list)
    log.info("loading bio next:")
    time_data_bio = process_ctd.load_all_ctd_files(ssscc_bio_list)
    btl_data_bio = process_bottle.load_all_btl_files(ssscc_bio_list)

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)
    process_ctd.apply_pressure_offset(time_data_bio)
    process_ctd.apply_pressure_offset(btl_data_bio)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)
    # bio only goes to 1000m, will not get good depth so skip this

    # calibrate temperature against reference
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)
    log.info("applying T coefs to bio cast")
    apply_T_coefs(ssscc_bio_list, time_data_bio, btl_data_bio)

    # calibrate temperature against reference
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)
    log.info("applying C coefs to bio cast")
    apply_C_coefs(ssscc_bio_list, time_data_bio, btl_data_bio)

    # calculate params needs for oxy/rinko calibration
    # TODO: move density matching to prepare_oxy
    oxy_fitting.prepare_oxy(btl_data_all, time_data_all, ssscc_list)

    # calibrate oxygen against reference
    oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    # skip applying to SBE43 for now...
    rinko.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    log.info("applying Rinko coefs to bio cast")
    apply_Rinko_coefs(ssscc_bio_list, time_data_bio, btl_data_bio)

    #####
    # Step 3: export data
    #####
    # merge core and bio casts
    btl_data_merged = pd.concat([btl_data_all, btl_data_bio])
    time_data_merged = pd.concat([time_data_all, time_data_bio])
    ssscc_list_merged = ssscc_list + ssscc_bio_list

    # export files for making cruise report figs
    process_bottle.export_report_data(btl_data_merged)

    # export to Exchange format
    # TODO: clean this up more
    process_ctd.export_ct1(time_data_merged, ssscc_list_merged)
    process_bottle.export_hy1(btl_data_merged)

    # run: ctd_to_bottle.py


def apply_T_coefs(ssscc_list, time_df, btl_df):
    T1_coefs = pd.read_csv("data/logs/fit_coef_t1.csv")
    T2_coefs = pd.read_csv("data/logs/fit_coef_t2.csv")
    T1_coefs.index = T1_coefs["SSSCC"].astype(str).str.slice(stop=3)
    T2_coefs.index = T2_coefs["SSSCC"].astype(str).str.slice(stop=3)
    T1_dict, T2_dict = T1_coefs.to_dict("index"), T2_coefs.to_dict("index")
    for ssscc in ssscc_list:
        time_rows = time_df["SSSCC"] == ssscc
        btl_rows = btl_df["SSSCC"] == ssscc

        # build coef tuples
        T1_t_coefs = tuple(T1_dict[ssscc[:3]][x] for x in ["c0", "ct1", "ct2"])
        T1_p_coefs = tuple(T1_dict[ssscc[:3]][x] for x in ["cp1", "cp2"])
        T2_t_coefs = tuple(T2_dict[ssscc[:3]][x] for x in ["c0", "ct1", "ct2"])
        T2_p_coefs = tuple(T2_dict[ssscc[:3]][x] for x in ["cp1", "cp2"])

        # time
        time_df.loc[time_rows, cfg.column["t1"]] = fit_ctd.apply_polyfit(
            time_df.loc[time_rows, cfg.column["t1"]],
            T1_t_coefs,
            (time_df.loc[time_rows, cfg.column["p"]], T1_p_coefs),
        )
        time_df.loc[time_rows, cfg.column["t1"]] = fit_ctd.apply_polyfit(
            time_df.loc[time_rows, cfg.column["t1"]],
            T2_t_coefs,
            (time_df.loc[time_rows, cfg.column["p"]], T2_p_coefs),
        )

        # bottle
        btl_df.loc[btl_rows, cfg.column["t1"]] = fit_ctd.apply_polyfit(
            btl_df.loc[btl_rows, cfg.column["t1"]],
            T1_t_coefs,
            (btl_df.loc[btl_rows, cfg.column["p"]], T1_p_coefs),
        )
        btl_df.loc[btl_rows, cfg.column["t1"]] - fit_ctd.apply_polyfit(
            btl_df.loc[btl_rows, cfg.column["t1"]],
            T2_t_coefs,
            (btl_df.loc[btl_rows, cfg.column["p"]], T2_p_coefs),
        )

    # flag data
    time_df["CTDTMP_FLAG_W"] = 2
    btl_df["CTDTMP_FLAG_W"] = 2


def apply_C_coefs(ssscc_list, time_df, btl_df):
    C1_coefs = pd.read_csv("data/logs/fit_coef_c1.csv")
    C2_coefs = pd.read_csv("data/logs/fit_coef_c2.csv")
    C1_coefs.index = C1_coefs["SSSCC"].astype(str).str.slice(stop=3)
    C2_coefs.index = C2_coefs["SSSCC"].astype(str).str.slice(stop=3)
    C1_dict, C2_dict = C1_coefs.to_dict("index"), C2_coefs.to_dict("index")
    for ssscc in ssscc_list:
        time_rows = time_df["SSSCC"] == ssscc
        btl_rows = btl_df["SSSCC"] == ssscc

        # build coef tuples
        C1_c_coefs = tuple(C1_dict[ssscc[:3]][x] for x in ["c0", "cc1", "cc2"])
        C1_p_coefs = tuple(C1_dict[ssscc[:3]][x] for x in ["cp1", "cp2"])
        C1_t_coefs = tuple(C1_dict[ssscc[:3]][x] for x in ["ct1", "ct2"])
        C2_c_coefs = tuple(C2_dict[ssscc[:3]][x] for x in ["c0", "ct1", "ct2"])
        C2_p_coefs = tuple(C2_dict[ssscc[:3]][x] for x in ["cp1", "cp2"])
        C2_t_coefs = tuple(C2_dict[ssscc[:3]][x] for x in ["ct1", "ct2"])

        # time
        time_df.loc[time_rows, cfg.column["c1"]] = fit_ctd.apply_polyfit(
            time_df.loc[time_rows, cfg.column["c1"]],
            C1_c_coefs,
            (time_df.loc[time_rows, cfg.column["p"]], C1_p_coefs),
            (time_df.loc[time_rows, cfg.column["t1"]], C1_t_coefs),
        )
        time_df.loc[time_rows, cfg.column["c2"]] = fit_ctd.apply_polyfit(
            time_df.loc[time_rows, cfg.column["c2"]],
            C2_c_coefs,
            (time_df.loc[time_rows, cfg.column["p"]], C2_p_coefs),
            (time_df.loc[time_rows, cfg.column["t1"]], C2_t_coefs),
        )

        # bottle
        btl_df.loc[btl_rows, cfg.column["c1"]] = fit_ctd.apply_polyfit(
            btl_df.loc[btl_rows, cfg.column["c1"]],
            C1_c_coefs,
            (btl_df.loc[btl_rows, cfg.column["p"]], C1_p_coefs),
            (time_df.loc[time_rows, cfg.column["t1"]], C1_t_coefs),
        )
        btl_df.loc[btl_rows, cfg.column["c2"]] - fit_ctd.apply_polyfit(
            btl_df.loc[btl_rows, cfg.column["c2"]],
            C2_c_coefs,
            (btl_df.loc[btl_rows, cfg.column["p"]], C2_p_coefs),
            (time_df.loc[time_rows, cfg.column["t1"]], C2_t_coefs),
        )

    # recalculate salinity with calibrated C/T
    time_df[cfg.column["sal"]] = gsw.SP_from_C(
        time_df[cfg.column["c1"]],
        time_df[cfg.column["t1"]],
        time_df[cfg.column["p"]],
    )
    btl_df[cfg.column["sal"]] = gsw.SP_from_C(
        btl_df[cfg.column["c1"]],
        btl_df[cfg.column["t1"]],
        btl_df[cfg.column["p"]],
    )

    time_df["CTDSAL_FLAG_W"] = 2
    btl_df["CTDSAL_FLAG_W"] = 2


def apply_Rinko_coefs(ssscc_list, time_df, btl_df):
    R_coefs = pd.read_csv("data/logs/rinko_coefs.csv")
    R_coefs.index = R_coefs["Unnamed: 0"].astype(str).str.slice(stop=3)
    R_coefs.drop("Unnamed: 0", axis=1, inplace=True)
    R_dict = R_coefs.to_dict("index")

    # calculate oxygen solubility
    btl_df["SA"] = gsw.SA_from_SP(
        btl_df[cfg.column["sal"]],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["lon"]],
        btl_df[cfg.column["lat"]],
    )
    time_df["SA"] = gsw.SA_from_SP(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )
    btl_df["CT"] = gsw.CT_from_t(
        btl_df["SA"],
        btl_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p"]],
    )
    time_df["CT"] = gsw.CT_from_t(
        time_df["SA"],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
    )
    btl_df["OS"] = gsw.O2sol(
        btl_df["SA"],
        btl_df["CT"],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["lon"]],
        btl_df[cfg.column["lat"]],
    )
    time_df["OS"] = gsw.O2sol(
        time_df["SA"],
        time_df["CT"],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )

    for ssscc in ssscc_list:
        time_rows = time_df["SSSCC"] == ssscc
        btl_rows = btl_df["SSSCC"] == ssscc

        # build coef tuples
        rinko_coefs_ssscc = tuple(x for x in R_dict[ssscc[:3]].values())

        btl_df.loc[btl_rows, "CTDRINKO"] = rinko._Uchida_DO_eq(
            rinko_coefs_ssscc,
            (
                btl_df.loc[btl_rows, cfg.column["rinko_oxy"]],
                btl_df.loc[btl_rows, cfg.column["p"]],
                btl_df.loc[btl_rows, cfg.column["t1"]],
                btl_df.loc[btl_rows, cfg.column["sal"]],
                btl_df.loc[btl_rows, "OS"],
            ),
        )
        time_df.loc[time_rows, "CTDRINKO"] = rinko._Uchida_DO_eq(
            rinko_coefs_ssscc,
            (
                time_df.loc[time_rows, cfg.column["rinko_oxy"]],
                time_df.loc[time_rows, cfg.column["p"]],
                time_df.loc[time_rows, cfg.column["t1"]],
                time_df.loc[time_rows, cfg.column["sal"]],
                time_df.loc[time_rows, "OS"],
            ),
        )

    # flag
    time_df["CTDRINKO_FLAG_W"] = 2
    btl_df["CTDRINKO_FLAG_W"] = 2


if __name__ == "__main__":
    odf_process_all_p02()

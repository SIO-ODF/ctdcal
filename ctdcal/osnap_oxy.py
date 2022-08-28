"""
A modification of the ODF oxy routine for oxygen data that is already in units of umol/kg.

* Reads in CTD data and converts units to umol/kg
* Loads Winkler data and appends SSSCC and bottle number to correct locations in CTD bottle file
    * Does some quick flagging
* Execute calibrate_oxy to acquire coefficients

Developed for OSNAP32, 2022 with the intention of creating a "generic" oxy loader that can be merged into 
ctdcal's other modules.
"""
from ctdcal import get_ctdcal_config, oxy_fitting, flagging
import pandas as pd
import numpy as np
import gsw

pd.options.mode.chained_assignment = None
cfg = get_ctdcal_config()


def ctd_oxy_converter(btl_df, time_df):
    """Convert SBE43's output into umol/kg"""
    #   Extract the import bits of load_all_btl_files

    btl_df["SA"] = gsw.SA_from_SP(
        btl_df[cfg.column["sal"]],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["lon"]],
        btl_df[cfg.column["lat"]],
    )
    btl_df["CT"] = gsw.CT_from_t(
        btl_df["SA"],
        btl_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p"]],
    )
    time_df["SA"] = gsw.SA_from_SP(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )
    time_df["CT"] = gsw.CT_from_t(
        time_df["SA"],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
    )

    # calculate sigma
    btl_df["sigma_btl"] = gsw.sigma0(btl_df["SA"], btl_df["CT"])
    time_df["sigma_btl"] = gsw.sigma0(time_df["SA"], time_df["CT"])

    # Calculate oxygen solubility in µmol/kg
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
    # Convert CTDOXY units
    btl_df["CTDOXY"] = oxy_fitting.oxy_ml_to_umolkg(
        btl_df["CTDOXY1"], btl_df["sigma_btl"]
    )
    time_df["CTDOXY"] = oxy_fitting.oxy_ml_to_umolkg(
        time_df["CTDOXY1"], time_df["sigma_btl"]
    )

    return btl_df, time_df


def osnap_oxy_main(
    btl_data_all,
    time_data_all,
    ssscc_list,
    filepath=cfg.dirs["oxygen"] + "Winkler.xlsx",
):
    """
    Objective of skipping "prepare_oxy". Meg's titrations are imported in load_all_btl_files
    Get oxygen values converted, apply some preliminary flags, and apply the fit.
    * Modified versions of ctd_oxy_converter, calibrate_oxy
    * Return the dataframes to the workspace
    """
    # oxy_data = pd.read_excel(filepath, sheet_name="Aaron")
    # oxy_data["SSSCC"] = (
    #     oxy_data["Station/Cast"].astype(int).astype(str).str.zfill(3)
    # )  #   Drop floats and then to string

    # btl_data_all, time_data_all = ctd_oxy_converter(btl_data_all, time_data_all)

    btl_data_all["OXYGEN"] = btl_data_all["OxygenValue"]  #   Meg's reference titrations
    btl_data_all["OXYGEN_FLAG_W"] = flagging.nan_values(
        btl_data_all[cfg.column["refO"]]
    )

    #   Due to troubleshooting difficulties, calibrate_oxy was modified rather than getting a new function in here.
    btl_data_fit, time_data_fit = oxy_fitting.calibrate_oxy(
        btl_data_all, time_data_all, ssscc_list
    )

    return btl_data_fit, time_data_fit

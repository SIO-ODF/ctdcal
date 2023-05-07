"""
Script for processing stations, reusing a set of other coefficients instead of generating new ones
Good for use should a cast:
* Not have refT, refC, refO data
* Requests be made to reuse another station for fitting
* Have bottle closures or discrete data that should otherwise be disregarded

Requires an edited sssc to assign coefficients

Written 2023 - DMB cruise
"""

import gsw
import pandas as pd
import numpy as np
import logging
from pathlib import Path
from datetime import datetime

from ctdcal import (
    convert,
    fit_ctd,
    get_ctdcal_config,
    io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
    flagging,
)

log = logging.getLogger(__name__)
cfg = get_ctdcal_config()


def odf_process_reuse():
    #   Import SSSCC list of stations, including which to pull from which
    ssscc_main_file = Path("data/ssscc_reuse.csv")
    ssscc_list_matches = pd.read_csv(ssscc_main_file, sep=",", dtype=str)

    print(ssscc_list_matches)

    #   Process HEX as normal
    convert.hex_to_ctd(ssscc_list_matches["raw"])

    convert.make_time_files(ssscc_list_matches["raw"])
    #   Check if this works when no bottles have been fired
    convert.make_btl_mean(ssscc_list_matches["raw"])

    #   If bottles *were* fired, run process_bottle.process_reft(ssscc_list_matches["raw"])
    process_bottle.process_reft(ssscc_list_matches["raw"])

    #   Skip temperature, salt, oxygen files
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list_matches["raw"])
    btl_data_all = process_bottle.load_all_btl_files(ssscc_list_matches["raw"])

    #   Copies for figures
    time_data_old = time_data_all.copy()
    btl_data_old = btl_data_all.copy()

    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)
    process_ctd.make_depth_log(time_data_all)

    #   Now start loading coefficients
    ssscc_coef_t1 = Path("data/logs/fit_coef_t1.csv")
    ssscc_coef_t1 = pd.read_csv(ssscc_coef_t1, dtype="str").squeeze(axis=1)
    ssscc_coef_t2 = Path("data/logs/fit_coef_t2.csv")
    ssscc_coef_t2 = pd.read_csv(ssscc_coef_t2, dtype="str").squeeze(axis=1)
    ssscc_coef_c1 = Path("data/logs/fit_coef_c1.csv")
    ssscc_coef_c1 = pd.read_csv(ssscc_coef_c1, dtype="str").squeeze(axis=1)
    ssscc_coef_c2 = Path("data/logs/fit_coef_c2.csv")
    ssscc_coef_c2 = pd.read_csv(ssscc_coef_c2, dtype="str").squeeze(axis=1)
    ssscc_coef_43 = Path("data/logs/sbe43_coefs.csv")
    ssscc_coef_43 = pd.read_csv(ssscc_coef_43, dtype="str", index_col=0).squeeze(axis=1)
    ssscc_coef_ri = Path("data/logs/rinko_coefs.csv")
    ssscc_coef_ri = pd.read_csv(ssscc_coef_ri, dtype="str", index_col=0).squeeze(axis=1)

    #   Use these to grab individual stations
    # ssscc_coef_ri.loc[ssscc_list_matches["from_ssscc"]]
    # ssscc_coef_ri.loc[ssscc_list_matches["from_ssscc"].iloc[0]]   #   0=i

    btl_df = btl_data_all.copy()  #   Just in case
    time_df = time_data_all.copy()
    #   For-loop through each ssscc
    #   For-loop of each t, c, sensor
    for ssscc in ssscc_list_matches["from_ssscc"]:
        t1 = ssscc_coef_t1.loc[ssscc_coef_t1["SSSCC"] == ssscc].astype(float)
        t2 = ssscc_coef_t2.loc[ssscc_coef_t2["SSSCC"] == ssscc].astype(float)
        c1 = ssscc_coef_c1.loc[ssscc_coef_c1["SSSCC"] == ssscc].astype(float)
        c2 = ssscc_coef_c2.loc[ssscc_coef_c2["SSSCC"] == ssscc].astype(float)

        if ssscc == "05001":
            print("Fish sucked up on 05101 - using secondary line")
            #   Why does this give NaN? Overwrite line 1 I guess...
            # time_df.loc[
            #     time_df["SSSCC"] == "05101",
            #     ["CTDTMP1", "CTDTMP2", "CTDCOND1", "CTDCOND2"],
            # ] = time_df.loc[
            #     time_df["SSSCC"] == "05101",
            #     ["CTDTMP2", "CTDTMP1", "CTDCOND2", "CTDCOND1"],
            # ]
            # btl_df.loc[
            #     btl_df["SSSCC"] == "05101",
            #     ["CTDTMP1", "CTDTMP2", "CTDCOND1", "CTDCOND2"],
            # ] = btl_df.loc[
            #     btl_df["SSSCC"] == "05101",
            #     ["CTDTMP2", "CTDTMP1", "CTDCOND2", "CTDCOND1"],
            # ]
            time_df.loc[time_df.SSSCC == "05101", "CTDTMP1"] = time_df.loc[
                time_df.SSSCC == "05101", "CTDTMP2"
            ]
            time_df.loc[time_df.SSSCC == "05101", "CTDCOND1"] = time_df.loc[
                time_df.SSSCC == "05101", "CTDCOND2"
            ]
            btl_df.loc[btl_df.SSSCC == "05101", "CTDTMP1"] = btl_df.loc[
                btl_df.SSSCC == "05101", "CTDTMP2"
            ]
            btl_df.loc[btl_df.SSSCC == "05101", "CTDCOND1"] = btl_df.loc[
                btl_df.SSSCC == "05101", "CTDCOND2"
            ]

        #   Temperature, order is ct1, ct2, then c0, cp1, cp2
        for tN in ["t1", "t2"]:
            #   Extract coefficients to tuple for polyfit
            if tN == "t1":
                #   Array, because polyfit doesn't take pandas series
                ts = np.array((t1.c0, t1.ct1, t1.ct2))
                ps = np.array((t1.cp1, t1.cp2))
            elif tN == "t2":
                ts = np.array((t2.c0, t2.ct1, t2.ct2))
                ps = np.array((t2.cp1, t2.cp2))
            btl_df[cfg.column[tN]] = fit_ctd.apply_polyfit(
                btl_df[cfg.column[tN]],
                ts,
                (btl_df[cfg.column["p"]], ps),
            )
            time_df[cfg.column[tN]] = fit_ctd.apply_polyfit(
                time_df[cfg.column[tN]],
                ts,
                (time_df[cfg.column["p"]], ps),
            )
        #   Conductivity, order is c0, cc1, cc2, then ct1, ct2, then cp1, cp2
        for cN, tN in zip(["c1", "c2"], ["t1", "t2"]):
            if cN == "c1":
                cs = np.array((c1.c0, c1.cc1, c1.cc2))
                ts = np.array((c1.ct1, c1.ct2))
                ps = np.array((c1.cp1, c1.cp2))
                # cs = (c1.c0, c1.cc1)
                # ts = c1.ct1
                # ps = c1.cp1
            elif cN == "c2":
                cs = np.array((c2.c0, c2.cc1, c2.cc2))
                ts = np.array((c2.ct1, c2.ct2))
                ps = np.array((c2.cp1, c2.cp2))
            btl_df[cfg.column[cN]] = fit_ctd.apply_polyfit(
                btl_df[cfg.column[cN]],
                cs,
                (btl_df[cfg.column["p"]], ps),
                (btl_df[cfg.column[tN]], ts),
            )
            time_df[cfg.column[cN]] = fit_ctd.apply_polyfit(
                time_df[cfg.column[cN]],
                cs,
                (time_df[cfg.column["p"]], ps),
                (time_df[cfg.column[tN]], ts),
            )
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

    #   Now use the fit data to get oxygen prepared
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
    btl_df["OXYGEN_FLAG_W"] = flagging.nan_values(np.array(btl_df["CTDOXYVOLTS"]))
    time_df["dv_dt"] = oxy_fitting.calculate_dV_dt(
        time_df["CTDOXYVOLTS"], time_df["scan_datetime"]
    )
    btl_df["dv_dt"] = time_df[
        "dv_dt"
    ].mean()  #   Grab a non-zero value from the downcasts
    #   TODO: extract from time files where a bottle is fired?
    #   Now loop through the oxygen routines for each station

    for ssscc in ssscc_list_matches["from_ssscc"]:
        c43 = sum(ssscc_coef_43.loc[[ssscc]].values.astype(float).tolist(), [])
        # sbe43_dict = {"00201": (5.0667e-1, -4.9031e-1, 1.4596e0, 3.4871e-4, 3.9292e-2)}
        btl_df["CTDOXY"] = oxy_fitting._PMEL_oxy_eq(
            c43,
            (
                btl_df[cfg.column["oxyvolts"]],
                btl_df[cfg.column["p"]],
                btl_df["CTDTMP1"],
                btl_df["dv_dt"],
                btl_df["OS"],
            ),
        )
        time_df["CTDOXY"] = oxy_fitting._PMEL_oxy_eq(
            c43,
            (
                time_df["CTDOXYVOLTS"],
                time_df[cfg.column["p"]],
                time_df[cfg.column["t1"]],
                time_df["dv_dt"],
                time_df["OS"],
            ),
        )
        #   Now for the RINKO
        ri = sum(ssscc_coef_ri.loc[[ssscc]].values.astype(float).tolist(), [])
        btl_df["CTDRINKO"] = rinko._Uchida_DO_eq(
            ri,
            (
                btl_df[cfg.column["rinko_oxy"]],
                btl_df[cfg.column["p"]],
                btl_df[cfg.column["t1"]],
                btl_df[cfg.column["sal"]],
                btl_df["OS"],
            ),
        )
        time_df["CTDRINKO"] = rinko._Uchida_DO_eq(
            ri,
            (
                time_df[cfg.column["rinko_oxy"]],
                time_df[cfg.column["p"]],
                time_df[cfg.column["t1"]],
                time_df[cfg.column["sal"]],
                time_df["OS"],
            ),
        )
    #   Assign flags and prep for writeout
    btl_df["STNNBR"] = [int(x[0:3]) for x in btl_df["SSSCC"]]
    btl_df["SALNTY"] = -999
    btl_df["SALNTY_FLAG_W"] = 9
    btl_df["CTDSAL_FLAG_W"] = 2
    btl_df["OXYGEN"] = -999
    btl_df["BTLCOND"] = -999
    btl_df["REFTMP"] = -999
    btl_df["REFTMP_FLAG_W"] = 9
    #   Flag temperature sensors to one another. If they don't agree, flag it
    btl_df["CTDTMP1_FLAG_W"] = flagging.by_residual(
        np.array(btl_df["CTDTMP1"]),
        np.array(btl_df["CTDTMP2"]),
        np.array(btl_df["CTDPRS"]),
    )
    btl_df["CTDTMP2_FLAG_W"] = flagging.by_residual(
        np.array(btl_df["CTDTMP1"]),
        np.array(btl_df["CTDTMP1"]),
        np.array(btl_df["CTDPRS"]),
    )
    #   Do conductivity flagging to other sensor. If they don't agree, flag it
    btl_df["CTDCOND1_FLAG_W"] = flagging.by_residual(
        np.array(btl_df["CTDCOND1"]),
        np.array(btl_df["CTDCOND2"]),
        np.array(btl_df["CTDPRS"]),
    )
    btl_df["CTDCOND2_FLAG_W"] = flagging.by_residual(
        np.array(btl_df["CTDCOND2"]),
        np.array(btl_df["CTDCOND1"]),
        np.array(btl_df["CTDPRS"]),
    )
    #   Flag oxygen relative to RINKO (which generally performs much better)
    btl_df["CTDOXY_FLAG_W"] = flagging.by_percent_diff(
        np.array(btl_df["CTDOXY"]), np.array(btl_df["CTDRINKO"])
    )
    btl_df["CTDRINKO_FLAG_W"] = 2

    #   Now for the continuous data
    time_df["CTDRINKO_FLAG_W"] = 2
    time_df["CTDTMP_FLAG_W"] = flagging.by_residual(
        time_df["CTDTMP1"], time_df["CTDTMP2"], time_df["CTDPRS"]
    )
    time_df["CTDSAL_FLAG_W"] = flagging.by_residual(
        time_df["CTDCOND1"], time_df["CTDCOND2"], time_df["CTDPRS"]
    )
    time_df["CTDOXY_FLAG_W"] = flagging.by_percent_diff(
        time_df["CTDOXY"], time_df["CTDRINKO"]
    )

    #   Now we have everything, so try to export the data
    export_report_data_reuse(btl_df)

    #   Same old ct1 writeout
    process_ctd.export_ct1(time_df, ssscc_list_matches["raw"])

    export_hy1_reuse(btl_df)


def export_report_data_reuse(df):
    """
    Copy of the other function in process_bottle,
    reworked to not overwrite other cruise report data
    """

    df["CTDPRS"] = df["CTDPRS"].round(1)
    cruise_report_cols = [
        "STNNBR",
        "CTDPRS",
        "CTDTMP1",
        "CTDTMP1_FLAG_W",
        "CTDTMP2",
        "CTDTMP2_FLAG_W",
        "REFTMP",
        "CTDCOND1",
        "CTDCOND1_FLAG_W",
        "CTDCOND2",
        "CTDCOND2_FLAG_W",
        "BTLCOND",
        "CTDSAL",
        "CTDSAL_FLAG_W",
        "SALNTY",
        "CTDOXY",
        "CTDOXY_FLAG_W",
        "CTDRINKO",
        "CTDRINKO_FLAG_W",
        "OXYGEN",
    ]

    df[cruise_report_cols].to_csv("data/report_data_reuse.csv", index=False)

    return


def export_hy1_reuse(df, out_dir=cfg.dirs["pressure"], org="ODF"):
    """
    Cloned from export_hy1 in process_bottle,
    reworked to not overwrite other files
    """
    log.info("Exporting bottle file")
    btl_data = df.copy()
    now = datetime.now()
    file_datetime = now.strftime("%Y%m%d")

    btl_columns = {
        "EXPOCODE": "",
        "SECT_ID": "",
        "STNNBR": "",
        "CASTNO": "",
        "SAMPNO": "",
        "BTLNBR": "",
        "BTLNBR_FLAG_W": "",
        "DATE": "",
        "TIME": "",
        "LATITUDE": "",
        "LONGITUDE": "",
        "DEPTH": "METERS",
        "CTDPRS": "DBAR",
        "CTDTMP": "ITS-90",
        "CTDSAL": "PSS-78",
        "CTDSAL_FLAG_W": "",
        "SALNTY": "PSS-78",
        "SALNTY_FLAG_W": "",
        "CTDOXY": "UMOL/KG",
        "CTDOXY_FLAG_W": "",
        "REFTMP": "ITS-90",
        "REFTMP_FLAG_W": "",
    }

    # rename outputs as defined in user_settings.yaml
    for param, attrs in cfg.ctd_outputs.items():
        if param not in btl_data.columns:
            btl_data.rename(columns={attrs["sensor"]: param}, inplace=True)

    btl_data["EXPOCODE"] = cfg.expocode
    btl_data["SECT_ID"] = cfg.section_id
    btl_data["STNNBR"] = [int(x[0:3]) for x in btl_data["SSSCC"]]
    btl_data["CASTNO"] = [int(x[3:]) for x in btl_data["SSSCC"]]
    btl_data["SAMPNO"] = btl_data["btl_fire_num"].astype(int)
    btl_data = process_bottle.add_btlnbr_cols(btl_data, btl_num_col="btl_fire_num")

    # sort by decreasing sample number (increasing pressure) and reindex
    btl_data = btl_data.sort_values(
        by=["STNNBR", "SAMPNO"], ascending=[True, False], ignore_index=True
    )

    # switch oxygen primary sensor to rinko
    print("Using RINKO as CTDOXY in hy1")
    btl_data["CTDOXY"] = btl_data.loc[:, "CTDRINKO"]
    btl_data["CTDOXY_FLAG_W"] = btl_data.loc[:, "CTDRINKO_FLAG_W"]

    # add depth
    depth_df = pd.read_csv(
        cfg.dirs["logs"] + "depth_log.csv", dtype={"SSSCC": str}, na_values=-999
    ).dropna()
    manual_depth_df = pd.read_csv(
        cfg.dirs["logs"] + "manual_depth_log.csv", dtype={"SSSCC": str}
    )
    full_depth_df = pd.concat([depth_df, manual_depth_df])
    full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)
    btl_data["DEPTH"] = -999
    for index, row in full_depth_df.iterrows():
        btl_data.loc[btl_data["SSSCC"] == row["SSSCC"], "DEPTH"] = int(row["DEPTH"])

    btl_data["REFTMP_FLAG_W"] = flagging.nan_values(
        btl_data["REFTMP_FLAG_W"], old_flags=btl_data["REFTMP_FLAG_W"]
    )
    btl_data = btl_data.where(~btl_data.isnull(), -999)

    # check columns
    try:
        btl_data[btl_columns.keys()]
    except KeyError as err:
        log.info("Column names not configured properly... attempting to correct")
        bad_cols = err.args[0].split("'")[1::2]  # every other str is a column name
        for col in bad_cols:
            if col.endswith("FLAG_W"):
                log.warning(col + " missing, flagging with 9s")
                btl_data[col] = 9
            else:
                log.warning(col + " missing, filling with -999s")
                btl_data[col] = -999

    btl_data = btl_data[btl_columns.keys()]
    time_stamp = file_datetime + org

    with open(out_dir + cfg.expocode + "_reuse_hy1.csv", mode="w+") as f:
        f.write("BOTTLE, %s\n" % (time_stamp))
        f.write(",".join(btl_columns.keys()) + "\n")
        f.write(",".join(btl_columns.values()) + "\n")
        btl_data.to_csv(f, header=False, index=False)
        f.write("\n" + "END_DATA")

    main_df = merge_hy1()

    with open(out_dir + cfg.expocode + "_combo_hy1.csv", mode="w+") as f:
        f.write("BOTTLE, %s\n" % (time_stamp))
        f.write(",".join(btl_columns.keys()) + "\n")
        f.write(",".join(btl_columns.values()) + "\n")
        main_df.to_csv(f, header=False, index=False)
        f.write("\n" + "END_DATA")

    return


def merge_hy1(
    main_file=Path("data/pressure/33RR20230409_hy1.csv"),
    reuse_file=Path("data/pressure/33RR20230409_reuse_hy1.csv"),
):
    """Merges two hy1 files, returning the combined dataframe"""

    main_df = io.load_exchange_btl(main_file)
    reuse_df = io.load_exchange_btl(reuse_file)

    main_df = main_df.append(reuse_df)
    main_df = main_df.sort_values(
        by=["STNNBR", "CASTNO", "SAMPNO"],
        ascending=[True, True, False],
        ignore_index=True,
    )

    return main_df


if __name__ == "__main__":
    odf_process_reuse()

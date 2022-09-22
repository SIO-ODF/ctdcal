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
import numpy as np
import gsw

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
    convert.make_time_files(
        ssscc_list, group, ssscc_list
    )  #   0914 want *all* upcasts fit

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
    # time_data_all["CTDTMP_FLAG_W"] = flagging.by_residual(
    #     time_data_all.CTDTMP1, time_data_all.CTDTMP2, time_data_all.CTDPRS
    # )
    # btl_data_all["CTDTMP_FLAG_W"] = flagging.by_residual(
    #     btl_data_all.CTDTMP1, btl_data_all.CTDTMP2, btl_data_all.CTDPRS
    # )  #   These columns should go in the output file.
    time_data_all["CTDTMP_FLAG_W"] = 1  #   ANo reference. By WOCE def.
    btl_data_all["CTDTMP_FLAG_W"] = 1

    #   Calibrate conductivity
    btl_data_all, time_data_all = fit_ctd.calibrate_cond(btl_data_all, time_data_all)
    skip_list = [
        "168",
        "169",
        "189",
        "191",
        "192",
        "194",
        "195",
        "200",
    ]  #   A list of stations where salts were not taken whatsoever. No ref, so use 1.
    btl_data_all.loc[btl_data_all.SSSCC.isin(skip_list), "CTDSAL_FLAG_W"] = 1
    time_data_all.loc[time_data_all.SSSCC.isin(skip_list), "CTDSAL_FLAG_W"] = 1

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
            "CTDTMP2",
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

    print("All 24 Hz downcast fit data exported.")

    #   UPCASTS
    #   From Jupyter with the purpose of:
    #   1) Extracting select SSSCCs and running upcasts through processing
    #   2) Providing upcast fits without overwriting figures or logs
    #   3) Visualizing data step by step for troubleshooting
    #   4) Merging final data products (MicroCATs want 24 Hz upcast + downcast, which neither
    #   a normal ctdcal run nor the upcast-only script could provide)
    #   5) Demonstrating to curious science party members
    #   Brought into whoi_process_all.py as the routine was refined and tested after a few days.
    #   testing.ipynb is bloated, so easier to document here.
    print("\n*** Starting up on the upcasts... ***\n")

    #   Concat all the microcat data into a single list from the upcast files
    #   -> Upcast is spit out when the time files are made from the converted.pkl
    df_list = []
    for ssscc in ssscc_list:
        time_file = (
            cfg.dirs["time"] + ssscc + "_time_upcast.pkl"
        )  #   Only real change made
        time_data = pd.read_pickle(time_file)
        time_data["SSSCC"] = str(ssscc)
        time_data["dv_dt"] = oxy_fitting.calculate_dV_dt(
            time_data["CTDOXYVOLTS"], time_data["scan_datetime"]
        )
        df_list.append(time_data)
        # print("** Finished TIME data station: " + ssscc + " **")
    time_data_mc = pd.concat(df_list, axis=0, sort=False)
    time_data_mc["master_index"] = range(len(time_data_mc))

    #   The bottle file is more or less identical but we'll need the btl data for fitting regardless...
    df_data_all = pd.DataFrame()
    filepath = cfg.dirs["oxygen"] + "Winkler.xlsx"
    oxy_data_all = pd.read_excel(filepath, sheet_name="Aaron")
    #   Reformat SSSCC
    oxy_data_all["SSSCC_oxy"] = (
        oxy_data_all["Station/Cast"].astype(int).astype(str).str.zfill(3)
    )
    for ssscc in ssscc_list:
        btl_file = cfg.dirs["bottle"] + ssscc + "_btl_mean.pkl"
        btl_data = process_bottle._load_btl_data(btl_file, None)
        ### load REFC data
        refc_file = cfg.dirs["salt"] + ssscc + "_salts.csv"
        try:
            refc_data = process_bottle._load_salt_data(refc_file, index_name="SAMPNO")
        except FileNotFoundError:
            log.warning(
                "Missing (or misnamed) REFC Data Station: "
                + ssscc
                + "...filling with NaNs"
            )
            refc_data = pd.DataFrame(
                index=btl_data.index,
                columns=["CRavg", "BathTEMP", "BTLCOND"],
                dtype=float,
            )
            refc_data["SAMPNO_SALT"] = btl_data["btl_fire_num"].astype(int)
        oxy_data = (
            oxy_data_all.loc[oxy_data_all.SSSCC_oxy == ssscc]
            .sort_values("Niskin")
            .reset_index(drop=True)
        )
        btl_data = pd.merge(
            btl_data,
            refc_data,
            left_on="btl_fire_num",
            right_on="SAMPNO_SALT",
            how="outer",
        )
        btl_data = pd.merge(
            btl_data,
            oxy_data,
            left_on="btl_fire_num",
            right_on="Niskin",
            how="outer",
        )
        btl_data = process_bottle._add_btl_bottom_data(btl_data, ssscc)
        # Merge cast into df_data_all
        try:
            df_data_all = pd.concat([df_data_all, btl_data], sort=False)
        except AssertionError:
            raise AssertionError(
                "Columns of " + ssscc + " do not match those of previous columns"
            )
    # Drop duplicated columns generated by concatenation
    btl_data_mc = df_data_all.loc[:, ~df_data_all.columns.duplicated()]
    btl_data_mc["master_index"] = range(len(btl_data_mc))

    #   Temp flagging and ctdcond prep (reset to 1 later as to not ignore any data)
    time_data_mc["CTDTMP_FLAG_W"] = flagging.by_residual(
        time_data_mc.CTDTMP1, time_data_mc.CTDTMP2, time_data_mc.CTDPRS
    )
    btl_data_mc["CTDTMP_FLAG_W"] = flagging.by_residual(
        btl_data_mc.CTDTMP1, btl_data_mc.CTDTMP2, btl_data_mc.CTDPRS
    )
    btl_df = btl_data_mc.copy()
    time_df = time_data_mc.copy()
    btl_df[cfg.column["refC"]] = convert.CR_to_cond(
        btl_df["CRavg"],
        btl_df["BathTEMP"],
        btl_df[cfg.column["t1"]],
        btl_df[cfg.column["p"]],
    )
    #   Skip most of the flagging code
    btl_df["SALNTY_FLAG_W"] = flagging.nan_values(btl_df["SALNTY"])
    #   Preparing to fit
    fit_yaml = fit_ctd.load_fit_yaml()

    #   CONDUCTIVITY FITTING
    for cN, tN in zip(["c1", "c2"], ["t1", "t2"]):
        C_flag, C_fit_coefs = pd.DataFrame(), pd.DataFrame()
        ssscc_sublist = ssscc_list  #   Acting as one loop through
        # 0) grab ssscc chunk to fit
        btl_rows = btl_df["SSSCC"].isin(ssscc_sublist).values
        good_rows = btl_rows & (btl_df["SALNTY_FLAG_W"] == 2)
        time_rows = time_df["SSSCC"].isin(ssscc_sublist).values

        # 1) plot pre-fit residual
        f_stem = "ssscc_c2"  #   Choosing c2 and not c1 because c2 is a first order fit of P, T, and C rather than just P
        ctd_plots._intermediate_residual_plot(
            btl_df.loc[btl_rows, cfg.column["refC"]]
            - btl_df.loc[btl_rows, cfg.column[cN]],
            btl_df.loc[btl_rows, cfg.column["p"]],
            btl_df.loc[btl_rows, "SSSCC"],
            xlabel=f"{cN.upper()} Residual (mS/cm)",
            f_out=f"{cfg.fig_dirs[cN]}upcasts/residual_{f_stem}_prefit_upcasts.pdf",
        )
        # 2) prepare data for fitting
        # NOTE: df_bad will be overwritten during post-fit data flagging
        # but is left here for future debugging (if necessary)
        df_good, df_bad = fit_ctd._prepare_fit_data(
            btl_df[good_rows],
            cfg.column[cN],
            cfg.column["refC"],
            zRange=fit_yaml[cN][f_stem]["zRange"],
        )
        ctd_plots._intermediate_residual_plot(
            df_good["Diff"],
            df_good[cfg.column["p"]],
            df_good["SSSCC"],
            xlabel=f"{cN.upper()} Residual (mS/cm)",
            f_out=f"{cfg.fig_dirs[cN]}upcasts/residual_{f_stem}_fit_data_upcasts.pdf",
        )
        # 3) calculate fit coefs
        # TODO: truncate coefs (10 digits? look at historical data)
        P_order = fit_yaml[cN][f_stem]["P_order"]
        T_order = fit_yaml[cN][f_stem]["T_order"]
        C_order = fit_yaml[cN][f_stem]["C_order"]
        coef_dict = fit_ctd.multivariate_fit(
            df_good["Diff"],
            (df_good[cfg.column["p"]], P_order),
            (df_good[cfg.column[tN]], T_order),
            (df_good[cfg.column[cN]], C_order),
            coef_names=["cp", "ct", "cc"],
        )

        # 4) apply fit
        P_coefs = tuple(coef_dict[f"cp{n}"] for n in np.arange(1, P_order + 1))
        T_coefs = tuple(coef_dict[f"ct{n}"] for n in np.arange(1, T_order + 1))
        C_coefs = tuple(coef_dict[f"cc{n}"] for n in np.arange(1, C_order + 1))
        btl_df.loc[btl_rows, cfg.column[cN]] = fit_ctd.apply_polyfit(
            btl_df.loc[btl_rows, cfg.column[cN]],
            (coef_dict["c0"],) + C_coefs,
            (btl_df.loc[btl_rows, cfg.column["p"]], P_coefs),
            (btl_df.loc[btl_rows, cfg.column[tN]], T_coefs),
        )
        time_df.loc[time_rows, cfg.column[cN]] = fit_ctd.apply_polyfit(
            time_df.loc[time_rows, cfg.column[cN]],
            (coef_dict["c0"],) + C_coefs,
            (time_df.loc[time_rows, cfg.column["p"]], P_coefs),
            (time_df.loc[time_rows, cfg.column[tN]], T_coefs),
        )
        # 4.5) flag CTDCOND and make residual plots
        df_ques, df_bad = fit_ctd._flag_btl_data(
            btl_df[btl_rows],
            param=cfg.column[cN],
            ref=cfg.column["refC"],
            f_out=f"{cfg.fig_dirs[cN]}upcasts/residual_{f_stem}_upcast.pdf",
        )

        # 5) handle quality flags
        C_flag = pd.concat([C_flag, df_bad, df_ques])

        # 6) handle fit params
        coef_df = pd.DataFrame()
        coef_df["SSSCC"] = ssscc_sublist
        coef_names = ["cp2", "cp1", "ct2", "ct1", "cc2", "cc1", "c0"]
        coef_df[coef_names] = 0.0
        for k, v in coef_dict.items():
            coef_df[k] = v

        C_fit_coefs = pd.concat([C_fit_coefs, coef_df])

        # one more fig with all cuts
        ctd_plots._intermediate_residual_plot(
            btl_df[cfg.column["refC"]] - btl_df[cfg.column[cN]],
            btl_df[cfg.column["p"]],
            btl_df["SSSCC"],
            xlabel=f"{cN.upper()} Residual (mS/cm)",
            show_thresh=True,
            f_out=f"{cfg.fig_dirs[cN]}upcasts/residual_all_postfit_upcast.pdf",
        )

        C_fit_coefs.to_csv(
            cfg.dirs["logs"] + f"fit_coef_{cN}_upcasts.csv", index=False
        )  #   Write coeffs
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
    # flag salinity data
    time_df[cfg.column["sal"] + "_FLAG_W"] = 2
    btl_df[cfg.column["sal"] + "_FLAG_W"] = flagging.by_residual(
        btl_df[cfg.column["sal"]],
        btl_df["SALNTY"],
        btl_df[cfg.column["p"]],
    )
    bad_rows = btl_df["SALNTY_FLAG_W"].isin([3, 4])
    btl_df.loc[bad_rows, cfg.column["sal"] + "_FLAG_W"] = 2  # bad salts not used for QC
    btl_df[cfg.column["sal"] + "_FLAG_W"] = flagging.nan_values(
        btl_df[cfg.column["sal"]], old_flags=btl_df[cfg.column["sal"] + "_FLAG_W"]
    )

    #   Adjust flags for select stations where no salt samples were taken
    skip_list = [
        "189",
        "191",
        "192",
        "194",
        "195",
        "200",
    ]
    print(f"Marking {skip_list} as flag 1.")
    btl_df.loc[btl_df.SSSCC.isin(skip_list), [cfg.column["sal"] + "_FLAG_W"]] = 1
    time_df.loc[time_df.SSSCC.isin(skip_list), [cfg.column["sal"] + "_FLAG_W"]] = 1

    #   OXY FITTING
    #   We know that the fitting routine does a poor job on the upcast without using the downcast coeffs.
    btl_df["OXYGEN"] = btl_df[
        "OxygenValue"
    ]  #   Meg's reference titrations have dif name
    btl_df["OXYGEN_FLAG_W"] = flagging.nan_values(btl_df[cfg.column["refO"]])
    btl_df, time_df = osnap_oxy.ctd_oxy_converter(btl_df, time_df)  #   Makes CTDOXY
    f_out = f"{cfg.fig_dirs['ox']}upcasts/sbe43_residual_all_prefit_upcast.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df["OXYGEN"] - btl_df["CTDOXY"],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 20),
    )

    all_sbe43_merged = pd.DataFrame()
    all_sbe43_fit = pd.DataFrame()
    btl_df["dv_dt"] = np.nan
    #   Fresh copies of the data (if using breakpoints)
    time_df_prefit = time_df.copy()
    btl_df_prefit = btl_df.copy()

    #   Boston University wants to use global coeffs, so assign the ox0 coeffs
    filename = cfg.dirs["logs"] + "sbe43_coefs.csv"
    sbe43_coefs = pd.read_csv(filename, index_col="Unnamed: 0")
    #   Ignore adjustments to the SBE43 dict in density matching (redundant, as btl extracts
    #   come from this, but need SA, OS assigned on indices)
    for ssscc in ssscc_list:
        time_data = time_df[time_df["SSSCC"] == ssscc].copy()
        btl_data = btl_df[btl_df["SSSCC"] == ssscc].copy()
        # can't calibrate without bottle oxygen ("OXYGEN")
        if (btl_data["OXYGEN_FLAG_W"] == 9).all():
            # sbe43_dict2[ssscc] = np.full(5, np.nan)
            print(ssscc + " skipped, all oxy data is NaN")
            #   Can't calibrate, but let's try to apply the default coeffs
        sbe43_merged = oxy_fitting.match_sigmas(
            btl_data[cfg.column["p"]],
            btl_data[cfg.column["refO"]],
            btl_data["CTDTMP1"],
            btl_data["SA"],
            time_data["OS"],
            time_data[cfg.column["p"]],
            time_data[cfg.column["t1"]],
            time_data["SA"],
            time_data[cfg.column["oxyvolts"]],
            time_data["scan_datetime"],
        )
        sbe43_merged = sbe43_merged.reindex(btl_data.index)  # add nan rows back in
        btl_df.loc[
            btl_df["SSSCC"] == ssscc, ["CTDOXYVOLTS", "dv_dt", "OS"]
        ] = sbe43_merged[["CTDOXYVOLTS", "dv_dt", "OS"]]
        sbe43_merged["SSSCC"] = ssscc
        all_sbe43_merged = pd.concat([all_sbe43_merged, sbe43_merged])
        print(ssscc + " density matching done")

    #   Apply the ox0 fit now
    time_df["CTDOXY"] = np.nan
    for ssscc in ssscc_list:
        btl_rows = (btl_df["SSSCC"] == ssscc).values
        time_rows = (time_df["SSSCC"] == ssscc).values
        btl_df.loc[btl_rows, "CTDOXY"] = oxy_fitting._PMEL_oxy_eq(
            sbe43_coefs.loc["ox0"].to_list(),
            (
                btl_df.loc[btl_rows, cfg.column["oxyvolts"]],
                btl_df.loc[btl_rows, cfg.column["p"]],
                btl_df.loc[btl_rows, cfg.column["t1"]],
                btl_df.loc[btl_rows, "dv_dt"],
                btl_df.loc[btl_rows, "OS"],
            ),
        )
        print(ssscc + " btl data fitting done")
        time_df.loc[time_rows, "CTDOXY"] = oxy_fitting._PMEL_oxy_eq(
            sbe43_coefs.loc["ox0"].to_list(),
            (
                time_df.loc[time_rows, cfg.column["oxyvolts"]],
                time_df.loc[time_rows, cfg.column["p"]],
                time_df.loc[time_rows, cfg.column["t1"]],
                time_df.loc[time_rows, "dv_dt"],
                time_df.loc[time_rows, "OS"],
            ),
        )
        print(ssscc + " time data fitting done")
        #   Check for asymptotic fitting from underconstrained curves at the surface
        check_region = time_df.loc[time_rows]
        if (np.std(check_region.CTDOXY.iloc[0:19]) > 10) & (
            (check_region.CTDOXY.iloc[0] - check_region.CTDOXY.iloc[1]) < 0
        ):
            idx = np.diff(np.diff(check_region.CTDOXY.iloc[0:29])).argmax() + 2
            fill_region = check_region.iloc[0:idx]
            coefs = [1.18855453, fill_region.iloc[-1]["CTDOXY"]]
            fn = np.poly1d(coefs)
            check_region.CTDOXY.iloc[0 : len(fill_region)] = np.flip(
                fn(fill_region.CTDPRS)
            )
            time_df.loc[time_rows, "CTDOXY"] = check_region.CTDOXY
    time_df[
        "CTDOXY_FLAG_W"
    ] = 2  # TODO: actual flagging of some kind? Oxy flags should match T flags, as oxy is highly dependent on T?
    btl_df["CTDOXY_FLAG_W"] = flagging.by_percent_diff(
        btl_df["CTDOXY"], btl_df["OXYGEN"], percent_thresh=1
    )

    f_out = f"{cfg.fig_dirs['ox']}upcasts/sbe43_residual_all_postfit_upcast.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df["OXYGEN"] - btl_df["CTDOXY"],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
        ylim=(4000, 0),
    )
    f_out = f"{cfg.fig_dirs['ox']}upcasts/sbe43_residual_all_postfit_upcast_flag2.pdf"
    flag2 = btl_df["CTDOXY_FLAG_W"] == 2
    ctd_plots._intermediate_residual_plot(
        btl_df.loc[flag2, "OXYGEN"] - btl_df.loc[flag2, "CTDOXY"],
        btl_df.loc[flag2, "CTDPRS"],
        btl_df.loc[flag2, "SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
        ylim=(4000, 0),
    )

    #   EXPORTING
    depth_df = pd.read_csv(
        cfg.dirs["logs"] + "depth_log.csv", dtype={"SSSCC": str}, na_values=-999
    ).dropna()
    manual_depth_df = pd.read_csv(
        cfg.dirs["logs"] + "manual_depth_log.csv", dtype={"SSSCC": str}
    )
    full_depth_df = pd.concat([depth_df, manual_depth_df])
    full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)
    btl_df["DEPTH"] = -999
    for index, row in full_depth_df.iterrows():
        btl_df.loc[btl_df["SSSCC"] == row["SSSCC"], "DEPTH"] = int(row["DEPTH"])
    btl_df["DateTime"] = btl_df.nmea_datetime
    time_df["DateTime"] = time_df.nmea_datetime
    btl_df = process_bottle.add_btlnbr_cols(btl_df, btl_num_col="btl_fire_num")
    btl_df = btl_df.rename(
        columns={"SALNTY": "BTL_SAL", "OXYGEN": "BTL_OXY", "CTDOXY": "CTDOXY_SIO"}
    )
    #   The upcast bottle file won't replace the downcast one, but write it out for comparison anyway
    #   Change output directory and filename
    outfile = cfg.dirs["pressure"] + "upcasts/bottle_data_upcasts"
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
        "CTDOXY_SIO": "UMOL/KG",
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
        btl_df[btl_cols.keys()].to_csv(f, header=False, index=False)
    save_btl = btl_df[btl_cols.keys()].to_xarray()
    save_btl.attrs = btl_cols
    save_btl.attrs["description"] = "OSNAP32 Bottle data"
    save_btl.attrs[
        "no_sample_stations"
    ] = skip_list  #   Provide the list of stations where no bottles were sampled
    save_btl.to_netcdf(path=outfile + ".nc")

    time_cols = {
        "SSSCC": "Station",
        "DateTime": "",
        "GPSLAT": "Dec Degrees",
        "GPSLON": "Dec Degrees",
        "CTDPRS": "DBAR",
        "CTDTMP1": "ITS-90",
        "CTDTMP2": "ITS-90",
        "CTDTMP_FLAG_W": "",
        "CTDCOND1": "mS/cm",
        "CTDCOND2": "mS/cm",
        "CTDSAL": "PSS-78",
        "CTDSAL_FLAG_W": "",
        "CTDOXY_SIO": "UMOL/KG",
        "CTDOXY_FLAG_W": "",
        "CTDOXY1": "ML/L",
        "CTDOXYVOLTS": "0-5VDC",
        "ALT": "M",
        "CTDFLUOR": "mg/m^3",
        "TURBIDITY": "0-5VDC",
        "CTDXMISS": "0-5VDC",
        "FLUOR_CDOM": "0-5VDC",
    }
    #   Clarify CTDOXY as SIO-calculated for Boston University
    time_df = time_df.rename(columns={"CTDOXY": "CTDOXY_SIO"})
    for ssscc in ssscc_list:
        print("Writing upcast/downcast/MicroCAT time files for", ssscc)
        upcast_out = (
            cfg.dirs["pressure"] + "upcasts/" + ssscc + "_profile_upcast_2db.nc"
        )
        downcast_out = (
            cfg.dirs["pressure"] + "upcasts/" + ssscc + "_profile_downcast_2db.nc"
        )
        downcast = (
            cfg.dirs["pressure"] + ssscc + "_profile.nc"
        )  #   Load the downcast for the station
        downcast = xr.open_dataset(downcast).to_dataframe()
        downcast = downcast.rename(columns={"CTDOXY": "CTDOXY_SIO"})
        time_ssscc = time_df.loc[
            time_df.SSSCC == ssscc
        ]  #   Slice the cumulative upcast time dataframe
        # time_ssscc["SSSCC"] = ssscc
        if ssscc in microcat_list:  #   Affix the downcast & export at 24 Hz
            time_out_full = (
                cfg.dirs["pressure"] + "upcasts/" + ssscc + "_profile_complete.nc"
            )
            full_cast = pd.concat([downcast, time_ssscc], ignore_index=True)
            full_cast = full_cast[time_cols.keys()].to_xarray()
            full_cast.attrs = time_cols
            full_cast.attrs[
                "description"
            ] = "OSNAP32 Full 24 Hz data of station downcast and upcast"
            full_cast.to_netcdf(path=time_out_full)
        #   Now create a 2 db version on the downcast and upcast, then resave them

        up_2db = process_ctd.pressure_sequence(time_ssscc, direction="up")
        up_2db["SSSCC"] = ssscc
        down_2db = process_ctd.pressure_sequence(downcast)
        down_2db["SSSCC"] = ssscc

        up_2db = up_2db[time_cols.keys()].to_xarray()
        down_2db = down_2db[time_cols.keys()].to_xarray()
        up_2db.attrs = time_cols
        down_2db.attrs = time_cols
        if ssscc in skip_list:
            up_2db.attrs["cond_warn"] = "No salinity samples taken"
            down_2db.attrs["cond_warn"] = "No salinity samples taken"

        up_2db.attrs["description"] = "OSNAP32 2 decibar average of station upcast"
        up_2db.to_netcdf(path=upcast_out)

        up_2db.attrs["description"] = "OSNAP32 2 decibar average of station downcast"
        down_2db.to_netcdf(path=downcast_out)

    print("All data fit and exported:")
    print("2 db upcast and downcasts")
    print("24 Hz MicroCat stations")
    print("24 Hz downcasts (on request)")

    full_depth_df = pd.concat([depth_df, manual_depth_df])
    full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)
    full_depth_df = full_depth_df.sort_values("SSSCC")
    import datetime

    cast_filename = cfg.dirs["logs"] + "cast_details.csv"
    cast_dat = pd.read_csv(cast_filename, dtype={"SSSCC": str})
    cast_dat["Time"] = cast_dat.end_time.apply(
        lambda x: datetime.datetime.fromtimestamp(x)
    )

    if len(cast_dat) == len(full_depth_df):
        cast_dat["Depth"] = full_depth_df.DEPTH.values
        cast_dat[["Lat", "Lon"]] = np.nan
        for i in range(0, len(cast_dat)):
            cast_dat["Lat"].iloc[i] = convert.deg_to_dms(
                cast_dat.latitude.iloc[i], pretty_print="latitude"
            )
            cast_dat["Lat"].iloc[i] = (
                cast_dat["Lat"].iloc[i][:-5] + cast_dat["Lat"].iloc[i][-3:]
            )  #   Remove extra zeros
            cast_dat["Lon"].iloc[i] = convert.deg_to_dms(
                -cast_dat.longitude.iloc[i], pretty_print="longitude"
            )
            cast_dat["Lon"].iloc[i] = (
                cast_dat["Lon"].iloc[i][:-5] + cast_dat["Lon"].iloc[i][-3:]
            )  #   Remove extra zeros
        #   Final formatting. Change depth to integer
        cast_dat.loc[cast_dat.Depth.isnull(), "Depth"] = -999
        cast_dat["Depth"] = cast_dat["Depth"].astype(int)
        file_out = cfg.dirs["logs"] + "cast_data_formatted.csv"
        cast_dat[["SSSCC", "Time", "Lat", "Lon", "max_pressure", "Depth"]].to_csv(
            file_out, index=False
        )
    else:
        print(
            "Cast data and depth log are not the same size. Could not write out formatted cast table."
        )


if __name__ == "__main__":
    whoi_process_all()

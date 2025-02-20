"""
A module for fitting CTD discrete/bottle and continuous data.
"""
import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from ctdcal.fitting.fit_common import get_node, NodeNotFoundError, multivariate_fit, apply_polyfit
from ctdcal import get_ctdcal_config
from ctdcal import process_ctd as process_ctd
from ctdcal.processors.functions_salt import CR_to_cond
from ctdcal.plotting.plot_fit import _intermediate_residual_plot
from ctdcal.fitting.fit_legacy import load_fit_yaml
from ctdcal.flagging.flag_common import _flag_btl_data, outliers, nan_values, by_residual

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def _prepare_fit_data(df, param, ref_param, zRange=None):
    """Remove non-finite data, trim to desired zRange, and remove extreme outliers"""

    good_data = df[np.isfinite(df[ref_param]) & np.isfinite(df[param])].copy()

    if zRange is not None:
        zMin, zMax = zRange.split(":")
        good_data = good_data[
            (good_data["CTDPRS"] > int(zMin)) & (good_data["CTDPRS"] < int(zMax))
        ]

    good_data["Diff"] = good_data[ref_param] - good_data[param]
    good_data["Flag"] = outliers(good_data["Diff"], n_sigma2=4)

    df_good = good_data[good_data["Flag"] == 2].copy()
    df_bad = good_data[good_data["Flag"] == 4].copy()

    return df_good, df_bad


def calibrate_temp(btl_df, time_df):
    """
    Least-squares fit CTD temperature data against reference data.

    Parameters
    -----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data

    Returns
    --------

    """
    log.info("Calibrating temperature")
    ssscc_subsets = sorted(Path(cfg.dirs["ssscc"]).glob("ssscc_t*.csv"))
    if not ssscc_subsets:  # if no t-segments exists, write one from full list
        log.debug(
            "No CTDTMP grouping file found... creating ssscc_t1.csv with all casts"
        )
        if not Path(cfg.dirs["ssscc"]).exists():
            Path(cfg.dirs["ssscc"]).mkdir()
        ssscc_list = process_ctd.get_ssscc_list()
        ssscc_subsets = [Path(cfg.dirs["ssscc"] + "ssscc_t1.csv")]
        pd.Series(ssscc_list).to_csv(ssscc_subsets[0], header=None, index=False)

    fit_yaml = load_fit_yaml()  # load fit polynomial order
    for tN in ["t1", "t2"]:
        T_flag, T_fit_coefs = pd.DataFrame(), pd.DataFrame()
        for f in ssscc_subsets:
            # 0) load ssscc subset to be fit together
            ssscc_sublist = (
                pd.read_csv(f, header=None, dtype="str").squeeze(axis=1).to_list()
            )
            btl_rows = btl_df["SSSCC"].isin(ssscc_sublist).values
            good_rows = btl_rows & (btl_df["REFTMP_FLAG_W"] == 2)
            time_rows = time_df["SSSCC"].isin(ssscc_sublist).values

            # 1) plot pre-fit residual
            f_stem = f.stem  # get "ssscc_t*" from path
            _intermediate_residual_plot(
                btl_df.loc[btl_rows, cfg.column["refT"]]
                - btl_df.loc[btl_rows, cfg.column[tN]],
                btl_df.loc[btl_rows, cfg.column["p"]],
                btl_df.loc[btl_rows, "SSSCC"],
                xlabel=f"{tN.upper()} Residual (T90 C)",
                f_out=f"{cfg.fig_dirs[tN]}residual_{f_stem}_prefit.pdf",
            )

            # 2) prepare data for fitting
            # NOTE: df_bad will be overwritten during post-fit data flagging but is
            # left here for future debugging (if necessary)
            df_good, df_bad = _prepare_fit_data(
                btl_df[good_rows],
                cfg.column[tN],
                cfg.column["refT"],
                zRange=fit_yaml[tN][f_stem]["zRange"],
            )
            _intermediate_residual_plot(
                df_good["Diff"],
                df_good[cfg.column["p"]],
                df_good["SSSCC"],
                xlabel=f"{tN.upper()} Residual (T90 C)",
                f_out=f"{cfg.fig_dirs[tN]}residual_{f_stem}_fit_data.pdf",
            )

            # 3) calculate fit coefs
            P_order = fit_yaml[tN][f_stem]["P_order"]
            T_order = fit_yaml[tN][f_stem]["T_order"]
            coef_dict = multivariate_fit(
                df_good["Diff"],
                (df_good[cfg.column["p"]], P_order),
                (df_good[cfg.column[tN]], T_order),
                coef_names=["cp", "ct"],
            )

            # 4) apply fit
            P_coefs = tuple(coef_dict[f"cp{n}"] for n in np.arange(1, P_order + 1))
            T_coefs = tuple(coef_dict[f"ct{n}"] for n in np.arange(1, T_order + 1))
            btl_df.loc[btl_rows, cfg.column[tN]] = apply_polyfit(
                btl_df.loc[btl_rows, cfg.column[tN]],
                (coef_dict["c0"],) + T_coefs,
                (btl_df.loc[btl_rows, cfg.column["p"]], P_coefs),
            )
            time_df.loc[time_rows, cfg.column[tN]] = apply_polyfit(
                time_df.loc[time_rows, cfg.column[tN]],
                (coef_dict["c0"],) + T_coefs,
                (time_df.loc[time_rows, cfg.column["p"]], P_coefs),
            )

            # 4.5) flag CTDTMP and make residual plots
            df_ques, df_bad = _flag_btl_data(
                btl_df[btl_rows],
                param=cfg.column[tN],
                ref=cfg.column["refT"],
                f_out=f"{cfg.fig_dirs[tN]}residual_{f_stem}.pdf",
            )

            # 5) handle quality flags
            T_flag = pd.concat([T_flag, df_bad, df_ques])

            # 6) handle fit params
            coef_df = pd.DataFrame()
            coef_df["SSSCC"] = ssscc_sublist
            coef_names = ["cp2", "cp1", "ct2", "ct1", "c0"]
            coef_df[coef_names] = 0.0
            for k, v in coef_dict.items():
                coef_df[k] = v

            T_fit_coefs = pd.concat([T_fit_coefs, coef_df])

        # one more fig with all cuts
        _intermediate_residual_plot(
            btl_df[cfg.column["refT"]] - btl_df[cfg.column[tN]],
            btl_df[cfg.column["p"]],
            btl_df["SSSCC"],
            xlabel=f"{tN.upper()} Residual (T90 C)",
            show_thresh=True,
            f_out=f"{cfg.fig_dirs[tN]}residual_all_postfit.pdf",
        )

        # export temp quality flags
        T_flag.sort_index().to_csv(f"{cfg.dirs['logs']}qual_flag_{tN}.csv", index=False)

        # export temp fit params (formated to 5 sig figs, scientific notation)
        T_fit_coefs[coef_names] = T_fit_coefs[coef_names].applymap(
            lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
        )
        T_fit_coefs.to_csv(cfg.dirs["logs"] + f"fit_coef_{tN}.csv", index=False)

    # flag temperature data
    time_df["CTDTMP_FLAG_W"] = 2

    return True


def calibrate_cond(btl_df, time_df, user_cfg, ref_node):
    """
    Least-squares fit CTD conductivity data against bottle salts.

    Parameters
    -----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data
    user_cfg : Munch object
        Dictionary of user configuration parameters
    ref_node : str
        Name of reference parameter

    Returns
    --------

    """
    log.info("Calibrating conductivity")
    # calculate BTLCOND values from autosal data
    btl_df[cfg.column["refC"]] = CR_to_cond(
        btl_df["CRavg"],
        btl_df["BathTEMP"],
        btl_df[cfg.column["t1"]],
        btl_df[cfg.column["p"]],
    )

    # merge in handcoded salt flags
    flag_file = Path(user_cfg.datadir, 'flag', user_cfg.bottleflags_man)
    salt_flags_manual = None
    if flag_file.exists():
        try:
            salt_flags_manual = get_node(flag_file, ref_node)
        except NodeNotFoundError:
            log.info("No previously flagged values for %s found in flag file." % ref_node)
    else:
        log.info("No pre-existing flag file found.")

    if salt_flags_manual is not None:
        log.info("Merging previously flagged values for %s." % ref_node)
        salt_flags_manual_df = pd.DataFrame.from_dict(salt_flags_manual)
        salt_flags_manual_df = salt_flags_manual_df.rename(
            columns={"cast_id": "SSSCC", "bottle_num": "btl_fire_num", "value": "SALNTY_FLAG_W"}
        ).drop(columns=["notes"])
        btl_df = btl_df.merge(salt_flags_manual_df, on=["SSSCC", "btl_fire_num"], how="left")
        btl_df["SALNTY_FLAG_W"] = nan_values(
            btl_df["SALNTY"], old_flags=btl_df["SALNTY_FLAG_W"]
        )
    else:
        btl_df["SALNTY_FLAG_W"] = nan_values(btl_df["SALNTY"])

    ssscc_subsets = sorted(Path(cfg.dirs["ssscc"]).glob("ssscc_c*.csv"))
    if not ssscc_subsets:  # if no c-segments exists, write one from full list
        log.debug(
            "No CTDCOND grouping file found... creating ssscc_c1.csv with all casts"
        )
        ssscc_list = process_ctd.get_ssscc_list()
        ssscc_subsets = [Path(cfg.dirs["ssscc"] + "ssscc_c1.csv")]
        pd.Series(ssscc_list).to_csv(ssscc_subsets[0], header=None, index=False)

    fit_yaml = load_fit_yaml()  # load fit polynomial order
    for cN, tN in zip(["c1", "c2"], ["t1", "t2"]):
        C_flag, C_fit_coefs = pd.DataFrame(), pd.DataFrame()
        for f in ssscc_subsets:
            # 0) grab ssscc chunk to fit
            ssscc_sublist = (
                pd.read_csv(f, header=None, dtype="str").squeeze(axis=1).to_list()
            )
            btl_rows = btl_df["SSSCC"].isin(ssscc_sublist).values
            good_rows = btl_rows & (btl_df["SALNTY_FLAG_W"] == 2)
            time_rows = time_df["SSSCC"].isin(ssscc_sublist).values

            # 1) plot pre-fit residual
            f_stem = f.stem  # get "ssscc_c*" from path
            _intermediate_residual_plot(
                btl_df.loc[btl_rows, cfg.column["refC"]]
                - btl_df.loc[btl_rows, cfg.column[cN]],
                btl_df.loc[btl_rows, cfg.column["p"]],
                btl_df.loc[btl_rows, "SSSCC"],
                xlabel=f"{cN.upper()} Residual (mS/cm)",
                f_out=f"{cfg.fig_dirs[cN]}residual_{f_stem}_prefit.pdf",
            )

            # 2) prepare data for fitting
            # NOTE: df_bad will be overwritten during post-fit data flagging
            # but is left here for future debugging (if necessary)
            df_good, df_bad = _prepare_fit_data(
                btl_df[good_rows],
                cfg.column[cN],
                cfg.column["refC"],
                zRange=fit_yaml[cN][f_stem]["zRange"],
            )
            _intermediate_residual_plot(
                df_good["Diff"],
                df_good[cfg.column["p"]],
                df_good["SSSCC"],
                xlabel=f"{cN.upper()} Residual (mS/cm)",
                f_out=f"{cfg.fig_dirs[cN]}residual_{f_stem}_fit_data.pdf",
            )

            # 3) calculate fit coefs
            P_order = fit_yaml[cN][f_stem]["P_order"]
            T_order = fit_yaml[cN][f_stem]["T_order"]
            C_order = fit_yaml[cN][f_stem]["C_order"]
            coef_dict = multivariate_fit(
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
            btl_df.loc[btl_rows, cfg.column[cN]] = apply_polyfit(
                btl_df.loc[btl_rows, cfg.column[cN]],
                (coef_dict["c0"],) + C_coefs,
                (btl_df.loc[btl_rows, cfg.column["p"]], P_coefs),
                (btl_df.loc[btl_rows, cfg.column[tN]], T_coefs),
            )
            time_df.loc[time_rows, cfg.column[cN]] = apply_polyfit(
                time_df.loc[time_rows, cfg.column[cN]],
                (coef_dict["c0"],) + C_coefs,
                (time_df.loc[time_rows, cfg.column["p"]], P_coefs),
                (time_df.loc[time_rows, cfg.column[tN]], T_coefs),
            )

            # 4.5) flag CTDCOND and make residual plots
            df_ques, df_bad = _flag_btl_data(
                btl_df[btl_rows],
                param=cfg.column[cN],
                ref=cfg.column["refC"],
                f_out=f"{cfg.fig_dirs[cN]}residual_{f_stem}.pdf",
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
        _intermediate_residual_plot(
            btl_df[cfg.column["refC"]] - btl_df[cfg.column[cN]],
            btl_df[cfg.column["p"]],
            btl_df["SSSCC"],
            xlabel=f"{cN.upper()} Residual (mS/cm)",
            show_thresh=True,
            f_out=f"{cfg.fig_dirs[cN]}residual_all_postfit.pdf",
        )

        # export cond quality flags
        C_flag.sort_index().to_csv(f"{cfg.dirs['logs']}qual_flag_{cN}.csv", index=False)

        # export cond fit params
        C_fit_coefs[coef_names] = C_fit_coefs[coef_names].applymap(
            lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
        )
        C_fit_coefs.to_csv(cfg.dirs["logs"] + f"fit_coef_{cN}.csv", index=False)

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
    btl_df[cfg.column["sal"] + "_FLAG_W"] = by_residual(
        btl_df[cfg.column["sal"]],
        btl_df["SALNTY"],
        btl_df[cfg.column["p"]],
    )
    bad_rows = btl_df["SALNTY_FLAG_W"].isin([3, 4])
    btl_df.loc[bad_rows, cfg.column["sal"] + "_FLAG_W"] = 2  # bad salts not used for QC
    btl_df[cfg.column["sal"] + "_FLAG_W"] = nan_values(
        btl_df[cfg.column["sal"]], old_flags=btl_df[cfg.column["sal"] + "_FLAG_W"]
    )

    return btl_df, time_df

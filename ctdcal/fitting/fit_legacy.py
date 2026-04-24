import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import scipy
import yaml

from ctdcal import get_ctdcal_config
from ctdcal.common import get_ssscc_list
from ctdcal.fitting.fit_common import multivariate_fit, apply_polyfit, get_node, NodeNotFoundError
from ctdcal.fitting.fit_ctd import _prepare_fit_data
from ctdcal.fitting.fit_oxy import calculate_weights, PMEL_oxy_weighted_residual
from ctdcal.flagging.flag_common import _flag_btl_data, nan_values, by_residual, by_percent_diff
from ctdcal.plotting.plot_fit import _intermediate_residual_plot
from ctdcal.processors.functions_oxy import oxy_ml_to_umolkg, calculate_dV_dt
from ctdcal.processors.proc_oxy_ctd import _get_sbe_coef, _PMEL_oxy_eq
from ctdcal.processors.proc_oxy_odf import calculate_bottle_oxygen

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def load_fit_yaml(fname=f"{cfg.dirs['logs']}fit_coefs.yaml", to_object=False):
    """Load polynomial fit order information from .yaml file."""

    if not Path(fname).exists():
        log.warning("Warning: Coefficients fit order YAML does not exist. Generating from scratch...")
        generate_yaml()

    with open(fname, "r") as f:
        ymlfile = yaml.safe_load(f)

    if to_object:
        return type("ymlfile", (object,), ymlfile)
    else:
        return ymlfile


def generate_yaml(fname="fit_coefs.yaml", outdir=cfg.dirs['logs']):
    """
    Create a default coeff. yaml file.
    """
    data = {
        't1': {
            'ssscc_t1': {
                'P_order': 1,
                'T_order': 0,
                'zRange': "1000:6000"
            }
        },
        'c1': {
            'ssscc_c1': {
                'P_order': 1,
                'T_order': 0,
                'C_order': 0,
                'zRange': "1000:6000"
            }
        },
        't2': {
            'ssscc_t1': {
                'P_order': 1,
                'T_order': 0,
                'zRange': "1000:6000"
            }
        },
        'c2': {
            'ssscc_c1': {
                'P_order': 1,
                'T_order': 0,
                'C_order': 0,
                'zRange': "1000:6000"
            }
        }
    }

    # Write the data to a YAML file
    with open(outdir + fname, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)


def write_fit_yaml():
    """For future use with automated fitting routine(s).
    i.e., iterate to find best fit parameters, save to file"""
    pass

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
    log.warning("Use of fit_legacy.calibrate_temp() is deprecated. Use fit_ctd.calibrate_temperature() instead.")
    ssscc_subsets = sorted(Path(cfg.dirs["ssscc"]).glob("ssscc_t*.csv"))
    if not ssscc_subsets:  # if no t-segments exists, write one from full list
        log.debug(
            "No CTDTMP grouping file found... creating ssscc_t1.csv with all casts"
        )
        if not Path(cfg.dirs["ssscc"]).exists():
            Path(cfg.dirs["ssscc"]).mkdir()
        ssscc_list = get_ssscc_list()
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
    log.warning("Use of fit_legacy.calibrate_cond() is deprecated. Use fit_ctd.calibrate_conductivity() instead.")
    # log.info("Calibrating conductivity")
    # # calculate BTLCOND values from autosal data
    # btl_df[cfg.column["refC"]] = CR_to_cond(
    #     btl_df["CRavg"],
    #     btl_df["BathTEMP"],
    #     btl_df[cfg.column["t1"]],
    #     btl_df[cfg.column["p"]],
    # )

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
        ssscc_list = get_ssscc_list()
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

def prepare_oxy(btl_df, time_df, ssscc_list, user_cfg, ref_node):
    """
    Calculate oxygen-related variables needed for calibration:
    sigma, oxygen solubility (OS), and bottle oxygen

    Parameters
    ----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data
    ssscc_list : list of str
        List of stations to process
    user_cfg : Munch object
        Munch dictionary of user-defined parameters
    ref_node : str
        Name of reference parameter

    Returns
    -------

    """
    log.warning("Use of fit_legacy.prepare_oxy() is deprecated. Use proc_oxy_ctd.prepare_oxy() instead.")
    # Calculate SA and CT
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
    btl_df["CTDOXY"] = oxy_ml_to_umolkg(btl_df["CTDOXY1"], btl_df["sigma_btl"])
    # Calculate bottle oxygen
    # TODO: this part should be in proc_oxy_btl module
    btl_df[cfg.column["refO"]] = calculate_bottle_oxygen(
        ssscc_list,
        btl_df["SSSCC"],
        btl_df["TITR_VOL"],
        btl_df["TITR_TEMP"],
        btl_df["FLASKNO"],
    )
    btl_df[cfg.column["refO"]] = oxy_ml_to_umolkg(
        btl_df[cfg.column["refO"]], btl_df["sigma_btl"]
    )
    btl_df["OXYGEN_FLAG_W"] = nan_values(btl_df[cfg.column["refO"]])

    # Load manual OXYGEN flags
    flag_file = Path(user_cfg.datadir, "flag", user_cfg.bottleflags_man)
    oxy_flags_manual = None
    if flag_file.exists():
        try:
            oxy_flags_manual = get_node(flag_file, ref_node)
        except NodeNotFoundError:
            log.info(
                "No previously flagged values for %s found in flag file." % ref_node
            )
    else:
        log.info("No pre-existing flag file found.")

    if oxy_flags_manual is not None:
        log.info("Merging previously flagged values for %s." % ref_node)
        oxy_flags_manual_df = pd.DataFrame.from_dict(oxy_flags_manual)
        oxy_flags_manual_df = oxy_flags_manual_df.rename(
            columns={
                "cast_id": "SSSCC",
                "bottle_num": "btl_fire_num",
                "value": "OXYGEN_FLAG_W",
            }
        )
        btl_df.set_index(["SSSCC", "btl_fire_num"], inplace=True)
        btl_df.update(oxy_flags_manual_df.set_index(["SSSCC", "btl_fire_num"]))
        btl_df.reset_index(inplace=True)

    return True


##
## Oxygen fitting from legacy fit_oxy.py
##

def sbe43_oxy_fit(merged_df, sbe_coef0=None, f_suffix=None):
    """
    Fit weighted oxygen data following match_sigmas with the option for initial coefficients.
    """

    # Plot data to be fit together
    f_out = f"{cfg.fig_dirs['ox']}sbe43_residual{f_suffix}_prefit.pdf"
    _intermediate_residual_plot(
        merged_df["REFOXY"] - merged_df["CTDOXY"],
        merged_df["CTDPRS"],
        merged_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
    )

    bad_df = pd.DataFrame()  # initialize DF for questionable values

    if sbe_coef0 is None:
        sbe_coef0 = _get_sbe_coef()  # load initial coefficient guess

    # Curve fit (weighted)
    weights = calculate_weights(merged_df["CTDPRS"])
    fit_vars = ["CTDOXYVOLTS", "CTDPRS", "CTDTMP", "dv_dt", "OS"]
    fit_data = tuple(merged_df[v] for v in fit_vars)
    res = scipy.optimize.minimize(
        PMEL_oxy_weighted_residual,
        x0=sbe_coef0,
        args=(weights, fit_data, merged_df["REFOXY"]),
        bounds=[(None, None), (None, None), (0, None), (None, None), (None, None)],
    )

    cfw_coefs = res.x
    merged_df["CTDOXY"] = _PMEL_oxy_eq(cfw_coefs, fit_data)
    merged_df["residual"] = merged_df["REFOXY"] - merged_df["CTDOXY"]
    cutoff = 2.8 * np.std(merged_df["residual"])
    thrown_values = merged_df[np.abs(merged_df["residual"]) > cutoff]
    bad_df = pd.concat([bad_df, thrown_values])
    merged_df = merged_df[np.abs(merged_df["residual"]) <= cutoff].copy()

    while not thrown_values.empty:  # runs as long as there are thrown_values

        p0 = tuple(cfw_coefs)  # initialize coefficients with previous results
        weights = calculate_weights(merged_df["CTDPRS"])
        fit_data = tuple(merged_df[v] for v in fit_vars)  # merged_df changes each loop
        res = scipy.optimize.minimize(
            PMEL_oxy_weighted_residual,
            x0=p0,
            args=(weights, fit_data, merged_df["REFOXY"]),
            bounds=[(None, None), (None, None), (0, None), (None, None), (None, None)],
        )

        cfw_coefs = res.x
        merged_df["CTDOXY"] = _PMEL_oxy_eq(cfw_coefs, fit_data)
        merged_df["residual"] = merged_df["REFOXY"] - merged_df["CTDOXY"]
        cutoff = 2.8 * np.std(merged_df["residual"])
        thrown_values = merged_df[np.abs(merged_df["residual"]) > cutoff]
        bad_df = pd.concat([bad_df, thrown_values])
        merged_df = merged_df[np.abs(merged_df["residual"]) <= cutoff].copy()

    # intermediate plots to diagnose data chunks goodness
    if f_suffix is not None:
        f_out = f"{cfg.fig_dirs['ox']}sbe43_residual{f_suffix}.pdf"
        _intermediate_residual_plot(
            merged_df["residual"],
            merged_df["CTDPRS"],
            merged_df["SSSCC"],
            xlabel="CTDOXY Residual (umol/kg)",
            f_out=f_out,
            xlim=(-10, 10),
        )

    merged_df["CTDOXY_FLAG_W"] = 2
    bad_df["CTDOXY_FLAG_W"] = 3
    df = pd.concat([merged_df, bad_df])

    return cfw_coefs, df


def match_sigmas(
    btl_prs,
    btl_oxy,
    btl_tmp,
    btl_SA,
    btl_idx,
    ctd_os,
    ctd_prs,
    ctd_tmp,
    ctd_SA,
    ctd_oxyvolts,
    ctd_time,
):
    """
    Density match time/btl oxy dataframes between up/downcasts.
    """

    # Construct Dataframe from bottle and ctd values for merging
    btl_data = pd.DataFrame(
        data={"CTDPRS": btl_prs, "REFOXY": btl_oxy, "CTDTMP": btl_tmp, "SA": btl_SA},
        index=btl_idx,
    )
    time_data = pd.DataFrame(
        data={
            "CTDPRS": ctd_prs,
            "OS": ctd_os,
            "CTDTMP": ctd_tmp,
            "SA": ctd_SA,
            "CTDOXYVOLTS": ctd_oxyvolts,
            "CTDTIME": ctd_time,
        }
    )
    time_data["dv_dt"] = calculate_dV_dt(time_data["CTDOXYVOLTS"], time_data["CTDTIME"])

    # Merge DF
    merged_df = pd.DataFrame(
        columns=["CTDPRS", "CTDOXYVOLTS", "CTDTMP", "dv_dt", "OS"], dtype=float
    )
    merged_df["REFOXY"] = btl_data["REFOXY"].copy()

    # calculate sigma referenced to multiple depths
    for idx, p_ref in enumerate([0, 1000, 2000, 3000, 4000, 5000, 6000]):

        # pandas 1.2.1 ufunc issue workaround
        btl_inputs = np.broadcast_arrays(
            btl_data["SA"], btl_data["CTDTMP"], btl_data["CTDPRS"], p_ref
        )
        time_inputs = np.broadcast_arrays(
            time_data["SA"], time_data["CTDTMP"], time_data["CTDPRS"], p_ref
        )

        btl_data[f"sigma{idx}"] = (
            gsw.pot_rho_t_exact(*btl_inputs)
            - 1000  # subtract 1000 to get potential density *anomaly*
        ) + 1e-8 * np.random.standard_normal(btl_data["SA"].size)
        time_data[f"sigma{idx}"] = (
            gsw.pot_rho_t_exact(*time_inputs)
            - 1000  # subtract 1000 to get potential density *anomaly*
        ) + 1e-8 * np.random.standard_normal(time_data["SA"].size)
        rows = (btl_data["CTDPRS"] > (p_ref - 500)) & (
            btl_data["CTDPRS"] < (p_ref + 500)
        )

        time_sigma_sorted = time_data[f"sigma{idx}"].sort_values().to_numpy()
        sigma_min = np.min(
            [np.min(btl_data.loc[rows, f"sigma{idx}"]), np.min(time_sigma_sorted)]
        )
        sigma_max = np.max(
            [np.max(btl_data.loc[rows, f"sigma{idx}"]), np.max(time_sigma_sorted)]
        )
        time_sigma_sorted = np.insert(time_sigma_sorted, 0, sigma_min - 1e-4)
        time_sigma_sorted = np.append(time_sigma_sorted, sigma_max + 1e-4)
        cols = ["CTDPRS", "CTDOXYVOLTS", "CTDTMP", "dv_dt", "OS"]
        inds = np.concatenate(([0], np.arange(0, len(time_data)), [len(time_data) - 1]))
        for col in cols:
            merged_df.loc[rows, col] = np.interp(
                btl_data.loc[rows, f"sigma{idx}"],
                time_sigma_sorted,
                time_data[col].iloc[inds],
            )

    # Apply coef and calculate CTDOXY
    sbe_coef0 = _get_sbe_coef()  # initial coefficient guess
    merged_df["CTDOXY"] = _PMEL_oxy_eq(
        sbe_coef0,
        (
            merged_df["CTDOXYVOLTS"],
            merged_df["CTDPRS"],
            merged_df["CTDTMP"],
            merged_df["dv_dt"],
            merged_df["OS"],
        ),
    )

    return merged_df


def calibrate_oxy(btl_df, time_df, ssscc_list):
    """
    Non-linear least squares fit chemical sensor oxygen against bottle oxygen.

    Parameters
    ----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data
    ssscc_list : list of str
        List of stations to process

    Returns
    -------

    """
    log.info("Calibrating oxygen (SBE43)")
    # Plot all pre fit data
    f_out = f"{cfg.fig_dirs['ox']}sbe43_residual_all_prefit.pdf"
    _intermediate_residual_plot(
        btl_df["OXYGEN"] - btl_df["CTDOXY"],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
    )
    # Prep vars, dfs, etc.
    all_sbe43_merged = pd.DataFrame()
    sbe43_dict = {}
    all_sbe43_fit = pd.DataFrame()

    btl_df.set_index("master_index", inplace=True)

    btl_df["dv_dt"] = np.nan  # initialize column
    # Density match time/btl oxy dataframes
    for ssscc in ssscc_list:
        time_data = time_df[time_df["SSSCC"] == ssscc].copy()
        btl_data = btl_df[btl_df["SSSCC"] == ssscc].copy()
        # can't calibrate without bottle oxygen ("OXYGEN")
        if (btl_data["OXYGEN_FLAG_W"] == 9).all():
            sbe43_dict[ssscc] = np.full(5, np.nan)
            log.warning(ssscc + " skipped, all oxy data is NaN")
            continue
        sbe43_merged = match_sigmas(
            btl_data[cfg.column["p"]],
            btl_data[cfg.column["refO"]],
            btl_data["CTDTMP1"],
            btl_data["SA"],
            btl_data.index,
            time_data["OS"],
            time_data[cfg.column["p"]],
            time_data[cfg.column["t1"]],
            time_data["SA"],
            time_data[cfg.column["oxyvolts"]],
            time_data["scan_datetime"],
        )
        sbe43_merged = sbe43_merged.reindex(btl_data.index)  # add nan rows back in
        btl_df.loc[btl_df["SSSCC"] == ssscc, ["CTDOXYVOLTS", "dv_dt", "OS"]] = (
            sbe43_merged[["CTDOXYVOLTS", "dv_dt", "OS"]]
        )
        sbe43_merged["SSSCC"] = ssscc
        all_sbe43_merged = pd.concat([all_sbe43_merged, sbe43_merged])
        log.info(ssscc + " density matching done")

    # Only fit using OXYGEN flagged good (2)
    all_sbe43_merged = all_sbe43_merged.loc[btl_df["OXYGEN_FLAG_W"] == 2].copy()

    # Fit ALL oxygen stations together to get initial coefficient guess
    (sbe_coef0, _) = sbe43_oxy_fit(all_sbe43_merged, f_suffix="_ox0")
    sbe43_dict["ox0"] = sbe_coef0

    # Fit each cast individually
    for ssscc in ssscc_list:
        sbe_coef, sbe_df = sbe43_oxy_fit(
            all_sbe43_merged.loc[all_sbe43_merged["SSSCC"] == ssscc].copy(),
            sbe_coef0=sbe_coef0,
            f_suffix=f"_{ssscc}",
        )
        # build coef dictionary
        if ssscc not in sbe43_dict.keys():  # don't overwrite NaN'd stations
            sbe43_dict[ssscc] = sbe_coef
        # all non-NaN oxygen data with flags
        all_sbe43_fit = pd.concat([all_sbe43_fit, sbe_df])

    # apply coefs
    time_df["CTDOXY"] = np.nan
    for ssscc in ssscc_list:
        if np.isnan(sbe43_dict[ssscc]).all():
            log.warning(
                f"{ssscc} missing oxy data, leaving nan values and flagging as 9"
            )
            time_df.loc[time_df["SSSCC"] == ssscc, "CTDOXY_FLAG_W"] = 9
            continue
        btl_rows = (btl_df["SSSCC"] == ssscc).values
        time_rows = (time_df["SSSCC"] == ssscc).values
        btl_df.loc[btl_rows, "CTDOXY"] = _PMEL_oxy_eq(
            sbe43_dict[ssscc],
            (
                btl_df.loc[btl_rows, cfg.column["oxyvolts"]],
                btl_df.loc[btl_rows, cfg.column["p"]],
                btl_df.loc[btl_rows, cfg.column["t1"]],
                btl_df.loc[btl_rows, "dv_dt"],
                btl_df.loc[btl_rows, "OS"],
            ),
        )
        log.info(ssscc + " btl data fitting done")
        time_df.loc[time_rows, "CTDOXY"] = _PMEL_oxy_eq(
            sbe43_dict[ssscc],
            (
                time_df.loc[time_rows, cfg.column["oxyvolts"]],
                time_df.loc[time_rows, cfg.column["p"]],
                time_df.loc[time_rows, cfg.column["t1"]],
                time_df.loc[time_rows, "dv_dt"],
                time_df.loc[time_rows, "OS"],
            ),
        )
        log.info(ssscc + " time data fitting done")

    # flag CTDOXY with more than 1% difference
    time_df["CTDOXY_FLAG_W"] = 2
    btl_df["CTDOXY_FLAG_W"] = by_percent_diff(
        btl_df["CTDOXY"], btl_df["OXYGEN"], percent_thresh=1
    )

    # Plot all post fit data
    f_out = f"{cfg.fig_dirs['ox']}sbe43_residual_all_postfit.pdf"
    _intermediate_residual_plot(
        btl_df["OXYGEN"] - btl_df["CTDOXY"],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
    )
    f_out = f"{cfg.fig_dirs['ox']}sbe43_residual_all_postfit_flag2.pdf"
    flag2 = btl_df["CTDOXY_FLAG_W"] == 2
    _intermediate_residual_plot(
        btl_df.loc[flag2, "OXYGEN"] - btl_df.loc[flag2, "CTDOXY"],
        btl_df.loc[flag2, "CTDPRS"],
        btl_df.loc[flag2, "SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
    )

    # export fitting coefs
    sbe43_coefs = pd.DataFrame.from_dict(
        sbe43_dict, orient="index", columns=["Soc", "Voffset", "Tau20", "Tcorr", "E"]
    ).applymap(lambda x: np.format_float_scientific(x, precision=4, exp_digits=1))
    sbe43_coefs.to_csv(cfg.dirs["logs"] + "sbe43_coefs.csv")

    return True

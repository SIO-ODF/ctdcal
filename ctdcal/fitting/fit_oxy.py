"""
Functions for fitting CTD oxy data to bottle oxy data.
"""
import logging

import gsw
import numpy as np
import pandas as pd
import scipy

from ctdcal import get_ctdcal_config
from ctdcal.flagging.flag_common import by_percent_diff
from ctdcal.plotting.plot_fit import _intermediate_residual_plot
from ctdcal.processors.functions_oxy import calculate_dV_dt
from ctdcal.processors.proc_oxy_ctd import _get_sbe_coef, _PMEL_oxy_eq

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def calculate_weights(pressure):
    """
    Calculate weights (as a function of pressure) for weighted least squares fitting.
    Deep measurements are weighted higher than shallow.

    Parameters
    ----------
    pressure : array-like
        Pressure values of oxygen measurements [dbar]

    Returns
    -------
    weights : array-like
        Weight factor for each pressure value
    """
    epsilon = 1e-5  # small offset to avoid interpolation issues

    # define piecewise weight function dependent on pressure
    p_bins = [
        0,
        100,
        100 + epsilon,
        300,
        300 + epsilon,
        500,
        500 + epsilon,
        1200,
        1200 + epsilon,
        2000,
        2000 + epsilon,
        7000,
    ]
    w_bins = [20, 20, 25, 25, 50, 50, 100, 100, 200, 200, 500, 500]
    wgt = scipy.interpolate.interp1d(p_bins, w_bins)

    weights = wgt(pressure)  # get weights from piecewise function

    return weights


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


def PMEL_oxy_weighted_residual(coefs, weights, inputs, refoxy, L_norm=2):
    """
    Do a weighted oxygen residual fit using PMEL's SBE43 equation.
    """
    # (abstracted from PMEL code oxygen_cal_ml.m)
    # unweighted L2: sum((ref - oxy)^2)  # if weighted fails
    # unweighted L4: sum((ref - oxy)^4)  # unsure of use case
    # unweighted L1: sum(abs(ref - oxy))  # very far from ideal
    # anything else? genericize with integer "norm" function input?

    residuals = np.sum(
        (weights * (refoxy - _PMEL_oxy_eq(coefs, inputs)) ** 2)
    ) / np.sum(weights**2)

    return residuals


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

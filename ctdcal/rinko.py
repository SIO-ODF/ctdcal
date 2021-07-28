from collections import namedtuple
from pathlib import Path

import logging
import numpy as np
import pandas as pd
import scipy

from . import (
    ctd_plots,
    get_ctdcal_config,
    flagging,
    process_ctd,
    oxy_fitting,
)

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

RinkoO2Cal = namedtuple("RinkoO2Cal", [*"ABCDEFGH"])
RinkoTMPCal = namedtuple("RinkoTMPCal", [*"ABCD"])


def rinko_DO(p_prime, G, H):
    """
    Calculates the dissolved oxygen percentage.
    """

    DO = G + H * p_prime

    return DO


def rinko_p_prime(N, t, A, B, C, D, E, F, G, H):
    """
    Per RinkoIII manual: 'The film sensing the water is affect by environment
    temperature and pressure at the depth where it is deployed. Based on experiments,
    an empirical algorithm as following is used to correct data dissolved oxygen.'

    Parameters
    ----------
    N : array-like
        Raw instrument output
    t : array-like
        Temperature [degC]
    A-H : float
        Calibration parameters
    """
    p_prime = A / (1 + D * (t - 25)) + B / ((N - F) * (1 + D * (t - 25)) + C + F)

    return p_prime


def correct_pressure(P, d, E):
    """
    Parameters
    ----------
    P : array-like
        Temperature-corrected DO [%]
    d : array-like
        Pressure [MPa]
    E : float
        Manufacturer calibration coefficient

    Returns
    -------
    P_d : array-like
        Temperature- and pressure-corrected DO [%]
    """

    # TODO: check d range to make sure it's MPa
    # what is the dbar ~ MPa?

    P_d = P * (1 + E * d)

    return P_d


def salinity_correction(DO_c, T, S):
    """
    Oxygen optode is not able to detect salinity, so a correction is applied to
    account for the effect of salt on oxygen concentration. See Uchida (2010) in
    GO-SHIP manual (pg. 6, eq. 9) for more info.

    Parameters
    ----------
    DO_c : array-like
        Pressure-corrected dissolved oxygen
    T : array-like
        Calibrated CTD temperature
    S : array-like
        Calibrated CTD salinity

    Returns
    -------
    DO_sc : array-like
        Pressure- and salinity-corrected dissolved oxygen
    """
    # solubility coefficients from Benson and Krause (1984),
    # as recommended by Garcia and Gordon (1992)
    B0 = -6.24523e-3
    B1 = -7.37614e-3
    B2 = -1.03410e-2
    B3 = -8.17083e-3
    C0 = -4.88682e-7

    # "scaled temperature"
    T_scaled = np.log((298.15 - T) / (273.15 + T))

    # correction equation
    DO_sc = DO_c * np.exp(
        S * (B0 + (B1 * T_scaled) + (B2 * T_scaled ** 2) + (B3 * T_scaled ** 3))
        + C0 * S ** 2
    )

    return DO_sc


def _Uchida_DO_eq(coefs, inputs):
    """
    See Uchida et. al (2008) for more info:
    https://doi.org/10.1175/2008JTECHO549.1
    and Uchida et. al (2010) - GO-SHIP manual

    Parameters
    ----------
    coefs : tuple
        (c0, c1, c2, d0, d1, d2, cp)
    inputs : tuple
        (raw voltage, pressure, temperature, salinity, oxygen solubility)
    """
    c0, c1, c2, d0, d1, d2, cp = coefs
    V_r, P, T, S, o2_sol = inputs

    K_sv = c0 + (c1 * T) + (c2 * T ** 2)  # Stern-Volmer constant (Tengberg et al. 2006)
    V0 = (1 + d0 * T)  # voltage at zero oxygen (Uchida 2010, eq. 10)
    Vc = (d1 + d2 * V_r)  # raw voltage (Uchida 2010, eq. 10)
    o2_sat = ((V0 / Vc) - 1) / K_sv  # oxygen saturation [%] (Uchida 2010, eq. 6)

    DO = o2_sat * o2_sol  # dissolved oxygen concentration
    DO_c = DO * (1 + cp * P / 1000) ** (1 / 3)  # pressure compensated DO
    DO_sc = salinity_correction(DO_c, T, S)  # salinity + pressure compensated DO

    return DO_sc


def oxy_weighted_residual(coefs, weights, inputs, refoxy, L_norm=2):
    # TODO: optionally include other residual types
    # (abstracted from PMEL code oxygen_cal_ml.m)
    # unweighted L2: sum((ref - oxy)^2)  # if weighted fails
    # unweighted L4: sum((ref - oxy)^4)  # unsure of use case
    # unweighted L1: sum(abs(ref - oxy))  # very far from ideal
    # anything else? genericize with integer "norm" function input?

    residuals = np.sum(
        (weights * (refoxy - _Uchida_DO_eq(coefs, inputs)) ** 2)
    ) / np.sum(weights ** 2)

    return residuals


def calibrate_oxy(btl_df, time_df, ssscc_list):
    """
    Non-linear least squares fit oxygen optode against bottle oxygen.

    Note: optode data that were obtained during bottle stops can be used for calibration
    instead of density matching the downcast (see Uchida 2010, pg. 7).

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

    log.info("Calibrating oxygen (RINKO)")

    # initialize coef df
    coefs_df = pd.DataFrame(columns=["c0", "c1", "c2", "d0", "d1", "d2", "cp"])

    # Only fit using OXYGEN flagged good (2)
    good_data = btl_df[btl_df["OXYGEN_FLAG_W"] == 2].copy()

    # Fit ALL oxygen stations together to get initial coefficient guess
    (rinko_coefs0, _) = rinko_oxy_fit(good_data, f_suffix="_r0")
    coefs_df.loc["r0"] = rinko_coefs0  # log for comparison

    # fit station groups, like T/C fitting (ssscc_r1, _r2, etc.)
    ssscc_subsets = sorted(Path(cfg.dirs["ssscc"]).glob("ssscc_r*.csv"))
    if not ssscc_subsets:  # if no r-segments exists, write one from full list
        log.debug(
            "No CTDRINKO grouping file found... creating ssscc_r1.csv with all casts"
        )
        if not Path(cfg.dirs["ssscc"]).exists():
            Path(cfg.dirs["ssscc"]).mkdir()
        ssscc_list = process_ctd.get_ssscc_list()
        ssscc_subsets = [Path(cfg.dirs["ssscc"] + "ssscc_r1.csv")]
        pd.Series(ssscc_list).to_csv(ssscc_subsets[0], header=None, index=False)
    for f in ssscc_subsets:
        ssscc_sublist = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()
        f_stem = f.stem
        (rinko_coefs_group, _) = rinko_oxy_fit(
            good_data.loc[good_data["SSSCC"].isin(ssscc_sublist)].copy(),
            rinko_coef0=rinko_coefs0,
            f_suffix=f"_{f_stem}",
        )
        coefs_df.loc[f_stem.split("_")[1]] = rinko_coefs_group  # log for comparison

        # deal with time dependent coefs by further fitting individual casts
        # NOTE (4/9/21): tried adding time drift term unsuccessfully
        # Uchida (2010) says fitting individual stations is the same (even preferred?)
        for ssscc in ssscc_sublist:
            (rinko_coefs_ssscc, _) = rinko_oxy_fit(
                good_data.loc[good_data["SSSCC"] == ssscc].copy(),
                rinko_coef0=rinko_coefs_group,
                f_suffix=f"_{ssscc}",
            )

            # check mean/stdev to see if new fit is better or worse
            btl_rows = btl_df["SSSCC"] == ssscc
            time_rows = time_df["SSSCC"] == ssscc
            group_resid = (
                _Uchida_DO_eq(
                    rinko_coefs_group,
                    (
                        btl_df.loc[btl_rows, cfg.column["rinko_oxy"]],
                        btl_df.loc[btl_rows, cfg.column["p"]],
                        btl_df.loc[btl_rows, cfg.column["t1"]],
                        btl_df.loc[btl_rows, cfg.column["sal"]],
                        btl_df.loc[btl_rows, "OS"],
                    ),
                )
                - btl_df.loc[btl_rows, "OXYGEN"]
            )
            ssscc_resid = (
                _Uchida_DO_eq(
                    rinko_coefs_ssscc,
                    (
                        btl_df.loc[btl_rows, cfg.column["rinko_oxy"]],
                        btl_df.loc[btl_rows, cfg.column["p"]],
                        btl_df.loc[btl_rows, cfg.column["t1"]],
                        btl_df.loc[btl_rows, cfg.column["sal"]],
                        btl_df.loc[btl_rows, "OS"],
                    ),
                )
                - btl_df.loc[btl_rows, "OXYGEN"]
            )
            worse_mean = np.abs(ssscc_resid.mean()) > np.abs(group_resid.mean())
            worse_stdev = ssscc_resid.std() > group_resid.std()
            if worse_mean and worse_stdev:
                log.info(
                    f"{ssscc} fit parameters worse than {f_stem} group – reverting back"
                )
                rinko_coefs_ssscc = rinko_coefs_group

            # apply coefficients
            btl_df.loc[btl_rows, "CTDRINKO"] = _Uchida_DO_eq(
                rinko_coefs_ssscc,
                (
                    btl_df.loc[btl_rows, cfg.column["rinko_oxy"]],
                    btl_df.loc[btl_rows, cfg.column["p"]],
                    btl_df.loc[btl_rows, cfg.column["t1"]],
                    btl_df.loc[btl_rows, cfg.column["sal"]],
                    btl_df.loc[btl_rows, "OS"],
                ),
            )
            time_df.loc[time_rows, "CTDRINKO"] = _Uchida_DO_eq(
                rinko_coefs_ssscc,
                (
                    time_df.loc[time_rows, cfg.column["rinko_oxy"]],
                    time_df.loc[time_rows, cfg.column["p"]],
                    time_df.loc[time_rows, cfg.column["t1"]],
                    time_df.loc[time_rows, cfg.column["sal"]],
                    time_df.loc[time_rows, "OS"],
                ),
            )

            # save coefficients to dataframe
            coefs_df.loc[ssscc] = rinko_coefs_ssscc

    # flag CTDRINKO with more than 1% difference
    time_df["CTDRINKO_FLAG_W"] = 2  # TODO: actual flagging of some kind?
    btl_df["CTDRINKO_FLAG_W"] = flagging.by_percent_diff(
        btl_df["CTDRINKO"], btl_df["OXYGEN"], percent_thresh=1
    )

    # Plot all post fit data
    f_out = f"{cfg.fig_dirs['rinko']}rinko_residual_all_postfit.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df["OXYGEN"] - btl_df["CTDRINKO"],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDRINKO Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
    )
    f_out = f"{cfg.fig_dirs['rinko']}rinko_residual_all_postfit_flag2.pdf"
    flag2 = btl_df["CTDRINKO_FLAG_W"] == 2
    ctd_plots._intermediate_residual_plot(
        btl_df.loc[flag2, "OXYGEN"] - btl_df.loc[flag2, "CTDRINKO"],
        btl_df.loc[flag2, "CTDPRS"],
        btl_df.loc[flag2, "SSSCC"],
        xlabel="CTDRINKO Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
    )

    # export fitting coefs
    coefs_df.applymap(
        lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
    ).to_csv(cfg.dirs["logs"] + "rinko_coefs.csv")

    return True


def rinko_oxy_fit(
    btl_df,
    rinko_coef0=(
        1.89890,
        1.71137e-2,
        1.59838e-4,
        -1.07941e-3,
        -1.23152e-1,
        3.06114e-1,
        4.50828e-2,
    ),
    f_suffix=None,
):
    """
    Iteratively fit Rinko DO data against bottle oxygen.

    Default coefficients come from an old cruise report:
    https://cchdo.ucsd.edu/data/2362/p09_49RY20100706do.txt
    (there's probably a better way – are there physical meanings?)
    """
    # TODO: this should all get turned into a reusable/semi-general function.
    # It's used 2x here and 2x in oxy_fitting.sbe43_oxy_fit()
    # something like this:
    #
    # coefs, thrown_values = iter_oxy_fit(inputs, _Uchida_DO_eq)
    # while not thrown_values.empty:
    #     coefs, thrown_values = iter_oxy_fit(inputs, _Uchida_DO_eq)

    # # Plot data to be fit together
    # f_out = f"{cfg.fig_dirs['ox']}rinko_residual{f_suffix}_prefit.pdf"
    # ctd_plots._intermediate_residual_plot(
    #     merged_df["REFOXY"] - merged_df[cfg.column["rinko_oxy"]],
    #     merged_df["CTDPRS"],
    #     merged_df["SSSCC"],
    #     xlabel="CTDRINKO Residual (umol/kg)",
    #     f_out=f_out,
    #     xlim=(-10, 10),
    # )

    bad_df = pd.DataFrame()
    weights = oxy_fitting.calculate_weights(btl_df["CTDPRS"])
    fit_data = (
        btl_df[cfg.column["rinko_oxy"]],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["t1"]],
        btl_df[cfg.column["sal"]],
        btl_df["OS"],
    )
    # bounds need to be specified for all, can't just do cp...
    coef_bounds = [
        (None, None),  # c0,
        (None, None),  # c1,
        (None, None),  # c2,
        (None, None),  # d0,
        (None, None),  # d1,
        (None, None),  # d2,
        (0, 0.2),  # cp, pressure compensation
    ]
    res = scipy.optimize.minimize(
        oxy_weighted_residual,
        x0=rinko_coef0,
        args=(weights, fit_data, btl_df[cfg.column["refO"]]),
        bounds=coef_bounds,
    )

    cfw_coefs = res.x
    btl_df["RINKO_OXY"] = _Uchida_DO_eq(cfw_coefs, fit_data)
    btl_df["residual"] = btl_df[cfg.column["refO"]] - btl_df["RINKO_OXY"]

    cutoff = 2.8 * np.std(btl_df["residual"])
    thrown_values = btl_df[np.abs(btl_df["residual"]) > cutoff]
    bad_df = pd.concat([bad_df, thrown_values])
    btl_df = btl_df[np.abs(btl_df["residual"]) <= cutoff].copy()

    while not thrown_values.empty:

        p0 = tuple(cfw_coefs)
        weights = oxy_fitting.calculate_weights(btl_df["CTDPRS"])
        fit_data = (
            btl_df[cfg.column["rinko_oxy"]],
            btl_df[cfg.column["p"]],
            btl_df[cfg.column["t1"]],
            btl_df[cfg.column["sal"]],
            btl_df["OS"],
        )
        res = scipy.optimize.minimize(
            oxy_weighted_residual,
            x0=p0,
            args=(weights, fit_data, btl_df[cfg.column["refO"]]),
            bounds=coef_bounds,
        )
        cfw_coefs = res.x
        btl_df["RINKO_OXY"] = _Uchida_DO_eq(cfw_coefs, fit_data)
        btl_df["residual"] = btl_df[cfg.column["refO"]] - btl_df["RINKO_OXY"]
        cutoff = 2.8 * np.std(btl_df["residual"])
        thrown_values = btl_df[np.abs(btl_df["residual"]) > cutoff]
        bad_df = pd.concat([bad_df, thrown_values])
        btl_df = btl_df[np.abs(btl_df["residual"]) <= cutoff].copy()

    # intermediate plots to diagnose data chunks goodness
    # TODO: implement into bokeh/flask dashboard
    if f_suffix is not None:
        f_out = f"{cfg.fig_dirs['rinko']}rinko_residual{f_suffix}.pdf"
        ctd_plots._intermediate_residual_plot(
            btl_df["residual"],
            btl_df["CTDPRS"],
            btl_df["SSSCC"],
            xlabel="CTDRINKO Residual (umol/kg)",
            f_out=f_out,
            xlim=(-10, 10),
        )

    btl_df["CTDRINKO_FLAG_W"] = 2
    bad_df["CTDRINKO_FLAG_W"] = 3
    df = pd.concat([btl_df, bad_df])

    return cfw_coefs, df

    # from . import ctd_plots
    # diff = all_rinko_merged["REFOXY"] - all_rinko_merged["RINKO_OXY"]
    # rows = all_rinko_merged["SSSCC"] == "00101"
    # plt.scatter(all_rinko_merged["REFOXY"], all_rinko_merged["CTDPRS"])
    # plt.scatter(all_rinko_merged["RINKO_OXY"], all_rinko_merged["CTDPRS"])
    # ctd_plots._intermediate_residual_plot(diff, all_rinko_merged["CTDPRS"], all_rinko_merged["SSSCC"], xlim=(-10,10), f_out="rinko_test_residual.pdf")


### Everything below is old code (not necessarily bad)



def rinko_o2_cal_parameters(**kwargs):
    """
    Calibration coefficients for the oxygen calculations (from RinkoIII manual).
    """
    A = kwargs.get("A", -4.524084e1)
    B = kwargs.get("B",  1.449377e2)
    C = kwargs.get("C", -3.051590e-1)
    D = kwargs.get("D",  1.065300e-2)
    E = kwargs.get("E",  4.000000e-3)
    F = kwargs.get("F",  6.250000e-5)
    G = kwargs.get("G",  0.000000e0)
    H = kwargs.get("H",  1.000000e0)
    return RinkoO2Cal(A,B,C,D,E,F,G,H)

def rinko_temperature_cal_parameters(**kwargs):
    A = kwargs.get("A", -5.305905e0)
    B = kwargs.get("B",  1.666857e1)
    C = kwargs.get("C", -2.142681e0)
    D = kwargs.get("D",  4.582805e-1)
    return RinkoTMPCal(A,B,C,D)

def rinko_temperature(v, tmp_cal:RinkoTMPCal):
    if type(tmp_cal) is not RinkoTMPCal:
        raise ValueError("tmp_cal must be of type RinkoTMPCal")

    A, B, C, D = tmp_cal
    return A + B*v + C*v**2 + D*v**3

def rinko_pprime_aro_cav(v, t, o2_cal:RinkoO2Cal):
    """
    Calculates Rinko P' of the equation P = G + H * P'
    where P is DO physical value IN PERCENT [%]
    """
    A, B, C, D, E, F, G, H = o2_cal

    term_1_denominator = 1 + D*(t-25) + F*(t-25)**2
    term_1 = A/term_1_denominator

    term_2_denominator = v * (1 + D*(t-25) + F*(t-25)**2) + C
    term_2 = B/term_2_denominator

    return term_1 + term_2

def rinko_saturation(pprime, o2_cal:RinkoO2Cal):
    """
        Calculates Rinko P of the equation P = G + H * P'
        where P is DO physical value IN PERCENT [%]
    """

    A, B, C, D, E, F, G, H = o2_cal

    return G + H * pprime

def rinko_correct_for_pressure(p, d, o2_cal:RinkoO2Cal):
    """Note that the pressure term, d, must be in MPa

    1 decibar = 0.01 Mpa
    """
    A, B, C, D, E, F, G, H = o2_cal

    return p*(1 + E*d)

def rinko_saturation(df, film="B", model="ARO-CAV", **kwargs):
    pass

def rinko_oxy_eq(press, temp, oxyvo, os, o2_cal:RinkoO2Cal):

    #Calculate pprime

    pprime = rinko_pprime_aro_cav(oxyvo,temp,o2_cal)

    # Calculate P (DO physical value in %)

    #p = rinko_saturation(pprime, o2_cal)

    # Correct for pressure * d is pressure in Mpa *

    d = press * 0.01

    p_corr = rinko_correct_for_pressure(pprime,d,o2_cal)

    # Divide by 100 to get percents in the form of 0.xx

    p_corr = p_corr / 100

    # Multiply by OS to get DO (os can be in either ml/l or umol/kg)

    DO = p_corr * os

    return DO

def rinko_curve_fit_eq(X, a, b, c, d, e, f, g, h):
    """
    Same as rinko_oxy_eq, but in a form that is more suitible for scipy's curve fit routine
    X contains pressure, temperature, voltage, and OS (the normal arguments for rinko_oxy_eq)
    """

    press, temp, oxyvo, os = X
    o2_cal = RinkoO2Cal(a, b, c, d, e, f, g, h)

    #Calculate pprime

    pprime = rinko_pprime_aro_cav(oxyvo,temp,o2_cal)

    # Calculate P (DO physical value in %)

    #p = rinko_saturation(pprime, o2_cal)

    # Correct for pressure * d is pressure in Mpa *

    d = press * 0.01

    p_corr = rinko_correct_for_pressure(pprime,d,o2_cal)

    # Divide by 100 to get percents in the form of 0.xx

    p_corr = p_corr / 100

    # Multiply by OS to get DO (os can be in either ml/l or umol/kg)

    DO = p_corr * os

    return DO

def rinko_oxygen_cal(o2_cal,pressure,temp,oxyvolts,os,ref_oxy,switch):

    """"

    Rinko oxygen fitting routine using the equation:
        Calculates Rinko P of the equation P = G + H * P'
        where P is DO physical value IN PERCENT [%]

    Fits 7 Coef:
    coef[0] = A
    coef[1] = B
    coef[2] = C
    coef[3] = D
    coef[4] = E
    coef[5] = F
    coef[6] = G


    """

    rinko_oxy = rinko_oxy_eq(pressure, temp, oxyvolts, os, o2_cal)

    #Weight Determination
    if switch == 1:

        weights = oxy_fitting.calculate_weights(pressure)

        resid = ((weights * (ref_oxy - rinko_oxy))**2) / (np.sum(weights)**2) #Original way (np.sum(weights)**2)

    elif switch == 2:
        #L2 Norm
        resid = (ref_oxy - rinko_oxy)**2

    elif switch == 3:
        #ODF residuals

        resid = np.sqrt(((ref_oxy - rinko_oxy)**2) / (np.std(rinko_oxy)**2))

    elif switch == 4:
        # Weighted ODF residuals

        weights = oxy_fitting.calculate_weights(pressure)
        resid = np.sqrt(weights * ((ref_oxy - rinko_oxy)**2) / (np.sum(weights)**2))#(np.std(ctd_oxy_mlL)**2))

    elif switch == 5:

        weights = oxy_fitting.calculate_weights(pressure)

        resid = ((weights * (ref_oxy - rinko_oxy))**2) / (np.sum(weights**2))


    return resid

def rinko_weighted_residual(coefs,weights,inputs,refoxy):
    a, b, c, d, e, f, g, h = coefs
    return np.sum((weights*(refoxy-rinko_curve_fit_eq(inputs, a, b, c, d, e, f, g, h))**2))/np.sum(weights**2)

def match_sigmas(btl_prs, btl_oxy, btl_sigma, ctd_sigma, ctd_os, ctd_prs, ctd_tmp, rinkovolts, btl_ssscc=None):

    # Construct Dataframe from bottle and ctd values for merging
    if 'btl_ssscc' in locals():
        btl_dict = {'CTDPRS_rinko_btl':btl_prs, 'REFOXY_rinko':btl_oxy, 'sigma_rinko_btl':btl_sigma, 'SSSCC_rinko':btl_ssscc}
    else:
        btl_dict = {'CTDPRS_rinko_btl':btl_prs, 'REFOXY_rinko':btl_oxy, 'sigma_rinko_btl':btl_sigma}
    btl_data = pd.DataFrame(btl_dict)
    time_dict = {'CTDPRS_rinko_ctd':ctd_prs, 'sigma_rinko_ctd':ctd_sigma, 'OS_rinko_ctd':ctd_os, 'CTDTMP_rinko':ctd_tmp, 'CTDRINKOVOLTS':rinkovolts}
    time_data = pd.DataFrame(time_dict)

    # Sort DataFrames by sigma0
    time_data.sort_values('sigma_rinko_ctd', inplace=True)
    btl_data.sort_values('sigma_rinko_btl', inplace=True)
    btl_data.dropna(subset=['REFOXY_rinko'], inplace=True)

    # Merge DF
    merged_df = pd.merge_asof(btl_data, time_data, left_on='sigma_rinko_btl', right_on='sigma_rinko_ctd', direction='nearest', suffixes=['_btl','_ctd'])

    # Apply coef and calculate CTDRINKO
    rinko_coef0 = rinko_o2_cal_parameters()
    merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],rinko_coef0)

    return merged_df

    
def rinko_oxygen_fit(merged_df, rinko_coef0=None, f_out=None):

    # Create DF for good and questionable values

    bad_df = pd.DataFrame()
    good_df = pd.DataFrame()

    if rinko_coef0 is None:
        # Load initial coefficient guess
        rinko_coef0 = rinko_o2_cal_parameters()

    p0 = [x for x in rinko_coef0]

    # Curve fit (weighted)
    weights = oxy_fitting.calculate_weights(merged_df['CTDPRS_rinko_ctd'])
    cfw_coefs = scipy.optimize.fmin(
        rinko_weighted_residual,
        x0=p0,
        args=(
            weights,(
                merged_df['CTDPRS_rinko_ctd'],
                merged_df['CTDTMP_rinko'],
                merged_df['CTDRINKOVOLTS'], 
                merged_df['OS_rinko_ctd']
            ),
            merged_df['REFOXY_rinko'],
        ),
        maxfun=10000,
    )
    merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coefs)
    
    merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
    stdres = np.std(merged_df['res_rinko'])
    cutoff = stdres * 2.8

    thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    bad_df = pd.concat([bad_df, thrown_values])
    merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

    while not thrown_values.empty:  # runs as long as there are thrown_values

        p0 = cfw_coefs
        weights = oxy_fitting.calculate_weights(merged_df['CTDPRS_rinko_ctd'])
        cfw_coefs = scipy.optimize.fmin(
            rinko_weighted_residual,
            x0=p0,
            args=(
                weights,(
                    merged_df['CTDPRS_rinko_ctd'],
                    merged_df['CTDTMP_rinko'],
                    merged_df['CTDRINKOVOLTS'], 
                    merged_df['OS_rinko_ctd']
                ),
                merged_df['REFOXY_rinko'],
            ),
            maxfun=10000,
        )
        merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coefs)

        merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
        stdres = np.std(merged_df['res_rinko'])
        cutoff = stdres * 2.8
        thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
        print(len(thrown_values))
        print(p0)
        bad_df = pd.concat([bad_df, thrown_values])
        merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

    # try:
    #     cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq,
    #                                               (merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']),
    #                                               merged_df['REFOXY_rinko'], coef0, sigma=weights, absolute_sigma=False, maxfev=50000)

    #     merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)
    #     merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
    #     stdres = np.std(merged_df['res_rinko'])
    #     cutoff = stdres * 2.8

    #     thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #     bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #     bad_df = pd.concat([bad_df, bad_values])
    #     merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]


    #     while not thrown_values.empty:

    #         p0 = cfw_coef
    #         # weights = 1/(np.power(merged_df['CTDPRS_rinko_ctd'],1/2))
    #         weights = 1/(oxy_fitting.calculate_weights(merged_df['CTDPRS_rinko_ctd']))
    #         cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq, (merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']), merged_df['REFOXY_rinko'], p0, sigma=weights, absolute_sigma=False, maxfev=50000)
    #         merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)
    #         merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
    #         stdres = np.std(merged_df['res_rinko'])
    #         cutoff = stdres * 2.8
    #         thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #         bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #         bad_df = pd.concat([bad_df, bad_values])
    #         merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

    # except RuntimeError:

    #         try:#Nested try/except could be better
    #             print('Weighted curve fitting failed for RINKO...using Unweighted Fitting')
    #             cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq, (merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']), merged_df['REFOXY_rinko'], coef0)
    #             merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)

    #             merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
    #             stdres = np.std(merged_df['res_rinko'])
    #             cutoff = stdres * 2.8

    #             thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #             bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #             bad_df = pd.concat([bad_df, bad_values])
    #             merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

    #             while not thrown_values.empty:

    #                 p0 = cfw_coef
    #                 cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq, (merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']), merged_df['REFOXY_rinko'], coef0)
    #                 merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_rinko_ctd'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)

    #                 merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
    #                 stdres = np.std(merged_df['res_rinko'])
    #                 cutoff = stdres * 2.8

    #                 thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #                 bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #                 bad_df = pd.concat([bad_df, bad_values])
    #                 merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

    #         except:
    #             print('Logging RINKO coef...')
    #             cfw_coef = coef0
    #             merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
    #             stdres = np.std(merged_df['res_rinko'])
    #             cutoff = stdres * 2.8
    #             thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #             bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
    #             merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]


    #coef_dict[station] = cfw_coef
    good_df = pd.concat([good_df, merged_df])
    good_df['CTDRINKO_FLAG_W'] = 2
    bad_df = pd.concat([bad_df, bad_values])
    bad_df['CTDRINKO_FLAG_W'] = 3
    df = pd.concat([good_df,bad_df])
    df.sort_values(by='CTDPRS_rinko_btl',ascending=False,inplace=True)
    oxy_df = df.copy()

    return cfw_coef, df



if __name__ == "__main__":
    print(rinko_temperature(1.565323565, rinko_temperature_cal_parameters()))
    print(rinko_correct_for_pressure(
    rinko_pprime_aro_cav(1.5421, 0.442, rinko_o2_cal_parameters()),
    2/100,
    rinko_o2_cal_parameters()
    ))

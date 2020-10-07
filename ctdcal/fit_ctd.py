#!/usr/bin/env python
from pathlib import Path

import config as cfg
import gsw
import numpy as np
import pandas as pd
import scipy
from scipy.ndimage.interpolation import shift

import ctdcal.ctd_plots as ctd_plots
import ctdcal.process_ctd as process_ctd
import ctdcal.sbe_equations_dict as sbe_eq

# TODO: clean up oxygen fitting things, they should be in ctdcal.oxy_fitting

M = 31.9988   #Molecular weight of O2
R = 831.432    #/* Gas constant, X 10 J Kmole-1 K-1 */
D = 1.42905481 #  /* density O2 g/l @ 0 C */


# Find nearest value to argument in array
# Return the index of that value
"""code_pruning: only used in old processing script odf_calibrate_ctd.py"""
def find_isopycnals(p_btl_col, t_btl_col, sal_btl_col, dov_btl_col, lat_btl_col, lon_btl_col, btl_data, p_col, t_col, sal_col, dov_col, lat_col, lon_col, time_data):
    """find_iscopycnals

    p_btl_col:   Pressure column for bottle data
    t_btl_col:   Temperature column for bottle data
    sal_btl_col: Salinity column for bottle data
    dov_btl_col: Oxygen voltage column for bottle data
    lat_btl_col: Latitude bottle column for bottle data
    lon_btl_col: Longitude bottle column for bottle data
    btl_data:    Bottle data ndarray
    p_col:       Pressure column for bottle data
    t_col:       Temperature column for bottle data
    sal_col:     Salinity column for bottle data
    dov_col:     Oxygen voltage column for bottle data
    lat_col:     Latitude column for bottle data
    lon_col:     Longitude column for bottle data
    time_data:   Time data ndarray

    """

    time_sigma = []
    CT = gsw.CT_from_t(time_data[sal_col],time_data[t_col],time_data[p_col])
    SA = gsw.SA_from_SP(time_data[sal_col],time_data[p_col],time_data[lon_col],time_data[lat_col])
    time_sigma = gsw.sigma0(SA,CT)

    # Pressure-bounded isopycnal search.
    # Based on maximum decent rate and package size.
    CT = gsw.CT_from_t(btl_data[sal_btl_col],btl_data[t_btl_col],btl_data[p_btl_col])
    SA = gsw.SA_from_SP(btl_data[sal_btl_col],btl_data[p_btl_col],btl_data[lon_btl_col],btl_data[lat_btl_col])
    btl_sigma = gsw.sigma0(SA,CT)
    for i in range(0,len(btl_data[p_btl_col])):
        #CT = gsw.CT_from_t(btl_data[sal_btl_col][i],btl_data[t_btl_col][i],btl_data[p_btl_col][i])
        #SA = gsw.SA_from_SP(btl_data[sal_btl_col][i],btl_data[p_btl_col][i],btl_data[lon_btl_col][i],btl_data[lat_btl_col][i])
        #btl_sigma = gsw.sigma0(SA,CT)
        p_indx = find_nearest(time_data[p_col], btl_data[p_btl_col][i])
        indx = find_nearest(time_sigma[p_indx-int(24*1.5):p_indx+int(24*1.5)], btl_sigma[i])

        #print('Bottle:')
        #print('Sigma: '+str(btl_sigma))
        #print('Pres: '+str(btl_data[p_btl_col][i]))
        #print('Temp: '+str(btl_data[t_btl_col][i]))
        #print('Salt: '+str(btl_data[sal_btl_col][i]))
        #print('Pressure: '+str(p_indx)+' '+str(indx+p_indx))
        #print('Sigma: '+str(time_sigma[indx+p_indx]))
        #print('Pres: '+str(time_data[p_col][indx+p_indx]))
        #print('Temp: '+str(time_data[t_col][indx+p_indx]))
        #print('Salt: '+str(time_data[sal_col][indx+p_indx]))
        '''code_pruning: code is identical except for input passed, could be turned to function (but not necessary)'''
        if indx+p_indx > len(time_sigma):
            btl_data[t_btl_col][i] = time_data[t_col][len(time_data)-1]
            btl_data[sal_btl_col][i] = time_data[sal_col][len(time_data)-1]
            btl_data[dov_btl_col][i] = time_data[dov_col][len(time_data)-1]
        else:
            btl_data[t_btl_col][i] = time_data[t_col][indx+p_indx]
            btl_data[sal_btl_col][i] = time_data[sal_col][indx+p_indx]
            btl_data[dov_btl_col][i] = time_data[dov_col][indx+p_indx]

    return btl_data


# Find nearest value to argument in array
# Return the index of that value
"""code_pruning: only used in find_isopycnals"""
def find_nearest(yarr, val):
    """find_nearest

    """
    indx = (np.abs(yarr-val)).argmin()
    return indx

"""code_pruning: only called by outdated oxy fitting routine.. marked for removal"""
# Residual calculation
def residualO2(calib, o2pl, P, K, T, S, V):
    """residual weighted difference of dissolved oxygen bottle data
       vs dissolved oxygen CTD data.

    This conversion is included for least squares fitting routine.

    calib is a dict holding Soc, Voffset, Tau20, A, B, C, E
    The following are single or list/tuple:
    calib is a list of oxy_dict coefficients to be optimized
    o2pl is dissolved oxygen winkler titrated data
    P is pressure in decibars
    K is temperature in Kelvin
    T is temperature in Celcius
    S is Practical Salinity Units
    V is Voltage from instrument
    """
    weight = []
    ctd_o2pl = sbe_eq.oxy_dict(calib, P, K, T, S, V)
    sig = np.std(ctd_o2pl)

    # Least sq residual
    for i in range(0, len(o2pl)):
        if o2pl[i] > 0:
            weight.append(scipy.sqrt((o2pl[i] - sbe_eq.oxy_dict(calib, P[i], K[i], T[i], S[i], V[i]))**2/sig**2))
    return weight


def _conductivity_polyfit(cond,temp,press,coef):

    fitted_cond = cond + (
        coef[0] * (press ** 2)
        + coef[1] * press
        + coef[2] * (temp ** 2)
        + coef[3] * temp
        + coef[4] * (cond ** 2)
        + coef[5] * cond
        + coef[6]
    )
    fitted_cond = fitted_cond.round(4)
     #fitted_sal =  gsw.SP_from_C(fitted_cond, temp, press)
    return fitted_cond#, fitted_sal


def cell_therm_mass_corr(temp, cond, sample_int=1/24, alpha=0.03, beta=1/7):
    """Correct conductivity signal for effects of cell thermal mass.

    Parameters
    ----------
    temp : array-like
        CTD temperature [degC]
    cond : array-like
        CTD conductivity [mS/cm]
    sample_int : float, optional
        CTD sample interval [seconds]
    alpha : float, optional
        Thermal anomaly amplitude
    beta : float, optional
        Thermal anomaly time constant

    Returns
    -------
    cond_corr : array-like
        Corrected CTD conductivity [mS/cm]

    Notes
    -----
    See Sea-Bird Seasoft V2 manual (Section 6, page 93) for equation information.
    Default alpha/beta values taken from Seasoft manual (page 92).
    c.f. "Thermal Inertia of Conductivity Cells: Theory" (Lueck 1990) for more info
    https://doi.org/10.1175/1520-0426(1990)007<0741:TIOCCT>2.0.CO;2
    """
    a = 2 * alpha / (sample_int * beta + 2)
    b = 1 - (2 * a / alpha)
    dc_dT = 0.1 * (1 + 0.006 * (temp - 20))
    dT = np.insert(np.diff(temp), 0, 0)  # forward diff reduces len by 1

    def calculate_CTM(b, CTM_0, a, dc_dT, dT):
        """Return CTM in units of [S/m]"""
        CTM = -1.0 * b * CTM_0 + a * (dc_dT) * dT
        return CTM

    CTM = calculate_CTM(b, 0, a, dc_dT, dT)
    CTM = calculate_CTM(b, shift(CTM, 1, order=0), a, dc_dT, dT)
    CTM = np.nan_to_num(CTM) * 10.0  # [S/m] to [mS/cm]
    cond_corr = cond + CTM

    return cond_corr


def _flag_btl_data(
    df, param=None, ref=None, thresh=[0.002, 0.005, 0.010, 0.020], f_out=None,
):
    """
    Flag CTD "btl" data against reference measurement (e.g. SBE35, bottle salts).

    Parameters
    ----------
    df : DataFrame,
        DataFrame containing btl data
    param : str
        Name of parameter to calibrate (e.g. "CTDCOND1", "CTDTMP2")
    ref : str
        Name of reference parameter to calibrate against (e.g. "BTLCOND", "T90")
    thresh : list of float, optional
        Maximum acceptable residual for each pressure range
    f_out : str, optional
        Path and filename to save residual vs. pressure plots

    Returns
    -------
    df_ques : DataFrame
        Data flagged as questionable (flag 3s)
    df_bad : DataFrame
        Data flagged as bad (flag 4s)

    """
    # TODO: thresh should probably be put in config/cast-by-cast config
    prs = cfg.column["p_btl"]

    # Remove extreme outliers and code bad
    df = df.reset_index(drop=True)
    df, df_bad = _wild_edit(
        df[param], df[ref], df[prs], df["SSSCC"], df["btl_fire_num"]
    )

    # Find values that are above the threshold and code questionable
    df.loc[(df[prs] > 2000) & (df["Diff"].abs() > thresh[0]), "Flag"] = 3
    df.loc[
        (df[prs] <= 2000) & (df[prs] > 1000) & (df["Diff"].abs() > thresh[1]), "Flag",
    ] = 3
    df.loc[
        (df[prs] <= 1000) & (df[prs] > 500) & (df["Diff"].abs() > thresh[2]), "Flag"
    ] = 3
    df.loc[(df[prs] <= 500) & (df["Diff"].abs() > thresh[3]), "Flag"] = 3
    df_good = df[df["Flag"] == 2].copy()
    df_ques = df[df["Flag"] == 3].copy()

    if f_out is not None:
        if param == cfg.column["t1_btl"]:
            xlabel = "T1 Residual (T90 C)"
        elif param == cfg.column["t2_btl"]:
            xlabel = "T2 Residual (T90 C)"
        elif param == cfg.column["c1_btl"]:
            xlabel = "C1 Residual (mS/cm)"
        elif param == cfg.column["c2_btl"]:
            xlabel = "C2 Residual (mS/cm)"
        f_out = f_out.split(".pdf")[0] + "_postfit.pdf"
        ctd_plots._intermediate_residual_plot(
            df["Diff"],
            df[prs],
            df["SSSCC"],
            show_thresh=True,
            xlabel=xlabel,
            f_out=f_out,
        )
        f_out = f_out.split(".pdf")[0] + "_flag2.pdf"
        ctd_plots._intermediate_residual_plot(
            df_good["Diff"],
            df_good[prs],
            df_good["SSSCC"],
            show_thresh=True,
            xlabel=xlabel,
            f_out=f_out,
        )

    return df_ques, df_bad


def _prepare_fit_data(df, param, ref_param, zRange=None):
    """Remove non-finite data, trim to desired zRange, and remove extreme outliers"""

    good_data = df[np.isfinite(df[ref_param]) & np.isfinite(df[param])].copy()
    if zRange is not None:
        zMin, zMax = zRange.split(":")
        good_data = good_data[
            (good_data["CTDPRS"] > int(zMin)) & (good_data["CTDPRS"] < int(zMax))
        ]
    df_good, df_bad = _wild_edit(
        good_data[param],
        good_data[ref_param],
        good_data["CTDPRS"],
        good_data["SSSCC"],
        good_data["btl_fire_num"],
        n_sigma=5,
    )

    return df_good, df_bad


def _wild_edit(param, ref_param, prs, ssscc, btl_num, n_sigma=10):
    """Calculate residual then find extreme outliers and flag as bad (code 4)"""

    diff = ref_param - param
    df = pd.concat([ssscc, btl_num, param, ref_param, prs], axis=1)
    df["Diff"] = ref_param - param
    outliers = (df["Diff"] - df["Diff"].mean()).abs() > (n_sigma * df["Diff"].std())
    df_good = df[~outliers].copy()
    df_good["Flag"] = 2
    df_bad = df[outliers].copy()
    df_bad["Flag"] = 4

    return df_good, df_bad


def _temperature_polyfit(temp,press,coef):
    
    fitted_temp = (
        temp
        + coef[0] * (press ** 2)
        + coef[1] * press
        + coef[2] * (temp ** 2)
        + coef[3] * temp
        + coef[4]
    )
    
    return fitted_temp.round(4)


def _get_T_coefs(df, T_col=None, P_order=2, T_order=2, zRange=None, f_stem=None):

    if T_col is None:
        print("Parameter invalid, specify what temp sensor is being calibrated")
        return
    P_col = cfg.column["p_btl"]

    # remove non-finite data and extreme outliers and trim to fit zRange
    df_good, df_bad = _prepare_fit_data(df, T_col, cfg.column["reft"], zRange)

    # plot data which will be used in fit (for debugging purposes)
    if f_stem is not None:
        if T_col == cfg.column["t1_btl"]:
            xlabel = "T1 Residual (T90 C)"
            f_out = f"{cfg.directory['t1_fit_figs']}residual_{f_stem}_fit_data.pdf"
        elif T_col == cfg.column["t2_btl"]:
            xlabel = "T2 Residual (T90 C)"
            f_out = f"{cfg.directory['t2_fit_figs']}residual_{f_stem}_fit_data.pdf"
        ctd_plots._intermediate_residual_plot(
            df_good["Diff"],
            df_good[cfg.column["p_btl"]],
            df_good["SSSCC"],
            xlabel=xlabel,
            f_out=f_out,
        )

    # Toggle columns based on desired polyfit order
    # (i.e. don't calculate 2nd order term if only doing 1st order fit)
    order_list = [[0, 0], [0, 1], [1, 1]]
    P_fit = order_list[P_order]
    T_fit = order_list[T_order]

    # Calculate coefficients using linear algebra.
    #
    # Columns are [P^2, P, T^2, T, 1] and give associated coefs for:
    # T_fit = c0*P^2 + c1*P + c2*T^2 + c3*T + c4
    fit_matrix = np.vstack(
        [
            P_fit[0] * df_good[cfg.column["p_btl"]] ** 2,
            P_fit[1] * df_good[cfg.column["p_btl"]],
            T_fit[0] * df_good[T_col] ** 2,
            T_fit[1] * df_good[T_col],
            np.ones(len(df_good[T_col])),
        ]
    )
    coefs = np.linalg.lstsq(fit_matrix.T, df_good["Diff"], rcond=None)[0]

    # Column of zeros can sometimes return a non-zero value (machine precision),
    # so force uncalculated fit terms to be truly zero
    coefs = coefs * np.concatenate((P_fit, T_fit, [1]))

    return coefs, df_bad


def calibrate_temp(btl_df, time_df):
    # TODO: break off parts of this to useful functions for all vars (C/T/O)
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
    print("Calibrating temperature")
    ssscc_subsets = sorted(Path(cfg.directory["ssscc"]).glob('ssscc_t*.csv'))
    if not ssscc_subsets:  # if no t-segments exists, write one from full list
        print("No CTDTMP grouping file found... creating ssscc_t1.csv with all casts")
        if not Path(cfg.directory["ssscc"]).exists():
            Path(cfg.directory["ssscc"]).mkdir()
        ssscc_list = process_ctd.get_ssscc_list()
        ssscc_subsets = [Path(cfg.directory["ssscc"] + "ssscc_t1.csv")]
        pd.Series(ssscc_list).to_csv(ssscc_subsets[0], header=None, index=False)
    qual_flag_t1 = pd.DataFrame()
    qual_flag_t2 = pd.DataFrame()
    coef_t1_all = pd.DataFrame()
    coef_t2_all = pd.DataFrame()

    for f in ssscc_subsets:
        # 0) load ssscc subset to be fit together
        ssscc_sublist = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()
        btl_rows = btl_df["SSSCC"].isin(ssscc_sublist).values
        time_rows = time_df["SSSCC"].isin(ssscc_sublist).values

        # 1) plot pre-fit residual
        f_stem = f.stem  # get "ssscc_t*" from path
        ctd_plots._intermediate_residual_plot(
            btl_df.loc[btl_rows, cfg.column["reft"]]
            - btl_df.loc[btl_rows, cfg.column["t1_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            btl_df.loc[btl_rows, "SSSCC"],
            xlabel="T1 Residual (T90 C)",
            f_out=f"{cfg.directory['t1_fit_figs']}residual_{f_stem}_prefit.pdf",
        )
        ctd_plots._intermediate_residual_plot(
            btl_df.loc[btl_rows, cfg.column["reft"]]
            - btl_df.loc[btl_rows, cfg.column["t2_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            btl_df.loc[btl_rows, "SSSCC"],
            xlabel="T2 Residual (T90 C)",
            f_out=f"{cfg.directory['t2_fit_figs']}residual_{f_stem}_prefit.pdf",
        )

        # TODO: allow for cast-by-cast T_order/P_order/zRange
        # TODO: truncate coefs (10 digits? look at historical data)
        # 2 & 3) calculate fit params
        # NOTE: df_bad_c1/2 will be overwritten during post-fit data flagging
        # but are left here for future debugging (if necessary)
        coef_t1, df_bad_t1 = _get_T_coefs(
            btl_df[btl_rows],
            T_col=cfg.column["t1_btl"],
            P_order=1,
            T_order=0,
            zRange="1000:6000",
            f_stem=f_stem,
        )
        coef_t2, df_bad_t2 = _get_T_coefs(
            btl_df[btl_rows],
            T_col=cfg.column["t2_btl"],
            P_order=1,
            T_order=0,
            zRange="1000:6000",
            f_stem=f_stem,
        )

        # 4) apply fit
        btl_df.loc[btl_rows, cfg.column["t1_btl"]] = _temperature_polyfit(
            btl_df.loc[btl_rows, cfg.column["t1_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            coef_t1,
        )
        btl_df.loc[btl_rows, cfg.column["t2_btl"]] = _temperature_polyfit(
            btl_df.loc[btl_rows, cfg.column["t2_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            coef_t2,
        )
        time_df.loc[time_rows, cfg.column["t1"]] = _temperature_polyfit(
            time_df.loc[time_rows, cfg.column["t1"]],
            time_df.loc[time_rows, cfg.column["p"]],
            coef_t1,
        )
        time_df.loc[time_rows, cfg.column["t2"]] = _temperature_polyfit(
            time_df.loc[time_rows, cfg.column["t2"]],
            time_df.loc[time_rows, cfg.column["p"]],
            coef_t2,
        )

        # 4.5) flag CTDTMP and make residual plots
        df_ques_t1, df_bad_t1 = _flag_btl_data(
            btl_df[btl_rows],
            param=cfg.column["t1_btl"],
            ref=cfg.column["reft"],
            f_out=f"{cfg.directory['t1_fit_figs']}residual_{f_stem}.pdf",
        )
        df_ques_t2, df_bad_t2 = _flag_btl_data(
            btl_df[btl_rows],
            param=cfg.column["t2_btl"],
            ref=cfg.column["reft"],
            f_out=f"{cfg.directory['t2_fit_figs']}residual_{f_stem}.pdf",
        )

        # 5) handle quality flags
        qual_flag_t1 = pd.concat([qual_flag_t1, df_bad_t1, df_ques_t1])
        qual_flag_t2 = pd.concat([qual_flag_t2, df_bad_t2, df_ques_t2])

        # 6) handle fit params
        coef_t1_df = pd.DataFrame()
        coef_t1_df["SSSCC"] = ssscc_sublist
        coef_t2_df = coef_t1_df.copy()
        coef_names = ["cp2", "cp1", "ct2", "ct1", "c0"]
        for idx, coef_name in enumerate(coef_names):
            coef_t1_df[coef_name] = coef_t1[idx]
            coef_t2_df[coef_name] = coef_t2[idx]

        coef_t1_all = pd.concat([coef_t1_all, coef_t1_df])
        coef_t2_all = pd.concat([coef_t2_all, coef_t2_df])

    # one more fig with all cuts
    ctd_plots._intermediate_residual_plot(
        btl_df[cfg.column["reft"]]
        - btl_df[cfg.column["t1_btl"]],
        btl_df[cfg.column["p_btl"]],
        btl_df["SSSCC"],
        xlabel="T1 Residual (T90 C)",
        show_thresh=True,
        f_out=f"{cfg.directory['t1_fit_figs']}residual_all_postfit.pdf",
    )
    ctd_plots._intermediate_residual_plot(
        btl_df[cfg.column["reft"]]
        - btl_df[cfg.column["t2_btl"]],
        btl_df[cfg.column["p_btl"]],
        btl_df["SSSCC"],
        xlabel="T2 Residual (T90 C)",
        show_thresh=True,
        f_out=f"{cfg.directory['t2_fit_figs']}residual_all_postfit.pdf",
    )

    # export temp quality flags
    qual_flag_t1.sort_index().to_csv(
        cfg.directory["logs"] + "qual_flag_t1.csv", index=False,
    )
    qual_flag_t2.sort_index().to_csv(
        cfg.directory["logs"] + "qual_flag_t2.csv", index=False
    )

    # export temp fit params
    coef_t1_all.to_csv(cfg.directory["logs"] + "fit_coef_t1.csv", index=False)
    coef_t2_all.to_csv(cfg.directory["logs"] + "fit_coef_t2.csv", index=False)

    return True


def _get_C_coefs(
    df, C_col=None, P_order=2, T_order=2, C_order=2, zRange=None, f_stem=None
):

    if C_col is None:
        print("Parameter invalid, specify what cond sensor is being calibrated")
        return
    elif C_col == cfg.column["c1_btl"]:
        T_col = cfg.column["t1_btl"]
    elif C_col == cfg.column["c2_btl"]:
        T_col = cfg.column["t2_btl"]
    P_col = cfg.column["p_btl"]

    df = df.reset_index().copy()
    # remove non-finite data and extreme outliers and trim to fit zRange
    df_good, df_bad = _prepare_fit_data(df, C_col, cfg.column["refc"], zRange)

    # add CTDTMP column
    df_good[T_col] = df.loc[df_good.index, T_col]

    # plot data which will be used in fit (for debugging purposes)
    if f_stem is not None:
        if C_col == cfg.column["c1_btl"]:
            xlabel = "C1 Residual (mS/cm)"
            f_out = f"{cfg.directory['c1_fit_figs']}residual_{f_stem}_fit_data.pdf"
        elif C_col == cfg.column["c2_btl"]:
            xlabel = "C2 Residual (mS/cm)"
            f_out = f"{cfg.directory['c2_fit_figs']}residual_{f_stem}_fit_data.pdf"
        ctd_plots._intermediate_residual_plot(
            df_good["Diff"],
            df_good[cfg.column["p_btl"]],
            df_good["SSSCC"],
            xlabel=xlabel,
            f_out=f_out,
        )

    # Toggle columns based on desired polyfit order
    # (i.e. don't calculate 2nd order term if only doing 1st order fit)
    order_list = [[0, 0], [0, 1], [1, 1]]
    P_fit = order_list[P_order]
    T_fit = order_list[T_order]
    C_fit = order_list[C_order]

    # Calculate coefficients using linear algebra.
    #
    # Columns are [P^2, P, T^2, T, C^2, C, 1] and give associated coefs for:
    # C_fit = c0*P^2 + c1*P + c2*T^2 + c3*T + c4*C^2 + c5*C + c6
    fit_matrix = np.vstack(
        [
            P_fit[0] * df_good[cfg.column["p_btl"]] ** 2,
            P_fit[1] * df_good[cfg.column["p_btl"]],
            T_fit[0] * df_good[T_col] ** 2,
            T_fit[1] * df_good[T_col],
            C_fit[0] * df_good[C_col] ** 2,
            C_fit[1] * df_good[C_col],
            np.ones(len(df_good[C_col])),
        ]
    )
    coefs = np.linalg.lstsq(fit_matrix.T, df_good["Diff"], rcond=None)[0]

    # Column of zeros can sometimes return a non-zero value (machine precision),
    # so force uncalculated fit terms to be truly zero
    coefs = coefs * np.concatenate((P_fit, T_fit, C_fit, [1]))

    return coefs, df_bad


def calibrate_cond(btl_df, time_df):
    # TODO: break off parts of this to useful functions for all vars (C/T/O)
    # TODO: salt subset lists aren't loading in increasing order:
    # (still functions properly but the fit_coef_c#.csv is confusing as a result)
    """
    Least-squares fit CTD conductivity data against bottle salts.

    Parameters
    -----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data

    Returns
    --------

    """
    print("Calibrating conductivity")
    # calculate BTLCOND values from autosal data
    btl_df[cfg.column["refc"]] = CR_to_cond(
        btl_df["CRavg"],
        btl_df["BathTEMP"],
        btl_df[cfg.column["t1_btl"]],
        btl_df[cfg.column["p_btl"]],
    )

    # merge in handcoded salt flags
    # TODO: make salt flagger move .csv somewhere else? or just always have it
    # somewhere else and read it from that location (e.g. in data/scratch_folder/salts)
    salt_file = "tools/salt_flags_handcoded.csv"  # abstract to config.py
    if Path(salt_file).exists():
        handcoded_salts = pd.read_csv(
            salt_file, dtype={"SSSCC": str, "salinity_flag": int}
        )
        handcoded_salts = handcoded_salts.rename(
            columns={"SAMPNO": "btl_fire_num", "salinity_flag": "SALNTY_FLAG_W"}
        ).drop(columns=["diff", "Comments"])
        btl_df = btl_df.merge(
            handcoded_salts, on=["SSSCC", "btl_fire_num"], how="left"
        )
        btl_df.loc[btl_df["SALNTY"].isnull(), "SALNTY_FLAG_W"] = 9
        btl_df["SALNTY_FLAG_W"] = btl_df["SALNTY_FLAG_W"].fillna(
            2, downcast="infer"  # fill remaining NaNs with 2s and cast to dtype int
        )
    else:
        btl_df["SALNTY_FLAG_W"] = 2

    ssscc_subsets = sorted(Path(cfg.directory["ssscc"]).glob('ssscc_c*.csv'))
    if not ssscc_subsets:  # if no c-segments exists, write one from full list
        print("No CTDCOND grouping file found... creating ssscc_c1.csv with all casts")
        ssscc_list = process_ctd.get_ssscc_list()
        ssscc_subsets = [Path(cfg.directory["ssscc"] + "ssscc_c1.csv")]
        pd.Series(ssscc_list).to_csv(ssscc_subsets[0], header=None, index=False)
    qual_flag_c1 = pd.DataFrame()
    qual_flag_c2 = pd.DataFrame()
    coef_c1_all = pd.DataFrame()
    coef_c2_all = pd.DataFrame()

    for f in ssscc_subsets:
        # 0) grab ssscc chunk to fit
        ssscc_sublist = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()
        btl_rows = (btl_df["SSSCC"].isin(ssscc_sublist).values) & (
            btl_df["SALNTY_FLAG_W"] == 2  # only use salts flagged good (e.g. 2)
        )
        time_rows = time_df["SSSCC"].isin(ssscc_sublist).values

        # 1) plot pre-fit residual
        f_stem = f.stem  # get "ssscc_c*" from path
        ctd_plots._intermediate_residual_plot(
            btl_df.loc[btl_rows, cfg.column["refc"]]
            - btl_df.loc[btl_rows, cfg.column["c1_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            btl_df.loc[btl_rows, "SSSCC"],
            xlabel="C1 Residual (mS/cm)",
            f_out=f"{cfg.directory['c1_fit_figs']}residual_{f_stem}_prefit.pdf",
        )
        ctd_plots._intermediate_residual_plot(
            btl_df.loc[btl_rows, cfg.column["refc"]]
            - btl_df.loc[btl_rows, cfg.column["c2_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            btl_df.loc[btl_rows, "SSSCC"],
            xlabel="C2 Residual (mS/cm)",
            f_out=f"{cfg.directory['c2_fit_figs']}residual_{f_stem}_prefit.pdf",
        )

        # TODO: allow for cast-by-cast T_order/P_order/zRange
        # TODO: truncate coefs (10 digits? look at historical data)
        # 2 & 3) calculate fit params
        # NOTE: df_bad_c1/2 will be overwritten during post-fit data flagging
        # but are left here for future debugging (if necessary)
        coef_c1, df_bad_c1 = _get_C_coefs(
            btl_df[btl_rows],
            C_col=cfg.column["c1_btl"],
            P_order=1,
            T_order=0,
            C_order=0,
            zRange="1000:5000",
            f_stem=f_stem,
        )
        coef_c2, df_bad_c2 = _get_C_coefs(
            btl_df[btl_rows],
            C_col=cfg.column["c2_btl"],
            P_order=1,
            T_order=0,
            C_order=0,
            zRange="1000:5000",
            f_stem=f_stem,
        )

        # 4) apply fit
        btl_df.loc[btl_rows, cfg.column["c1_btl"]] = _conductivity_polyfit(
            btl_df.loc[btl_rows, cfg.column["c1_btl"]],
            btl_df.loc[btl_rows, cfg.column["t1_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            coef_c1,
        )
        btl_df.loc[btl_rows, cfg.column["c2_btl"]] = _conductivity_polyfit(
            btl_df.loc[btl_rows, cfg.column["c2_btl"]],
            btl_df.loc[btl_rows, cfg.column["t2_btl"]],
            btl_df.loc[btl_rows, cfg.column["p_btl"]],
            coef_c2,
        )
        time_df.loc[time_rows, cfg.column["c1"]] = _conductivity_polyfit(
            time_df.loc[time_rows, cfg.column["c1"]],
            time_df.loc[time_rows, cfg.column["t1"]],
            time_df.loc[time_rows, cfg.column["p"]],
            coef_c1,
        )
        time_df.loc[time_rows, cfg.column["c2"]] = _conductivity_polyfit(
            time_df.loc[time_rows, cfg.column["c2"]],
            time_df.loc[time_rows, cfg.column["t2"]],
            time_df.loc[time_rows, cfg.column["p"]],
            coef_c2,
        )

        # 4.5) flag CTDCOND and make residual plots
        df_ques_c1, df_bad_c1 = _flag_btl_data(
            btl_df[btl_rows],
            param=cfg.column["c1_btl"],
            ref=cfg.column["refc"],
            f_out=f"{cfg.directory['c1_fit_figs']}residual_{f_stem}.pdf",
        )
        df_ques_c2, df_bad_c2 = _flag_btl_data(
            btl_df[btl_rows],
            param=cfg.column["c2_btl"],
            ref=cfg.column["refc"],
            f_out=f"{cfg.directory['c2_fit_figs']}residual_{f_stem}.pdf",
        )

        # 5) handle quality flags
        qual_flag_c1 = pd.concat([qual_flag_c1, df_bad_c1, df_ques_c1])
        qual_flag_c2 = pd.concat([qual_flag_c2, df_bad_c2, df_ques_c2])

        # 6) handle fit params
        coef_c1_df = pd.DataFrame()
        coef_c1_df["SSSCC"] = ssscc_sublist
        coef_c2_df = coef_c1_df.copy()
        coef_names = ["cp2", "cp1", "ct2", "ct1", "cc2", "cc1", "c0"]
        for idx, coef_name in enumerate(coef_names):
            coef_c1_df[coef_name] = coef_c1[idx]
            coef_c2_df[coef_name] = coef_c2[idx]

        coef_c1_all = pd.concat([coef_c1_all, coef_c1_df])
        coef_c2_all = pd.concat([coef_c2_all, coef_c2_df])

    # one more fig with all cuts
    ctd_plots._intermediate_residual_plot(
        btl_df[cfg.column["refc"]]
        - btl_df[cfg.column["c1_btl"]],
        btl_df[cfg.column["p_btl"]],
        btl_df["SSSCC"],
        xlabel="C1 Residual (mS/cm)",
        show_thresh=True,
        f_out=f"{cfg.directory['c1_fit_figs']}residual_all_postfit.pdf",
    )
    ctd_plots._intermediate_residual_plot(
        btl_df[cfg.column["refc"]]
        - btl_df[cfg.column["c2_btl"]],
        btl_df[cfg.column["p_btl"]],
        btl_df["SSSCC"],
        xlabel="C2 Residual (mS/cm)",
        show_thresh=True,
        f_out=f"{cfg.directory['c2_fit_figs']}residual_all_postfit.pdf",
    )

    # export cond quality flags
    qual_flag_c1.sort_index().to_csv(
        cfg.directory["logs"] + "qual_flag_c1.csv", index=False,
    )
    qual_flag_c2.sort_index().to_csv(
        cfg.directory["logs"] + "qual_flag_c2.csv", index=False,
    )

    # export cond fit params
    coef_c1_all.to_csv(cfg.directory["logs"] + "fit_coef_c1.csv", index=False)
    coef_c2_all.to_csv(cfg.directory["logs"] + "fit_coef_c2.csv", index=False)

    # recalculate salinity with calibrated C/T
    time_df[cfg.column["sal"]] = gsw.SP_from_C(
        time_df[cfg.column["c1"]],
        time_df[cfg.column["t1"]],
        time_df[cfg.column["p"]],
    )
    btl_df[cfg.column["sal"]] = gsw.SP_from_C(
        btl_df[cfg.column["c1_btl"]],
        btl_df[cfg.column["t1_btl"]],
        btl_df[cfg.column["p_btl"]],
    )

    return btl_df, time_df
        

#def load_qual(path):
#    comment_dict = {}
#    with open(path) as f:
#        reader = csv.reader(f)
#        # ignore first line
#        next(reader)
#        for line in reader:
#            sta = int(line[0])
#            cast = int(line[1])
#            bottle = int(line[2])
#            param = line[3]
#            flag = int(line[4])
#            comment = line[5]
#
#            comment_dict[(sta, cast, bottle, param)] = [flag, comment]
#
#    return comment_dict
#
#
#
#salts = requests.get("http://go-ship.rrevelle.sio.ucsd.edu/api/salt").json()
#def o2_calc(path, o2_payload, thio_ns):

    
def CR_to_cond(cr,bath_t,ref_t,btl_p):

    """
    Convert AutoSal double conductivity ratio (CR) to conductivity using
    GSW conversion routines.

    Parameters
    ----------
    cr : array-like
        Double conductivity ratio from AutoSal, unitless
    bath_t : array-like
        AutoSal water bath temperature, degrees C
    ref_t : array-like
        CTD temperature at bottle stop, degrees C
    btl_p : array-like
        CTD pressure at bottle stop, dbar

    Returns
    -------
    cond : array-like
        Converted reference conductivity, mS/cm

    """
    
    salinity = gsw.SP_salinometer((cr / 2.0),bath_t)
    cond = gsw.C_from_SP(salinity,ref_t,btl_p)  
    
    # ignore RunTimeWarning from (np.nan <= 1)
    with np.errstate(invalid="ignore"):
        cond[cond <= 1] = np.nan
    
    return cond


'''code_pruning: looks like not used. marked for removal
def write_calib_coef(ssscc,coef,param):
    """ Write coef to csv
    
    
    """
    df = pd.DataFrame()
    df['SSSCC'] = ssscc  
    
    if param == 'T':

        df['coef_0'] = coef[0]
        df['coef_1'] = coef[1]
        df['coef_2'] = coef[2]
        df['coef_3'] = coef[3]
        df['coef_4'] = coef[4]
        df['coef_5'] = coef[5]
        df['coef_6'] = coef[6]
        
    if param == 'C':
       
        df['coef_0'] = coef[0]
        df['coef_1'] = coef[1]
        df['coef_2'] = coef[2]
        df['coef_3'] = coef[3]
        df['coef_4'] = coef[4]
        

    return df

'''
  
##
#            key = (station, cast, bottle, "o2")
#            if key in qual:
#                flag, comment = qual[key]
#                row_dict["o2_qual"] = flag
#                row_dict["comment"] = comment
#
#            o2_payload.append(row_dict)
#
#
#            # stoichiometric relation between mols of thio and mols of o2
#            #print(station, cast, bottle, o2ml, o2kg)
#            print("{:>2}, {:8.1f}".format(bottle, o2kg))
#
#        thio_ns.append([station, thio_n])
#
#o2_payload = []
#thio_ns = []
#for root, dirs, files in os.walk("/Volumes/public/O2Backup/O2/"):
#    for file in files:
#        if file.startswith("9") or file.startswith("8") or file.startswith("."):
#            continue
#        path = os.path.join(root, file)
#        o2_calc(path, o2_payload, thio_ns)
#    break
#
#with open("o2kg.json", "w") as f:
#    json.dump(o2_payload,f, indent=2)
#with open("thios.json", "w") as f:
#    json.dump(thio_ns, f)

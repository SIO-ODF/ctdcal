#!/usr/bin/env python
from pathlib import Path

import gsw
import logging
import numpy as np
import pandas as pd
import scipy
from scipy.ndimage.interpolation import shift
import yaml

from . import convert as convert
from . import ctd_plots as ctd_plots
from . import flagging as flagging
from . import get_ctdcal_config
from . import process_ctd as process_ctd

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def load_fit_yaml(fname=f"{cfg.dirs['logs']}fit_coefs.yaml", to_object=False):
    """Load polynomail fit order information from .yaml file."""

    with open(fname, "r") as f:
        ymlfile = yaml.safe_load(f)

    if to_object:
        return type("ymlfile", (object,), ymlfile)
    else:
        return ymlfile


def write_fit_yaml():
    """For future use with automated fitting routine(s).
    i.e., iterate to find best fit parameters, save to file"""
    pass


def _conductivity_polyfit(cond, temp, press, coef):

    fitted_cond = cond + (
        coef[0] * (press ** 2)
        + coef[1] * press
        + coef[2] * (temp ** 2)
        + coef[3] * temp
        + coef[4] * (cond ** 2)
        + coef[5] * cond
        + coef[6]
    )

    return fitted_cond


def cell_therm_mass_corr(temp, cond, sample_int=1 / 24, alpha=0.03, beta=1 / 7):
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
    df,
    param=None,
    ref=None,
    thresh=[0.002, 0.005, 0.010, 0.020],
    f_out=None,
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
    prs = cfg.column["p"]

    # Remove extreme outliers and code bad
    df = df.reset_index(drop=True)
    df["Diff"] = df[ref] - df[param]
    df["Flag"] = flagging.outliers(df["Diff"])

    # Find values that are above the threshold and code questionable
    df["Flag"] = flagging.by_residual(df[param], df[ref], df[prs], old_flags=df["Flag"])
    df_good = df[df["Flag"] == 2].copy()
    df_ques = df[df["Flag"] == 3].copy()
    df_bad = df[df["Flag"] == 4].copy()

    if f_out is not None:
        if param == cfg.column["t1"]:
            xlabel = "T1 Residual (T90 C)"
        elif param == cfg.column["t2"]:
            xlabel = "T2 Residual (T90 C)"
        elif param == cfg.column["c1"]:
            xlabel = "C1 Residual (mS/cm)"
        elif param == cfg.column["c2"]:
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

    good_data["Diff"] = good_data[ref_param] - good_data[param]
    good_data["Flag"] = flagging.outliers(good_data["Diff"], n_sigma2=4)

    df_good = good_data[good_data["Flag"] == 2].copy()
    df_bad = good_data[good_data["Flag"] == 4].copy()

    return df_good, df_bad


def _temperature_polyfit(temp, press, coef):

    fitted_temp = (
        temp
        + coef[0] * (press ** 2)
        + coef[1] * press
        + coef[2] * (temp ** 2)
        + coef[3] * temp
        + coef[4]
    )

    return fitted_temp


def _get_T_coefs(df, T_col=None, P_order=2, T_order=2, zRange=None, f_stem=None):

    if T_col is None:
        print("Parameter invalid, specify what temp sensor is being calibrated")
        return

    # remove non-finite data and extreme outliers and trim to fit zRange
    df_good, df_bad = _prepare_fit_data(df, T_col, cfg.column["refT"], zRange)

    # plot data which will be used in fit (for debugging purposes)
    if f_stem is not None:
        if T_col == cfg.column["t1"]:
            xlabel = "T1 Residual (T90 C)"
            f_out = f"{cfg.fig_dirs['t1']}residual_{f_stem}_fit_data.pdf"
        elif T_col == cfg.column["t2"]:
            xlabel = "T2 Residual (T90 C)"
            f_out = f"{cfg.fig_dirs['t2']}residual_{f_stem}_fit_data.pdf"
        ctd_plots._intermediate_residual_plot(
            df_good["Diff"],
            df_good[cfg.column["p"]],
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
            P_fit[0] * df_good[cfg.column["p"]] ** 2,
            P_fit[1] * df_good[cfg.column["p"]],
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


def multivariate_fit(y, *args, coef_names=None, const_name="c0"):
    """
    Least-squares fit data using multiple dependent variables. Dependent variables must
    be provided in tuple pairs of (data, order) as positional arguments.

    If coef_names are defined, coefficients will be returned as a dict. Otherwise,
    coefficients are return as an array in the order of the dependent variables, sorted
    by decreasing powers.

    Parameters
    ----------
    y : array-like
        Indepedent variable to be fit
    args : tuple
        Pairs of dependent variable data and fit order (i.e., (data, order))
    coef_names : list-like, optional
        Base names for coefficients (i.e., "a" for 2nd order yields ["a2", "a1"])
    const_name : str, optional
        Name for constant offset term

    Returns
    -------
    coefs : array-like
        Least-squares fit coefficients in decreasing powers

    Examples
    --------
    Behavior when coef_names is None:

    >>> z = [1, 4, 9]
    >>> x = [1, 3, 5]
    >>> y = [1, 2, 3]
    >>> multivariate_fit(z, (x, 2), (y, 1))
    array([0.25, 0.375, 0.25, 0.125])  # [c1, c2, c3, c4]

    where z = (c1 * x ** 2) + (c2 * x) + (c3 * y) + c4

    Behavior when coef_names is given:

    >>> z = [1, 4, 9]
    >>> x = [1, 3, 5]
    >>> y = [1, 2, 3]
    >>> multivariate_fit(z, (x, 2), (y, 1), coef_names=["a", "b"])
    {"a2": 0.25, "a1": 0.375, "b1": 0.25, "c0": 0.125}

    where z = (a2 * x ** 2) + (a1 * x) + (b1 * y) + c0
    """
    to_dict = True
    if coef_names is None:
        to_dict = False
        coef_names = [""] * len(args)  # needed to prevent zip() error
    elif len(args) != len(coef_names):
        raise ValueError("coef_names must be the same length as dependent variables")

    # iteratively build fit matrix
    rows, names = [], []
    for arg, coef_root in zip(args, coef_names):
        if type(arg) is not tuple:
            raise TypeError(f"Positional args must be tuples, not {type(arg)}")

        series, order = arg
        for n in np.arange(1, order + 1)[::-1]:
            rows.append(series ** n)  # n is np.int64 so series will cast to np.ndarray
            names.append(f"{coef_root}{n}")

    # add constant offset term
    rows.append(np.ones(len(y)))
    names.append(const_name)

    fit_matrix = np.vstack(rows)
    coefs = np.linalg.lstsq(fit_matrix.T, y, rcond=None)[0]

    return dict(zip(names, coefs)) if to_dict else coefs


def apply_polyfit(y, y_coefs, *args):
    """
    Apply a polynomial correction to series of data. Coefficients should be provided in
    increasing order (i.e., a0, a1, a2 for y_fit = y + a2 * y ** 2 + a1 * y + a0)

    For the independent variables (y), coefficients start from the zero-th order (i.e.,
    constant offset). For dependent variables (args), coefficients start from the first
    order (i.e., linear term).

    Parameters
    ----------
    y : array-like
        Independent variable data to be corrected
    y_coefs : tuple of float
        Independent variable fit coefficients (i.e., (data, coef0, ..., coefN))
    args : tuple of (array-like, (float, float, ...))
        Dependent variable data and fit coefficients (i.e., (data, coef1, ..., coefN))

    Returns
    -------
    fitted_y : array-like
        Independent variable data with polynomial fit correction applied
    """
    fitted_y = np.copy(y)
    for n, coef in enumerate(y_coefs):
        fitted_y += coef * y ** n

    for arg in args:
        if type(arg) is not tuple:
            raise TypeError(f"Positional args must be tuples, not {type(arg)}")

        series, coefs = arg
        for n, coef in enumerate(coefs):
            fitted_y += coef * series ** (n + 1)

    return fitted_y


def calibrate_temp(btl_df, time_df):
    # TODO: make this an ODF script in ctdcal/scripts?
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
            ssscc_sublist = pd.read_csv(f, header=None, dtype="str").squeeze().to_list()
            btl_rows = btl_df["SSSCC"].isin(ssscc_sublist).values
            good_rows = btl_rows & (btl_df["REFTMP_FLAG_W"] == 2)
            time_rows = time_df["SSSCC"].isin(ssscc_sublist).values

            # 1) plot pre-fit residual
            f_stem = f.stem  # get "ssscc_t*" from path
            ctd_plots._intermediate_residual_plot(
                btl_df.loc[btl_rows, cfg.column["refT"]]
                - btl_df.loc[btl_rows, cfg.column[tN]],
                btl_df.loc[btl_rows, cfg.column["p"]],
                btl_df.loc[btl_rows, "SSSCC"],
                xlabel=f"{tN} Residual (T90 C)",
                f_out=f"{cfg.fig_dirs[tN]}residual_{f_stem}_prefit.pdf",
            )

            # 2) prepare data for fitting
            # NOTE: df_bad will be overwritten during post-fit data flagging but are
            # left here for future debugging (if necessary)
            df_good, df_bad = _prepare_fit_data(
                btl_df[good_rows],
                cfg.column[tN],
                cfg.column["refT"],
                zRange=fit_yaml[tN][f_stem]["zRange"],
            )

            # 3) calculate fit coefs
            # TODO: truncate coefs (10 digits? look at historical data)
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
                f_out=f"{cfg.fig_dirs['t1']}residual_{f_stem}.pdf",
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
        ctd_plots._intermediate_residual_plot(
            btl_df[cfg.column["refT"]] - btl_df[cfg.column[tN]],
            btl_df[cfg.column["p"]],
            btl_df["SSSCC"],
            xlabel=f"{tN.upper()} Residual (T90 C)",
            show_thresh=True,
            f_out=f"{cfg.fig_dirs[tN]}residual_all_postfit.pdf",
        )

        # export temp quality flags
        # TODO: make these flags useful/less cluttered
        T_flag.sort_index().to_csv(f"{cfg.dirs['logs']}qual_flag_{tN}.csv", index=False)

        # export temp fit params (formated to 5 sig figs, scientific notation)
        T_fit_coefs[coef_names] = T_fit_coefs[coef_names].applymap(
            lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
        )
        T_fit_coefs.to_csv(cfg.dirs["logs"] + f"fit_coef_{tN}.csv", index=False)

    # flag temperature data
    # TODO: CTDTMP_FLAG_W historically not included in hy1 file... should it be?
    time_df["CTDTMP_FLAG_W"] = 2  # TODO: flag w/ REFT somehow? discrete vs continuous

    return True


def _get_C_coefs(
    df, C_col=None, P_order=2, T_order=2, C_order=2, zRange=None, f_stem=None
):

    if C_col is None:
        print("Parameter invalid, specify what cond sensor is being calibrated")
        return
    elif C_col == cfg.column["c1"]:
        T_col = cfg.column["t1"]
    elif C_col == cfg.column["c2"]:
        T_col = cfg.column["t2"]

    df = df.reset_index().copy()
    # remove non-finite data and extreme outliers and trim to fit zRange
    df_good, df_bad = _prepare_fit_data(df, C_col, cfg.column["refC"], zRange)

    # add CTDTMP column
    df_good[T_col] = df.loc[df_good.index, T_col]

    # plot data which will be used in fit (for debugging purposes)
    if f_stem is not None:
        if C_col == cfg.column["c1"]:
            xlabel = "C1 Residual (mS/cm)"
            f_out = f"{cfg.fig_dirs['c1']}residual_{f_stem}_fit_data.pdf"
        elif C_col == cfg.column["c2"]:
            xlabel = "C2 Residual (mS/cm)"
            f_out = f"{cfg.fig_dirs['c2']}residual_{f_stem}_fit_data.pdf"
        ctd_plots._intermediate_residual_plot(
            df_good["Diff"],
            df_good[cfg.column["p"]],
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
            P_fit[0] * df_good[cfg.column["p"]] ** 2,
            P_fit[1] * df_good[cfg.column["p"]],
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
    log.info("Calibrating conductivity")
    # calculate BTLCOND values from autosal data
    btl_df[cfg.column["refC"]] = convert.CR_to_cond(
        btl_df["CRavg"],
        btl_df["BathTEMP"],
        btl_df[cfg.column["t1"]],
        btl_df[cfg.column["p"]],
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
        btl_df = btl_df.merge(handcoded_salts, on=["SSSCC", "btl_fire_num"], how="left")
        # TODO: may be easier to try using flagging._merge_flags()?
        btl_df["SALNTY_FLAG_W"] = flagging.nan_values(
            btl_df["SALNTY"], old_flags=btl_df["SALNTY_FLAG_W"]
        )
    else:
        btl_df["SALNTY_FLAG_W"] = flagging.nan_values(btl_df["SALNTY"])

    ssscc_subsets = sorted(Path(cfg.dirs["ssscc"]).glob("ssscc_c*.csv"))
    if not ssscc_subsets:  # if no c-segments exists, write one from full list
        log.debug(
            "No CTDCOND grouping file found... creating ssscc_c1.csv with all casts"
        )
        ssscc_list = process_ctd.get_ssscc_list()
        ssscc_subsets = [Path(cfg.dirs["ssscc"] + "ssscc_c1.csv")]
        pd.Series(ssscc_list).to_csv(ssscc_subsets[0], header=None, index=False)
    qual_flag_c1 = pd.DataFrame()
    qual_flag_c2 = pd.DataFrame()
    coef_c1_all = pd.DataFrame()
    coef_c2_all = pd.DataFrame()

    fit_yaml = load_fit_yaml(to_object=True)  # load fit polynomial order

    for f in ssscc_subsets:
        # 0) grab ssscc chunk to fit
        ssscc_sublist = pd.read_csv(f, header=None, dtype="str", squeeze=True).to_list()
        btl_rows = btl_df["SSSCC"].isin(ssscc_sublist).values
        good_rows = btl_rows & (btl_df["SALNTY_FLAG_W"] == 2)
        time_rows = time_df["SSSCC"].isin(ssscc_sublist).values

        # 1) plot pre-fit residual
        f_stem = f.stem  # get "ssscc_c*" from path
        ctd_plots._intermediate_residual_plot(
            btl_df.loc[btl_rows, cfg.column["refC"]]
            - btl_df.loc[btl_rows, cfg.column["c1"]],
            btl_df.loc[btl_rows, cfg.column["p"]],
            btl_df.loc[btl_rows, "SSSCC"],
            xlabel="C1 Residual (mS/cm)",
            f_out=f"{cfg.fig_dirs['c1']}residual_{f_stem}_prefit.pdf",
        )
        ctd_plots._intermediate_residual_plot(
            btl_df.loc[btl_rows, cfg.column["refC"]]
            - btl_df.loc[btl_rows, cfg.column["c2"]],
            btl_df.loc[btl_rows, cfg.column["p"]],
            btl_df.loc[btl_rows, "SSSCC"],
            xlabel="C2 Residual (mS/cm)",
            f_out=f"{cfg.fig_dirs['c2']}residual_{f_stem}_prefit.pdf",
        )

        # TODO: allow for cast-by-cast T_order/P_order/zRange
        # TODO: truncate coefs (10 digits? look at historical data)
        # 2 & 3) calculate fit params
        # NOTE: df_bad_c1/2 will be overwritten during post-fit data flagging
        # but are left here for future debugging (if necessary)
        coef_c1, df_bad_c1 = _get_C_coefs(
            btl_df[good_rows],
            C_col=cfg.column["c1"],
            P_order=fit_yaml.fit_orders1[f_stem]["P_order"],
            T_order=fit_yaml.fit_orders1[f_stem]["T_order"],
            C_order=fit_yaml.fit_orders1[f_stem]["C_order"],
            zRange=fit_yaml.fit_orders1[f_stem]["zRange"],
            f_stem=f_stem,
        )
        coef_c2, df_bad_c2 = _get_C_coefs(
            btl_df[good_rows],
            C_col=cfg.column["c2"],
            P_order=fit_yaml.fit_orders2[f_stem]["P_order"],
            T_order=fit_yaml.fit_orders2[f_stem]["T_order"],
            C_order=fit_yaml.fit_orders2[f_stem]["C_order"],
            zRange=fit_yaml.fit_orders2[f_stem]["zRange"],
            f_stem=f_stem,
        )

        # 4) apply fit
        btl_df.loc[btl_rows, cfg.column["c1"]] = _conductivity_polyfit(
            btl_df.loc[btl_rows, cfg.column["c1"]],
            btl_df.loc[btl_rows, cfg.column["t1"]],
            btl_df.loc[btl_rows, cfg.column["p"]],
            coef_c1,
        )
        btl_df.loc[btl_rows, cfg.column["c2"]] = _conductivity_polyfit(
            btl_df.loc[btl_rows, cfg.column["c2"]],
            btl_df.loc[btl_rows, cfg.column["t2"]],
            btl_df.loc[btl_rows, cfg.column["p"]],
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
            param=cfg.column["c1"],
            ref=cfg.column["refC"],
            f_out=f"{cfg.fig_dirs['c1']}residual_{f_stem}.pdf",
        )
        df_ques_c2, df_bad_c2 = _flag_btl_data(
            btl_df[btl_rows],
            param=cfg.column["c2"],
            ref=cfg.column["refC"],
            f_out=f"{cfg.fig_dirs['c2']}residual_{f_stem}.pdf",
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
        btl_df[cfg.column["refC"]] - btl_df[cfg.column["c1"]],
        btl_df[cfg.column["p"]],
        btl_df["SSSCC"],
        xlabel="C1 Residual (mS/cm)",
        show_thresh=True,
        f_out=f"{cfg.fig_dirs['c1']}residual_all_postfit.pdf",
    )
    ctd_plots._intermediate_residual_plot(
        btl_df[cfg.column["refC"]] - btl_df[cfg.column["c2"]],
        btl_df[cfg.column["p"]],
        btl_df["SSSCC"],
        xlabel="C2 Residual (mS/cm)",
        show_thresh=True,
        f_out=f"{cfg.fig_dirs['c2']}residual_all_postfit.pdf",
    )

    # export cond quality flags
    qual_flag_c1.sort_index().to_csv(
        cfg.dirs["logs"] + "qual_flag_c1.csv",
        index=False,
    )
    qual_flag_c2.sort_index().to_csv(
        cfg.dirs["logs"] + "qual_flag_c2.csv",
        index=False,
    )

    # export cond fit params
    coef_c1_all[coef_names] = coef_c1_all[coef_names].applymap(
        lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
    )
    coef_c2_all[coef_names] = coef_c2_all[coef_names].applymap(
        lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
    )
    coef_c1_all.to_csv(cfg.dirs["logs"] + "fit_coef_c1.csv", index=False)
    coef_c2_all.to_csv(cfg.dirs["logs"] + "fit_coef_c2.csv", index=False)

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
    # TODO: flag time using handcoded salts somehow? discrete vs continuous
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

    return btl_df, time_df

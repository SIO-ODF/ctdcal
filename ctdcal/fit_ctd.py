"""
A module for fitting CTD discrete/bottle and continuous data.
"""

#!/usr/bin/env python
import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import yaml
from scipy.ndimage import shift

from ctdcal.fitting.common import get_node, NodeNotFoundError
from . import convert as convert
from . import ctd_plots as ctd_plots
from . import flagging as flagging
from . import get_ctdcal_config
from . import process_ctd as process_ctd
from .common import validate_dir

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def load_fit_yaml(fname, to_object=False):
    """Load polynomial fit order information from .yaml file."""

    if not Path(fname).exists():
        log.warning("Warning: Coefficients fit order YAML does not exist. Generating from scratch...")
        generate_yaml(fname)

    with open(fname, "r") as f:
        ymlfile = yaml.safe_load(f)

    if to_object:
        return type("ymlfile", (object,), ymlfile)
    else:
        return ymlfile

def generate_yaml(fname):
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
    with open(fname, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)

def write_fit_yaml():
    """For future use with automated fitting routine(s).
    i.e., iterate to find best fit parameters, save to file"""
    pass


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
        if param == 'CTDTMP1':
            xlabel = "T1 Residual (T90 C)"
        elif param == 'CTDTMP2':
            xlabel = "T2 Residual (T90 C)"
        elif param == "CTDCOND1":
            xlabel = "C1 Residual (mS/cm)"
        elif param == 'CTDCOND2':
            xlabel = "C2 Residual (mS/cm)"
        f_out = Path(str(f_out).replace(".pdf", "_postfit.pdf"))
        # f_out = f_out.split(".pdf")[0] + "_postfit.pdf"
        ctd_plots._intermediate_residual_plot(
            df["Diff"],
            df[prs],
            df["SSSCC"],
            show_thresh=True,
            xlabel=xlabel,
            f_out=f_out,
        )
        f_out = Path(str(f_out).replace(".pdf", "_flag2.pdf"))
        # f_out = f_out.split(".pdf")[0] + "_flag2.pdf"
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
        raise ValueError(
            "length of coef_names must match the number of dependent variables"
        )

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
        Independent variable fit coefficients (i.e., (coef0, ..., coefN))
    args : tuple of (array-like, (float, float, ...))
        Dependent variable data and fit coefficients (i.e., (data, (coef1, ..., coefN)))

    Returns
    -------
    fitted_y : array-like
        Independent variable data with polynomial fit correction applied

    Examples
    --------
    Behavior without additional args:

    >>> y = [2, 4, 6]
    >>> apply_polyfit(y, (1, 2, 3))  # y0 = 1; y1 = 2; y2 = 3
    array([ 19.,  61., 127.])

    where fitted_y = y + y0 + (y1 * y) + (y2 * y ** 2)

    Behavior with additional args:

    >>> y = [2, 4, 6]
    >>> x = [1, 2, 3]
    >>> apply_polyfit(y, (1,), (x, (2, 3)))  # y0 = 1; x1 = 2; x2 = 3
    array([ 8., 21., 40.])

    where fitted_y = y + y0 + (x1 * x) + (x2 * x ** 2)
    """
    fitted_y = np.copy(y).astype(float)
    for n, coef in enumerate(y_coefs):
        fitted_y += coef * np.power(y, n)

    for arg in args:
        if type(arg) is not tuple:
            raise TypeError(f"Positional args must be tuples, not {type(arg)}")

        series, coefs = arg
        for n, coef in enumerate(coefs):
            fitted_y += coef * np.power(series, n + 1)

    return fitted_y


def calibrate_temp(btl_df, time_df, datadir, inst, cast_list):
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
    fit_dir = validate_dir(Path(datadir, 'fit', inst), create=True)
    fit_groups_dir = Path(datadir, 'fit_groups', inst)
    validate_dir(fit_groups_dir, create=True)
    ssscc_subsets = sorted(fit_groups_dir.glob("ssscc_t*.csv"))
    if not ssscc_subsets:  # if no t-segments exists, write one from full list
        log.debug(
            "No CTDTMP grouping file found... creating ssscc_t1.csv with all casts"
        )
        ssscc_subsets = [Path(fit_groups_dir, "ssscc_t1.csv")]
        pd.Series(cast_list).to_csv(ssscc_subsets[0], header=None, index=False)

    fit_yaml = load_fit_yaml(Path(fit_groups_dir, 'fit_coeffs.yml'))  # load fit polynomial order
    for tN in ['CTDTMP1', 'CTDTMP2']:
        fig_dir = validate_dir(Path(datadir, 'fig', 'inst', tN), create=True)
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
            ctd_plots._intermediate_residual_plot(
                btl_df.loc[btl_rows, "REFTMP"]
                - btl_df.loc[btl_rows, tN],
                btl_df.loc[btl_rows, "CTDPRS"],
                btl_df.loc[btl_rows, "SSSCC"],
                xlabel=f"{tN.upper()} Residual (T90 C)",
                f_out=Path(fig_dir, f"{tN}residual_{f_stem}_prefit.pdf"),
            )

            # 2) prepare data for fitting
            # NOTE: df_bad will be overwritten during post-fit data flagging but is
            # left here for future debugging (if necessary)
            df_good, df_bad = _prepare_fit_data(
                btl_df[good_rows],
                tN,
                "REFTMP",
                zRange=fit_yaml[tN][f_stem]["zRange"],
            )
            ctd_plots._intermediate_residual_plot(
                df_good["Diff"],
                df_good["CTDPRS"],
                df_good["SSSCC"],
                xlabel=f"{tN.upper()} Residual (T90 C)",
                f_out=Path(fig_dir, f"{tN}residual_{f_stem}_fit_data.pdf"),
            )

            # 3) calculate fit coefs
            P_order = fit_yaml[tN][f_stem]["P_order"]
            T_order = fit_yaml[tN][f_stem]["T_order"]
            coef_dict = multivariate_fit(
                df_good["Diff"],
                (df_good["CTDPRS"], P_order),
                (df_good[tN], T_order),
                coef_names=["cp", "ct"],
            )

            # 4) apply fit
            P_coefs = tuple(coef_dict[f"cp{n}"] for n in np.arange(1, P_order + 1))
            T_coefs = tuple(coef_dict[f"ct{n}"] for n in np.arange(1, T_order + 1))
            btl_df.loc[btl_rows, tN] = apply_polyfit(
                btl_df.loc[btl_rows, tN],
                (coef_dict["c0"],) + T_coefs,
                (btl_df.loc[btl_rows, "CTDPRS"], P_coefs),
            )
            time_df.loc[time_rows, tN] = apply_polyfit(
                time_df.loc[time_rows, tN],
                (coef_dict["c0"],) + T_coefs,
                (time_df.loc[time_rows, "CTDPRS"], P_coefs),
            )

            # 4.5) flag CTDTMP and make residual plots
            df_ques, df_bad = _flag_btl_data(
                btl_df[btl_rows],
                param=tN,
                ref="REFTMP",
                f_out=Path(fig_dir, f"{tN}residual_{f_stem}.pdf"),
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
            btl_df["REFTMP"] - btl_df[tN],
            btl_df["CTDPRS"],
            btl_df["SSSCC"],
            xlabel=f"{tN.upper()} Residual (T90 C)",
            show_thresh=True,
            f_out=Path(fig_dir, f"{tN}residual_all_postfit.pdf"),
        )

        # export temp quality flags
        T_flag.sort_index().to_csv(Path(fit_dir, f"qual_flag_{tN}.csv"), index=False)

        # export temp fit params (formated to 5 sig figs, scientific notation)
        T_fit_coefs[coef_names] = T_fit_coefs[coef_names].applymap(
            lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
        )
        T_fit_coefs.to_csv(Path(fit_dir, f"fit_coef_{tN}.csv"), index=False)

    # flag temperature data
    time_df["CTDTMP_FLAG_W"] = 2

    return True


def calibrate_cond(btl_df, time_df, datadir, inst, ref, cast_list, bottleflags_man):
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
    fit_dir = validate_dir(Path(datadir, 'fit', inst), create=True)
    fit_groups_dir = Path(datadir, 'fit_groups', inst)
    validate_dir(fit_groups_dir, create=True)
    ssscc_subsets = sorted(fit_groups_dir.glob("ssscc_c*.csv"))
    # calculate BTLCOND values from lab salinity
    if not ssscc_subsets:  # if no c-segments exists, write one from full list
        log.debug(
                "No CTDCOND grouping file found... creating ssscc_c1.csv with all casts"
        )
        ssscc_subsets = [Path(fit_groups_dir, "ssscc_c1.csv")]
        pd.Series(cast_list).to_csv(ssscc_subsets[0], header=None, index=False)

    btl_df["BTLCOND"] = convert.sal_to_cond(
        btl_df["SALNTY"],
        btl_df["CTDTMP1"],
        btl_df["CTDPRS"],
    )

    # merge in handcoded salt flags
    flag_file = Path(datadir, 'flag', bottleflags_man)
    salt_flags_manual = None
    if flag_file.exists():
        try:
            salt_flags_manual = get_node(flag_file, ref)
        except NodeNotFoundError:
            log.info("No previously flagged values for %s found in flag file." % ref)
    else:
        log.info("No pre-existing flag file found.")

    if salt_flags_manual is not None:
        log.info("Merging previously flagged values for %s." % ref)
        salt_flags_manual_df = pd.DataFrame.from_dict(salt_flags_manual)
        salt_flags_manual_df = salt_flags_manual_df.rename(
            columns={"cast_id": "SSSCC", "bottle_num": "btl_fire_num", "value": "SALNTY_FLAG_W"}
        ).drop(columns=["notes"])
        btl_df = btl_df.merge(salt_flags_manual_df, on=["SSSCC", "btl_fire_num"], how="left")
        btl_df["SALNTY_FLAG_W"] = flagging.nan_values(
            btl_df["SALNTY"], old_flags=btl_df["SALNTY_FLAG_W"]
        )
    else:
        btl_df["SALNTY_FLAG_W"] = flagging.nan_values(btl_df["SALNTY"])

    fit_yaml = load_fit_yaml(Path(fit_groups_dir, 'fit_coeffs.yml'))  # load fit polynomial order
    for cN, tN in zip(['CTDCOND1', 'CTDCOND2'], ['CTDTMP1', 'CTDTMP2']):
        fig_dir = validate_dir(Path(datadir, 'fig', 'inst', cN), create=True)
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
            ctd_plots._intermediate_residual_plot(
                btl_df.loc[btl_rows, 'BTLCOND']
                - btl_df.loc[btl_rows, cN],
                btl_df.loc[btl_rows, 'CTDPRS'],
                btl_df.loc[btl_rows, "SSSCC"],
                xlabel=f"{cN.upper()} Residual (mS/cm)",
                f_out=Path(fig_dir, f"{cN}residual_{f_stem}_prefit.pdf"),
            )

            # 2) prepare data for fitting
            # NOTE: df_bad will be overwritten during post-fit data flagging
            # but is left here for future debugging (if necessary)
            df_good, df_bad = _prepare_fit_data(
                btl_df[good_rows],
                cN,
                'BTLCOND',
                zRange=fit_yaml[cN][f_stem]["zRange"],
            )
            ctd_plots._intermediate_residual_plot(
                df_good["Diff"],
                df_good['CTDPRS'],
                df_good["SSSCC"],
                xlabel=f"{cN.upper()} Residual (mS/cm)",
                f_out=Path(fig_dir, f"{cN}residual_{f_stem}_fit_data.pdf"),
            )

            # 3) calculate fit coefs
            P_order = fit_yaml[cN][f_stem]["P_order"]
            T_order = fit_yaml[cN][f_stem]["T_order"]
            C_order = fit_yaml[cN][f_stem]["C_order"]
            coef_dict = multivariate_fit(
                df_good["Diff"],
                (df_good['CTDPRS'], P_order),
                (df_good[tN], T_order),
                (df_good[cN], C_order),
                coef_names=["cp", "ct", "cc"],
            )

            # 4) apply fit
            P_coefs = tuple(coef_dict[f"cp{n}"] for n in np.arange(1, P_order + 1))
            T_coefs = tuple(coef_dict[f"ct{n}"] for n in np.arange(1, T_order + 1))
            C_coefs = tuple(coef_dict[f"cc{n}"] for n in np.arange(1, C_order + 1))
            btl_df.loc[btl_rows, cN] = apply_polyfit(
                btl_df.loc[btl_rows, cN],
                (coef_dict["c0"],) + C_coefs,
                (btl_df.loc[btl_rows, 'CTDPRS'], P_coefs),
                (btl_df.loc[btl_rows, tN], T_coefs),
            )
            time_df.loc[time_rows, cN] = apply_polyfit(
                time_df.loc[time_rows, cN],
                (coef_dict["c0"],) + C_coefs,
                (time_df.loc[time_rows, 'CTDPRS'], P_coefs),
                (time_df.loc[time_rows, tN], T_coefs),
            )

            # 4.5) flag CTDCOND and make residual plots
            df_ques, df_bad = _flag_btl_data(
                btl_df[btl_rows],
                param=cN,
                ref='BTLCOND',
                f_out=Path(fig_dir, f"{cN}residual_{f_stem}.pdf"),
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
            btl_df['BTLCOND'] - btl_df[cN],
            btl_df['CTDPRS'],
            btl_df["SSSCC"],
            xlabel=f"{cN.upper()} Residual (mS/cm)",
            show_thresh=True,
            f_out=Path(fig_dir, f"{cN}residual_all_postfit.pdf"),
        )

        # export cond quality flags
        C_flag.sort_index().to_csv(Path(fit_dir, f"qual_flag_{cN}.csv"), index=False)

        # export cond fit params
        C_fit_coefs[coef_names] = C_fit_coefs[coef_names].applymap(
            lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
        )
        C_fit_coefs.to_csv(Path(fit_dir, f"fit_coef_{cN}.csv"), index=False)

    # recalculate salinity with calibrated C/T
    time_df['CTDSAL'] = gsw.SP_from_C(
        time_df['CTDCOND1'],
        time_df['CTDTMP1'],
        time_df['CTDPRS'],
    )
    btl_df['CTDSAL'] = gsw.SP_from_C(
        btl_df['CTDCOND1'],
        btl_df['CTDTMP1'],
        btl_df['CTDPRS'],
    )

    # flag salinity data
    time_df['CTDSAL_FLAG_W'] = 2
    btl_df['CTDSAL_FLAG_W'] = flagging.by_residual(
        btl_df['CTDSAL'],
        btl_df["SALNTY"],
        btl_df['CTDPRS'],
    )
    bad_rows = btl_df["SALNTY_FLAG_W"].isin([3, 4])
    btl_df.loc[bad_rows, 'CTDSAL_FLAG_W'] = 2  # bad salts not used for QC
    btl_df["CTDSAL_FLAG_W"] = flagging.nan_values(
        btl_df['CTDSAL'], old_flags=btl_df["CTDSAL_FLAG_W"]
    )

    return btl_df, time_df

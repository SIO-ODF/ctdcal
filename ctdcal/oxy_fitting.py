"""
Module for processing oxygen from CTD and bottle samples.
"""

import csv
import logging
import xml.etree.cElementTree as ET
from collections import OrderedDict
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import scipy

from . import ctd_plots as ctd_plots
from . import flagging as flagging
from . import get_ctdcal_config
from . import process_ctd as process_ctd
from . import sbe_reader as sbe_rd

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def load_winkler_oxy(oxy_file):
    """
    Load Winkler oxygen titration data file.

    Parameters
    ----------
    oxy_file : str or Path
        Path to oxygen file

    Returns
    -------
    df : DataFrame
        Oxygen data
    params : list of str
        List of oxygen parameters used in titration
    """

    with open(oxy_file, newline="") as f:
        oxyF = csv.reader(
            f, delimiter=" ", quoting=csv.QUOTE_NONE, skipinitialspace="True"
        )
        oxy_array = []
        for row in oxyF:
            if len(row) > 9:
                row = row[:9]
            oxy_array.append(row)

    # TODO turn params into a dict with useful labels
    params = oxy_array.pop(0)  # save file header info for later (Winkler values)
    cols = OrderedDict(
        [
            ("STNNO_OXY", int),
            ("CASTNO_OXY", int),
            ("BOTTLENO_OXY", int),
            ("FLASKNO", int),
            ("TITR_VOL", float),
            ("TITR_TEMP", float),
            ("DRAW_TEMP", float),
            ("TITR_TIME", int),
            ("END_VOLTS", float),
        ]
    )

    df = pd.DataFrame(oxy_array, columns=cols.keys()).astype(cols)
    df = df[df["BOTTLENO_OXY"] != 99]  # remove "Dummy Data"
    df = df[df["TITR_VOL"] > 0]  # remove "ABORTED DATA"
    df = df.sort_values("BOTTLENO_OXY").reset_index(drop=True)
    df["FLASKNO"] = df["FLASKNO"].astype(str)

    return df, params


def load_flasks(flask_file=cfg.dirs["oxygen"] + "o2flasks.vol", comment="#"):
    """
    Load oxygen flask information from .vol file.

    Parameters
    ----------
    flask_file : str or Path, optional
        Path to flask file
    comment : str, optional
        Identifier signifying line is a comment and should be skipped

    Returns
    -------
    flasks : DataFrame
        Flask numbers and volumes
    """
    with open(flask_file, "r") as f:
        flasks = []
        for line in f:
            is_comment = line.strip().startswith(comment)
            if ("Volume" in line) or is_comment:
                continue
            num, vol = line.strip().split()[:2]  # only need first two cols (#, volume)
            flasks.append([str(num), float(vol)])

    flasks = pd.DataFrame(flasks, columns=["FLASKNO", "FLASK_VOL"])

    return flasks


def correct_flask_vol(flask_vol, t=20.0, glass="borosilicate"):
    """
    Correct flask volume for changes from thermal expansion of glass.

    Parameters
    ----------
    flask_vol : array-like
        Flask volumes at standard temperature (20C)
    t : float, optional
        New temperature to calculate volume
    glass : str, optional
        Type of glass ("borosilicate" or "soft)

    Returns
    -------
    corrected_vol : array-like
        Flask volumes are new temperature

    Notes
    -----
    Flask volume equation from 2007 Best Practices for Ocean CO2 Measurements,
    SOP 13 - Gravimetric calibration of volume contained using water
    """
    alpha = {  # thermal expansion coefficient
        "borosilicate": 1.0e-5,
        "soft": 2.5e-3,
    }
    if glass not in alpha.keys():
        raise KeyError(f"Glass type not found, must be one of {list(alpha.keys())}")
    standard_t = 20.0
    corrected_vol = flask_vol * (1.0 + alpha[glass] * (t - standard_t))

    return corrected_vol


def gather_oxy_params(oxy_file):
    """
    Collect Winkler oxygen measurement parameters from LabVIEW data file headers.

    Parameters
    ----------
    oxy_file : str or Path
        Path to oxygen file

    Returns
    -------
    df : DataFrame
        Oxygen measurement parameters
    """
    with open(oxy_file, newline="") as f:
        header = f.readline()

    param_list = header.split()[:6]
    params = pd.DataFrame(param_list, dtype=float).transpose()
    params.columns = ["V_std", "V_blank", "N_KIO3", "V_KIO3", "T_KIO3", "T_thio"]

    return params


def calculate_bottle_oxygen(ssscc_list, ssscc_col, titr_vol, titr_temp, flask_nums):
    """
    Wrapper function for collecting parameters and calculating oxygen values from
    Winkler titrations.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to process
    ssscc_col : array-like
        Station/cast for each sample taken
    titr_vol : array-like
        Titration volume [mL]
    titr_temp : array-like
        Temperature of titration [degC]
    flask_nums : array-like
        Oxygen flask used for each sample

    Returns
    -------
    oxy_mL_L : array-like
        Oxygen concentration [mL/L]

    Notes
    -----
    Titration equation comes from WHP Operations and Methods, Culberson (1991):
    https://cchdo.github.io/hdo-assets/documentation/manuals/pdf/91_1/culber2.pdf

    """
    params = pd.DataFrame()
    for ssscc in ssscc_list:
        df = gather_oxy_params(cfg.dirs["oxygen"] + ssscc)
        df["SSSCC"] = ssscc
        params = pd.concat([params, df])

    # get flask volumes and merge with titration parameters
    flask_df = load_flasks()  # TODO: volume correction from thermal expansion?
    volumes = pd.merge(flask_nums, flask_df, how="left")["FLASK_VOL"].values
    params = pd.merge(ssscc_col, params, how="left")

    # find 20degC equivalents
    rho_20C = gsw.rho_t_exact(0, 20, 0)
    rho_T_KIO3 = gsw.rho_t_exact(0, params["T_KIO3"], 0)
    N_KIO3_20C = params["N_KIO3"] * (rho_T_KIO3 / rho_20C)

    # TODO: does KIO3 volume get corrected? what is the recorded value?
    # V_KIO3_20C = correct_flask_vol(params["V_KIO3"], t=params["T_KIO3"])

    # calculate O2 concentration (in mL/L)
    E = 5598  # stoichiometric relationship between thio_n and DO
    DO_reg = 0.0017  # correction for oxygen added by reagents
    V_reg = 2.0  # volume of reagents (mL)

    oxy_mL_L = (
        (
            ((titr_vol.values - params["V_blank"]) * params["V_KIO3"] * N_KIO3_20C * E)
            / (params["V_std"] - params["V_blank"])
            - 1000 * DO_reg
        )
    ) / (volumes - V_reg)

    return oxy_mL_L.values


def hysteresis_correction(oxygen, pressure, H1=-0.033, H2=5000, H3=1450, freq=24):
    """
    Remove hysteresis effects from oxygen concentration values.

    Oxygen hysteresis can be corrected before conversion from volts to oxygen
    concentration, see equations_sbe.sbe43_hysteresis_voltage()

    # TODO: should this just be a wrapper that calls sbe43_hysteresis_voltage()?

    Parameters
    ----------
    oxygen : array-like
        Oxygen concentration values
    pressure : array-like
        CTD pressure values (dbar)
    H1 : scalar, optional
        Amplitude of hysteresis correction function (range: -0.02 to -0.05)
    H2 : scalar, optional
        Function constant or curvature function for hysteresis
    H3 : scalar, optional
        Time constant for hysteresis (seconds) (range: 1200 to 2000)
    freq : scalar, optional
        CTD sampling frequency (Hz)

    Returns
    -------
    oxy_corrected : array-like
        Hysteresis-corrected oxygen concentration values (with same units as input)

    Notes
    -----
    See Application Note 64-3 for more information.
    """
    # TODO: vectorize (if possible), will probably require matrix inversion
    dt = 1 / freq
    D = 1 + H1 * (np.exp(pressure / H2) - 1)
    C = np.exp(-1 * dt / H3)

    oxy_corrected = np.zeros(oxygen.shape)
    oxy_corrected[0] = oxygen[0]
    for i in np.arange(1, len(oxygen)):
        oxy_corrected[i] = (
            oxygen[i] + (oxy_corrected[i - 1] * C * D[i]) - (oxygen[i - 1] * C)
        ) / D[i]

    return oxy_corrected


def oxy_ml_to_umolkg(oxy_mL_L, sigma0):
    """Convert dissolved oxygen from units of mL/L to micromol/kg.

    Parameters
    ----------
    oxy_mL_L : array-like
        Dissolved oxygen in units of [mL/L]
    sigma0 : array-like
        Potential density anomaly (i.e. sigma - 1000) referenced to 0 dbar [kg/m^3]

    Returns
    -------
    oxy_umol_kg : array-like
        Dissolved oxygen in units of [umol/kg]

    Notes
    -----
    Conversion value 44660 is exact for oxygen gas and derived from the ideal gas law.
    (c.f. Sea-Bird Application Note 64, pg. 6)
    """

    oxy_umol_kg = oxy_mL_L * 44660 / (sigma0 + 1000)

    return oxy_umol_kg


def oxy_umolkg_to_ml(oxy_umol_kg, sigma0):
    """Convert dissolved oxygen from units of micromol/kg to mL/L.

    Parameters
    ----------
    oxy_umol_kg : array-like
        Dissolved oxygen in units of [umol/kg]
    sigma0 : array-like
        Potential density anomaly (i.e. sigma - 1000) referenced to 0 dbar [kg/m^3]

    Returns
    -------
    oxy_mL_L : array-like
        Dissolved oxygen in units of [mL/L]

    Notes
    -----
    Conversion value 44660 is exact for oxygen gas and derived from the ideal gas law.
    (c.f. Sea-Bird Application Note 64, pg. 6)
    """

    oxy_mL_L = oxy_umol_kg * (sigma0 + 1000) / 44660

    return oxy_mL_L


def calculate_dV_dt(oxy_volts, time, nan_replace=True):
    """
    Calculate the time derivative of oxygen voltage.

    Parameters
    ----------
    oxy_volts : array-like
        Oxygen sensor voltage output
    time : array-like
        Time from oxygen sensor (must be same length as oxy_volts)
    nan_replace : bool, optional
        Replace nans in time derivative with the mean value

    Returns
    -------
    dV_dt : array-like
        Time derivative of oxygen voltage
    """
    # TODO: experiment with dt, filtering
    # Uchida (2008): dV/dt "estimated by linear fits over 2 second intervals"
    # should dt just be 1 / freq? i.e. 1/24 Hz

    dV = np.diff(oxy_volts)  # central differences shorten vectors by 1
    dt = np.diff(time)
    dt[dt == 0] = np.median(dt[dt > 0])  # replace with median to avoid dividing by zero

    dV_dt = dV / dt
    dV_dt = np.insert(dV_dt, 0, 0)  # add zero in front to match original length
    dV_dt[np.isinf(dV_dt)] = np.nan  # this check is probably unnecessary

    if nan_replace:
        dV_dt = np.nan_to_num(dV_dt, nan=np.nanmean(dV_dt))

    # TODO: should we do some kind of filtering? e.g.:
    # (PMEL does this calculation on binned data already so filtering is not the same)
    # a = 1
    # windowsize = 5
    # b = (1 / windowsize) * np.ones(windowsize)
    # filtered_dvdt = scipy.signal.filtfilt(b, a, dv_dt)

    return dV_dt  # filtered_dvdt


def _get_sbe_coef(idx=0):
    """
    Get SBE oxygen coefficients from raw .xmlcon files.
    Defaults to using first station in ssscc.csv file.

    Returns the following tuple of coefficients: Soc, offset, Tau20, Tcor, E
    """
    # TODO: does scipy's minimize function needs a tuple? can this be improved further?

    station = process_ctd.get_ssscc_list()[idx]
    xmlfile = cfg.dirs["raw"] + station + ".XMLCON"

    tree = ET.parse(xmlfile)
    root_eq0 = tree.find(".//CalibrationCoefficients[@equation='0']")  # Owens-Millard
    root_eq1 = tree.find(".//CalibrationCoefficients[@equation='1']")  # SBE equation

    coefs = {c.tag: float(c.text) for c in root_eq1}
    coefs["Tcor"] = float(root_eq0.find("Tcor").text)  # only coef needed from eq0
    keep_keys = ["Soc", "offset", "Tau20", "Tcor", "E"]

    return tuple(coefs[key] for key in keep_keys)


def calculate_weights(pressure):
    """
    Calculate weights (as a function of pressure) for weighted least squares fitting.
    Deep measurements are weighted higher than shallow.

    Parameters
    ----------
    presssure : array-like
        Pressure values of oxygen measurements [dbar]

    Returns
    -------
    weights : array-like
        Weight factor for each pressure value
    """
    # TODO: automatic weight calculation rather than hardcoded (machine learning?)

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


"""code_pruning: should this be here or in equations_sbe? somewhere else?"""


def _PMEL_oxy_eq(coefs, inputs, cc=[1.92634e-4, -4.64803e-2]):
    """
    Modified oxygen equation for SBE 43 used by NOAA/PMEL
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E
    """
    Soc, Voff, Tau20, Tcorr, E = coefs
    oxyvolts, pressure, temp, dvdt, os = inputs
    o2 = (
        Soc
        * (
            oxyvolts
            + Voff
            + Tau20 * np.exp(cc[0] * pressure + cc[1] * (temp - 20)) * dvdt
        )
        * os
        * np.exp(Tcorr * temp)
        * np.exp((E * pressure) / (temp + 273.15))
    )

    return o2


def PMEL_oxy_weighted_residual(coefs, weights, inputs, refoxy, L_norm=2):
    # TODO: optionally include other residual types
    # (abstracted from PMEL code oxygen_cal_ml.m)
    # unweighted L2: sum((ref - oxy)^2)  # if weighted fails
    # unweighted L4: sum((ref - oxy)^4)  # unsure of use case
    # unweighted L1: sum(abs(ref - oxy))  # very far from ideal
    # anything else? genericize with integer "norm" function input?

    residuals = np.sum(
        (weights * (refoxy - _PMEL_oxy_eq(coefs, inputs)) ** 2)
    ) / np.sum(weights ** 2)

    return residuals


def match_sigmas(
    btl_prs,
    btl_oxy,
    btl_tmp,
    btl_SA,
    ctd_os,
    ctd_prs,
    ctd_tmp,
    ctd_SA,
    ctd_oxyvolts,
    ctd_time,
):

    # Construct Dataframe from bottle and ctd values for merging
    btl_data = pd.DataFrame(
        data={"CTDPRS": btl_prs, "REFOXY": btl_oxy, "CTDTMP": btl_tmp, "SA": btl_SA}
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
        # TODO: can this be vectorized?
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


def sbe43_oxy_fit(merged_df, sbe_coef0=None, f_suffix=None):

    # Plot data to be fit together
    f_out = f"{cfg.fig_dirs['ox']}sbe43_residual{f_suffix}_prefit.pdf"
    ctd_plots._intermediate_residual_plot(
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
        # TODO: get some kind of logging in here in case things go awry
        # e.g. count of thrown values, start/final stdev, failing to converge, etc.
        bad_df = pd.concat([bad_df, thrown_values])
        merged_df = merged_df[np.abs(merged_df["residual"]) <= cutoff].copy()

    # intermediate plots to diagnose data chunks goodness
    # TODO: implement into bokeh/flask dashboard
    if f_suffix is not None:
        f_out = f"{cfg.fig_dirs['ox']}sbe43_residual{f_suffix}.pdf"
        ctd_plots._intermediate_residual_plot(
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


def prepare_oxy(btl_df, time_df, ssscc_list):
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

    Returns
    -------

    """
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
    btl_df["OXYGEN_FLAG_W"] = flagging.nan_values(btl_df[cfg.column["refO"]])
    # Load manual OXYGEN flags
    if Path("data/oxygen/manual_oxy_flags.csv").exists():
        manual_flags = pd.read_csv(
            "data/oxygen/manual_oxy_flags.csv", dtype={"SSSCC": str}
        )
        for _, flags in manual_flags.iterrows():
            df_row = (btl_df["SSSCC"] == flags["SSSCC"]) & (
                btl_df["btl_fire_num"] == flags["SAMPNO"]
            )
            btl_df.loc[df_row, "OXYGEN_FLAG_W"] = flags["Flag"]

    return True


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
    ctd_plots._intermediate_residual_plot(
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
        log.info(ssscc + " density matching done")

    # Only fit using OXYGEN flagged good (2)
    all_sbe43_merged = all_sbe43_merged[btl_df["OXYGEN_FLAG_W"] == 2].copy()

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

    # TODO: save outlier data from fits?
    # TODO: secondary oxygen flagging step (instead of just taking outliers from fit routine)

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
    time_df["CTDOXY_FLAG_W"] = 2  # TODO: actual flagging of some kind?
    btl_df["CTDOXY_FLAG_W"] = flagging.by_percent_diff(
        btl_df["CTDOXY"], btl_df["OXYGEN"], percent_thresh=1
    )

    # Plot all post fit data
    f_out = f"{cfg.fig_dirs['ox']}sbe43_residual_all_postfit.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df["OXYGEN"] - btl_df["CTDOXY"],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10, 10),
    )
    f_out = f"{cfg.fig_dirs['ox']}sbe43_residual_all_postfit_flag2.pdf"
    flag2 = btl_df["CTDOXY_FLAG_W"] == 2
    ctd_plots._intermediate_residual_plot(
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

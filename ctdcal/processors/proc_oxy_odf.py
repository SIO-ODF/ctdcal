"""
Module for parsing and processing oxygen from ODF LabVIEW oxygen
titration software
"""
import csv
import logging
from collections import OrderedDict

import gsw
import numpy as np
import pandas as pd

from ctdcal import get_ctdcal_config

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def load_winkler_oxy(oxy_file):
    """
    Load ODF Winkler oxygen raw titration data file.

    Files are generated from ODF LabVIEW software, generating a file named
    with the current run's first SSSCC.

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

    ssscc = oxy_file.stem
    with open(oxy_file, newline="") as f:
        oxyF = csv.reader(
            f, delimiter=" ", quoting=csv.QUOTE_NONE, skipinitialspace=True
        )
        oxy_array = []
        for row in oxyF:
            if row[0].startswith("#"):
                continue
            elif "ABORT" in row:
                log.warning(f"Aborted value found in oxy file for {ssscc} Skipping.")
                continue
            elif len(row) > 9:
                row = row[:9]
            oxy_array.append(row)

    params = oxy_array.pop(0)  # save file header info for later (Winkler values)
    if params in oxy_array:
        #   Check for duplicate standardizations. Check with analyst to see which standard to keep.
        log.warning(f"Raw Winkler contains duplicate standardization in {ssscc}")
    elif any(len(row[0]) > 4 for row in oxy_array):
        #   Check for cases where the ODF software is writing out a standard line (Winkler writes out 4 chars for station)
        log.warning(f"Raw Winkler contains abnormal row in {ssscc}")

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

    df = pd.DataFrame(oxy_array, columns=cols.keys()).astype(
        cols
    )  #   Force dtypes and build dataframe
    if df["BOTTLENO_OXY"].value_counts().get(99, 0) != 1:
        #   Multiple "dummy" data present
        log.warning(f"Multiple dummy lines in raw titration file {ssscc}")
    df = df[df["BOTTLENO_OXY"] != 99]  # remove "Dummy Data"

    if not df[df["TITR_VOL"] <= 0].empty:
        #   Check if titration volumes are all positive
        log.warning(f"Non-positive entries found for titration volume in {ssscc}")
    if any(df.duplicated()):
        #   Check for duplicates
        log.warning(f"Raw Winkler contains duplicate values in {ssscc}")
    if len(df["CASTNO_OXY"].value_counts()) != 1:
        #   Check if there are multiple casts in the file
        log.warning(f"Multiple casts reported in {ssscc}")
    elif len(df["STNNO_OXY"].value_counts()) != 1:
        #   Check if there are multiple stations in the file
        log.warning(f"Multiple stations reported in {ssscc}")
    if any((df["END_VOLTS"] < 0) | df["END_VOLTS"] > 5):
        #   Scale from 0-5 V
        log.warning(f"Titration file has erroneous voltage reported in {ssscc}")

    df = df[df["TITR_VOL"] > 0]  # remove "ABORTED DATA" or otherwise bad points
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
    titr_columns = ["V_std", "V_blank", "N_KIO3", "V_KIO3", "T_KIO3", "T_thio"]
    try:
        with open(oxy_file, newline="") as f:
            header = f.readline()

        param_list = header.split()[:6]
        params = pd.DataFrame(param_list, dtype=float).transpose()
        params.columns = titr_columns

    except FileNotFoundError:
        # fitting data before titration file has been received from oxygen analyst
        log.info(f"Failed to load {oxy_file} titration file, filling with NaNs")
        params = pd.DataFrame(np.nan, index=[0], columns=titr_columns)

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
    flask_df = load_flasks()
    volumes = pd.merge(flask_nums, flask_df, how="left")["FLASK_VOL"].values
    params = pd.merge(ssscc_col, params, how="left")

    # find 20degC equivalents
    rho_20C = gsw.rho_t_exact(0, 20, 0)
    rho_T_KIO3 = gsw.rho_t_exact(0, params["T_KIO3"], 0)
    N_KIO3_20C = params["N_KIO3"] * (rho_T_KIO3 / rho_20C)

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


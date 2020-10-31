"""
Module for processing oxygen from CTD and bottle samples.
"""

import csv
from pathlib import Path

from collections import OrderedDict
import config as cfg
import gsw
import numpy as np
import pandas as pd
import scipy

import ctdcal.ctd_plots as ctd_plots
import ctdcal.process_ctd as process_ctd
import ctdcal.sbe_reader as sbe_rd


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


def load_flasks(flask_file=cfg.directory["oxy"] + "o2flasks.vol", comment="#"):
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
        for l in f:
            is_comment = l.strip().startswith(comment)
            if ("Volume" in l) or is_comment:
                continue
            num, vol = l.strip().split()[:2]  # only need first two cols (# and volume)
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
        "borosilicate": 0.00001,
        "soft": 0.000025,
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
    params.columns = ["TITR", "BLANK", "KIO3_N", "KIO3_V", "KIO3_T", "THIO_T"]

    return params


def calculate_thio_norm(params, ssscc):
    """
    Calculate normality of thiosulfate used in Winkler oxygen titrations.

    Parameters
    ----------
    params : DataFrame
        Calibration parameters used for each station
    ssscc : Series
        Column of station/cast values FOR EACH measurement taken
        (i.e. length of ssscc should equal the total number of oxygen titrations)

    Returns
    -------
    thio_n : array-like
        Thiosulfate normality
    """
    # TODO: arraylike inputs instead of df? is this a function worth exposing?
    params = pd.merge(ssscc, params, how="left")

    rho_stp = gsw.rho_t_exact(0, 20, 0)
    rho_kio = gsw.rho_t_exact(0, params["KIO3_T"], 0)
    rho_thio = gsw.rho_t_exact(0, params["THIO_T"], 0)

    kio3_v_20c = params["KIO3_V"] * (rho_kio / rho_stp)
    thio_v_20c = params["TITR"] * (rho_thio / rho_stp) - params["BLANK"]

    thio_n = kio3_v_20c * params["KIO3_N"] / thio_v_20c

    return thio_n.values


def calculate_20C_vol(titr_vol, titr_temp):
    """
    Calculate the 20degC equivalent titration volume.

    Parameters
    ----------
    titr_vol : array-like
        Titration volume
    titr_temp : array-like
        Temperature of titration

    Returns
    -------
    titr_vol_20c : array-like
        Titration volume equivalent at 20degC
    """
    # TODO: is this func useful beyond titrations?
    rho_titr = gsw.rho_t_exact(0, titr_temp, 0)
    rho_stp = gsw.rho_t_exact(0, 20, 0)
    titr_vol_20c = titr_vol * rho_titr / rho_stp

    return titr_vol_20c


def calculate_bottle_oxygen(ssscc_list, ssscc_col, titr_vol, titr_temp, flask_nums):
    """
    Calculates oxygen values from Winkler titrations

    Parameters
    ----------

    ssscc_list :array-like
                list containing all of the stations to be calculated


    """
    params = pd.DataFrame()
    for ssscc in ssscc_list:
        df = gather_oxy_params(cfg.directory["oxy"] + ssscc)
        df["SSSCC"] = ssscc
        params = pd.concat([params, df])

    thio_n = calculate_thio_norm(params, ssscc_col)
    titr_20C = calculate_20C_vol(titr_vol, titr_temp)

    flask_df = load_flasks()
    volumes = pd.merge(flask_nums, flask_df, how="left")["FLASK_VOL"].values
    params = pd.merge(ssscc_col, params, how="left")

    oxy_mlL = oxygen_eq(titr_20C, params["BLANK"].values, thio_n, volumes)

    return oxy_mlL

def flag_winkler_oxygen(oxygen):
    flag = pd.Series(oxygen).copy()
    flag.loc[flag.notnull()] = 2
    flag.loc[flag.isnull()] = 9
    flag = flag.astype(int)
    return flag


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


def oxygen_eq(titr, blank, thio_n, flask_vol):
    E = 5.598  # L O2 equivalent (?)
    DO_rgts = 0.0017  # correction for oxygen added by reagents
    V_rgts = 2e-3  # volume of reagents (L)
    KIO3_V = 10.0  # volume of KIO3 standard (mL)
    KIO3_N = 0.01  # normality of KIO3 standard (N)
    oxyMl_L = ((titr - blank) * KIO3_V * KIO3_N * E) / ((flask_vol * 1e-3) - V_rgts) - DO_rgts
    breakpoint()

    # TODO: where does this eq come from? what are the magic numbers?
    oxy_mlL = (((titr - blank) * thio_n * 5.598 - 0.0017) / ((flask_vol - 2.0) * 0.001))

    return oxy_mlL

def oxy_from_ctd_eq(coef, oxyvolts, pressure, temp, dvdt, os, cc=[1.92634e-4,-4.64803e-2]):
    """
    Modified oxygen equation for SBE 43 used by NOAA/PMEL

    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E
    """


   # ox=c(1)*(doxvo+c(2)+c(3)*exp(cc(1)*dpres+cc(2)*dtepr).*dsdvdt).*os.*exp(c(4)*dtepr).*exp(c(5)*dpres./(273.15+dtepr));
    oxy_mlL = coef[0] * (oxyvolts + coef[1] + coef[2] * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt) * os \
            * np.exp(coef[3] * temp) \
            * np.exp((coef[4] * pressure) \
            / (temp + 273.15))


    return oxy_mlL

def SB_oxy_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc=[1.92634e-4,-4.64803e-2]):
    """
    Oxygen (ml/l) = Soc * (V + Voffset) * (1.0 + A * T + B * T^2 + C * T^3 ) * OxSol(T,S) * exp(E * P / K)
    Original equation from calib sheet dated 2014:

    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = A
    coef[3] = B
    coef[4] = C
    coef[5] = E
    coef[6] = Tau20

    """

    oxygen = (coef[0]) * ((oxyvolts + coef[1] + (coef[6] * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt))
                 * (1.0 + coef[2] * temp + coef[3] * np.power(temp,2) + coef[4] * np.power(temp,3) )
                 * os
                 * np.exp((coef[5] * pressure) / (temp + 273.15)))

    return oxygen

def get_SB_coef(hexfile,xmlfile):

    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1,
                        'ranking': 5, 'column': 'CTDOXYVOLTS',
                        'sensor_info': rawConfig['Sensors'][i]}


    Soc = oxy_meta['sensor_info']['Soc'] # Oxygen slope
    Voff = oxy_meta['sensor_info']['offset'] # Sensor output offset voltage
    Tau20 = oxy_meta['sensor_info']['Tau20'] #Sensor time constant at 20 deg C and 1 Atm
    Tcorr = oxy_meta['sensor_info']['Tcor'] # Temperature correction
    E  = oxy_meta['sensor_info']['E'] #Compensation Coef for pressure effect on membrane permeability (Atkinson et al)
    A = oxy_meta['sensor_info']['A'] # Compensation Coef for temp effect on mem perm
    B = oxy_meta['sensor_info']['B'] # Compensation Coef for temp effect on mem perm
    C = oxy_meta['sensor_info']['C'] # Compensation Coef for temp effect on mem perm
    D = [oxy_meta['sensor_info']['D0'],oxy_meta['sensor_info']['D1'],oxy_meta['sensor_info']['D2']] # Press effect on time constant Usually not fitted for.

    #coef = {'Soc':Soc,'Voff':Voff,'Tau20':Tau20,'Tcorr':Tcorr,'E':E,'A':A,'B':B,'C':C,'D':D}

    # Make into an array in order to do keast squares fitting easier
#        coef0s:
#    coef[0] = Soc
#    coef[1] = Voffset
#    coef[2] = Tau20
#    coef[3] = Tcorr
#    coef[4] = E
#
#    cc[0] = D1
#    cc[1] = D2
    coef = [Soc,Voff,A,B,C,E,Tau20]

    return coef


def sigma_from_CTD(sal, temp, press, lon, lat, ref=0):
    """
    Calculates potential density from CTD parameters at various reference pressures

    Parameters
    ----------

    sal :array-like
         Salinity in PSU (PSS-78)

    temp :array_like
          In-situ temperature in deg C

    press :array_like
           Pressure in dbar

    lon :array_like
         longitude in decimal degrees

    lat :array_like
         latitute in decimal degrees

    ref :int
         reference pressure point for caluclateing sigma0 (in ref * 1000 dBar ref[0-4])

    Returns
    -------

    simga :array-like
           Potential density calculated at a reference pressure of ref * 1000 dBar

    """

    CT = gsw.CT_from_t(sal, temp, press)

    SA = gsw.SA_from_SP(sal, press, lon, lat)

    # Reference pressure in ref*1000 dBars
    if ref == 0:
        sigma = gsw.sigma0(SA, CT)
    elif ref == 1:
        sigma = gsw.sigma1(SA, CT)
    elif ref == 2:
        sigma = gsw.sigma2(SA, CT)
    elif ref == 3:
        sigma = gsw.sigma3(SA, CT)
    elif ref == 4:
        sigma = gsw.sigma4(SA, CT)

    return sigma

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

def calculate_dVdT(oxyvolts, time):
#
     doxyv = np.diff(oxyvolts)


     dt = np.diff(time)

     # Replace 0 values with median value
     m = np.median(dt[dt > 0])
     dt[dt == 0] = m


     dv_dt = doxyv/dt
     dv_dt = np.insert(dv_dt,0,0)
     dv_dt_mean = np.mean(dv_dt[(dv_dt != np.inf) & (dv_dt !=-np.inf)])
     dv_dt[(dv_dt == np.inf) | (dv_dt == -np.inf)] = np.NaN
     #dv_dt = dv_dt.replace([np.inf,-np.inf],np.NaN)
     #dv_dt = dv_dt.fillna(dv_dt_mean)

     np.nan_to_num(dv_dt_mean)

     #dv_dt_conv = np.convolve(dv_dt,[0.5,0.5])

#     a = 1
#     windowsize = 5
#     b = (1/windowsize)*np.ones(windowsize)
#     #filtered_dvdt = scipy.signal.filtfilt(b,a,dv_dt_conv)
#     filtered_dvdt = scipy.signal.filtfilt(b,a,dv_dt)


     return dv_dt#filtered_dvdt


def _get_sbe_coef(idx=0):
    """
    Get SBE oxygen coefficients from raw data files.
    Defaults to using first station in ssscc.csv file.
    """
    ssscc_list = process_ctd.get_ssscc_list()
    station = ssscc_list[idx]
    hexfile = cfg.directory["raw"] + station + ".hex"
    xmlfile = cfg.directory["raw"] + station + ".XMLCON"

    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1,
                        'ranking': 5, 'column': 'CTDOXYVOLTS',
                        'sensor_info': rawConfig['Sensors'][i]}


    Soc = oxy_meta['sensor_info']['Soc'] # Oxygen slope
    Voff = oxy_meta['sensor_info']['offset'] # Sensor output offset voltage
    Tau20 = oxy_meta['sensor_info']['Tau20'] #Sensor time constant at 20 deg C and 1 Atm
    Tcorr = oxy_meta['sensor_info']['Tcor'] # Temperature correction
    E  = oxy_meta['sensor_info']['E'] #Compensation Coef for pressure effect on membrane permeability (Atkinson et al)
    A = oxy_meta['sensor_info']['A'] # Compensation Coef for temp effect on mem perm
    B = oxy_meta['sensor_info']['B'] # Compensation Coef for temp effect on mem perm
    C = oxy_meta['sensor_info']['C'] # Compensation Coef for temp effect on mem perm
    D = [oxy_meta['sensor_info']['D0'],oxy_meta['sensor_info']['D1'],oxy_meta['sensor_info']['D2']] # Press effect on time constant Usually not fitted for.

    coef = [Soc,Voff,Tau20,Tcorr,E]#,A,B,C,D]

    return coef

def calculate_weights(pressure):

    eps=1e-5
    wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
    wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]#[20,20,25,25,50,50,100,100,200,200,500,500]
    wgt = scipy.interpolate.interp1d(wrow1,wrow2)

    weights = wgt(pressure)

    return weights


def oxy_equation(X, Soc, Voffset, A, B, C, E, Tau20):
    # eq (3) in Uchida CTD manual
    cc=[1.92634e-4,-4.64803e-2]
    oxyvolts, pressure, temp, dvdt, os = X

    oxygen = (Soc * ((oxyvolts + Voffset + (Tau20 * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt))
                 * (1.0 + A * temp + B * np.power(temp,2) + C * np.power(temp,3) )
                 * os
                 * np.exp((E * pressure) / (temp + 273.15))))

    return oxygen


def _PMEL_oxy_eq(coefs,inputs,cc=[1.92634e-4,-4.64803e-2]):
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
    o2 = Soc * (oxyvolts + Voff + Tau20 * np.exp(cc[0] * pressure + cc[1] * (temp - 20)) * dvdt) * os \
            * np.exp(Tcorr * temp) \
            * np.exp((E * pressure) / (temp + 273.15))

    return o2

# TODO: optionally include other residual types
# (abstracted from PMEL code oxygen_cal_ml.m)
# unweighted L2: sum((ref - oxy)^2)  # if weighted fails
# unweighted L4: sum((ref - oxy)^4)  # unsure of use case
# unweighted L1: sum(abs(ref - oxy))  # very far from ideal
# anything else? genericize with integer "norm" function input?
def PMEL_oxy_weighted_residual(coefs,weights,inputs,refoxy):
    return np.sum((weights*(refoxy-_PMEL_oxy_eq(coefs, inputs))**2))/np.sum(weights**2)

def match_sigmas(btl_prs, btl_oxy, btl_tmp, btl_SA, ctd_os, ctd_prs, ctd_tmp, ctd_SA, ctd_oxyvolts, ctd_time):

    # Construct Dataframe from bottle and ctd values for merging
    btl_data = pd.DataFrame(data={
        "CTDPRS": btl_prs,
        "REFOXY": btl_oxy,
        "CTDTMP": btl_tmp,
        "SA": btl_SA,
    })
    time_data = pd.DataFrame(data={
        "CTDPRS": ctd_prs,
        "OS": ctd_os,
        "CTDTMP": ctd_tmp,
        "SA": ctd_SA,
        "CTDOXYVOLTS": ctd_oxyvolts,
        "CTDTIME": ctd_time,
    })
    time_data["dv_dt"] = calculate_dVdT(time_data["CTDOXYVOLTS"], time_data["CTDTIME"])

    # Merge DF
    merged_df = pd.DataFrame(
        columns=["CTDPRS", "CTDOXYVOLTS", "CTDTMP", "dv_dt", "OS"], dtype=float
    )
    merged_df["REFOXY"] = btl_data["REFOXY"].copy()

    # calculate sigma referenced to multiple depths
    for idx, p_ref in enumerate([0, 1000, 2000, 3000, 4000, 5000, 6000]):
        btl_data[f"sigma{idx}"] = (
            gsw.pot_rho_t_exact(
                btl_data["SA"],
                btl_data["CTDTMP"],
                btl_data["CTDPRS"],
                p_ref,
            )
            - 1000  # subtract 1000 to get potential density *anomaly*
        ) + 1e-8*np.random.standard_normal(btl_data["SA"].size)
        time_data[f"sigma{idx}"] = (
            gsw.pot_rho_t_exact(
                time_data["SA"],
                time_data["CTDTMP"],
                time_data["CTDPRS"],
                p_ref,
            )
            - 1000  # subtract 1000 to get potential density *anomaly*
        ) + 1e-8*np.random.standard_normal(time_data["SA"].size)
        rows = (btl_data["CTDPRS"] > (p_ref - 500)) & (btl_data["CTDPRS"] < (p_ref + 500))
        time_sigma_sorted = time_data[f"sigma{idx}"].sort_values().to_numpy()
        sigma_min = np.min([np.min(btl_data.loc[rows, f"sigma{idx}"]), np.min(time_sigma_sorted)])
        sigma_max = np.max([np.max(btl_data.loc[rows, f"sigma{idx}"]), np.max(time_sigma_sorted)])
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
    sbe_coef0 = _get_sbe_coef() # initial coefficient guess
    merged_df['CTDOXY'] = _PMEL_oxy_eq(sbe_coef0, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']))

    return merged_df


def sbe43_oxy_fit(merged_df, sbe_coef0=None, f_suffix=None):

    # Plot data to be fit together
    f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual{f_suffix}_prefit.pdf"
    ctd_plots._intermediate_residual_plot(
        merged_df['REFOXY'] - merged_df['CTDOXY'],
        merged_df["CTDPRS"],
        merged_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10,10)
    )

    # Create DF for good and questionable values
    bad_df = pd.DataFrame()
    good_df = pd.DataFrame()

    if sbe_coef0 is None:
        sbe_coef0 = _get_sbe_coef()  # load initial coefficient guess

    p0 = sbe_coef0[0], sbe_coef0[1], sbe_coef0[2], sbe_coef0[3], sbe_coef0[4]
    
    # Curve fit (weighted)
    weights = calculate_weights(merged_df['CTDPRS'])
    res = scipy.optimize.minimize(PMEL_oxy_weighted_residual,x0=p0,args=(weights, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']), merged_df['REFOXY']), bounds=[(None,None),(None,None),(0,None),(None,None),(None,None)])
    cfw_coefs = res.x
    merged_df['CTDOXY'] = _PMEL_oxy_eq(cfw_coefs, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']))        

    merged_df['residual'] = merged_df['REFOXY'] - merged_df['CTDOXY']
    stdres = np.std(merged_df['residual'])
    cutoff = stdres * 2.8

    thrown_values = merged_df[np.abs(merged_df['residual']) > cutoff]
    bad_df = pd.concat([bad_df, thrown_values])
    merged_df = merged_df[np.abs(merged_df['residual']) <= cutoff]

    while not thrown_values.empty: # runs as long as there are thrown_values

        p0 = cfw_coefs[0], cfw_coefs[1], cfw_coefs[2], cfw_coefs[3], cfw_coefs[4]
        weights = calculate_weights(merged_df['CTDPRS'])
        res = scipy.optimize.minimize(PMEL_oxy_weighted_residual,x0=p0,args=(weights, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']), merged_df['REFOXY']), bounds=[(None,None),(None,None),(0,None),(None,None),(None,None)])
        cfw_coefs = res.x
        merged_df['CTDOXY'] = _PMEL_oxy_eq(cfw_coefs, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']))

        merged_df['residual'] = merged_df['REFOXY'] - merged_df['CTDOXY']
        stdres = np.std(merged_df['residual'])
        cutoff = stdres * 2.8
        thrown_values = merged_df[np.abs(merged_df['residual']) > cutoff]
        print(len(thrown_values))
        print(p0)
        bad_df = pd.concat([bad_df, thrown_values])
        merged_df = merged_df[np.abs(merged_df['residual']) <= cutoff]

    # intermediate plots to diagnose data chunks goodness
    # TODO: implement into bokeh/flask dashboard
    if f_suffix is not None:
        f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual{f_suffix}.pdf"
        ctd_plots._intermediate_residual_plot(
            merged_df["residual"],
            merged_df["CTDPRS"],
            merged_df["SSSCC"],
            xlabel="CTDOXY Residual (umol/kg)",
            f_out=f_out,
            xlim=(-10,10)
        )

    # good_df = pd.concat([good_df, merged_df])
    merged_df['CTDOXY_FLAG_W'] = 2
    # bad_df = pd.concat([bad_df, bad_values])
    bad_df['CTDOXY_FLAG_W'] = 3
    df = pd.concat([merged_df,bad_df])

    # df['SSSCC_int'] = df['SSSCC_sbe43'].astype(int)
    # df.sort_values(by=['SSSCC_int','btl_fire_num'],ascending=[True,True],inplace=True)
    # df.drop('SSSCC_int', axis=1, inplace=True)

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
        btl_df[cfg.column["p_btl"]],
        btl_df[cfg.column["lon_btl"]],
        btl_df[cfg.column["lat_btl"]],
    )
    btl_df["CT"] = gsw.CT_from_t(
        btl_df["SA"],
        btl_df[cfg.column["t1_btl"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p_btl"]],
    )
    time_df["SA"] = gsw.SA_from_SP(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon_btl"]],
        time_df[cfg.column["lat_btl"]],
    )
    time_df["CT"] = gsw.CT_from_t(
        time_df["SA"],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
    )
    # calculate sigma
    btl_df["sigma_btl"] = sigma_from_CTD(
        btl_df[cfg.column["sal"]],
        btl_df[cfg.column["t1_btl"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p_btl"]],
        btl_df[cfg.column["lon_btl"]],
        btl_df[cfg.column["lat_btl"]],
    )
    time_df["sigma_ctd"] = sigma_from_CTD(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon_btl"]],
        time_df[cfg.column["lat_btl"]],
    )
    # Calculate oxygen solubility in µmol/kg
    btl_df["OS"] = gsw.O2sol(
        btl_df["SA"],
        btl_df["CT"],
        btl_df[cfg.column["p_btl"]],
        btl_df[cfg.column["lon_btl"]],
        btl_df[cfg.column["lat_btl"]],
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
    btl_df[cfg.column["oxy_btl"]] = calculate_bottle_oxygen(
        ssscc_list,
        btl_df["SSSCC"],
        btl_df["TITR_VOL"],
        btl_df["TITR_TEMP"],
        btl_df["FLASKNO"],
    )
    btl_df[cfg.column["oxy_btl"]] = oxy_ml_to_umolkg(
        btl_df[cfg.column["oxy_btl"]], btl_df["sigma_btl"]
    )
    btl_df["OXYGEN_FLAG_W"] = flag_winkler_oxygen(btl_df[cfg.column["oxy_btl"]])
    # Load manual OXYGEN flags
    if Path("data/oxygen/manual_oxy_flags.csv").exists():
        manual_flags = pd.read_csv(
            "data/oxygen/manual_oxy_flags.csv", dtype={"SSSCC": str}
        )
        for _, flags in manual_flags.iterrows():
            df_row = (btl_df["SSSCC"] == flags["SSSCC"]) & (
                btl_df["btl_fire_num"] == flags["Bottle"]
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
    # Plot all pre fit data
    f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual_all_prefit.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df['OXYGEN'] - btl_df['CTDOXY'],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10,10)
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
            print(ssscc + " skipped, all oxy data is NaN")
            continue
        sbe43_merged = match_sigmas(
            btl_data[cfg.column["p_btl"]],
            btl_data[cfg.column["oxy_btl"]],
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
        btl_df.loc[btl_df["SSSCC"] == ssscc, ["CTDOXYVOLTS","dv_dt","OS"]] = sbe43_merged[["CTDOXYVOLTS","dv_dt","OS"]]
        sbe43_merged["SSSCC"] = ssscc
        all_sbe43_merged = pd.concat([all_sbe43_merged, sbe43_merged])
        print(ssscc + " density matching done")

    # Only fit using OXYGEN flagged good (2)
    all_sbe43_merged = all_sbe43_merged[btl_df["OXYGEN_FLAG_W"] == 2].copy()

    # Fit ALL oxygen stations together to get initial coefficient guess
    (sbe_coef0, _) = sbe43_oxy_fit(all_sbe43_merged, f_suffix="_ox0")

    # Fit each cast individually
    for ssscc in ssscc_list:
        sbe_coef, sbe_df = sbe43_oxy_fit(
            all_sbe43_merged.loc[all_sbe43_merged["SSSCC"] == ssscc],
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
            print(ssscc + " missing oxy data, leaving nan values and flagging as 9")
            time_df.loc[time_df["SSSCC"] == ssscc, "CTDOXY_FLAG_W"] = 9
            time_df.loc[time_df["SSSCC"] == ssscc, "RINKO_FLAG_W"] = 9
            continue
        btl_rows = (btl_df["SSSCC"] == ssscc).values
        time_rows = (time_df["SSSCC"] == ssscc).values
        btl_df.loc[btl_rows, "CTDOXY"] = _PMEL_oxy_eq(
            sbe43_dict[ssscc],
            (
                btl_df.loc[btl_rows, cfg.column["oxyvolts"]],
                btl_df.loc[btl_rows, cfg.column["p_btl"]],
                btl_df.loc[btl_rows, cfg.column["t1_btl"]],
                btl_df.loc[btl_rows, "dv_dt"],
                btl_df.loc[btl_rows, "OS"],
            ),
        )
        print(ssscc + " btl data fitting done")
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
        print(ssscc + " time data fitting done")

    # TODO: flag oxy data here? compare w/ T/C routines

    # Plot all post fit data
    f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual_all_postfit.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df['OXYGEN'] - btl_df['CTDOXY'],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10,10)
    )

    # export fitting coefs
    sbe43_coefs = pd.DataFrame.from_dict(
        sbe43_dict, orient="index", columns=["Soc", "Voffset", "Tau20", "Tcorr", "E"]
    )
    sbe43_coefs.to_csv(cfg.directory["logs"] + "sbe43_coefs.csv")
    
    return True

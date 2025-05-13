"""
Oxygen functions for processing sensor data, unit conversions and derived variables.
"""
from collections import namedtuple

import gsw
import numpy as np

from ctdcal.processors.functions_ctd import _check_coefs, _check_volts


def sbe43(volts, p, t, c, coefs, lat=0.0, lon=0.0, decimals=4):
    # NOTE: lat/lon = 0 is not "acceptable" for GSW, come up with something else?
    """
    SBE equation for converting SBE43 engineering units to oxygen (ml/l).
    SensorID: 38

    Parameters
    ----------
    volts : array-like
        Raw voltage
    p : array-like
        Converted pressure (dbar)
    t : array-like
        Converted temperature (Celsius)
    c : array-like
        Converted conductivity (mS/cm)
    coefs : dict
        Dictionary of calibration coefficients (Soc, offset, Tau20, A, B, C, E)
    lat : array-like, optional
        Latitude (decimal degrees north)
    lon : array-like, optional
        Longitude (decimal degrees)

    Returns
    -------
    oxy_ml_l : array-like
        Converted oxygen (mL/L)
    """
    _check_coefs(coefs, ["Soc", "offset", "Tau20", "A", "B", "C", "E"])
    volts = _check_volts(volts)
    t_Kelvin = np.array(t) + 273.15

    SP = gsw.SP_from_C(c, t, p)
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    sigma0 = gsw.sigma0(SA, CT)
    o2sol = gsw.O2sol(SA, CT, p, lon, lat)  # umol/kg
    o2sol_ml_l = oxy_umolkg_to_ml(o2sol, sigma0)  # equation expects mL/L

    # NOTE: lat/lon always required to get o2sol (and need SA/CT for sigma0 anyway)
    # the above is equivalent to:
    # pt = gsw.pt0_from_t(SA, t, p)
    # o2sol = gsw.O2sol_SP_pt(s, pt)

    oxy_ml_l = (
        coefs["Soc"]
        * (volts + coefs["offset"])
        * (
            1.0
            + coefs["A"] * np.array(t)
            + coefs["B"] * np.power(t, 2)
            + coefs["C"] * np.power(t, 3)
        )
        * o2sol_ml_l
        * np.exp(coefs["E"] * np.array(p) / t_Kelvin)
    )
    return np.around(oxy_ml_l, decimals)


def sbe43_hysteresis_voltage(volts, p, coefs, sample_freq=24):
    """
    SBE equation for removing hysteresis from raw voltage values. This function must
    be run before the sbe43 conversion function above.

    Oxygen hysteresis can be corrected after conversion from volts to oxygen
    concentration, see oxy_fitting.hysteresis_correction()

    Parameters
    ----------
    volts : array-like
        Raw voltage
    p : array-like
        CTD pressure values (dbar)
    coefs : dict
        Dictionary of calibration coefficients (H1, H2, H3, offset)
    sample_freq : scalar, optional
        CTD sampling frequency (Hz)

    Returns
    -------
    volts_corrected : array-like
        Hysteresis-corrected voltage

    Notes
    -----
    The hysteresis algorithm is backward-looking so scan 0 must be skipped (as no
    information is available before the first scan).

    See Application Note 64-3 for more information.
    """
    _check_coefs(coefs, ["H1", "H2", "H3", "offset"])
    volts = _check_volts(volts)

    dt = 1 / sample_freq
    D = 1 + coefs["H1"] * (np.exp(np.array(p) / coefs["H2"]) - 1)
    C = np.exp(-1 * dt / coefs["H3"])

    oxy_volts = volts + coefs["offset"]
    oxy_volts_new = np.zeros(oxy_volts.shape)
    oxy_volts_new[0] = oxy_volts[0]
    for i in np.arange(1, len(oxy_volts)):
        oxy_volts_new[i] = (
            (oxy_volts[i] + (oxy_volts_new[i - 1] * C * D[i])) - (oxy_volts[i - 1] * C)
        ) / D[i]

    volts_corrected = oxy_volts_new - coefs["offset"]

    return volts_corrected


def hysteresis_correction(oxygen, pressure, H1=-0.033, H2=5000, H3=1450, freq=24):
    """
    Remove hysteresis effects from oxygen concentration values.

    Oxygen hysteresis can be corrected before conversion from volts to oxygen
    concentration, see equations_sbe.sbe43_hysteresis_voltage()

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

    # (PMEL does this calculation on binned data already so filtering is not the same)
    # a = 1
    # windowsize = 5
    # b = (1 / windowsize) * np.ones(windowsize)
    # filtered_dvdt = scipy.signal.filtfilt(b, a, dv_dt)

    return dV_dt  # filtered_dvdt


# Rinko functions
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
    """
    A weighted residual fit in a similar manner to that of the SBE43 method of oxy_fitting.py.
    """
    # (abstracted from PMEL code oxygen_cal_ml.m)
    # unweighted L2: sum((ref - oxy)^2)  # if weighted fails
    # unweighted L4: sum((ref - oxy)^4)  # unsure of use case
    # unweighted L1: sum(abs(ref - oxy))  # very far from ideal
    # anything else? genericize with integer "norm" function input?

    residuals = np.sum(
        (weights * (refoxy - _Uchida_DO_eq(coefs, inputs)) ** 2)
    ) / np.sum(weights ** 2)

    return residuals


def rinko_temperature(v, tmp_cal:RinkoTMPCal):
    """
    Calculate rinko temperature from voltage and calibration coeffieicnets.
    """
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

#   20240612 AJM commenting out method, possibly redundant
# def rinko_saturation(df, film="B", model="ARO-CAV", **kwargs):
#     """
#     Derive oxygen saturation for the RINKO. Unused as of 2021.
#     """
#     pass

def rinko_oxy_eq(press, temp, oxyvo, os, o2_cal:RinkoO2Cal):
    """
    Derive corrected RINKO dissolved oxygen.
    """
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


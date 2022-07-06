"""
A module for SBE conversion equations and related helper equations.

SBE internally references each sensor with an assigned SensorID for indexing. The
given sensor numbers are determined empirically from .XMLCON files in 2016-2017 and
not via an official document, and may change according to SBE wishes.
"""

import inspect
import logging

import gsw
import numpy as np

from ctdcal.oxy_fitting import oxy_umolkg_to_ml

log = logging.getLogger(__name__)


def _check_coefs(coefs_in, expected):
    """Compare function input coefs with expected"""
    missing_coefs = sorted(set(expected) - set(coefs_in))
    if missing_coefs != []:
        raise KeyError(f"Coefficient dictionary missing keys: {missing_coefs}")


def _check_freq(freq):
    """Convert to np.array, NaN out zeroes, convert to float if needed"""
    freq = np.array(freq)
    sensor = inspect.stack()[1].function  # name of function calling _check_freq

    if freq.dtype != float:  # can sometimes come in as object
        log.warning(f"Attempting to convert {freq.dtype} to float for {sensor}")
        freq = freq.astype(float)

    if 0 in freq:
        N_zeroes = (freq == 0).sum()
        log.warning(
            f"Found {N_zeroes} zero frequency readings in {sensor}, replacing with NaN"
        )
        freq[freq == 0] = np.nan

    return freq


def _check_volts(volts, v_min=0, v_max=5):
    """Convert to np.array, NaN out values outside of 0-5V, convert to float if needed"""
    volts = np.array(volts)
    sensor = inspect.stack()[1].function  # name of function calling _check_volts

    if volts.dtype != float:  # can sometimes come in as object
        log.warning(f"Attempting to convert {volts.dtype} to float for {sensor}")
        volts = volts.astype(float)

    if any(volts < v_min) or any(volts > v_max):
        log.warning(
            f"{sensor} has values outside of {v_min}-{v_max}V, replacing with NaN"
        )
        volts[volts < v_min] = np.nan
        volts[volts > v_max] = np.nan

    return volts


def sbe3(freq, coefs, decimals=4):
    """
    SBE equation for converting SBE3 frequency to temperature.
    SensorID: 55

    Parameters
    ----------
    freq : array-like
        Raw frequency (Hz)
    coefs : dict
        Dictionary of calibration coefficients (G, H, I, J, F0)

    Returns
    -------
    t_ITS90 : array-like
        Converted temperature (ITS-90)
    """
    _check_coefs(coefs, ["G", "H", "I", "J", "F0"])
    freq = _check_freq(freq)

    t_ITS90 = (
        1
        / (
            coefs["G"]
            + coefs["H"] * (np.log(coefs["F0"] / freq))
            + coefs["I"] * np.power((np.log(coefs["F0"] / freq)), 2)
            + coefs["J"] * np.power((np.log(coefs["F0"] / freq)), 3)
        )
        - 273.15
    )
    return np.around(t_ITS90, decimals)


def sbe4(freq, t, p, coefs, decimals=4):
    """
    SBE equation for converting SBE4 frequency to conductivity. This conversion
    is valid for both SBE4C (profiling) and SBE4M (mooring).
    SensorID: 3

    Parameters
    ----------
    freq : array-like
        Raw frequency (Hz)
    t : array-like
        Converted temperature (ITS-90 degrees C)
    p : array-like
        Converted pressure (dbar)
    coefs : dict
        Dictionary of calibration coefficients (G, H, I, J, CPcor, CTcor)

    Returns
    -------
    c_mS_cm : array-like
        Converted conductivity (mS/cm)
    """
    _check_coefs(coefs, ["G", "H", "I", "J", "CPcor", "CTcor"])
    freq_kHz = _check_freq(freq) * 1e-3  # equation expects kHz

    c_S_m = (
        coefs["G"]
        + coefs["H"] * np.power(freq_kHz, 2)
        + coefs["I"] * np.power(freq_kHz, 3)
        + coefs["J"] * np.power(freq_kHz, 4)
    ) / (10 * (1 + coefs["CTcor"] * np.array(t) + coefs["CPcor"] * np.array(p)))
    c_mS_cm = c_S_m * 10  # S/m to mS/cm

    return np.around(c_mS_cm, decimals)


def sbe9(freq, t_probe, coefs, decimals=4):
    """
    SBE/STS(?) equation for converting SBE9 frequency to pressure.
    SensorID: 45

    Parameters
    ----------
    freq : array-like
        Raw frequency (Hz)
    t_probe : array-like
        Raw integer measurement from the Digiquartz temperature probe
    coefs : dict
        Dictionary of calibration coefficients
        (T1, T2, T3, T4, T5, C1, C2, C3, D1, D2, AD590M, AD590B)

    Returns
    -------
    p_dbar : array-like
        Converted pressure (dbar)
    """
    _check_coefs(
        coefs,
        (
            ["T1", "T2", "T3", "T4", "T5"]
            + ["C1", "C2", "C3"]
            + ["D1", "D2"]
            + ["AD590M", "AD590B"]
        ),
    )
    freq_MHz = _check_freq(freq) * 1e-6  # equation expects MHz
    t_probe = (coefs["AD590M"] * np.array(t_probe).astype(int)) + coefs["AD590B"]

    T0 = (
        coefs["T1"]
        + coefs["T2"] * t_probe
        + coefs["T3"] * np.power(t_probe, 2)
        + coefs["T4"] * np.power(t_probe, 3)
    )
    w = 1 - T0 * T0 * freq_MHz * freq_MHz
    p_dbar = 0.6894759 * (
        (coefs["C1"] + coefs["C2"] * t_probe + coefs["C3"] * t_probe * t_probe)
        * w
        * (1 - (coefs["D1"] + coefs["D2"] * t_probe) * w)
        - 14.7
    )
    return np.around(p_dbar, decimals)


def sbe_altimeter(volts, coefs, decimals=1):
    """
    SBE equation for converting altimeter voltages to meters. This conversion
    is valid for altimeters integrated with any Sea-Bird CTD (e.g. 9+, 19, 25).
    Sensor ID: 0

    Parameters
    ----------
    volts : array-like
        Raw voltages
    coefs : dict
        Dictionary of calibration coefficients (ScaleFactor, Offset)

    Returns
    -------
    bottom_distance : array-like
        Distance from the altimeter to an object below it (meters)

    Notes
    -----
    Equation provdided by SBE in Application Note 95, page 1.

    While the SBE documentation refers to a Teledyne Benthos or Valeport altimeter,
    the equation works for all altimeters typically found in the wild.
    """
    _check_coefs(coefs, ["ScaleFactor", "Offset"])
    volts = _check_volts(volts)

    bottom_distance = (300 * volts / coefs["ScaleFactor"]) + coefs["Offset"]

    return np.around(bottom_distance, decimals)


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
    # TODO: is there any reason for this to output mL/L? if oxygen eq uses o2sol
    # in umol/kg, result is in umol/kg... which is what we use at the end anyway?

    _check_coefs(coefs, ["Soc", "offset", "Tau20", "A", "B", "C", "E"])
    volts = _check_volts(volts)
    t_Kelvin = np.array(t) + 273.15

    SP = gsw.SP_from_C(c, t, p)
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    sigma0 = gsw.sigma0(SA, CT)
    o2sol = gsw.O2sol(SA, CT, p, lon, lat)  # umol/kg
    o2sol_ml_l = oxy_umolkg_to_ml(o2sol, sigma0)  # equation expects mL/L (see TODO)

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
    # TODO: vectorize (if possible), will probably require matrix inversion
    # TODO: any NaNs will be propagated through entire timeseries: feature or bug?

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


def wetlabs_eco_fl(volts, coefs, decimals=4):
    """
    SBE equation for converting ECO-FL fluorometer voltage to concentration.
    SensorID: 20

    Parameters
    ----------
    volts : array-like
        Raw voltage
    coefs : dict
        Dictionary of calibration coefficients (ScaleFactor, DarkOutput/Vblank)

    Returns
    -------
    chl : array-like
        Converted chlorophyll concentration

    Notes
    -----
    Chlorophyll units depend on scale factor (e.g. ug/L-volt, ug/L-counts, ppb/volts),
    see Application Note 62 for more information.
    """
    volts = _check_volts(volts)

    if "DarkOutput" in coefs.keys():
        chl = coefs["ScaleFactor"] * (volts - coefs["DarkOutput"])
    elif "Vblank" in coefs.keys():  # from older calibration sheets
        chl = coefs["ScaleFactor"] * (volts - coefs["Vblank"])
    else:
        log.warning(
            "No dark cast coefficient ('DarkOutput' or 'Vblank'), returning voltage."
        )
        chl = volts

    return np.around(chl, decimals)


def wetlabs_cstar(volts, coefs, decimals=4):
    """
    SBE equation for converting C-Star transmissometer voltage to light transmission.
    SensorID: 71

    Parameters
    ----------
    volts : array-like
        Raw voltage
    coefs : dict
        Dictionary of calibration coefficients (M, B, PathLength)

    Returns
    -------
    xmiss : array-like
        Light transmission [%]
    c : array-like
        Beam attenuation coefficient

    Notes
    -----
    M and B can be recalculated in the field by measuring voltage in air (A1) and
    voltage with the path blocked (Y1):
        M = (Tw / (W0 - Y0)) * (A0 - Y0) / (A1 - Y1)
        B = -M * Y1
    where A0, Y0, and W0 are factory values.
    For transmission relative to water, set Tw = 100%.

    See Application Note 91 for more information.
    """
    _check_coefs(coefs, ["M", "B", "PathLength"])
    volts = _check_volts(volts)
    xmiss = (coefs["M"] * volts) + coefs["B"]  # xmiss as a percentage
    c = -(1 / coefs["PathLength"]) * np.log(xmiss * 100)  # needs xmiss as a decimal

    return np.around(xmiss, decimals), np.around(c, decimals)


def seapoint_fluor(volts, coefs, decimals=6):
    """
    Raw voltage supplied from fluorometer right now, after looking at xmlcon.
    The method will do nothing but spit out the exact values that came in.
    SensorID: 11

    Parameters
    ----------
    volts : array-like
        Raw voltage
    coefs : dict
        Dictionary of calibration coefficients (GainSetting, Offset)

    Returns
    -------
    fluoro : array-like
        Raw voltage

    Notes
    -----
    According to .xmlcon, GainSetting "is an array index, not the actual gain setting."
    """
    _check_coefs(coefs, ["GainSetting", "Offset"])
    # TODO: actual calibration/conversion/something?
    # TODO: move this to different module? edge case since it's the only Seapoint sensor
    volts = np.array(volts)
    fluoro = np.around(volts, decimals)

    return fluoro

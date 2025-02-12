"""
Conductivity, temperature and pressure functions for processing sensor data, unit
conversions and derived variables.
"""

import inspect
import logging

import numpy as np


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

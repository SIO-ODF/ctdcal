"""
Auxiliary sensor functions to convert data from original to engineering units,
engineering to science units, and to produce derived variables.
"""
import logging

import numpy as np

from ctdcal.processors.functions_ctd import _check_coefs, _check_volts


log = logging.getLogger(__name__)


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


# Fluorometers
#
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
    volts = np.array(volts)
    fluoro = np.around(volts, decimals)

    return fluoro


def sbe_flntu_chl(volts, coefs, decimals=4):
    """
    SBE equation for converting SeaBird fluorometer and nepholometric turbidity
    combo sensor's chlorophyll and CDOM fluorometer.
    SensorID: 19
    Example coefs: {'ScaleFactor':11, 'DarkVoltage':0.079}

    Paramters
    ----------
    volts : array-like
        Raw voltage
    coefs : dict
        Dictionary of calibration coefficients (ScaleFactor, Vblank)

    Returns
    -------
    chl : array-like
        Converted chlorophyll concentration in μg/l

    Notes:
    ------
    ScaleFactor is usually given in volts, whereas Vblank are in μg/l/V in
    analogue calculations.
    """
    _check_coefs(coefs, ["ScaleFactor", "Vblank"])
    volts = _check_volts(volts)
    chl = np.around((volts - coefs["Vblank"]) * coefs["ScaleFactor"], decimals)

    return chl


def sbe_flntu_ntu(volts, coefs, decimals=4):
    """
    SBE equation for converting SeaBird fluorometer and nepholometric turbidity
    combo sensor's turbidity channel.
    SensorID: 67
    Example coefs: {'ScaleFactor':5, 'DarkVoltage':0.05}

    Paramters
    ----------
    volts : array-like
        Raw voltage
    coefs : dict
        Dictionary of calibration coefficients (ScaleFactor, DarkVoltage)

    Returns
    -------
    turb : array-like
        Converted turbidity (units expressed in NTU)

    Notes:
    ------
    ScaleFactor is usually given in volts, whereas DarkVoltage are in NTU/V in
    analogue calculations.
    Future release may factor this into sbe_flntu_chl (same equation)
    """

    _check_coefs(coefs, ["ScaleFactor", "DarkVoltage"])
    volts = _check_volts(volts)
    turb = np.around((volts - coefs["DarkVoltage"]) * coefs["ScaleFactor"], decimals)

    return turb

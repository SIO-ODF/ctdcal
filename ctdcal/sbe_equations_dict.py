"""A module for SBE conversion equations and related helper equations.

Eventual goal is to convert all outputs to numpy arrays to make compatible with
gsw libraries, and remove written wrappers.

SBE internally references each sensor with an assigned SensorID for indexing. The
given sensor numbers are determined empirically from .XMLCON files in 2016-2017 and
not via an official document, and may change according to SBE wishes.

"""

'''code_pruning: this module really needs a new name, and all functions inside it as well.'''
'''code_pruning: Also all functions should be checked to work with vectorized versions and possibly remove single datapoint version.'''

import gsw
import numpy as np
from ctdcal.oxy_fitting import oxy_umolkg_to_ml


def temp_its90(coefs, freq):
    """
    SBE equation for converting engineering units to Celsius according to ITS-90.
    SensorID: 55

    Parameters
    ----------
    coefs : dict
        Dictionary of calibration coefficients (G, H, I, J, F0)
    freq : array-like
        Raw frequency (Hz)

    Returns
    -------
    t_ITS90 : array-like
        Converted temperature (ITS-90)
    """
    freq = np.array(freq)

    if freq.dtype != float:  # can sometimes come in as object
        freq = freq.astype(float)
        # TODO: (logger) e.g. "warning: converting {dtype} to float" or something

    if 0 in freq: # TODO: is this actually needed? what about for other conversion funcs?
        freq[freq == 0] = np.nan  # nan out zero frequencies
        # TODO: (logger) e.g. "warning: converting zero frequency found in {} to nan"

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
    return np.around(t_ITS90, 4)


def sbe9(coefs, freq, t_probe):
    """
    SBE/STS(?) equation for converting SBE9 frequency to pressure.
    SensorID: 45

    Parameters
    ----------
    coefs : dict
        Dictionary of calibration coefficients
        (T1, T2, T3, T4, T5, C1, C2, C3, D1, D2, AD590M, AD590B)
    freq : array-like
        Raw frequency (Hz)
    t_probe : array-like
        Raw integer measurement from the Digiquartz temperature probe

    Returns
    -------
    p_dbar : array-like
        Converted pressure (dbar)
    """
    freq = np.array(freq)
    t_probe = np.array(t_probe).astype(int)
    freq_MHz = freq * 1e-6  # equation expects MHz
    t_probe = (coefs["AD590M"] * t_probe) + coefs["AD590B"]
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
    return np.around(p_dbar, 4)


def sbe4c(coefs, freq, t, p):
    """
    SBE equation for converting SBE4C frequency to conductivity.
    SensorID: 3

    Parameters
    ----------
    coefs : dict
        Dictionary of calibration coefficients (G, H, I, J, CPcor, CTcor)
    freq : array-like
        Raw frequency (Hz)
    t : array-like
        Converted temperature (ITS-90 degrees C)
    p : array-like
        Converted pressure (dbar)

    Returns
    -------
    c_mS_cm : array-like
        Converted conductivity (mS/cm)
    """
    freq_kHz = freq * 1e-3  # equation expects kHz
    c_S_m = (
        coefs["G"]
        + coefs["H"] * np.power(freq_kHz, 2)
        + coefs["I"] * np.power(freq_kHz, 3)
        + coefs["J"] * np.power(freq_kHz, 4)
    ) / (10 * (1 + coefs["CTcor"] * t + coefs["CPcor"] * p))
    c_mS_cm = c_S_m * 10  # S/m to mS/cm

    return np.around(c_mS_cm, 5)


def sbe43(coefs, p, t, c, V, lat=0., lon=0.):
    """
    SBE equation for converting SBE43 engineering units to oxygen (ml/l).
    SensorID: 38

    Parameters
    ----------
    coefs : dict
        Dictionary of calibration coefficients (Soc, Voffset, Tau20, A, B, C, E)
    p : array-like
        Converted pressure (dbar)
    t : array-like
        Converted temperature (Celsius)
    c : array-like
        Converted conductivity (mS/cm)
    V : array-like
        Raw voltage
    lat : array-like, optional
        Latitude (decimal degrees north)
    lon : array-like, optional
        Longitude (decimal degrees)

    Returns
    -------
    oxygen : array-like
        Converted oxygen (mL/L)
    """
    # TODO: is there any reason for this to output mL/L? if oxygen eq uses o2sol
    # in umol/kg, result is in umol/kg... which is what we use at the end anyway?

    t_Kelvin = t + 273.15

    SP = gsw.SP_from_C(c, t, p)
    SA = gsw.SA_from_SP(SP, p, lat, lon)
    CT = gsw.CT_from_t(SA, t, p)
    sigma0 = gsw.sigma0(SA, CT)
    o2sol = gsw.O2sol(SA, CT, p, lon, lat)  # umol/kg
    o2sol_ml_l = oxy_umolkg_to_ml(o2sol, sigma0)  # equation expects mL/L (see TODO)

    # NOTE: lat/lon always required to get o2sol (and need SA/CT for sigma0 anyway)
    # the above is equivalent to:
    # pt = gsw.pt0_from_t(SA, t, p)
    # o2sol = gsw.O2sol_SP_pt(s, pt)

    oxygen = (
        coefs["Soc"]
        * (V + coefs["offset"])
        * (
            1.0
            + coefs["A"] * t
            + coefs["B"] * np.power(t, 2)
            + coefs["C"] * np.power(t, 3)
        )
        * o2sol_ml_l
        * np.exp(coefs["E"] * p / t_Kelvin)
    )
    return np.around(oxygen, 4)


def oxy_hysteresis_voltage(calib, voltage, scan_window=48):
    '''SBE equation for computing hysteresis from raw voltage.
    Must be run before oxy_dict.
    Because of looking backwards, must skip i = 0.

    Input:
    calib: a dict holding H1, H2, H3, Voffset
    voltage: a sequence of engineering voltages
    scan_window: an int for scans to skip between, default 48 scans at 24Hz OR 2 seconds
    '''
    output = []

    for i, x in enumerate(voltage):
        if i == 0:
            continue

        D = 1 + calib['H1']*(exp(P(i)/calib['H2']) - 1)
        C = exp(-1 * ())

    return output

def sbe43_hysteresis_voltage(calib, voltage):
    '''NOT TESTED NOT FINISHED

    SBE equation for computing hysteresis from raw voltage in relation to SBE43.
    Must be run before oxy_dict.
    Because of looking backwards, must skip i = 0.

    Parameters
    ----------
    calib : dict
        calib is a dict holding H1, H2, H3, Voffset
    voltage : array-like
        voltage is a sequence of engineering voltages

    Returns
    -------
    output : array-like
        output is the
    '''
    output = []

    for i, x in enumerate(voltage):
        if i == 0:
            continue

        D = 1 + calib['H1']*(exp(P(i)/calib['H2']) - 1)
        C = exp(-1 * ())

    return output


def wetlabs_flrtd_chl_dict(calib, counts):
    """Wetlabs

    UNFINISHED

    """
    chl = calib['scale_factor'] * (output - calib['darkcounts'])
    return chl


def wetlabs_transmissometer_cstar_dict(calib, signal):
    """Wetlabs CStar Transmissiometer.
    Equation from calib sheet for S/N#: CST-479-DR, Date: October 31, 2014
    SensorID: 71

    Inputs:
    calib is a dictionary of constants/coefficients
        calib['dark'] = voltage when beam is blocked. V_d on calib sheet
        calib['air'] = voltage with clear beam path in air. V_air on calib sheet
        calib['reference'] = voltage with beam path in clean water. V_ref on calib sheet
    signal: dict/single of signal voltage

    Relationship of transmittance (Tr) to beam attenuation coefficient (c), and pathlength (x, in meters): Tr = e^-ex
    beam attenuation coefficient is determined as: c = -1/x * ln (Tr)

    """

    #array mode
    try:
        tx = []
        for signal_x in signal:
            temp = (signal_x - calib['dark'])/(calib['reference'] - calib['dark'])
            tx.append(temp)
    #single mode
    except:
        tx = (signal - calib['dark'])/(calib['reference'] - calib['dark'])
    return tx


def benthos_psa916_dict(calib, signal):
    """Equation for determining altitude from a Benthos PSA-916 altimeter.
    Equation provided by SBE as AN95, or here: http://www.seabird.com/document/an95-setting-teledyne-benthos-altimeter-sea-bird-profiling-ctd
    Equation stated as: altimeter height = [300 * voltage / scale factor] + offset
    SensorID: 0

    Inputs:
    calib is a dictionary of coefficients
        calib['ScaleFactor']: scaling factor to be applied
        calib['Offset']: offset to be applied

    signal: signal voltage
    """

    #array mode
    try:
        altitude = []
        for signal_x in signal:
            temp = (300 * signal_x / calib['ScaleFactor']) + calib['Offset']
            altitude.append(temp)
    #single mode
    except:
        altitude = (300 * signal / calib['ScaleFactor']) + calib['Offset']
    return altitude


def fluoro_seapoint_dict(calib, signal):
    """
    Raw voltage supplied from fluorometer right now, after looking at xmlcon.
    The method will do nothing but spit out the exact values that came in.
    SensorID: 11

    Inputs:
    calib is a dictionary of coefficients(?)
        GainSetting: the gain applied. according to xmlcon,
            "<!-- The following is an array index, not the actual gain setting. -->"
        Offset: offset applied

    signal: signal voltage
    """
    try:
        fluoro = []
        for signal_x in signal:
            temp = signal_x
            fluoro.append(round(temp,6))
    #single mode
    except:
        fluoro = round(signal,6)
    return fluoro

def altimeter_voltage(calib, volts):
    """
    SBE Equation for converting voltages from an altimeter to meters.

    While the SBE documentation refers to a Teledyne Benthos altimeter, the equation
    works for all altimeters typically found in the wild.

    Sensor ID: 0

    Parameters
    ----------
    calib : dict
        calib is a dict holding SerialNumber, CalibrationDate, ScaleFactor, and Offset
    volts : array-like
        volts is the voltage from the altimeter

    Returns
    -------
    bottom_distance : array-like
        bottom_distance is the distance from the altimeter to an object below it, in meters.

    Equation provided by SBE as AN95, or here:
    http://www.seabird.com/document/an95-setting-teledyne-benthos-altimeter-sea-bird-profiling-ctd
    Equation stated as: altimeter height = [300 * voltage / scale factor] + offset
    """

    # The array might come in as type=object, which throws AttributeError. Maybe this should be in try/except?
    volts = volts.astype(float)

    bottom_distance = np.around((
                                (300 * volts / calib['ScaleFactor'])
                                + calib['Offset']
                                ),1)
    return bottom_distance

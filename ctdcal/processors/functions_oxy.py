"""
Oxygen functions for processing sensor data, unit conversions and derived variables.
"""
import gsw
import numpy as np

from ctdcal.oxy_fitting import oxy_umolkg_to_ml
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

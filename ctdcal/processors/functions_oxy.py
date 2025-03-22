"""
Oxygen functions for processing sensor data, unit conversions and derived variables.
"""
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

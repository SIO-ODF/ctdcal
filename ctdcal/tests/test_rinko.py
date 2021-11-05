import numpy as np
import pytest

from ctdcal import rinko


def test_salinity_correction():

    DO_c = np.array([100, 100, 0, 0, 0, np.nan])
    T = np.array([1, 1, 0, 0, np.nan, 0])
    C = np.array([1, 1, 0, np.nan, 0, 0])

    # check values are converted correctly
    corr = rinko.salinity_correction(DO_c, T, C)
    assert corr.shape == (6,)
    assert corr[0] == corr[1]
    assert corr[2] == 0
    assert all(np.isnan(corr[3:]))


def test_Uchida_DO_eq():

    coefs = (-1, 1.5, 2, 4.5, 5, 6, 7)
    inputs = (2.5, 1000, 2, 36.6, 90)

    # check behavior with scalars
    data = rinko._Uchida_DO_eq(coefs, inputs)
    assert np.round(data, 4) == -7.0004

    # check behavior with vectors
    vector_inputs = (np.full(20, input) for input in inputs)
    data = rinko._Uchida_DO_eq(coefs, vector_inputs)
    assert len(data) == 20
    assert all(np.round(data, 4) == -7.0004)

    # error if wrong number of coefs or inputs
    # note: need escape char since () is special in regex match
    with pytest.raises(ValueError, match=r"coefficients \(expected 7, got 6\)"):
        rinko._Uchida_DO_eq(coefs[:-1], inputs)
    with pytest.raises(ValueError, match=r"inputs \(expected 5, got 4\)"):
        rinko._Uchida_DO_eq(coefs, inputs[:-1])


### this section is specific to Rinko equations from JFE Advantech (not used by ODF)
def test_rinko_temperature():
    volts = np.array([2.0] * 9 + [np.nan])
    t_cal = rinko.T_cal(0.1, 0.2, 0.3, 0.4)

    # check values are converted correctly
    T = rinko.rinko_temperature(volts, t_cal)
    assert len(T) == 10
    assert all(T[:-1] == 4.9)
    assert np.isnan(T[-1])

    # check same behavior if calibration is a normal tuple
    T = rinko.rinko_temperature(volts, (0.1, 0.2, 0.3, 0.4))
    assert len(T) == 10
    assert all(T[:-1] == 4.9)
    assert np.isnan(T[-1])


def test_rinko_p_prime():
    volts = np.array([2.0] * 9 + [np.nan])
    t = np.full(volts.shape, 25)
    o2_cal = rinko.O2_cal(1, 2, 3, 4, 5, 6, 7, 8)

    # check values are converted correctly
    P_prime = rinko.rinko_p_prime(volts, t, o2_cal)
    assert len(P_prime) == 10
    assert all(P_prime[:-1] == 1.4)
    assert np.isnan(P_prime[-1])

    # check same behavior if calibration is a normal tuple
    P_prime = rinko.rinko_p_prime(volts, t, (1, 2, 3, 4, 5, 6, 7, 8))
    assert len(P_prime) == 10
    assert all(P_prime[:-1] == 1.4)
    assert np.isnan(P_prime[-1])


def test_rinko_saturation():
    P_prime = np.array([1.4] * 9 + [np.nan])
    o2_cal = rinko.O2_cal(1, 2, 3, 4, 5, 6, 7, 8)

    # check values are converted correctly
    P = rinko.rinko_saturation(P_prime, o2_cal)
    assert len(P) == 10
    assert all(P[:-1] == 18.2)
    assert np.isnan(P_prime[-1])

    # check same behavior if calibration is a normal tuple
    P = rinko.rinko_saturation(P_prime, (1, 2, 3, 4, 5, 6, 7, 8))
    assert len(P) == 10
    assert all(P[:-1] == 18.2)
    assert np.isnan(P_prime[-1])


def test_pressure_correction(caplog):
    P = np.array([18.2] * 9 + [np.nan])
    pressure = np.full(P.shape, 25)
    o2_cal = rinko.O2_cal(1, 2, 3, 4, 5, 6, 7, 8)

    # check values are converted correctly
    P_d = rinko.rinko_pressure_correction(P, pressure, o2_cal)
    assert len(P_d) == 10
    assert all(P_d[:-1] == 2293.2)
    assert np.isnan(P_d[-1])

    # check same behavior if calibration is a normal tuple
    P_d = rinko.rinko_pressure_correction(P, pressure, (1, 2, 3, 4, 5, 6, 7, 8))
    assert len(P_d) == 10
    assert all(P_d[:-1] == 2293.2)
    assert np.isnan(P_d[-1])

    # error (but continue) if pressures seem to be in dbar (instead of MPa)
    P_d_err = rinko.rinko_pressure_correction(P, np.full(P.shape, 80), o2_cal)
    assert "Found pressures >60 MPa" in caplog.messages[0]
    assert len(P_d_err) == 10
    assert all(P_d_err[:-1] == 7298.2)
    assert np.isnan(P_d_err[-1])


def test_rinko_DO(caplog):
    oxy_volts = np.array([2.0] * 9 + [np.nan])
    p = np.full(oxy_volts.shape, 2500)
    t = np.full(oxy_volts.shape, 25)
    OS = np.full(oxy_volts.shape, 100)
    o2_cal = rinko.O2_cal(1, 2, 3, 4, 5, 6, 7, 8)

    # check values are converted correctly
    DO = rinko.rinko_DO(p, t, oxy_volts, OS, o2_cal)
    assert len(DO) == 10
    assert all(DO[:-1] == 2293.2)
    assert np.isnan(DO[-1])

    # check same behavior if calibration is normal tuple
    DO = rinko.rinko_DO(p, t, oxy_volts, OS, (1, 2, 3, 4, 5, 6, 7, 8))
    assert len(DO) == 10
    assert all(DO[:-1] == 2293.2)
    assert np.isnan(DO[-1])

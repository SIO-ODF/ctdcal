import numpy as np
import pandas as pd
import pytest

from ctdcal import equations_sbe as eqs


def make_coefs(*args, value=0):
    return {x: value for x in args[0]}


def test_sbe3(caplog):
    freq = 99 * [0] + [1]
    coefs = ["G", "H", "I", "J", "F0"]

    # check values are converted correctly
    cnv = eqs.sbe3(freq, make_coefs(coefs, value=1))  # avoid divide by zero
    assert "int" in caplog.records[0].message  # check logging output
    assert "sbe3" in caplog.records[1].message  # check logging output
    assert all(np.isnan(cnv[:-1]))
    assert cnv[-1] == -272.15  # t_ITS90 = (1 / 1) - 273.15
    assert cnv.dtype == float

    # check there are no warnings if freq is float with no zeroes
    no_err = eqs.sbe3(np.ones(len(freq)), make_coefs(coefs, value=1))
    assert len(caplog.records) == 2
    assert all(no_err == -272.15)
    assert no_err.dtype == float

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="F0"):
        eqs.sbe3(freq, make_coefs(coefs[:-1]))


def test_sbe4(caplog):
    freq = 99 * [0] + [1]
    t = len(freq) * [0]
    p = len(freq) * [0]
    coefs = ["G", "H", "I", "J", "CPcor", "CTcor"]

    # check values are converted correctly
    cnv = eqs.sbe4(freq, t, p, make_coefs(coefs, value=1))  # avoid divide by zero
    assert "int" in caplog.records[0].message  # check logging output
    assert "sbe4" in caplog.records[1].message  # check logging output
    assert all(np.isnan(cnv[:-1]))
    assert cnv[-1] == 1
    assert cnv.dtype == float

    # check there are no warnings if freq is float with no zeroes
    no_err = eqs.sbe4(np.ones(len(freq)), t, p, make_coefs(coefs, value=1))
    assert len(caplog.records) == 2
    assert all(no_err == 1)
    assert no_err.dtype == float

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="CTcor"):
        eqs.sbe4(freq, t, p, make_coefs(coefs[:-1]))


def test_sbe9(caplog):
    freq = 99 * [0] + [1]
    t = len(freq) * [0]
    coefs = (
        ["T1", "T2", "T3", "T4", "T5"]
        + ["C1", "C2", "C3"]
        + ["D1", "D2"]
        + ["AD590M", "AD590B"]
    )

    # check values are converted correctly
    cnv = eqs.sbe9(freq, t, make_coefs(coefs))
    assert "int" in caplog.records[0].message  # check logging output
    assert "sbe9" in caplog.records[1].message  # check logging output
    assert all(np.isnan(cnv[:-1]))
    assert cnv[-1] == -10.1353
    assert cnv.dtype == float

    # check there are no warnings if freq is float with no zeroes
    no_err = eqs.sbe9(np.ones(len(freq)), t, make_coefs(coefs))
    assert len(caplog.records) == 2
    assert all(no_err == -10.1353)
    assert no_err.dtype == float

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="AD590B"):
        eqs.sbe9(freq, t, make_coefs(coefs[:-1]))


def test_sbe_altimeter():
    volts = 99 * [0] + [1]
    coefs = ["ScaleFactor", "Offset"]

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="Offset"):
        eqs.sbe_altimeter(volts, make_coefs(coefs[:-1]))


def test_sbe43():
    freq = 99 * [0] + [1]
    t = len(freq) * [0]
    p = len(freq) * [0]
    c = len(freq) * [0]
    coefs = ["Soc", "Voffset", "Tau20", "A", "B", "C", "E"]

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="E"):
        eqs.sbe43(freq, p, t, c, make_coefs(coefs[:-1]))


def test_sbe43_hysteresis_voltage():
    volts = 99 * [0] + [1]
    p = len(volts) * [0]
    coefs = ["H1", "H2", "H3", "offset"]

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="offset"):
        eqs.sbe43_hysteresis_voltage(volts, p, make_coefs(coefs[:-1]))


# def test_wetlabs_eco_fl():
#     eqs.wetlabs_eco_fl()


def test_wetlabs_cstar():
    volts = 99 * [0] + [1]
    coefs = ["M", "B", "PathLength"]

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="PathLength"):
        eqs.wetlabs_cstar(volts, make_coefs(coefs[:-1]))


def test_seapoint_fluor():
    volts = 99 * [0] + [1]
    coefs = ["GainSetting", "Offset"]

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="Offset"):
        eqs.seapoint_fluor(volts, make_coefs(coefs[:-1]))

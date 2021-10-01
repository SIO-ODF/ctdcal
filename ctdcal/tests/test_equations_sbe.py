import numpy as np
import pandas as pd
import pytest

from ctdcal import equations_sbe as eqs


def make_coefs(*args, value=0):
    return {x: value for x in args[0]}


def test_sbe3():
    freq = 99 * [0] + [1]
    coefs = ["G", "H", "I", "J", "F0"]

    cnv = eqs.sbe3(freq, make_coefs(coefs, value=1))

    # check values are converted correctly
    assert all(np.isnan(cnv[:-1]))
    assert cnv[-1] == -272.15  # t_ITS90 = (1 / 1) - 273.15
    assert cnv.dtype == np.float64

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="F0"):
        eqs.sbe3(freq, make_coefs(coefs[:-1]))


def test_sbe4():
    freq = 99 * [0] + [1]
    t = len(freq) * [0]
    p = len(freq) * [0]
    coefs = ["G", "H", "I", "J", "CPcor", "CTcor"]

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="CTcor"):
        eqs.sbe4(freq, t, p, make_coefs(coefs[:-1]))


def test_sbe9():
    freq = 99 * [0] + [1]
    t = len(freq) * [0]
    coefs = (
        ["T1", "T2", "T3", "T4", "T5"]
        + ["C1", "C2", "C3"]
        + ["D1", "D2"]
        + ["AD590M", "AD590B"]
    )

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

import numpy as np
import pytest

from ctdcal.processors import functions_ctd as ctd
from ctdcal.processors import functions_aux as aux
from ctdcal.processors import functions_oxy as oxy


def make_coefs(*args, value=0):
    return {x: value for x in args[0]}


def test_sbe3(caplog):
    freq = 99 * [0] + [1]
    coefs = ["G", "H", "I", "J", "F0"]

    # check values are converted correctly
    cnv = ctd.sbe3(freq, make_coefs(coefs, value=1))  # avoid divide by zero
    assert "int" in caplog.records[0].message  # check logging output
    assert "sbe3" in caplog.records[1].message  # check logging output
    assert all(np.isnan(cnv[:-1]))
    assert cnv[-1] == -272.15  # t_ITS90 = (1 / 1) - 273.15
    assert cnv.dtype == float

    # check there are no warnings if freq is float with no zeroes
    no_warn = ctd.sbe3(np.ones(len(freq)), make_coefs(coefs, value=1))
    assert len(caplog.records) == 2
    assert all(no_warn == -272.15)
    assert no_warn.dtype == float

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="F0"):
        ctd.sbe3(freq, make_coefs(coefs[:-1]))


def test_sbe4(caplog):
    freq = 99 * [0] + [1]
    t = len(freq) * [0]
    p = len(freq) * [0]
    coefs = ["G", "H", "I", "J", "CPcor", "CTcor"]

    # check values are converted correctly
    cnv = ctd.sbe4(freq, t, p, make_coefs(coefs, value=1))  # avoid divide by zero
    assert "int" in caplog.records[0].message  # check logging output
    assert "sbe4" in caplog.records[1].message  # check logging output
    assert all(np.isnan(cnv[:-1]))
    assert cnv[-1] == 1
    assert cnv.dtype == float

    # check there are no warnings if freq is float with no zeroes
    no_warn = ctd.sbe4(np.ones(len(freq)), t, p, make_coefs(coefs, value=1))
    assert len(caplog.records) == 2
    assert all(no_warn == 1)
    assert no_warn.dtype == float

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="CTcor"):
        ctd.sbe4(freq, t, p, make_coefs(coefs[:-1]))


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
    cnv = ctd.sbe9(freq, t, make_coefs(coefs))
    assert "int" in caplog.records[0].message  # check logging output
    assert "sbe9" in caplog.records[1].message  # check logging output
    assert all(np.isnan(cnv[:-1]))
    assert cnv[-1] == -10.1353
    assert cnv.dtype == float

    # check there are no warnings if freq is float with no zeroes
    no_warn = ctd.sbe9(np.ones(len(freq)), t, make_coefs(coefs))
    assert len(caplog.records) == 2
    assert all(no_warn == -10.1353)
    assert no_warn.dtype == float

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="AD590B"):
        ctd.sbe9(freq, t, make_coefs(coefs[:-1]))


def test_sbe_altimeter(caplog):
    volts = 99 * [0] + [1]
    coefs = ["ScaleFactor", "Offset"]

    # check values are converted correctly
    cnv = aux.sbe_altimeter(volts, make_coefs(coefs, value=1))  # avoid divide by zero
    assert "int" in caplog.records[0].message
    assert "sbe_altimeter" in caplog.records[0].message
    assert all(cnv[:-1] == 1)
    assert cnv[-1] == 301
    assert cnv.dtype == float

    # check there are no warnings if volts are float with no zeroes
    no_warn = aux.sbe_altimeter(np.ones(len(volts)), make_coefs(coefs, value=1))
    assert len(caplog.records) == 1
    assert all(no_warn == 301)
    assert no_warn.dtype == float

    # check voltages outside 0-5V are replaced with NaNs
    bad_range = aux.sbe_altimeter(np.arange(-1, 7), make_coefs(coefs, value=1))
    assert "sbe_altimeter" in caplog.records[-1].message
    assert "outside of 0-5V" in caplog.records[-1].message
    assert np.isnan(bad_range[0])
    assert np.isnan(bad_range[-1])
    assert not all(np.isnan(bad_range[1:-1]))

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="Offset"):
        aux.sbe_altimeter(volts, make_coefs(coefs[:-1]))


def test_sbe43(caplog):
    volts = 99 * [0] + [1]
    p = len(volts) * [0]
    t = len(volts) * [0]
    c = len(volts) * [0]
    coefs = ["Soc", "offset", "Tau20", "A", "B", "C", "E"]

    # check values are converted correctly
    cnv = oxy.sbe43(volts, p, t, c, make_coefs(coefs, value=1))
    assert "int" in caplog.records[0].message
    assert "sbe43" in caplog.records[0].message
    assert all(cnv[:-1] == 10.2314)
    assert cnv[-1] == 20.4628
    assert cnv.dtype == float

    # check there are no warnings if volts are float with no zeroes
    no_warn = oxy.sbe43(np.ones(len(volts)), p, t, c, make_coefs(coefs, value=1))
    assert len(caplog.records) == 1
    assert all(no_warn == 20.4628)
    assert no_warn.dtype == float

    # check voltages outside 0-5V are replaced with NaNs
    volts = np.arange(-1, 7)
    p, t, c = tuple(np.ones(len(volts)) for x in [1, 2, 3])  # definition shorthand
    bad_range = oxy.sbe43(volts, p, t, c, make_coefs(coefs, value=1))
    assert "sbe43" in caplog.records[-1].message
    assert "outside of 0-5V" in caplog.records[-1].message
    assert np.isnan(bad_range[0])
    assert np.isnan(bad_range[-1])
    assert not all(np.isnan(bad_range[1:-1]))

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="E"):
        oxy.sbe43(volts, p, t, c, make_coefs(coefs[:-1]))


def test_sbe43_hysteresis_voltage(caplog):
    volts = 99 * [0] + [1]
    p = len(volts) * [0]
    coefs = ["H1", "H2", "H3", "offset"]

    # check values are converted correctly
    cnv = oxy.sbe43_hysteresis_voltage(volts, p, make_coefs(coefs, value=1))
    assert "int" in caplog.records[0].message
    assert "sbe43_hysteresis_voltage" in caplog.records[0].message
    assert all(cnv[:-1] == 0)
    assert cnv[-1] == 1
    assert cnv.dtype == float

    # check there are no warnings if volts are float with no zeroes
    no_warn = oxy.sbe43_hysteresis_voltage(
        np.ones(len(volts)), p, make_coefs(coefs, value=1)
    )
    assert len(caplog.records) == 1
    assert all(no_warn == 1)
    assert no_warn.dtype == float

    # check voltages outside 0-5V are replaced with NaNs
    volts = np.arange(0, 7)
    p = np.ones(len(volts))
    bad_range = oxy.sbe43_hysteresis_voltage(volts, p, make_coefs(coefs, value=1))
    assert "sbe43_hysteresis_voltage" in caplog.records[-1].message
    assert "outside of 0-5V" in caplog.records[-1].message
    assert not all(np.isnan(bad_range[:-1]))
    assert np.isnan(bad_range[-1])

    # check NaNs are propagated through hysteresis correction
    volts = np.arange(-1, 7)
    p = np.ones(len(volts))
    bad_range = oxy.sbe43_hysteresis_voltage(volts, p, make_coefs(coefs, value=1))
    assert all(np.isnan(bad_range))

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="offset"):
        oxy.sbe43_hysteresis_voltage(volts, p, make_coefs(coefs[:-1]))


def test_wetlabs_eco_fl(caplog):
    volts = 99 * [0] + [1]
    coefs = ["ScaleFactor", "DarkOutput"]

    # check values are converted correctly
    cnv = aux.wetlabs_eco_fl(volts, make_coefs(coefs, value=1))
    assert "int" in caplog.records[0].message
    assert "wetlabs_eco_fl" in caplog.records[0].message
    assert all(cnv[:-1] == -1)
    assert cnv[-1] == 0
    assert cnv.dtype == float

    # check for same result with "Vblank"
    coefs = ["ScaleFactor", "Vblank"]
    cnv_Vblank = aux.wetlabs_eco_fl(volts, make_coefs(coefs, value=1))
    assert all(cnv == cnv_Vblank)

    # check there are no warnings if volts are float with no zeroes
    no_warn = aux.wetlabs_eco_fl(np.ones(len(volts)), make_coefs(coefs, value=1))
    assert len(caplog.records) == 2
    assert all(no_warn == 0)
    assert no_warn.dtype == float

    # check voltages outside 0-5V are replaced with NaNs
    bad_range = aux.wetlabs_eco_fl(np.arange(-1, 7), make_coefs(coefs, value=1))
    assert "wetlabs_eco_fl" in caplog.records[-1].message
    assert "outside of 0-5V" in caplog.records[-1].message
    assert np.isnan(bad_range[0])
    assert np.isnan(bad_range[-1])
    assert not all(np.isnan(bad_range[1:-1]))

    # return volts if "DarkOutput" or "Vblank" not in coefs
    no_cnv = aux.wetlabs_eco_fl(np.ones(len(volts)), make_coefs(["ScaleFactor"]))
    assert "returning voltage" in caplog.records[-1].message
    assert all(no_cnv == np.ones(len(volts)))


def test_wetlabs_cstar(caplog):
    volts = 99 * [0] + [1]
    coefs = ["M", "B", "PathLength"]

    # check values are converted correctly
    cnv = aux.wetlabs_cstar(volts, make_coefs(coefs, value=1))  # avoid divide by zero
    xmiss, c = cnv
    assert "int" in caplog.records[0].message
    assert "wetlabs_cstar" in caplog.records[0].message
    assert all(xmiss[:-1] == 1)
    assert xmiss[-1] == 2
    assert all(c[:-1] == -4.6052)
    assert c[-1] == -5.2983

    # check there are no warnings if volts are float with no zeroes
    no_warn = aux.wetlabs_cstar(np.ones(len(volts)), make_coefs(coefs, value=1))
    assert len(caplog.records) == 1
    assert all(no_warn[0] == 2)
    assert all(no_warn[1] == -5.2983)

    # check voltages outside 0-5V are replaced with NaNs
    bad_range = aux.wetlabs_cstar(np.arange(-1, 7), make_coefs(coefs, value=1))
    xmiss, c = bad_range
    assert "wetlabs_cstar" in caplog.records[-1].message
    assert "outside of 0-5V" in caplog.records[-1].message
    assert np.isnan(xmiss[0])
    assert np.isnan(xmiss[-1])
    assert not all(np.isnan(xmiss[1:-1]))
    assert np.isnan(c[0])
    assert np.isnan(c[-1])
    assert not all(np.isnan(c[1:-1]))

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="PathLength"):
        aux.wetlabs_cstar(volts, make_coefs(coefs[:-1]))


def test_seapoint_fluor():
    volts = 99 * [0] + [1]
    coefs = ["GainSetting", "Offset"]

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="Offset"):
        aux.seapoint_fluor(volts, make_coefs(coefs[:-1]))

def test_sbe_flntu_chl(caplog):
    volts = np.concatenate((np.zeros(10), np.array([1, 2, 3, 4, 5])))
    coefs = {'ScaleFactor':11, 'Vblank':0.079}

    chl = aux.sbe_flntu_chl(volts, coefs)

    #   Test hand-calculated values with example coefs
    assert all(chl[0:10] == -0.869)
    assert chl[-1] == 54.131 

    # check there are no warnings if volts are float with no zeroes
    no_warn = aux.sbe_flntu_chl(np.ones(len(volts)), make_coefs(coefs, value=1))
    assert len(caplog.records) == 0
    assert all(no_warn == 0)

    # check voltages outside 0-5V are replaced with NaNs
    bad_range = aux.sbe_flntu_chl(np.arange(-1, 7), make_coefs(coefs, value=1))
    assert "sbe_flntu_chl" in caplog.records[-1].message
    assert "outside of 0-5V" in caplog.records[-1].message
    assert np.isnan(bad_range[0])
    assert np.isnan(bad_range[-1])
    assert bad_range[1] == -1
    assert not all(np.isnan(bad_range[1:-1]))

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="Vblank"):
        aux.sbe_flntu_chl(volts, {"ScaleFactor":1, "dark_counts":1})
        assert "dictionary missing keys" in caplog.records[-1].message

def test_sbe_flntu_ntu(caplog):
    #   Practically identical to sbe_flntu_chl
    volts = np.concatenate((np.zeros(10), np.array([1, 2, 3, 4, 5])))
    coefs = {'ScaleFactor':5, 'DarkVoltage':0.050}

    ntu = aux.sbe_flntu_ntu(volts, coefs)

    #   Test hand-calculated values with example coefs
    assert all(ntu[0:10] == -0.25)
    assert ntu[-1] == 24.75

    # check there are no warnings if volts are float with no zeroes
    no_warn = aux.sbe_flntu_ntu(np.ones(len(volts)), make_coefs(coefs, value=1))
    assert len(caplog.records) == 0
    assert all(no_warn == 0)

    # check voltages outside 0-5V are replaced with NaNs
    bad_range = aux.sbe_flntu_ntu(np.arange(-1, 7), make_coefs(coefs, value=1))
    assert "sbe_flntu_ntu" in caplog.records[-1].message
    assert "outside of 0-5V" in caplog.records[-1].message
    assert np.isnan(bad_range[0])
    assert np.isnan(bad_range[-1])
    assert bad_range[1] == -1
    assert not all(np.isnan(bad_range[1:-1]))

    # error saying which keys are missing from coef dict
    with pytest.raises(KeyError, match="DarkVoltage"):
        aux.sbe_flntu_ntu(volts, {"ScaleFactor":1, "dark_counts":1})
        assert "dictionary missing keys" in caplog.records[-1].message
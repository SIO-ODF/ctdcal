import logging
from pathlib import Path

import numpy as np
import scipy

from ctdcal import oxy_fitting


def test_gather_oxy_params(caplog, tmp_path):
    # make fake file path
    f_path = tmp_path / "90909"
    # breakpoint()
    assert not Path(tmp_path / "90909").exists()
    with caplog.at_level(logging.INFO):
        oxy_params = oxy_fitting.gather_oxy_params(f_path)
        assert "Failed to load" in caplog.messages[0]
        assert oxy_params.isnull().values.all()

def test_calculate_weights():
    #   Set some quick db values
    pressure = np.array([50, 150, 250, 400, 800, 1500, 2500, 5000])
    p_bins = [
        0,
        100,
        100 + 1e-5, #   Epsilon = 1e-5
        300,
        300 + 1e-5,
        500,
        500 + 1e-5,
        1200,
        1200 + 1e-5,
        2000,
        2000 + 1e-5,
        7000,
    ]
    w_bins = [20, 20, 25, 25, 50, 50, 100, 100, 200, 200, 500, 500]
    wgt_manual = scipy.interpolate.interp1d(p_bins, w_bins)(pressure)   #   Pull weights out

    wgt = oxy_fitting.calculate_weights(pressure)

    assert np.array_equal(wgt, wgt_manual)
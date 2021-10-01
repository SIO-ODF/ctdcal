import numpy as np
import pandas as pd
import pytest

from ctdcal import ctd_plots


def test_residual_vs_pressure():
    param = pd.Series(data=[0, 1, 2], name="param")
    ref = pd.Series(data=[0, 1, 2], name="ref")
    prs = pd.Series(data=[0, 3000, 6000], name="prs")
    stn = pd.Series(data=[0, 1, 2], name="stn")

    # check titles, labels, limits
    axes = ctd_plots.residual_vs_pressure(param, ref, prs, stn)
    assert axes.get_title() == "ref-param vs. prs"
    assert axes.get_xlabel() == "Residual"
    assert axes.get_ylabel() == "Pressure (dbar)"
    assert axes.get_xlim() == (-0.02, 0.02)
    assert axes.get_ylim() == (6000, 0)

    # check titles, labels, limits (with deep setting)
    deep = ctd_plots.residual_vs_pressure(param, ref, prs, stn, deep=True)
    y_data = deep.collections[0].get_offsets().data[:, 1]
    assert deep.get_title() == "ref-param (>2000 dbar) vs. prs"
    assert deep.get_xlabel() == "Residual"
    assert deep.get_ylabel() == "Pressure (dbar)"
    assert deep.get_xlim() == (-0.02, 0.02)
    assert deep.get_ylim() == (6000, 0)
    assert all(y_data > 2000)

    # check behavior with unlabeled data inputs
    x = np.array([0, 0, 0])
    unlabeled = ctd_plots.residual_vs_pressure(x, x, x, x)
    assert unlabeled.get_title() == ""

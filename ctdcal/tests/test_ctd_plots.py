import numpy as np
import pandas as pd
import pytest

from ctdcal import ctd_plots


def test_residual_vs_pressure():
    param = pd.Series(data=[0, 1, 2], name="param")
    ref = pd.Series(data=[0, 1, 2], name="ref")
    prs = pd.Series(data=[0, 3000, 6000], name="prs")
    stn = pd.Series(data=[0, 1, 2], name="stn")

    # check titles, labels, limits, grid
    axes = ctd_plots.residual_vs_pressure(param, ref, prs, stn)
    grid_lines = axes.get_xgridlines() + axes.get_ygridlines()
    assert axes.get_title() == "ref-param vs. prs"
    assert axes.get_xlabel() == "Residual"
    assert axes.get_ylabel() == "Pressure (dbar)"
    assert axes.get_xlim() == (-0.02, 0.02)
    assert axes.get_ylim() == (6000, 0)
    assert not all(line.get_visible() for line in grid_lines)

    # check titles, labels, limits, grid (with deep setting)
    deep = ctd_plots.residual_vs_pressure(param, ref, prs, stn, deep=True)
    grid_lines = deep.get_xgridlines() + deep.get_ygridlines()
    y_data = deep.collections[0].get_offsets().data[:, 1]
    assert deep.get_title() == "ref-param (>2000 dbar) vs. prs"
    assert deep.get_xlabel() == "Residual"
    assert deep.get_ylabel() == "Pressure (dbar)"
    assert deep.get_xlim() == (-0.02, 0.02)
    assert deep.get_ylim() == (6000, 0)
    assert not all(line.get_visible() for line in grid_lines)
    assert all(y_data > 2000)

    # check settings are applied properly
    settings = ctd_plots.residual_vs_pressure(
        param,
        ref,
        prs,
        stn,
        xlim=(-1, 1),
        ylim=(2000, 1000),
        xlabel="x",
        ylabel="y",
        auto_title=False,
        grid=True,
    )
    grid_lines = settings.get_xgridlines() + settings.get_ygridlines()
    assert settings.get_xlim() == (-1, 1)
    assert settings.get_ylim() == (2000, 1000)
    assert settings.get_xlabel() == "x"
    assert settings.get_ylabel() == "y"
    assert settings.get_title() == ""
    assert all(line.get_visible() for line in grid_lines)

    # check behavior with unlabeled data inputs
    x = np.array([0, 0, 0])
    unlabeled = ctd_plots.residual_vs_pressure(x, x, x, x)
    assert unlabeled.get_title() == ""


def test_residual_vs_station():
    param = pd.Series(data=[0, 1, 2], name="param")
    ref = pd.Series(data=[0, 1, 2], name="ref")
    prs = pd.Series(data=[0, 3000, 6000], name="prs")
    stn = pd.Series(data=[0, 1, 2], name="stn")

    # check titles, labels, limits
    axes = ctd_plots.residual_vs_station(param, ref, prs, stn)
    grid_lines = axes.get_xgridlines() + axes.get_ygridlines()
    assert axes.get_title() == "ref-param vs. stn"
    assert axes.get_xlabel() == "Station Number"
    assert axes.get_ylabel() == "Residual"
    # xlim intentionally skipped as it's not currently a setting for this function
    assert axes.get_ylim() == (-0.02, 0.02)
    assert not all(line.get_visible() for line in grid_lines)

    # check titles, labels, limits (with deep setting)
    deep = ctd_plots.residual_vs_station(param, ref, prs, stn, deep=True)
    grid_lines = deep.get_xgridlines() + deep.get_ygridlines()
    z_data = deep.collections[0].get_array().data
    assert deep.get_title() == "ref-param (>2000 dbar) vs. stn"
    assert deep.get_xlabel() == "Station Number"
    assert deep.get_ylabel() == "Residual"
    # xlim intentionally skipped as it's not currently a setting for this function
    assert deep.get_ylim() == (-0.02, 0.02)
    assert not all(line.get_visible() for line in grid_lines)
    assert all(z_data > 2000)

    # check settings are applied properly
    settings = ctd_plots.residual_vs_station(
        param,
        ref,
        prs,
        stn,
        ylim=(-1, 1),
        xlabel="x",
        ylabel="y",
        grid=True,
    )
    grid_lines = settings.get_xgridlines() + settings.get_ygridlines()
    assert settings.get_ylim() == (-1, 1)
    assert settings.get_xlabel() == "x"
    assert settings.get_ylabel() == "y"
    assert all(line.get_visible() for line in grid_lines)

    # check behavior with unlabeled data inputs
    x = np.array([0, 0, 0])
    unlabeled = ctd_plots.residual_vs_station(x, x, x, x)
    assert unlabeled.get_title() == ""

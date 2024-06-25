import numpy as np
import pandas as pd
import pytest

from ctdcal import ctd_plots


def test_residual_vs_pressure(tmp_path):
    param = pd.Series(data=[0, 1, 2], name="param")
    ref = pd.Series(data=[0, 1, 2], name="ref")
    prs = pd.Series(data=[0, 3000, 6000], name="prs")
    stn = pd.Series(data=[0, 1, 2], name="stn")

    # check titles, labels, limits, grid
    axes = ctd_plots.residual_vs_pressure(param, ref, prs, stn=stn)
    grid_lines = axes.get_xgridlines() + axes.get_ygridlines()
    assert axes.get_title() == "ref-param vs. prs"
    assert axes.get_xlabel() == "Residual"
    assert axes.get_ylabel() == "Pressure (dbar)"
    assert axes.get_xlim() == (-0.02, 0.02)
    assert axes.get_ylim() == (6000, 0)
    assert axes.collections[0].colorbar is not None
    assert not all(line.get_visible() for line in grid_lines)

    # check titles, labels, limits, grid (with deep setting)
    deep = ctd_plots.residual_vs_pressure(param, ref, prs, stn=stn, deep=True)
    grid_lines = deep.get_xgridlines() + deep.get_ygridlines()
    y_data = deep.collections[0].get_offsets().data[:, 1]
    assert deep.get_title() == "ref-param (>2000 dbar) vs. prs"
    assert deep.get_xlabel() == "Residual"
    assert deep.get_ylabel() == "Pressure (dbar)"
    assert deep.get_xlim() == (-0.02, 0.02)
    assert deep.get_ylim() == (6000, 0)
    assert deep.collections[0].colorbar is not None
    assert not all(line.get_visible() for line in grid_lines)
    assert all(y_data > 2000)

    # check settings are applied properly
    settings = ctd_plots.residual_vs_pressure(
        param,
        ref,
        prs,
        stn=stn,
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
    assert settings.collections[0].colorbar is not None
    assert all(line.get_visible() for line in grid_lines)

    # check behavior with unlabeled data inputs and no stations
    x = np.array([0, 0, 0])
    unlabeled = ctd_plots.residual_vs_pressure(x, x, x)
    assert unlabeled.get_title() == ""
    assert unlabeled.collections[0].colorbar is None

    # check figure saving
    # TODO: pymark parameterize?
    # for ext in [".jpg", ".png", ".pdf"]:
    #     with (tmp_path / "figures" / f"fig1{ext}") as f_out:
    #         assert not f_out.exists()
    #         ctd_plots.residual_vs_pressure(x, x, x, stn=x, f_out=f_out)
    #         assert f_out.exists()

    for ext in [".jpg", ".png", ".pdf"]:
        output_path = tmp_path / "figures"
        output_path.mkdir(parents=True, exist_ok=True)
        file_path = output_path / f"fig1{ext}"

        assert not file_path.exists()
        ctd_plots.residual_vs_pressure(x, x, x, stn=x, f_out=file_path)
        assert file_path.exists()   # Verify that the file has been created


def test_residual_vs_station(tmp_path):
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

    # check figure saving
    # TODO: pymark parameterize?
    for ext in [".jpg", ".png", ".pdf"]:
        # with (tmp_path / "figures" / f"fig2{ext}") as f_out:
        #     assert not f_out.exists()
        #     ctd_plots.residual_vs_station(x, x, x, x, f_out=f_out)
        #     assert f_out.exists()
        output_path = tmp_path / "figures"
        output_path.mkdir(parents=True, exist_ok=True)
        file_path = output_path / f"fig2{ext}"

        assert not file_path.exists()
        ctd_plots.residual_vs_pressure(x, x, x, stn=x, f_out=file_path)
        assert file_path.exists()   # Verify that the file has been created


def test_intermediate_residual_plot(tmp_path):
    param = pd.Series(data=[0, 1, 2], name="param")
    ref = pd.Series(data=[0, 1, 2], name="ref")
    prs = pd.Series(data=[0, 3000, 6000], name="prs")
    stn = pd.Series(data=[0, 1, 2], name="stn")

    # check settings are applied properly
    settings = ctd_plots._intermediate_residual_plot(
        ref - param, prs, stn, xlim=(-1, 1), xlabel="x", show_thresh=True
    )
    grid_lines = settings.get_xgridlines() + settings.get_ygridlines()
    thresh_left = settings.lines[0].get_xydata()
    thresh_right = settings.lines[1].get_xydata()
    assert settings.get_xlim() == (-1, 1)
    assert settings.get_xlabel() == "x"
    assert settings.get_title() == "Mean: 0.0 / Stdev: 0.0"
    assert all(line.get_visible() for line in grid_lines)
    assert all(thresh_left[:, 0] == np.array([0.002, 0.005, 0.010, 0.020, 0.020]))
    assert all(thresh_right[:, 0] == -1 * np.array([0.002, 0.005, 0.010, 0.020, 0.020]))
    assert all(thresh_left[:, 1] == np.array([6000, 2000, 1000, 500, 0]))
    assert all(thresh_right[:, 1] == np.array([6000, 2000, 1000, 500, 0]))

    # check figure saving
    # TODO: pymark parameterize?
    x = np.array([0, 0, 0])
    for ext in [".jpg", ".png", ".pdf"]:
        # with (tmp_path / "figures" / f"fig3{ext}") as f_out:
        #     assert not f_out.exists()
        #     ctd_plots._intermediate_residual_plot(x, x, x, f_out=f_out)
        #     assert f_out.exists()
        output_path = tmp_path / "figures"
        output_path.mkdir(parents=True, exist_ok=True)
        file_path = output_path / f"fig3{ext}"

        assert not file_path.exists()
        ctd_plots.residual_vs_pressure(x, x, x, stn=x, f_out=file_path)
        assert file_path.exists()   # Verify that the file has been created

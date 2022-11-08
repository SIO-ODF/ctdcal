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
    for ext in [".jpg", ".png", ".pdf"]:
        with (tmp_path / "figures" / f"fig1{ext}") as f_out:
            assert not f_out.exists()
            ctd_plots.residual_vs_pressure(x, x, x, stn=x, f_out=f_out)
            assert f_out.exists()


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
        with (tmp_path / "figures" / f"fig2{ext}") as f_out:
            assert not f_out.exists()
            ctd_plots.residual_vs_station(x, x, x, x, f_out=f_out)
            assert f_out.exists()


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
        with (tmp_path / "figures" / f"fig3{ext}") as f_out:
            assert not f_out.exists()
            ctd_plots._intermediate_residual_plot(x, x, x, f_out=f_out)
            assert f_out.exists()


def test_haversine():
    #   Check that values are calculated correctly
    lon1 = 15.5
    lat1 = 0.0
    lon2 = 21.75
    lat2 = -6.35  #   Or 6.35 DD south
    x = ctd_plots.haversine(lon1, lat1, lon2, lat2)
    assert int(x) == 964
    assert type(x) == float


@pytest.fixture(scope="module")
def df():
    #   Fake data for dataframe tests
    df = pd.DataFrame(
        [
            [32.857, -117.263, 97.5, 2320.6, 34.6, 12.41, "001"],
            [32.857, -117.263, 504.6, 2320.6, 34.5, 4.46, "001"],
            [32.857, -117.263, 1126.3, 2320.6, 34.4, 1.81, "001"],
            [32.857, -117.263, 2303.5, 2320.6, 34.4, 1.1, "001"],
            [32.851, -117.271, 97.5, 2381.5, 34.7, 12.52, "002"],
            [32.851, -117.271, 504.6, 2381.5, 34.2, 4.44, "002"],
            [32.851, -117.271, 1126.3, 2381.5, 34.4, 1.80, "002"],
            [32.851, -117.271, 2303.5, 2381.5, 34.4, 1.1, "002"],
            [32.849, -117.288, 97.5, 2341.5, 34.9, 12.53, "003"],
            [32.849, -117.288, 504.6, 2341.5, 34.3, 4.47, "003"],
            [32.849, -117.288, 1126.3, 2341.5, 34.4, 1.85, "003"],
            [32.849, -117.288, 2303.5, 2341.5, 34.4, 1.1, "003"],
            [32.858, -117.305, 97.5, 3001.6, 34.9, 12.59, "004"],
            [32.858, -117.305, 504.6, 3001.6, 34.3, 4.42, "004"],
            [32.858, -117.305, 1126.3, 3001.6, 34.4, 1.83, "004"],
            [32.858, -117.305, 2890.5, 3001.6, 34.45, 0.7, "004"],
        ],
        columns=[
            "LATITUDE",
            "LONGITUDE",
            "CTDPRS",
            "DEPTH",
            "CTDSAL",
            "CTDTMP1",
            "SSSCC",
        ],
    )
    return df


def test_section_bottle_plot(df, tmp_path):
    fig = ctd_plots.section_bottle_plot(df)
    #   Check that figure limits were correctly assigned (haversine and depth rounding worked)
    assert fig.get_ylim() == (3500.0, 0.0)
    assert fig.get_xlim() == (0.0, 4.482320177048369)
    assert fig.get_ylabel() == "CTDPRS (dbar)"
    assert fig.get_xlabel() == "Section Distance (km)"
    #   TODO: Get colorbar information?

    # check figure saving
    for ext in [".jpg", ".png", ".pdf"]:
        with (tmp_path / "figures" / f"fig3{ext}") as f_out:
            assert not f_out.exists()
            ctd_plots.section_bottle_plot(df, f_out=f_out)
            assert f_out.exists()


def test_bottle_TS_plot(df, tmp_path):
    # #   TODO: Refactor bottle_TS_plot to return figure information
    # fig = ctd_plots.bottle_TS_plot(df)
    # assert fig.get_ylim() == (0.44099999999999984, 13.440999999999999)

    # check figure saving
    for ext in [".jpg", ".png", ".pdf"]:
        with (tmp_path / "figures" / f"fig3{ext}") as f_out:
            assert not f_out.exists()
            ctd_plots.bottle_TS_plot(df, f_out=f_out)
            assert f_out.exists()

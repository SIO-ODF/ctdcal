import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

def _btl_ql(filepath):
    """
    A bottle quick-loader for testing plotting routines on hy1 file
    Code originally by MK, data_qc.py
    """
    from ctdcal import io
    btl_data = io.load_exchange_btl(filepath).replace(-999, np.nan)
    btl_data["SSSCC"] = btl_data["STNNBR"].apply(lambda x: f"{x:03d}") + btl_data[
        "CASTNO"
    ].apply(lambda x: f"{x:02d}")
    return btl_data

def _btl_fl():
    """
    Load all of the bottle data for testing plotting routines prior to hy1 file generation
    """
    from ctdcal import process_bottle
    if not 'ssscc_list' in locals():
        from ctdcal import process_ctd
        ssscc_list = process_ctd.get_ssscc_list()
    btl_data = process_bottle.load_all_btl_files(ssscc_list)
    return btl_data

def _apply_default_fmt(xlim, ylim, xlabel, ylabel, title, grid, fontsize=12):
    plt.xlim(xlim)
    plt.xticks(rotation=45)
    plt.ylim(ylim)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.grid(grid)
    plt.tight_layout()


def _save_fig(ax, f_out):
    if f_out is not None:
        f_parent = Path(f_out).parent
        if not f_parent.exists():
            log.info(f"Parent folder '{f_parent.as_posix()}' doesn't exist... creating")
            Path(f_out).parent.mkdir(parents=True)
        plt.savefig(f_out)
        plt.close()
    else:
        return ax


def residual_vs_pressure(
    param,
    ref,
    prs,
    stn=None,
    xlim=(-0.02, 0.02),
    ylim=(6000, 0),
    xlabel="Residual",
    ylabel="Pressure (dbar)",
    auto_title=True,
    grid=False,
    deep=False,
    f_out=None,
):
    """
    Plot residuals (ref - param) as a function of pressure.

    Parameters
    ----------
    param : pd.Series or array-like
        Input variable
    ref : pd.Series or array-like
        Reference variable
    prs : pd.Series or array-like
        Pressure variable
    stn : pd.Series or array-like, optional
        Station variable
    xlim : tuple, optional
        Lower and upper x-limits
    ylim : tuple, optional
        Lower and upper y-limits
    xlabel : str, optional
        Label for the x-axis
    ylabel : str, optional
        Label for the y-axis
    auto_title : bool, optional
        Generate title from input attributes (iff dtype is pd.Series)
    grid : bool, optional
        Overlay grid on figure
    deep : bool, optional
        Whether to plot only values >2000 dbar
    f_out : path-like, optional
        Path to save figure

    Returns
    -------
    ax : matplotlib.axes, optional
        Formatted scatter plot (if f_out is not given)
    """
    diff = ref - param
    deep_text = " (>2000 dbar)" if deep else ""
    if deep:
        deep_rows = prs > 2000
        diff, prs, stn = diff[deep_rows], prs[deep_rows], stn[deep_rows]

    # initialize figure
    plt.figure(figsize=(7, 6))
    ax = plt.axes()

    # color scatter by stations if given
    if stn is not None:
        idx, uniques = pd.factorize(stn)  # find unique stations #s and index them
        sc = ax.scatter(diff, prs, c=idx, marker="+")
        cbar = plt.colorbar(sc, ax=ax, pad=0.1)  # set cbar ticks to station names
        tick_inds = cbar.get_ticks().astype(int)
        cbar.ax.yaxis.set_major_locator(ticker.FixedLocator(tick_inds))
        cbar.ax.set_yticklabels(uniques[tick_inds])
        cbar.ax.set_title("Station")
    else:
        sc = ax.scatter(diff, prs, marker="+")

    # formatting
    title = None
    if auto_title:
        try:
            title = f"{ref.name}-{param.name}{deep_text} vs. {prs.name}"
        except AttributeError:
            log.warning(
                "Failed to set title from variable names (requires dtype pd.Series)"
            )
            log.info("Set afterward using 'ax.set_title(\"title\")'")
    _apply_default_fmt(xlim, ylim, xlabel, ylabel, title, grid)

    # save to path or return axis
    return _save_fig(ax, f_out)


def residual_vs_station(
    param,
    ref,
    prs,
    stn,
    ylim=(-0.02, 0.02),
    xlabel="Station Number",
    ylabel="Residual",
    grid=False,
    deep=False,
    f_out=None,
):
    """
    Plot residuals (ref - param) as a function of station number.

    Parameters
    ----------
    param : pd.Series or array-like
        Input variable
    ref : pd.Series or array-like
        Reference variable
    prs : pd.Series or array-like
        Pressure variable
    stn : pd.Series or array-like, optional
        Station variable
    ylim : tuple, optional
        Lower and upper y-limits
    xlabel : str, optional
        Label for the x-axis
    ylabel : str, optional
        Label for the y-axis
    grid : bool, optional
        Overlay grid on figure
    deep : bool, optional
        Whether to plot only values >2000 dbar
    f_out : path-like, optional
        Path to save figure

    Returns
    -------
    ax : matplotlib.axes, optional
        Formatted scatter plot (if f_out is not given)
    """
    diff = ref - param
    deep_text = " (>2000 dbar)" if deep else ""
    if deep:
        deep_rows = prs > 2000
        diff, prs, stn = diff[deep_rows], prs[deep_rows], stn[deep_rows]

    plt.figure(figsize=(7, 6))
    ax = plt.axes()
    sc = ax.scatter(stn, diff, c=prs, marker="+")
    cbar = plt.colorbar(sc, pad=0.1)
    cbar.set_label("Pressure (dbar)")

    # formatting
    try:
        title = f"{ref.name}-{param.name}{deep_text} vs. {stn.name}"
    except AttributeError:
        title = None
        log.warning(
            "Failed to set title from variable names (requires dtype pd.Series)"
        )
        log.info('Set afterward using \'ax.set_title("title")`')
    _apply_default_fmt(None, ylim, xlabel, ylabel, title, grid)

    # save to path or return axis
    return _save_fig(ax, f_out)


def _intermediate_residual_plot(
    diff,
    prs,
    stn,
    xlim=(-0.02, 0.02),
    xlabel="Residual",
    show_thresh=False,
    f_out=None,
):
    """
    Internal function to make figures at intermediate processing stages for debugging.
    """
    ax = residual_vs_pressure(
        0, diff, prs, stn=stn, xlim=xlim, xlabel=xlabel, auto_title=False, grid=True
    )

    if show_thresh:
        thresh = np.array([0.002, 0.005, 0.010, 0.020])
        p_range = np.array([6000, 2000, 1000, 500])
        thresh = np.append(thresh, thresh[-1])  # this should still work fine even when
        p_range = np.append(p_range, 0)  # thresh/p_range are defined elsewhere
        plt.step(thresh, p_range, ":k")
        plt.step(-thresh, p_range, ":k")

    mean = np.round(np.nanmean(diff), 4)
    stdev = np.round(np.nanstd(diff), 4)
    ax.set_title(f"Mean: {mean} / Stdev: {stdev}")

    # save to path or return axis (primarily for testing)
    return _save_fig(ax, f_out)

def residual_3(x1, x2, xr, stn, f_out=None
):
    """
    Visualize Δ Τ1 and ref, T2 and ref, T1 and T2 (or cond, or SBE43 and RINKO vs winkler).
    Intended for 5 stations or more.
    x : Pandas Series
        Respective "CTDVARX", where xr is the reference parameter (CTDCOND1, CTDCOND2, SALNTY for example)
    stn: Pandas Series
        Station number, direct from bottle file (_btl_fl)
    """
    idx, uniques = pd.factorize(stn)

    delX1r = xr - x1    #   Δ = ref - meas, consistent with other methods in ctdcal
    delX2r = xr - x2
    delX12 = x2 - x1
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('Reference - Primary')
    ax.set_ylabel('Reference - Secondary')
    ax.set_zlabel('Primary - Secondary')
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlim(-0.5, 0.5)
    ax.set_zlim(-0.5, 0.5)
    sc = ax.scatter(delX1r, delX2r, delX12, c=idx, marker = '+')
    if len(uniques > 4):
        cbar = plt.colorbar(sc, ax=ax, pad=0.15)  # set cbar ticks to station names
        tick_locator = ticker.MaxNLocator(nbins=5)  #   Set 5 bins, 
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.set_yticklabels(uniques[0:-1:int(len(uniques)/5)])

    return _save_fig(ax, f_out)

def residual_bars(

):
    """
    Take the absolute sum of the residuals for the given params for given SSSCC
    """

def residual_hist(

):
    """
    Count the residuals within a specific range (A05 2020 example).
    """

def residual_depth_bar(

):
    """
    Bin data by depth of size X and stack bars of SSSCC for that depth
    """

def residual_btl_scatter(

):
    """
    Take the residuals for individual bottles of each cast and plot against depth (color code btl#/SSSCC)
    """

def TC_residual_corr(

):
    """
    Scatter of TX vs CX residuals, linear fit line w/ statistics
    """

def sensor_drift_plot(

):
    """
    Scatter of residuals from first of SSSCC sublist, and last of SSSCC sublist
    """

def autosal_offset(

):
    """
    For all SSSCC provided, do a scatter of the autosal offset applied to each run (individual files; A05 2020 example).
    """

def pressure_offset_vis(

):
    """
    For all SSSCC provided, do a scatter of the SBE9+ on-deck pre- vs post-cast pressures. Include dotted line as a "healthy" offset boundary.
    """

def btl_fits_plot(
):
    """
    Plot TX vs depth before and after fitting from the btl files. Include coefs in title.
    """

def interpolate_ssscc(

):
    """
    Plot param vs depth, calculating distance between SSSCC and interpolating between them.
    """

def ssscc_geo(

):
    """
    Plot station locations, with SSSCC and depth labels
    """

def plot_bio(

):
    """
    Depth plot of biologically-important sensors: TS, Fluor, Xmiss, Oxy
    """

def corr(

):
    """
    Create color-coded correlation matrix between CTD parameters
    """

def TS_route(

):
    """
    Plot temperature vs salinity, colored by depth, for individual SSSCCs.
    """


# TODO: more plots! what does ODV have?
# section plots (w/ PyDiva)
# parameter-parameter plots

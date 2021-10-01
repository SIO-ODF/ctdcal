import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def _apply_default_fmt(xlim, ylim, xlabel, ylabel, title):
    plt.xlim(xlim)
    plt.xticks(rotation=45)
    plt.ylim(ylim)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=12)
    plt.tight_layout()


def residual_vs_pressure(
    param,
    ref,
    prs,
    stn,
    xlim=(-0.02, 0.02),
    ylim=(6000, 0),
    xlabel="Residual",
    ylabel="Pressure (dbar)",
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
    deep_text = " (>2000 dbar) " if deep else " "
    if deep:
        deep_rows = prs > 2000
        diff, prs, stn = diff[deep_rows], prs[deep_rows], stn[deep_rows]

    # attempt to make title from pd.Series.name
    try:
        title = f"{ref.name}-{param.name}{deep_text}vs. {prs.name}"
    except AttributeError:
        title = None
        log.warning(
            "Failed to set title from variable names (requires dtype pd.Series)"
        )
        log.info('Set afterward using \'ax.set_title("title")`')

    # initialize figure
    plt.figure(figsize=(7, 6))
    ax = plt.axes()

    # color scatter by stations if given
    idx, uniques = pd.factorize(stn)  # find unique stations #s and index them
    if len(uniques) >= 1:
        sc = ax.scatter(diff, prs, c=idx, marker="+")
        cbar = plt.colorbar(sc, ax=ax, pad=0.1)  # set cbar ticks to station names
        tick_inds = cbar.get_ticks().astype(int)
        cbar.ax.yaxis.set_major_locator(ticker.FixedLocator(tick_inds))
        cbar.ax.set_yticklabels(uniques[tick_inds])
        cbar.set_label("Station Number")
    else:
        sc = ax.scatter(diff, prs, marker="+")

    # formatting
    _apply_default_fmt(xlim, ylim, xlabel, ylabel, title)

    # save to path or return axis
    if f_out is not None:
        if not Path(f_out).parent.exists():
            log.info(
                f"Path {Path(f_out).parent.as_posix()} does not exists... creating"
            )
            Path(f_out).parent.mkdir(parents=True)
        plt.savefig(f_out)
        plt.close()
    else:
        return ax


def residual_vs_station(
    param,
    ref,
    prs,
    stn,
    ylim=(-0.02, 0.02),
    xlabel="Station Number",
    ylabel="Residual",
    deep=False,
    f_out=None,
):

    title = f"{ref.name}-{param.name} vs. {stn.name}"
    diff = ref - param
    if deep:
        title = f"{ref.name}-{param.name} (>2000 dbar) vs. {stn.name}"
        deep_rows = prs > 2000
        diff = diff[deep_rows]
        prs = prs[deep_rows]
        stn = stn[deep_rows]

    plt.figure(figsize=(7, 6))
    plt.scatter(stn, diff, c=prs, marker="+")
    plt.xticks(rotation=45)
    plt.ylim(ylim)
    cbar = plt.colorbar(pad=0.1)
    cbar.set_label("Pressure (dbar)")
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=12)
    plt.tight_layout()
    if f_out is not None:
        if not Path(f_out).parent.exists():
            log.info(
                f"Path {Path(f_out).parent.as_posix()} does not exists... creating"
            )
            Path(f_out).parent.mkdir(parents=True)
        plt.savefig(f_out)
    plt.close()

    return True


def _intermediate_residual_plot(
    diff,
    prs,
    ssscc,
    xlim=(-0.02, 0.02),
    ylim=(6000, 0),
    xlabel="Residual",
    ylabel="CTDPRS",
    show_thresh=False,
    f_out=None,
):

    idx, uniques = ssscc.factorize()  # find unique SSSCC and index them

    plt.figure(figsize=(6, 6))
    plt.scatter(diff, prs, c=idx, marker="+")
    if show_thresh:
        # TODO: thresh should probably be put in config/cast-by-cast config
        thresh = np.array([0.002, 0.005, 0.010, 0.020])
        p_range = np.array([6000, 2000, 1000, 500])
        thresh = np.append(thresh, thresh[-1])  # this should still work fine even when
        p_range = np.append(p_range, 0)  # thresh/p_range are defined elsewhere
        plt.step(thresh, p_range, ":k")
        plt.step(-thresh, p_range, ":k")

    plt.xlim(xlim)
    plt.xticks(rotation=45)
    plt.ylim(ylim)
    cbar = plt.colorbar(pad=0.1)  # set cbar ticks to SSSCC names
    if not uniques.empty:
        tick_inds = cbar.get_ticks().astype(int)
        cbar.ax.yaxis.set_major_locator(ticker.FixedLocator(tick_inds))
        cbar.ax.set_yticklabels(uniques[tick_inds])
        mean = np.round(np.nanmean(diff), 4)
        stdev = np.round(np.nanstd(diff), 4)
        plt.title(f"Mean: {mean} / Stdev: {stdev}")
    cbar.ax.set_title("SSSCC")
    plt.grid()
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.tight_layout()
    if f_out is not None:
        if not Path(f_out).parent.exists():
            log.info(
                f"Path {Path(f_out).parent.as_posix()} does not exists... creating"
            )
            Path(f_out).parent.mkdir(parents=True)
        plt.savefig(f_out)
    plt.close()

    return True

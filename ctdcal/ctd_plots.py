import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


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


def two_element(x1, x2, y, ssscc=None, caption=None, f_out=None):

    plt.style.use("ggplot")

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey="all")

    ax1.plot(x1, y)
    ax1.plot(x2, y)
    ax1.invert_yaxis()
    plt.tight_layout()
    ax1.tick_params(labelrotation=45)
    plt.grid()

    ax2.plot(x2 - x1, y)
    plt.tight_layout()
    ax2.tick_params(labelrotation=45)
    plt.grid()

    plt.gcf().subplots_adjust(bottom=0.2, left=0.1)

    if type(x1) == pd.core.series.Series:
        ax1.legend([x1.name, x2.name])
        ax1.set_ylabel(y.name)
        ax1.set_xlabel(x1.name + ", " + x2.name)
        ax2.set_xlabel("Residual " + x2.name + "-" + x1.name)

    if ssscc is not None:
        fig.suptitle(ssscc, x=0.53, y=1)

    if caption is not None:
        plt.figtext(
            0.5, 0.01, caption, wrap=True, horizontalalignment="center", fontsize=12
        )
    return _save_fig(fig, f_out)


# TODO: more plots! what does ODV have?
# section plots (w/ PyDiva)
# parameter-parameter plots

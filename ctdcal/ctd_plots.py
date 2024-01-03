import logging
from pathlib import Path
from ctdcal import get_ctdcal_config

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)
cfg = get_ctdcal_config()

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
    ylim=(5000, 0),
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
        if tick_inds[-1] > len(uniques)-1:
            tick_inds = tick_inds[:-1]
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
        p_range = np.array([5000, 2000, 1000, 500])
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

def fit_comparison(pre, post, ref, ssscc, varlabel="CTDOXY", grid=False):
    """Take in a variable pre and post time data, as well as reference, and plot scatter"""
    plt.figure(figsize=(8, 5))
    ax = plt.axes()
    plt.scatter(pre.CTDOXY, pre.CTDPRS, label="Prefit")
    plt.scatter(post.CTDOXY, post.CTDPRS, label="Postfit")
    plt.scatter(ref.BTL_OXY, ref.CTDPRS, s=300, marker="*", label="Bottle Reference")
    plt.gca().invert_yaxis()

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        fancybox=True,
        shadow=True,
        ncol=3,
    )
    title_label = str(varlabel + " " + ssscc)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    _apply_default_fmt(xlim, ylim, varlabel, "CTDPRES (dbar)", title_label, grid)
    plt.xlabel(varlabel, horizontalalignment="right", x=1.0)
    logs = "logs"
    f_out = f"{cfg.dirs[logs]}postfit/{varlabel}_{ssscc}.png"
    # _save_fig(ax, f_out)
    plt.savefig(f_out, bbox_inches="tight")

def plot_TS(df, f_out=None):
    """
    Line plot of temperature and salinity for a single cast or series of casts.
    Modified from: https://oceanpython.org/2013/02/17/t-s-diagram/
    """
    import matplotlib.ticker as tick
    import gsw

    #   Demarcate temperature, salinty, and limits for density
    temp = df.CTDTMP1
    salt = df.CTDSAL
    smin = salt.min() - (0.01 * salt.min())
    smax = salt.max() + (0.01 * salt.max())
    tmin = temp.min() - (0.1 * temp.max())
    tmax = temp.max() + (0.1 * temp.max())
    #   Gridding for a contour
    xdim = round((smax - smin) / 0.1 + 1.0)
    ydim = round((tmax - tmin) + 1.0)
    #   Creating density vectors to fill the grid with densities
    dens = np.zeros((ydim, xdim))
    si = np.linspace(1, xdim - 1, xdim) * 0.1 + smin
    ti = np.linspace(1, ydim - 1, ydim) + tmin
    for j in range(0, int(ydim)):
        for i in range(0, int(xdim)):
            dens[j, i] = gsw.rho(si[i], ti[j], 0)
    dens = dens - 1000  #   Convert to sigma-theta

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot()
    #   Density contours
    cs = plt.contour(si, ti, dens, linestyles="dashed", colors="k")
    plt.clabel(cs, fontsize=12, inline=1, fmt="%0.1f")
    #   Scatter of sample values
    sc = ax.scatter(salt, temp, c=df.SSSCC.astype(int), marker="+")
    sc.set_clim(vmin=df.SSSCC.astype(int).min(), vmax=df.SSSCC.astype(int).max())
    ax.set_xlabel("Salinity")
    ax.set_ylabel("Temperature (ºC)")
    cbar = plt.colorbar(sc, format=tick.FormatStrFormatter("%.0f"))
    cbar.set_label("Cast Number")

    _save_fig(ax, f_out)
# TODO: more plots! what does ODV have?
# section plots (w/ PyDiva)
# parameter-parameter plots


def TCcoherence_plot(btl_df, outdir=cfg.dirs["figs"], ext=".pdf"):
    plt.figure(figsize=(7, 6))
    plt.scatter(
        btl_df["CTDTMP1"] - btl_df["CTDTMP2"],
        btl_df["CTDCOND1"] - btl_df["CTDCOND2"],
        c=btl_df["CTDPRS"],
        marker="+",
    )
    plt.xticks(rotation=45)
    cbar = plt.colorbar(pad=0.1)
    cbar.set_label("Pressure (dbar)")
    plt.xlim((-0.05, 0.05))
    plt.ylim((-0.05, 0.05))
    plt.xlabel("CTDTMP1-CTDTMP2 Residual (T90 C)", fontsize=12)
    plt.ylabel("CTDCOND1-CTDCOND2 Residual (mS/cm)", fontsize=12)
    plt.title("CTDCOND1-CTDCOND2 vs. CTDTMP1-CTDTMP2", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{outdir}c_t_coherence{ext}")
    plt.close()

def conductivity_overlap(ssscc, btl_df, time_df, btl_df2=None, time_df2=None, title_lead="", outdir = cfg.dirs["figs"], ext=".pdf"):
    """
    Plot the bottle conductivity against the continuous CTD
    data streams on the primary and secondary lines.
    Optionally, pass in a second set of dataframes to layer
    them and demonstrate post-fit adjustments.

    Parameters
    ----------
    ssscc : string
        Current SSSCC
    btl_df : pd.Series or array-like
        Bottle dataframe for a specific SSSCC
    time_df : pd.Series or array-like
        Continuous ct1 dataframe for a specific SSSCC
    btl_df2 : pd.Series or array-like, optional
        Second set of bottle data to overlap over the given
        SSSCC (postfit)
    time_df2 : pd.Series or array-like, optional
        Second set of continuous time data to overlap
        over the given SSSCC (postfit)
    title_lead : string, optional
        Title for the plot, definable elsewhere by configuration
    outdir : string, optional
        The path to write the figure
    ext : string, optional
        Filename extension when saving plot
    """
    #   Make BTLCOND before running this
    
    legend_entries = ["Salinometer Cond", "Extracted COND1", "Extracted COND2", 
                      "Post-Fit Extracted COND1", "Post-Fit Extracted COND2", 
                      "CTD COND1 Trace", "CTD COND2 Trace", "Post-Fit CTD COND1 Trace", 
                      "Post-Fit CTD COND2 Trace"]
    plt.figure(figsize=(12,7))
    
    if any(~btl_df.BTLCOND.isnull()):
        plt.scatter(btl_df.BTLCOND,btl_df.CTDPRS,marker="+",label=legend_entries[0])

    plt.scatter(btl_df.CTDCOND1,btl_df.CTDPRS,marker="+",label=legend_entries[1])
    plt.scatter(btl_df.CTDCOND2,btl_df.CTDPRS,marker="+",label=legend_entries[2])

    if btl_df2 is not None:
        plt.scatter(btl_df2.CTDCOND1,btl_df2.CTDPRS,marker="+",label=legend_entries[3])
        plt.scatter(btl_df2.CTDCOND2,btl_df2.CTDPRS,marker="+",label=legend_entries[4])
        ext = "_fit"+ext
    
    plt.plot(time_df.CTDCOND1,time_df.CTDPRS,label=legend_entries[5])
    plt.plot(time_df.CTDCOND2,time_df.CTDPRS,label=legend_entries[6])

    if time_df2 is not None:
        plt.plot(time_df2.CTDCOND1,time_df2.CTDPRS,label=legend_entries[7])
        plt.plot(time_df2.CTDCOND1,time_df2.CTDPRS,label=legend_entries[8])
        if not "fit" in ext:
            ext = "_fit" + ext

    plt.gca().invert_yaxis()
    plt.ylabel("CTDPRS")
    plt.xlabel("Conductivity")
    plt.title(title_lead)
    plt.grid()
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"{outdir}/{ssscc}/cond_overlay{ext}")
    plt.close()
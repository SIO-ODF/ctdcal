import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from ctdcal import get_ctdcal_config

cfg = get_ctdcal_config()
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

    if not "ssscc_list" in locals():
        from ctdcal import process_ctd

        ssscc_list = process_ctd.get_ssscc_list()
    btl_data = process_bottle.load_all_btl_files(ssscc_list)
    return btl_data


#   Not sure where to put this
def haversine(lat1, lon1, lat2, lon2, units="km"):
    """
    Calculate the great circle distance in kilometers between two points
    on the earth using the haversine formula.
    For use on small distances (<4 degrees)

    Parameters
    ----------
    lat1, lon1, lat2, lon2 : Float
        Latitude and longitude of coordinate points 1 and 2, respectively
        in decimal degrees
    units : String
        Abbreviation of desired distance units (km, mi, or nmi)

    Returns
    -------
    dout : Float
        Calculated distance between points with units specified by 'units'
    """
    from math import radians, cos, sin, asin, sqrt

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    #   Discern output
    if units == "mi":
        print("Desired output units of miles in haversine formula.")
        r = 3956
    elif units == "nmi":
        print("Desired output units of nautical miles in haversine formula.")
        r = 3440.1
    else:
        if units != "km":
            print(
                "Provided units:",
                units,
                "not recognized. Assuming kilomters output in haversine formula.",
            )
        r = 6371  #   Assume kilometers
    dout = c * r
    return dout


def residual_3(x1, x2, xr, stn, f_out=None):
    """
    Visualize Δ Τ1 and ref, T2 and ref, T1 and T2 (or cond, or SBE43 and RINKO vs winkler).
    Intended for 5 stations or more, reading from bottle (hy1) file. Extract using _btl_ql.
    x : Pandas Series
        Respective "CTDVARX", where xr is the reference parameter (CTDCOND1, CTDCOND2, SALNTY for example)
    stn: Pandas Series
        Station number, direct from bottle file (_btl_fl)
    """
    idx, uniques = pd.factorize(stn)

    delX1r = xr - x1  #   Δ = ref - meas, consistent with other methods in ctdcal
    delX2r = xr - x2
    delX12 = x2 - x1
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_xlabel("Reference - Primary")
    ax.set_ylabel("Reference - Secondary")
    ax.set_zlabel("Primary - Secondary")
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlim(-0.5, 0.5)
    ax.set_zlim(-0.5, 0.5)
    sc = ax.scatter(delX1r, delX2r, delX12, c=idx, marker="+")
    if len(uniques > 4):
        cbar = plt.colorbar(sc, ax=ax, pad=0.15)  # set cbar ticks to station names
        tick_locator = ticker.MaxNLocator(nbins=5)  #   Set 5 bins,
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.set_yticklabels(uniques[0 : -1 : int(len(uniques) / 5)])

    return _save_fig(ax, f_out)


def section_bottle_plot(df, varname="CTDSAL", f_out=None, interp=True, cmap="viridis"):
    """
    Scatter plot of a desired variable against section distance for a bottle file.

    Parameters
    ----------
    df : Pandas DataFrame
        Bottle DataFrame with latitude, longitude, CTD pressure, and other variables
    varname : String
        String of variable intended to be plotted. Defaults to CTDSAL.

    """
    from math import ceil

    depth_max = ceil(df.DEPTH.max() / 500) * 500
    dot_width = 2

    #   Create section distance column
    x = df.LATITUDE.diff()
    x.reset_index(inplace=True, drop=True)
    idx = x[x != 0.0].index.tolist()
    df["dist"] = np.nan
    for i in idx:
        if i == 0:
            df["dist"].iloc[i] = 0.0
            val_add = 0.0
        else:
            df["dist"].iloc[i] = (
                haversine(
                    df.LATITUDE.iloc[i - 1],
                    df.LONGITUDE.iloc[i - 1],
                    df.LATITUDE.iloc[i],
                    df.LONGITUDE.iloc[i],
                )
                + val_add
            )  #   Uses last assignment
        val_add = df["dist"].iloc[
            i
        ]  #       Define new value to add for the beginning of the next iteration
    df["dist"].fillna(method="ffill", inplace=True)  #   Forward fill the nans

    if interp:
        plt.set_cmap(cmap)
        plt.figure(figsize=(7, 6))
        fig, ax = plt.subplots()
        ax.tricontour(
            df.dist, df.CTDPRS, df[varname], 20
        )  #   Creates the irregularly-spaced bounds for the contour
        plt.fill_between(
            df["dist"],
            df["DEPTH"],
            depth_max,
            interpolate=True,
            color="gray",
            zorder=2,
        )
        cntr2 = ax.tricontourf(
            df.dist, df.CTDPRS, df[varname], 20
        )  #   Fills in the contour
        ax.scatter(df.dist, df.CTDPRS, s=dot_width, color="k", zorder=3)
        cbar = fig.colorbar(cntr2, ax=ax)
    else:
        plt.figure(figsize=(7, 6))
        ax = plt.axes()
        a = ax.scatter(df["dist"], df["CTDPRS"], c=df[varname])
        plt.colorbar(a, ax=ax)
        plt.fill_between(
            df["dist"],
            df["DEPTH"],
            depth_max,
            interpolate=True,
            color="gray",
        )

    plt.ylim(0, depth_max)
    plt.ylabel("CTDPRS (dbar)")
    plt.xlabel("Section Distance (km)")
    plt.tight_layout()
    plt.gca().invert_yaxis()

    return _save_fig(ax, f_out)


def osnap_suite(btl_prefit, btl_fit, time_prefit, time_fit):

    #   Section plots with key params for each hydrographic line
    #   Could write as 2D array and plot iteratively...
    line1 = [
        "002",
        "003",
        "003",
        "004",
        "005",
        "006",
        "007",
        "008",
        "009",
        "010",
        "011",
        "012",
        "013",
        "014",
        "015",
        "016",
        "017",
        "018",
        "019",
        "020",
        "021",
        "022",
        "023",
        "024",
        "025",
        "026",
        "027",
        "028",
        "029",
        "030",
        "031",
    ]
    plt_btl1 = btl_fit.loc[btl_fit.SSSCC.isin(line1)]
    section_bottle_plot(
        plt_btl1,
        varname="CTDTMP1",
        f_out=cfg.dirs["logs"] + "CTDTEMP1_to31.pdf",
        cmap="RdPu",
    )
    section_bottle_plot(
        plt_btl1, varname="CTDSAL", f_out=cfg.dirs["logs"] + "CTDSAL_to31.pdf"
    )
    section_bottle_plot(
        plt_btl1, varname="CTDOXY1", f_out=cfg.dirs["logs"] + "CTDOXY_to31.pdf"
    )
    osnap_plot_TS(plt_btl1, f_out=cfg.dirs["logs"] + "TS_to31.pdf")
    #   Second section
    line2 = [
        "052",
        "053",
        "054",
        "055",
        "056",
        "057",
        "058",
        "059",
        "060",
        "061",
        "062",
        "063",
        "064",
        "065",
    ]
    plt_btl2 = btl_fit.loc[btl_fit.SSSCC.isin(line2)]
    section_bottle_plot(
        plt_btl2,
        varname="CTDTMP1",
        f_out=cfg.dirs["logs"] + "CTDTEMP1_northbox.pdf",
        cmap="RdPu",
    )
    section_bottle_plot(
        plt_btl2, varname="CTDSAL", f_out=cfg.dirs["logs"] + "CTDSAL_northbox.pdf"
    )
    section_bottle_plot(
        plt_btl2, varname="CTDOXY1", f_out=cfg.dirs["logs"] + "CTDOXY_northbox.pdf"
    )
    osnap_plot_TS(plt_btl1, f_out=cfg.dirs["logs"] + "TS_northbox.pdf")
    osnap_plot_TS(btl_fit, f_out=cfg.dirs["logs"] + "TS.pdf")
    all_residuals(btl_fit)

    print("Starting the prefit-postfit comparisons...")
    for ssscc in btl_fit.SSSCC.unique():
        pre = time_prefit[["CTDPRS", "CTDOXY"]].loc[time_prefit.SSSCC == ssscc]
        post = time_fit[["CTDPRS", "CTDOXY"]].loc[time_fit.SSSCC == ssscc]
        ref = btl_fit[["CTDPRS", "OXYGEN"]].loc[btl_fit.SSSCC == ssscc]
        fit_comparison(pre, post, ref, ssscc)


def all_residuals(btl_df, outdir="data/logs/postfit/", ext=".pdf"):

    #   Temperature residuals
    param = "t1"
    residual_vs_pressure(
        btl_df[cfg.column[param]],
        btl_df[cfg.column["t2"]],
        btl_df["CTDPRS"],
        stn=btl_df["SSSCC"],
        xlabel=f"{cfg.column[param]} Residual (mS/cm)",
        f_out=f"{outdir}t2-{param}_vs_p{ext}",
    )
    residual_vs_station(
        btl_df[cfg.column[param]],
        btl_df[cfg.column["t2"]],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        ylabel=f"{cfg.column[param]} Residual (mS/cm)",
        f_out=f"{outdir}t2-{param}_vs_stn{ext}",
    )
    #   Conductivity residuals
    for param, ref in zip(["c1", "c2", "c2"], ["refC", "refC", "c1"]):
        residual_vs_pressure(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            stn=btl_df["SSSCC"],
            xlabel=f"{cfg.column[param]} Residual (mS/cm)",
            f_out=f"{outdir}{ref}-{param}_vs_p{ext}",
        )
        residual_vs_station(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            btl_df["SSSCC"],
            ylabel=f"{cfg.column[param]} Residual (mS/cm)",
            f_out=f"{outdir}{ref}-{param}_vs_stn{ext}",
        )
        residual_vs_station(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            btl_df["SSSCC"],
            ylabel=f"{cfg.column[param]} Residual (mS/cm)",
            deep=True,
            f_out=f"{outdir}{ref}-{param}_vs_stn_deep{ext}",
        )


def osnap_plot_sensor_comparison():
    """
    Line plot of primary-secondary sensor difference at a given station or series of stations.
    Inputting either the raw (.pkl) or calibrated (.extension_tbd) continuous downcast data.
    """
    pass


def fit_comparison(pre, post, ref, ssscc, varlabel="CTDOXY"):
    """Take in a variable pre and post time data, as well as reference, and plot scatter"""
    plt.figure(figsize=(7, 6))
    ax = plt.axes()
    a = plt.scatter(pre.CTDOXY, pre.CTDPRS)
    b = plt.scatter(post.CTDOXY, post.CTDPRS)
    c = plt.scatter(ref.OXYGEN, ref.CTDPRS, s=500, marker="*")
    plt.title(str(varlabel + " " + ssscc))
    plt.tight_layout()
    plt.gca().invert_yaxis()
    plt.ylabel("CTDPRES")
    plt.legend()
    logs = "logs"
    f_out = f"{cfg.dirs[logs]}postfit/{varlabel}_{ssscc}.png"
    _save_fig(ax, f_out)


def osnap_plot_sensor_comparison_3p():
    """
    3-panel version of "osnap_plot_sensor_comparison", reading in conductivity, temperature, and derived salinity in 3 panels.
    Inputting either the raw (.pkl) or calibrated (.extension_tbd) continuous downcast data.
    """


def osnap_plot_autosal_residuals():
    """
    Scatter plot of conductivity sensor data relative to AutoSal, with uncalibrated residuals in top panel
    and calibrated residuals in bottom panel. Imports both C1 and C2 sensor data.
    """
    pass


def osnap_plot_TS(df, f_out=None):
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
    sc.set_clim(vmin=1, vmax=df.SSSCC.astype(int).max())
    ax.set_xlabel("Salinity")
    ax.set_ylabel("Temperature (ºC)")
    cbar = plt.colorbar(sc, format=tick.FormatStrFormatter("%.0f"))
    cbar.set_label("Cast Number")

    _save_fig(ax, f_out)


def osnap_plot_rho():
    """
    Line plot of density and pressure for a single cast or series of casts.
    Inputting either the raw (.pkl) or calibrated (.extension_tbd) continuous downcast data.
    """
    pass


def osnap_plot_param_calibration():
    """
    Line plot of sensor reading vs depth, with calibrated/fit sensor data shown on same axis in alternative color.
    For OSNAP, if the parameter is temperature, do not plot calibration (as there will be none)
    Inputting either the raw (.pkl) or calibrated (.extension_tbd) continuous downcast data.
    """
    pass


def osnap_oxygen_residual():
    """
    A scatter plot of oxygen residuals, similar to that of "osnap_plot_autosal_residuals".
    Only one sensor is read in, and only one panel is necessary.
    """
    pass


# TODO: more plots! what does ODV have?
# section plots (w/ PyDiva)
# parameter-parameter plots

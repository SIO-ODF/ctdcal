"""
visualize cast data for single or groups of casts
"""
import gsw
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker

from ctdcal.plotting.plot_fit import _save_fig


def param_vs_param(param1, label1, param2, label2, f_out=None, stn=None, tsT=None, tsS=None):
    """
    Scatter plots two CTD parameters against each other with support for TS plotting.

    Support for bottle or time data.

    Most TS plots have temperature on the Y axis and salinity on the X axis.

    Parameters
    ----------
    param1 : array-like of float
        The first column/array for plotting on the X axis
    label1 : string
        Label for the X axis, attributed to param1
    param2 : array-like of float
        The second column/array for plotting on the Y-axis
    label2 : string
        Label for the Y axis, attributed to param2
    f_out : String, optional
        Path and filename to save TS plot
    stn : array-like, optional
        An array of station numbers or labels for
    tsT : array-like of float, optional
        Temperature array for use in density gradient calculations
    tsS : array-like of float, optional
        Salinity array for use in density gradient calculations
    """

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot()

    if stn is not None:
        stn = np.array(stn)  # Ensure stn is a NumPy array
        sc = ax.scatter(param1, param2, c=stn.astype(int), marker="+")
        sc.set_clim(vmin=stn.astype(int).min(), vmax=stn.astype(int).max())  # Color mapping based on station number
        cbar = plt.colorbar(sc, format=ticker.FormatStrFormatter("%.0f"))
        cbar.set_label("Station/Cast Label")
    else:
        sc = ax.scatter(param1, param2, color='black', marker="+")

    ax.set_xlabel(label1)
    ax.set_ylabel(label2)

    #   TS plotting - define density contours
    if tsT is not None:
        #   Demarcate temperature, salinty, and limits for density
        smin = tsS.min() - (0.01 * tsS.min())
        smax = tsS.max() + (0.01 * tsS.max())
        tmin = tsT.min() - (0.1 * tsT.max())
        tmax = tsT.max() + (0.1 * tsT.max())
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
        dens = dens - 1000  # Convert to sigma-theta
        #   Density contours formatting
        cs = plt.contour(si, ti, dens, linestyles="dashed", colors="k")
        plt.clabel(cs, fontsize=12, inline=1, fmt="%0.1f")

    return _save_fig(ax, f_out)  # Save the figure if f_out is given else return the scatter axis

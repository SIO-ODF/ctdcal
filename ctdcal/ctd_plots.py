import logging
from pathlib import Path

import cmocean
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

#####
# The following section is a little brittle due to hardcoded names, but we'll fix
# that later. Code copy pasted from jupyter notebook.
#####


def make_cruise_report_plots(df):
    """Create and output all plots"""
    btl_t1_residuals_pressure_plot(df)
    btl_t2_residuals_pressure_plot(df)
    t1_t2_residuals_pressure_plot(df)
    btl_t1_residuals_station_plot(df)
    btl_t2_residuals_station_plot(df)
    t1_t2_residuals_station_plot(df)
    btl_t1_residuals_station_deep_plot(df)
    btl_t2_residuals_station_deep_plot(df)
    t1_t2_residuals_station_deep_plot(df)
    btl_c1_residuals_pressure_plot(df)
    btl_c2_residuals_pressure_plot(df)
    c1_c2_residuals_pressure_plot(df)
    btl_c1_residuals_station_plot(df)
    btl_c2_residuals_station_plot(df)
    c1_c2_residuals_station_plot(df)
    btl_c1_residuals_station_deep_plot(df)
    btl_c2_residuals_station_deep_plot(df)
    c1_c2_residuals_station_deep_plot(df)
    c_t_coherence_plot(df)
    btl_c1_residuals_compare_plot(df)
    btl_c2_residuals_compare_plot(df)
    c1_c2_residuals_compare_plot(df)
    btl_c1_residuals_station_uncorrected_plot(df)
    btl_c2_residuals_station_uncorrected_plot(df)
    c1_c2_residuals_station_uncorrected_plot(df)
    btl_sal_pressure_plot(df)
    btl_sal_station_plot(df)
    btl_sal_station_deep_plot(df)
    btl_oxy_residuals_pressure_plot(df)
    btl_oxy_residuals_station_plot(df)
    btl_oxy_residuals_station_deep_plot(df)
    btl_oxy_residuals_temperature_plot(df)
    btl_oxy_residuals_station_temperature_plot(df)
    btl_oxy_residuals_station_deep_temperature_plot(df)
    btl_oxy_residuals_pressure_concentration_plot(df)
    btl_oxy_residuals_station_concentration_plot(df)
    return None

    #################################################################
    ##### Here lies the temperature plots, long may they rest.  #####
    #################################################################


def btl_t1_residuals_pressure_plot(reft_vals, t1_vals, press, stnno):

    reft_t1 = reft_vals - t1_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(reft_t1, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-0.02, 0.02)
    ax.set_title("REFTMP-CTDTMP1 vs CTDPRS")
    ax.set_xlabel("T1 Residual (T90 C)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/reftmp_t1_p.svg", format="svg")
    fig.savefig("./data/images/reftmp_t1_p.pdf", format="pdf")
    plt.close()
    return None


def btl_t2_residuals_pressure_plot(reft_vals, t2_vals, press, stnno):

    reft_t2 = reft_vals - t2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(reft_t2, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-0.02, 0.02)
    ax.set_title("REFTMP-CTDTMP2 vs CTDPRS")
    ax.set_xlabel("T2 Residual (T90 C)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/reftmp_t2_p.svg", format="svg")
    fig.savefig("./data/images/reftmp_t2_p.pdf", format="pdf")
    plt.close()
    return None


def t1_t2_residuals_pressure_plot(t1_vals, t2_vals, press, stnno):

    t1_t2 = t1_vals - t2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(t1_t2, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-0.02, 0.02)
    ax.set_title("CTDTMP1-CTDTMP2 vs CTDPRS")
    ax.set_xlabel("T1-T2 Residual (T90 C)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/t1_t2_p.svg", format="svg")
    fig.savefig("./data/images/t1_t2_p.pdf", format="pdf")
    plt.close()
    return None


def btl_t1_residuals_station_plot(reft_vals, t1_vals, press, stnno):

    reft_t1 = reft_vals - t1_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, reft_t1, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("REFTMP-CTDTMP1 vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("T1 Residual (T90 C)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/reftmp_t1_stn.svg", format="svg")
    fig.savefig("./data/images/reftmp_t1_stn.pdf", format="pdf")
    plt.close()
    return None


def btl_t2_residuals_station_plot(t1_vals, t2_vals, press, stnno):

    t1_t2 = t1_vals - t2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, t1_t2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("REFTMP-CTDTMP2 vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("T2 Residual (T90 C)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/reftmp_t2_stn.svg", format="svg")
    fig.savefig("./data/images/reftmp_t2_stn.pdf", format="pdf")
    plt.close()
    return None


def t1_t2_residuals_station_plot(t1_vals, t2_vals, press, stnno):

    t1_t2 = t1_vals - t2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, t1_t2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("CTDTMP1-CTDTMP2 vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("T1-T2 Residual (T90 C)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/t1_t2_stn.svg", format="svg")
    fig.savefig("./data/images/t1_t2_stn.pdf", format="pdf")
    plt.close()
    return None


def btl_t1_residuals_station_deep_plot(reft_vals, t1_vals, press, stnno):

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["REFT_T1"] = reft_vals - t1_vals
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["REFT_T1"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("REFTMP-CTDTMP1 (>2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("T1 Residual (T90 C)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/reftmp_t1_stn_deep.svg", format="svg")
    fig.savefig("./data/images/reftmp_t1_stn_deep.pdf", format="pdf")
    plt.close()
    return None


def btl_t2_residuals_station_deep_plot(reft_vals, t2_vals, press, stnno):

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["REFT_T2"] = reft_vals - t2_vals
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["REFT_T2"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("REFTMP-CTDTMP2 (>2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("T2 Residual (T90 C)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/reftmp_t2_stn_deep.svg", format="svg")
    fig.savefig("./data/images/reftmp_t2_stn_deep.pdf", format="pdf")
    plt.close()
    return None


def t1_t2_residuals_station_deep_plot(t1_vals, t2_vals, press, stnno):

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["T1_T2"] = t1_vals - t2_vals
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["T1_T2"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("CTDTMP1-CTDTMP2 (>2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("T1-T2 Residual (T90 C)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/t1_t2_stn_deep.svg", format="svg")
    fig.savefig("./data/images/t1_t2_stn_deep.pdf", format="pdf")
    plt.close()
    return None

    #################################################################
    ##### Here lies the conductivity plots, long may they rest. #####
    #################################################################


def btl_c1_residuals_pressure_plot(refc_vals, c1_vals, press, stnno):

    refc_c1 = refc_vals - c1_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(refc_c1, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-0.02, 0.02)
    ax.set_title("BTLCOND-CTDCOND1 vs CTDPRS")
    ax.set_xlabel("C1 Residual (mS/cm)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/btlcond_c1_p.svg", format="svg")
    fig.savefig("./data/images/btlcond_c1_p.pdf", format="pdf")
    plt.close()
    return None


def btl_c2_residuals_pressure_plot(refc_vals, c2_vals, press, stnno):

    refc_c2 = refc_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(refc_c2, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_title("BTLCOND-CTDCOND2 vs CTDPRS")
    ax.set_xlabel("C2 Residual (mS/cm)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/btlcond_c2_p.svg", format="svg")
    fig.savefig("./data/images/btlcond_c2_p.pdf", format="pdf")
    plt.close()
    return None


def c1_c2_residuals_pressure_plot(c1_vals, c2_vals, press, stnno):

    c1_c2 = c1_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(c1_c2, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-0.02, 0.02)
    ax.set_title("CTDCOND1-CTDCOND2 vs CTDPRS")
    ax.set_xlabel("C1-C2 Residual (mS/cm)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/c1_c2_p.svg", format="svg")
    fig.savefig("./data/images/c1_c2_p.pdf", format="pdf")
    plt.close()
    return None


def btl_c1_residuals_station_plot(refc_vals, c1_vals, press, stnno):

    refc_c1 = refc_vals - c1_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, refc_c1, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND-CTDCOND1 vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C1 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c1_stn.svg", format="svg")
    fig.savefig("./data/images/btlcond_c1_stn.pdf", format="pdf")
    plt.close()
    return None


def btl_c2_residuals_station_plot(refc_vals, c2_vals, press, stnno):

    refc_c2 = refc_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, refc_c2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND-CTDCOND2 vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c2_stn.svg", format="svg")
    fig.savefig("./data/images/btlcond_c2_stn.pdf", format="pdf")
    plt.close()
    return None


def c1_c2_residuals_station_plot(c1_vals, c2_vals, press, stnno):

    c1_c2 = c1_vals - c2_vals
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, c1_c2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("CTDCOND1-CTDCOND2 vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C1-C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/c1_c2_stn.svg", format="svg")
    fig.savefig("./data/images/c1_c2_stn.pdf", format="pdf")
    plt.close()
    return None


def btl_c1_residuals_station_deep_plot(refc_vals, c1_vals, press, stnno):

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["REFT_C1"] = refc_vals - c1_vals
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["REFT_C1"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND-CTDCOND1 (>2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C1 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c1_stn_deep.svg", format="svg")
    fig.savefig("./data/images/btlcond_c1_stn_deep.pdf", format="pdf")
    plt.close()
    return None


def btl_c2_residuals_station_deep_plot(refc_vals, c2_vals, press, stnno):

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["REFT_C2"] = refc_vals - c2_vals
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["REFT_C2"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND-CTDCOND2 (>2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c2_stn_deep.svg", format="svg")
    fig.savefig("./data/images/btlcond_c2_stn_deep.pdf", format="pdf")
    plt.close()
    return None


def c1_c2_residuals_station_deep_plot(c1_vals, c2_vals, press, stnno):

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["C1_C2"] = c1_vals - c2_vals
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["C1_C2"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("CTDCOND1-CTDCOND2 (>2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C1-C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/c1_c2_stn_deep.svg", format="svg")
    fig.savefig("./data/images/c1_c2_stn_deep.pdf", format="pdf")
    plt.close()
    return None


def c_t_coherence_plot(t1_vals, t2_vals, c1_vals, c2_vals, press):

    t1_t2 = t1_vals - t2_vals
    c1_c2 = c1_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(t1_t2, c1_c2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_xlim(-0.02, 0.02)
    ax.set_title("T1-T2 vs C1-C2")
    ax.set_xlabel("T1-T2 Residual (T90 C)")
    ax.set_ylabel("C1-C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/c_t_coherence_p.svg", format="svg")
    fig.savefig("./data/images/c_t_coherence_p.pdf", format="pdf")
    plt.close()
    return None


def btl_c1_residuals_compare_plot(refc_vals, c1_vals, press):

    refc_c1 = refc_vals - c1_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(refc_vals, refc_c1, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND vs BTLCOND-CTDCOND1")
    ax.set_xlabel("Reference Conductivity (mS/cm)")
    ax.set_ylabel("C1 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c1_compare.svg", format="svg")
    fig.savefig("./data/images/btlcond_c1_compare.pdf", format="pdf")
    plt.close()
    return None


def btl_c2_residuals_compare_plot(refc_vals, c2_vals, press):

    refc_c2 = refc_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(refc_vals, refc_c2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND vs BTLCOND-CTDCOND2")
    ax.set_xlabel("Reference Conductivity (mS/cm)")
    ax.set_ylabel("C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c2_compare.svg", format="svg")
    fig.savefig("./data/images/btlcond_c2_compare.pdf", format="pdf")
    plt.close()
    return None


def c1_c2_residuals_compare_plot(refc_vals, c1_vals, c2_vals, press):

    c1_c2 = c1_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(refc_vals, c1_c2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND vs CTDCOND1-CTDCOND2")
    ax.set_xlabel("Reference Conductivity (mS/cm)")
    ax.set_ylabel("C1-C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/c1_c2_compare.svg", format="svg")
    fig.savefig("./data/images/c1_c2_compare.pdf", format="pdf")
    plt.close()
    return None


def btl_c1_residuals_station_uncorrected_plot(refc_vals, c1_vals, press, stnno):

    refc_c1 = refc_vals - c1_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, refc_c1, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND-CTDCOND1 (Uncorrected) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C1 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c1_stn_uncorrected.svg", format="svg")
    fig.savefig("./data/images/btlcond_c1_stn_uncorrected.pdf", format="pdf")
    plt.close()
    return None


def btl_c2_residuals_station_uncorrected_plot(refc_vals, c2_vals, press, stnno):

    refc_c2 = refc_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, refc_c2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("BTLCOND-CTDCOND2 (Uncorrected) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlcond_c2_stn_uncorrected.svg", format="svg")
    fig.savefig("./data/images/btlcond_c2_stn_uncorrected.pdf", format="pdf")
    plt.close()
    return None


def c1_c2_residuals_station_uncorrected_plot(c1_vals, c2_vals, press, stnno):

    c1_c2 = c1_vals - c2_vals

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, c1_c2, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("CTDCOND1-CTDCOND2 (Uncorrected) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("C1-C2 Residual (mS/cm)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/c1_c2_stn_uncorrected.svg", format="svg")
    fig.savefig("./data/images/c1_c2_stn_uncorrected.pdf", format="pdf")
    plt.close()
    return None


def btl_sal_pressure_plot(btl_sal, ctd_sal, press, stnno):
    sal_res = btl_sal - ctd_sal
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(sal_res, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-0.02, 0.02)
    ax.set_title("SALNTY-CTDSAL vs CTDPRS")
    ax.set_xlabel("CTDSAL Residual (mPSU)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/btlsal_sal_p.svg", format="svg")
    fig.savefig("./data/images/btlsal_sal_p.pdf", format="pdf")
    plt.close()
    return None


def btl_sal_station_plot(btl_sal, ctd_sal, press, stnno):
    sal_res = btl_sal - ctd_sal
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, sal_res, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("SALNTY-CTDSAL vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("CTDSAL Residual (mPSU)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlsal_sal_stn.svg", format="svg")
    fig.savefig("./data/images/btlsal_sal_stn.pdf", format="pdf")
    plt.close()
    return None


def btl_sal_station_deep_plot(btl_sal, ctd_sal, press, stnno):
    sal_res = btl_sal - ctd_sal

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["BTL_SAL"] = sal_res
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["BTL_SAL"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-0.01, 0.01)
    ax.set_title("SALNTY-CTDSAL (>2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("CTDSAL Residual (mPSU)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btlsal_sal_stn_deep.svg", format="svg")
    fig.savefig("./data/images/btlsal_sal_stn_deep.pdf", format="pdf")
    plt.close()
    return None

    #################################################################
    ######## Here lies the oxygen plots, long may they rest. ########
    #################################################################


def btl_oxy_residuals_pressure_plot(ref_oxy, ctdoxy, press, stnno):

    btl_o = ref_oxy - ctdoxy

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(btl_o, -press, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY vs CTDPRS")
    ax.set_xlabel("CTDOXY Residual (umol/kg)")
    ax.set_ylabel("Pressure (dbar)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/btl_oxy_p.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_p.pdf", format="pdf")
    plt.close()
    return None


def btl_oxy_residuals_station_plot(ref_oxy, ctdoxy, press, stnno):

    btl_o = ref_oxy - ctdoxy

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, btl_o, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_ylim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("CTDOXY Residual (umol/kg)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btl_oxy_stn.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_stn.pdf", format="pdf")
    plt.close()
    return None


def btl_oxy_residuals_station_deep_plot(ref_oxy, ctdoxy, press, stnno):

    df = pd.DataFrame()
    df["CTDPRS"] = press
    df["BTL_O"] = ref_oxy - ctdoxy
    df["STNNBR"] = stnno

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["BTL_O"],
        marker="+",
        c=df_deep["CTDPRS"],
        cmap=plt.cm.viridis_r,
    )
    ax.set_ylim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY (> 2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("CTDOXY Residual (umol/kg)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btl_oxy_stn_deep.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_stn_deep.pdf", format="pdf")
    plt.close()
    return None


def btl_oxy_residuals_temperature_plot(ref_oxy, ctd_oxy, t1_vals, stnno):

    btl_o = ref_oxy - ctd_oxy

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(btl_o, t1_vals, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY vs CTDPRS")
    ax.set_xlabel("CTDOXY Residual (umol/kg)")
    ax.set_ylabel("Temperature (degrees C)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/btl_oxy_t.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_t.pdf", format="pdf")
    plt.close()
    return None


def btl_oxy_residuals_station_temperature_plot(ref_oxy, ctd_oxy, t1_vals, stnno):

    btl_o = ref_oxy - ctd_oxy

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(stnno, btl_o, marker="+", c=t1_vals, cmap=cmocean.cm.thermal)
    ax.set_ylim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("CTDOXY Residual (umol/kg)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Temperature (degrees C)")

    fig.savefig("./data/images/btl_oxy_stn_t.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_stn_t.pdf", format="pdf")
    plt.close()
    return None


def btl_oxy_residuals_station_deep_temperature_plot(
    ref_oxy, ctdoxy, t1_vals, press, stnno
):

    df = pd.DataFrame()
    df["CTDTMP1"] = t1_vals
    df["BTL_O"] = ref_oxy - ctdoxy
    df["STNNBR"] = stnno
    df["CTDPRS"] = press

    df_deep = df[df["CTDPRS"] > 2000]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(
        df_deep["STNNBR"],
        df_deep["BTL_O"],
        marker="+",
        c=df_deep["CTDTMP1"],
        cmap=cmocean.cm.thermal,
    )
    ax.set_ylim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY (> 2000 db) vs STNNBR")
    ax.set_xlabel("Station Number")
    ax.set_ylabel("CTDOXY Residual (umol/kg)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Temperature (degrees C)")

    fig.savefig("./data/images/btl_oxy_stn_deep_t.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_stn_deep_t.pdf", format="pdf")
    plt.close()
    return None


def btl_oxy_residuals_pressure_concentration_plot(ref_oxy, ctd_oxy, stnno):

    btl_o = ref_oxy - ctd_oxy

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(btl_o, ref_oxy, marker="+", c=stnno, cmap=plt.cm.tab20c_r)
    ax.set_xlim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY vs OXYGEN")
    ax.set_xlabel("CTDOXY Residual (umol/kg)")
    ax.set_ylabel("Dissolved Oxygen (umol/kg)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Station Number")

    fig.savefig("./data/images/btl_oxy_p_concentration.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_p_concentration.pdf", format="pdf")
    plt.close()
    return None


def btl_oxy_residuals_station_concentration_plot(ref_oxy, ctd_oxy, press, stnno):

    btl_o = ref_oxy - ctd_oxy

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    cm = ax.scatter(btl_o, ref_oxy, marker="+", c=press, cmap=plt.cm.viridis_r)
    ax.set_xlim(-10, 10)
    ax.set_title("OXYGEN-CTDOXY vs OXYGEN")
    ax.set_xlabel("CTDOXY Residual (umol/kg)")
    ax.set_ylabel("Dissolved Oxygen (umol/kg)")
    cbar = fig.colorbar(cm)
    cbar.set_label("Pressure (dbar)")

    fig.savefig("./data/images/btl_oxy_stn_concentration.svg", format="svg")
    fig.savefig("./data/images/btl_oxy_stn_concentration.pdf", format="pdf")
    plt.close()
    return None


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

    title = f"{ref.name}-{param.name} vs. {prs.name}"
    diff = ref - param
    if deep:
        title = f"{ref.name}-{param.name} (>2000 dbar) vs. {prs.name}"
        deep_rows = prs > 2000
        diff = diff[deep_rows]
        prs = prs[deep_rows]
        stn = stn[deep_rows]

    idx, uniques = stn.factorize()  # find unique stations #s and index them

    plt.figure(figsize=(7, 6))
    plt.scatter(diff, prs, c=idx, marker="+")
    plt.xlim(xlim)
    plt.xticks(rotation=45)
    plt.ylim(ylim)
    cbar = plt.colorbar(pad=0.1)  # set cbar ticks to station names
    if not uniques.empty:
        tick_inds = cbar.get_ticks().astype(int)
        cbar.ax.yaxis.set_major_locator(ticker.FixedLocator(tick_inds))
        cbar.ax.set_yticklabels(uniques[tick_inds])
    cbar.set_label("Station Number")
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

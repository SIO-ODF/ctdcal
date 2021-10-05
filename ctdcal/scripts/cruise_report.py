import logging

import numpy as np
import pandas as pd

from .. import ctd_plots, get_ctdcal_config

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

btl_file = "data/scratch_folder/report_data.csv"
btl_df = pd.read_csv(btl_file)

# TODO: these were in old plots module; worth adding?
# (BTLCOND - C1) vs. BTLCOND
# (BTLCOND - C2) vs. BTLCOND
# (C1 - C2) vs. BTLCOND
# (BTLCOND - C1_uncorrected) vs. STNNBR
# (BTLCOND - C2_uncorrected) vs. STNNBR
# (C1_uncorrected - C2_uncorrected) vs. STNNBR
# T1 vs. (OXYGEN - CTDOXY)
# (OXYGEN - CTDOXY) vs. STNNBR (c=T1)
# deep(OXYGEN - CTDOXY) vs. STNNBR (c=T1)
# OXYGEN vs. (OXYGEN - CTDOXY) (c=STNNBR)
# OXYGEN vs. (OXYGEN - CTDOXY) (c=CTDPRS)


def plot_residuals(outdir="data/report_figs/", ext=".pdf"):

    #################################################################
    ##### Here lies the temperature plots, long may they rest.  #####
    #################################################################
    log.info("Generating temperature residual plots")
    for param, ref in zip(["t1", "t2", "t2"], ["refT", "refT", "t1"]):
        ctd_plots.residual_vs_pressure(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            stn=btl_df["STNNBR"],
            xlabel=f"{cfg.column[param]} Residual (T90 C)",
            f_out=f"{outdir}{ref}-{param}_vs_p{ext}",
        )
        ctd_plots.residual_vs_station(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            btl_df["STNNBR"],
            ylabel=f"{cfg.column[param]} Residual (T90 C)",
            f_out=f"{outdir}{ref}-{param}_vs_stn{ext}",
        )
        ctd_plots.residual_vs_station(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            btl_df["STNNBR"],
            ylabel=f"{cfg.column[param]} Residual (T90 C)",
            deep=True,
            f_out=f"{outdir}{ref}-{param}_vs_stn_deep{ext}",
        )

    #################################################################
    ##### Here lies the conductivity plots, long may they rest. #####
    #################################################################
    log.info("Generating conductivity residual plots")
    for param, ref in zip(["c1", "c2", "c2"], ["refC", "refC", "c1"]):
        ctd_plots.residual_vs_pressure(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            stn=btl_df["STNNBR"],
            xlabel=f"{cfg.column[param]} Residual (mS/cm)",
            f_out=f"{outdir}{ref}-{param}_vs_p{ext}",
        )
        ctd_plots.residual_vs_station(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            btl_df["STNNBR"],
            ylabel=f"{cfg.column[param]} Residual (mS/cm)",
            f_out=f"{outdir}{ref}-{param}_vs_stn{ext}",
        )
        ctd_plots.residual_vs_station(
            btl_df[cfg.column[param]],
            btl_df[cfg.column[ref]],
            btl_df["CTDPRS"],
            btl_df["STNNBR"],
            ylabel=f"{cfg.column[param]} Residual (mS/cm)",
            deep=True,
            f_out=f"{outdir}{ref}-{param}_vs_stn_deep{ext}",
        )

    # coherence plot doesn't have its own function...
    import matplotlib.pyplot as plt

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

    # salinity plots
    log.info("Generating salinity residual plots")
    ctd_plots.residual_vs_pressure(
        btl_df["CTDSAL"],
        btl_df["SALNTY"],
        btl_df["CTDPRS"],
        stn=btl_df["STNNBR"],
        xlabel="CTDSAL Residual (PSU)",
        f_out=f"{outdir}btlsal-sal_vs_p{ext}",
    )
    ctd_plots.residual_vs_station(
        btl_df["CTDSAL"],
        btl_df["SALNTY"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        ylabel="CTDSAL Residual (PSU)",
        f_out=f"{outdir}btlsal-sal_vs_stn{ext}",
    )
    ctd_plots.residual_vs_station(
        btl_df["CTDSAL"],
        btl_df["SALNTY"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        ylabel="CTDSAL Residual (PSU)",
        deep=True,
        f_out=f"{outdir}btlsal-sal_vs_stn_deep{ext}",
    )
    #################################################################
    ######## Here lies the oxygen plots, long may they rest. ########
    #################################################################
    # SBE43 oxygen plots
    log.info("Generating oxygen (SBE43) residual plots")
    ctd_plots.residual_vs_pressure(
        btl_df["CTDOXY"],
        btl_df["OXYGEN"],
        btl_df["CTDPRS"],
        stn=btl_df["STNNBR"],
        xlim=(-10, 10),
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f"{outdir}oxy-43_vs_p{ext}",
    )
    ctd_plots.residual_vs_station(
        btl_df["CTDOXY"],
        btl_df["OXYGEN"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        ylim=(-10, 10),
        ylabel="CTDOXY Residual (umol/kg)",
        f_out=f"{outdir}oxy-43_vs_stn{ext}",
    )
    ctd_plots.residual_vs_station(
        btl_df["CTDOXY"],
        btl_df["OXYGEN"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        deep=True,
        ylim=(-10, 10),
        ylabel="CTDOXY Residual (umol/kg)",
        f_out=f"{outdir}oxy-43_vs_stn_deep{ext}",
    )

    # RINKO oxygen plots
    log.info("Generating oxygen (RINKO) residual plots")
    ctd_plots.residual_vs_pressure(
        btl_df["CTDRINKO"],
        btl_df["OXYGEN"],
        btl_df["CTDPRS"],
        stn=btl_df["STNNBR"],
        xlim=(-10, 10),
        xlabel="CTDRINKO Residual (umol/kg)",
        f_out=f"{outdir}oxy-rinko_vs_p{ext}",
    )
    ctd_plots.residual_vs_station(
        btl_df["CTDRINKO"],
        btl_df["OXYGEN"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        ylim=(-10, 10),
        ylabel="CTDRINKO Residual (umol/kg)",
        f_out=f"{outdir}oxy-rinko_vs_stn{ext}",
    )
    ctd_plots.residual_vs_station(
        btl_df["CTDRINKO"],
        btl_df["OXYGEN"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        deep=True,
        ylim=(-10, 10),
        ylabel="CTDRINKO Residual (umol/kg)",
        f_out=f"{outdir}oxy-rinko_vs_stn_deep{ext}",
    )


def pressure_offset():

    data = pd.read_csv("data/logs/ondeck_pressure.csv")
    print(f"Average deck pressure:\n{data.describe().loc[['min', 'max', 'mean']]}\n")
    print(
        f"Average offset:\n{(data['ondeck_end_p'] - data['ondeck_start_p']).describe()}\n"
    )


def fit_coefficients():
    """Write code to build table for fit groups?"""
    pass


def calculate_residuals():

    log.info("Calculating T/C/O residuals\n")

    low_grad_rows = (btl_df["CTDTMP1"] - btl_df["CTDTMP2"]).abs() < 0.002
    deep_rows = btl_df["CTDPRS"] > 2000

    T_flag2 = (btl_df["CTDTMP1_FLAG_W"] == 2) | (btl_df["CTDTMP2_FLAG_W"] == 2)
    C_flag2 = (btl_df["CTDCOND1_FLAG_W"] == 2) | (btl_df["CTDCOND2_FLAG_W"] == 2)
    S_flag2 = (btl_df["CTDSAL_FLAG_W"] == 2) | (btl_df["CTDSAL_FLAG_W"] == 2)
    O_flag2 = btl_df["CTDOXY_FLAG_W"] == 2
    R_flag2 = btl_df["CTDRINKO_FLAG_W"] == 2

    def fmt(x):
        return np.format_float_positional(x, precision=5)

    # Temperature
    low_grad = btl_df[low_grad_rows & T_flag2]
    deep = btl_df[deep_rows & T_flag2]
    print(f"REFT-T1 = {fmt((low_grad['REFTMP'] - low_grad['CTDTMP1']).std() * 2)}")
    print(f"REFT-T2 = {fmt((low_grad['REFTMP'] - low_grad['CTDTMP2']).std() * 2)}")
    print(f"T1-T2 = {fmt((low_grad['CTDTMP1'] - low_grad['CTDTMP2']).std() * 2)}")
    print(f"REFT-T1 (deep) = {fmt((deep['REFTMP'] - deep['CTDTMP1']).std() * 2)}")
    print(f"REFT-T2 (deep) = {fmt((deep['REFTMP'] - deep['CTDTMP2']).std() * 2)}")
    print(f"T1-T2 (deep) = {fmt((deep['CTDTMP1'] - deep['CTDTMP2']).std() * 2)}")
    print("")

    # Conductivity
    low_grad = btl_df[low_grad_rows & C_flag2]
    deep = btl_df[deep_rows & C_flag2]
    print(f"REFC-C1 = {fmt((low_grad['BTLCOND'] - low_grad['CTDCOND1']).std() * 2)}")
    print(f"REFC-C2 = {fmt((low_grad['BTLCOND'] - low_grad['CTDCOND2']).std() * 2)}")
    print(f"C1-C2 = {fmt((low_grad['CTDCOND1'] - low_grad['CTDCOND2']).std() * 2)}")
    print(f"REFC-C1 (deep) = {fmt((deep['BTLCOND'] - deep['CTDCOND1']).std() * 2)}")
    print(f"REFC-C2 (deep) = {fmt((deep['BTLCOND'] - deep['CTDCOND2']).std() * 2)}")
    print(f"C1-C2 (deep) = {fmt((deep['CTDCOND1'] - deep['CTDCOND2']).std() * 2)}")
    print("")

    # Salinity
    low_grad = btl_df[low_grad_rows & S_flag2]
    deep = btl_df[deep_rows & S_flag2]
    print(f"REFS-SAL = {fmt((low_grad['SALNTY'] - low_grad['CTDSAL']).std() * 2)}")
    print(f"REFS-SAL (deep) = {fmt((deep['SALNTY'] - deep['CTDSAL']).std() * 2)}")
    print("")

    # SBE43
    # low_grad = btl_df[low_grad_rows & O_flag2]
    flag2 = btl_df[O_flag2]
    deep = btl_df[deep_rows & O_flag2]
    print(f"OXY-SBE43 = {fmt((flag2['OXYGEN'] - flag2['CTDOXY']).std() * 2)}")
    print(f"OXY-SBE43 (deep) = {fmt((deep['OXYGEN'] - deep['CTDOXY']).std() * 2)}")
    print("")

    # RINKO
    # low_grad = btl_df[low_grad_rows & R_flag2]
    flag2 = btl_df[R_flag2]
    deep = btl_df[deep_rows & R_flag2]
    print(f"OXY-RINKO = {fmt((flag2['OXYGEN'] - flag2['CTDRINKO']).std() * 2)}")
    print(f"OXY-RINKO (deep) = {fmt((deep['OXYGEN'] - deep['CTDRINKO']).std() * 2)}")
    print("")


def cruise_report_residuals():
    """Do all the things for cruise report"""
    plot_residuals()
    pressure_offset()
    fit_coefficients()
    calculate_residuals()


# if __name__ == "__main__":
#     make_cruise_report_plots()

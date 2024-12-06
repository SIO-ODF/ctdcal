import logging
from pathlib import Path

import numpy as np
import pandas as pd

from ctdcal import ctd_plots
from ctdcal.common import load_user_config, validate_file, validate_dir

log = logging.getLogger(__name__)

USERCONFIG = '/Users/als026/data/ices/ices.yaml'
cfg = load_user_config(validate_file(USERCONFIG))
INST = 'ctd'

btl_file = Path(cfg.datadir, 'report', INST, "report_data.csv")
BTLDF = pd.read_csv(btl_file)


def plot_residuals(outdir, ext=".pdf"):

    #################################################################
    ##### Here lies the temperature plots, long may they rest.  #####
    #################################################################
    outdir = validate_dir(outdir, create=True)
    log.info("Generating temperature residual plots")
    for param, ref in zip(["CTDTMP1", "CTDTMP2", "CTDTMP2"], ["REFTMP", "REFTMP", "CTDTMP1"]):
        ctd_plots.residual_vs_pressure(
            BTLDF[param],
            BTLDF[ref],
            BTLDF["CTDPRS"],
            stn=BTLDF["CAST"],
            xlabel=f"{param} Residual (T90 C)",
            f_out=Path(outdir, f"{ref}-{param}_vs_p{ext}"),
        )
        ctd_plots.residual_vs_station(
            BTLDF[param],
            BTLDF[ref],
            BTLDF["CTDPRS"],
            BTLDF["CAST"],
            ylabel=f"{param} Residual (T90 C)",
            f_out=Path(outdir, f"{ref}-{param}_vs_stn{ext}"),
        )
        ctd_plots.residual_vs_station(
            BTLDF[param],
            BTLDF[ref],
            BTLDF["CTDPRS"],
            BTLDF["CAST"],
            ylabel=f"{param} Residual (T90 C)",
            deep=True,
            f_out=Path(outdir, f"{ref}-{param}_vs_stn_deep{ext}"),
        )

    #################################################################
    ##### Here lies the conductivity plots, long may they rest. #####
    #################################################################
    log.info("Generating conductivity residual plots")
    for param, ref in zip(["CTDCOND1", "CTDCOND2", "CTDCOND2"], ["BTLCOND", "BTLCOND", "CTDCOND1"]):
        ctd_plots.residual_vs_pressure(
            BTLDF[param],
            BTLDF[ref],
            BTLDF["CTDPRS"],
            stn=BTLDF["CAST"],
            xlabel=f"{param} Residual (mS/cm)",
            f_out=Path(outdir, f"{ref}-{param}_vs_p{ext}"),
        )
        ctd_plots.residual_vs_station(
            BTLDF[param],
            BTLDF[ref],
            BTLDF["CTDPRS"],
            BTLDF["CAST"],
            ylabel=f"{param} Residual (mS/cm)",
            f_out=Path(outdir, f"{ref}-{param}_vs_stn{ext}"),
        )
        ctd_plots.residual_vs_station(
            BTLDF[param],
            BTLDF[ref],
            BTLDF["CTDPRS"],
            BTLDF["CAST"],
            ylabel=f"{param} Residual (mS/cm)",
            deep=True,
            f_out=Path(outdir, f"{ref}-{param}_vs_stn_deep{ext}"),
        )

    # coherence plot doesn't have its own function...
    import matplotlib.pyplot as plt

    plt.figure(figsize=(7, 6))
    plt.scatter(
            BTLDF["CTDTMP1"] - BTLDF["CTDTMP2"],
            BTLDF["CTDCOND1"] - BTLDF["CTDCOND2"],
        c=BTLDF["CTDPRS"],
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
    plt.savefig(Path(outdir, f"c_t_coherence{ext}"))
    plt.close()

    # salinity plots
    log.info("Generating salinity residual plots")
    ctd_plots.residual_vs_pressure(
        BTLDF["CTDSAL"],
        BTLDF["SALNTY"],
        BTLDF["CTDPRS"],
        stn=BTLDF["CAST"],
        xlabel="CTDSAL Residual (PSU)",
        f_out=Path(outdir, f"btlsal-sal_vs_p{ext}"),
    )
    ctd_plots.residual_vs_station(
        BTLDF["CTDSAL"],
        BTLDF["SALNTY"],
        BTLDF["CTDPRS"],
        BTLDF["CAST"],
        ylabel="CTDSAL Residual (PSU)",
        f_out=Path(outdir, f"btlsal-sal_vs_stn{ext}"),
    )
    ctd_plots.residual_vs_station(
        BTLDF["CTDSAL"],
        BTLDF["SALNTY"],
        BTLDF["CTDPRS"],
        BTLDF["CAST"],
        ylabel="CTDSAL Residual (PSU)",
        deep=True,
        f_out=Path(outdir, f"btlsal-sal_vs_stn_deep{ext}"),
    )
    #################################################################
    ######## Here lies the oxygen plots, long may they rest. ########
    #################################################################
    # SBE43 oxygen plots
    log.info("Generating oxygen (SBE43) residual plots")
    ctd_plots.residual_vs_pressure(
        BTLDF["CTDOXY"],
        BTLDF["OXYGEN"],
        BTLDF["CTDPRS"],
        stn=BTLDF["CAST"],
        xlim=(-10, 10),
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=Path(outdir, f"oxy-43_vs_p{ext}"),
    )
    ctd_plots.residual_vs_station(
        BTLDF["CTDOXY"],
        BTLDF["OXYGEN"],
        BTLDF["CTDPRS"],
        BTLDF["CAST"],
        ylim=(-10, 10),
        ylabel="CTDOXY Residual (umol/kg)",
        f_out=Path(outdir, f"oxy-43_vs_stn{ext}"),
    )
    ctd_plots.residual_vs_station(
        BTLDF["CTDOXY"],
        BTLDF["OXYGEN"],
        BTLDF["CTDPRS"],
        BTLDF["CAST"],
        deep=True,
        ylim=(-10, 10),
        ylabel="CTDOXY Residual (umol/kg)",
        f_out=Path(outdir, f"oxy-43_vs_stn_deep{ext}"),
    )

    # # RINKO oxygen plots
    # log.info("Generating oxygen (RINKO) residual plots")
    # ctd_plots.residual_vs_pressure(
    #     BTLDF["CTDRINKO"],
    #     BTLDF["OXYGEN"],
    #     BTLDF["CTDPRS"],
    #     stn=BTLDF["CAST"],
    #     xlim=(-10, 10),
    #     xlabel="CTDRINKO Residual (umol/kg)",
    #     f_out=f"{outdir}oxy-rinko_vs_p{ext}",
    # )
    # ctd_plots.residual_vs_station(
    #     BTLDF["CTDRINKO"],
    #     BTLDF["OXYGEN"],
    #     BTLDF["CTDPRS"],
    #     BTLDF["CAST"],
    #     ylim=(-10, 10),
    #     ylabel="CTDRINKO Residual (umol/kg)",
    #     f_out=f"{outdir}oxy-rinko_vs_stn{ext}",
    # )
    # ctd_plots.residual_vs_station(
    #     BTLDF["CTDRINKO"],
    #     BTLDF["OXYGEN"],
    #     BTLDF["CTDPRS"],
    #     BTLDF["CAST"],
    #     deep=True,
    #     ylim=(-10, 10),
    #     ylabel="CTDRINKO Residual (umol/kg)",
    #     f_out=f"{outdir}oxy-rinko_vs_stn_deep{ext}",
    # )


def pressure_offset(infile):

    data = pd.read_csv(infile)
    print(f"Average deck pressure:\n{data.describe().loc[['min', 'max', 'mean']]}\n")
    print(
        f"Average offset:\n{(data['pressure_start'] - data['pressure_end']).describe()}\n"
    )


def fit_coefficients():
    """Write code to build table for fit groups?"""
    pass


def calculate_residuals():

    log.info("Calculating T/C/O residuals\n")

    low_grad_rows = (BTLDF["CTDTMP1"] - BTLDF["CTDTMP2"]).abs() < 0.002
    deep_rows = BTLDF["CTDPRS"] > 2000

    T_flag2 = (BTLDF["CTDTMP1_FLAG_W"] == 2) | (BTLDF["CTDTMP2_FLAG_W"] == 2)
    C_flag2 = (BTLDF["CTDCOND1_FLAG_W"] == 2) | (BTLDF["CTDCOND2_FLAG_W"] == 2)
    S_flag2 = (BTLDF["CTDSAL_FLAG_W"] == 2) | (BTLDF["CTDSAL_FLAG_W"] == 2)
    O_flag2 = BTLDF["CTDOXY_FLAG_W"] == 2
    # R_flag2 = BTLDF["CTDRINKO_FLAG_W"] == 2

    def fmt(x):
        return np.format_float_positional(x, precision=5)

    # Temperature
    low_grad = BTLDF[low_grad_rows & T_flag2]
    deep = BTLDF[deep_rows & T_flag2]
    print(f"REFT-T1 = {fmt((low_grad['REFTMP'] - low_grad['CTDTMP1']).std() * 2)}")
    print(f"REFT-T2 = {fmt((low_grad['REFTMP'] - low_grad['CTDTMP2']).std() * 2)}")
    print(f"T1-T2 = {fmt((low_grad['CTDTMP1'] - low_grad['CTDTMP2']).std() * 2)}")
    print(f"REFT-T1 (deep) = {fmt((deep['REFTMP'] - deep['CTDTMP1']).std() * 2)}")
    print(f"REFT-T2 (deep) = {fmt((deep['REFTMP'] - deep['CTDTMP2']).std() * 2)}")
    print(f"T1-T2 (deep) = {fmt((deep['CTDTMP1'] - deep['CTDTMP2']).std() * 2)}")
    print("")

    # Conductivity
    low_grad = BTLDF[low_grad_rows & C_flag2]
    deep = BTLDF[deep_rows & C_flag2]
    print(f"REFC-C1 = {fmt((low_grad['BTLCOND'] - low_grad['CTDCOND1']).std() * 2)}")
    print(f"REFC-C2 = {fmt((low_grad['BTLCOND'] - low_grad['CTDCOND2']).std() * 2)}")
    print(f"C1-C2 = {fmt((low_grad['CTDCOND1'] - low_grad['CTDCOND2']).std() * 2)}")
    print(f"REFC-C1 (deep) = {fmt((deep['BTLCOND'] - deep['CTDCOND1']).std() * 2)}")
    print(f"REFC-C2 (deep) = {fmt((deep['BTLCOND'] - deep['CTDCOND2']).std() * 2)}")
    print(f"C1-C2 (deep) = {fmt((deep['CTDCOND1'] - deep['CTDCOND2']).std() * 2)}")
    print("")

    # Salinity
    low_grad = BTLDF[low_grad_rows & S_flag2]
    deep = BTLDF[deep_rows & S_flag2]
    print(f"REFS-SAL = {fmt((low_grad['SALNTY'] - low_grad['CTDSAL']).std() * 2)}")
    print(f"REFS-SAL (deep) = {fmt((deep['SALNTY'] - deep['CTDSAL']).std() * 2)}")
    print("")

    # SBE43
    # low_grad = BTLDF[low_grad_rows & O_flag2]
    flag2 = BTLDF[O_flag2]
    deep = BTLDF[deep_rows & O_flag2]
    print(f"OXY-SBE43 = {fmt((flag2['OXYGEN'] - flag2['CTDOXY']).std() * 2)}")
    print(f"OXY-SBE43 (deep) = {fmt((deep['OXYGEN'] - deep['CTDOXY']).std() * 2)}")
    print("")

    # # RINKO
    # # low_grad = BTLDF[low_grad_rows & R_flag2]
    # flag2 = BTLDF[R_flag2]
    # deep = BTLDF[deep_rows & R_flag2]
    # print(f"OXY-RINKO = {fmt((flag2['OXYGEN'] - flag2['CTDRINKO']).std() * 2)}")
    # print(f"OXY-RINKO (deep) = {fmt((deep['OXYGEN'] - deep['CTDRINKO']).std() * 2)}")
    # print("")


def cruise_report_residuals():
    """Do all the things for cruise report"""
    plot_residuals(outdir=Path(cfg.datadir, 'report', INST))
    pressure_offset(infile=Path(cfg.datadir, 'logs', "ondeck_pressure.csv"))
    # fit_coefficients()
    calculate_residuals()


# if __name__ == "__main__":
#     make_cruise_report_plots()

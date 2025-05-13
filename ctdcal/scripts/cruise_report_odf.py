import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .. import get_ctdcal_config
from ctdcal.plotting.plot_fit import residual_vs_pressure, residual_vs_station
from ctdcal.common import validate_dir

# cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def plot_residuals(btl_df, outdir, ext=".pdf"):

    # set up outdir
    outdir = validate_dir(outdir, create=True)

    # set up params dict
    params = {
            'temperature': {
                    'param_col': ['CTDTMP1', 'CTDTMP2', 'CTDTMP2'],
                    'ref_col': ['REFTMP', 'REFTMP', 'CTDTMP1'],
                    'units': 'T90 C',
                    'lim': 0.02
            },
            'conductivity': {
                    'param_col': ['CTDCOND1', 'CTDCOND2', 'CTDCOND2'],
                    'ref_col': ['BTLCOND', 'BTLCOND', 'CTDCOND1'],
                    'units': 'mS/cm',
                    'lim': 0.02
            },
            'SBE43 oxygen': {
                    'param_col': ['CTDOXY'],
                    'ref_col': ['OXYGEN'],
                    'units': 'umol/kg',
                    'lim': 10
            },
            'Rinko oxygen': {
                    'param_col': ['CTDRINKO'],
                    'ref_col': ['OXYGEN'],
                    'units': 'umol/kg',
                    'lim': 10
            },
    }

    for param_name in params:
        param_dict = params[param_name]
        log.info("Generating %s residual plots" % param_name)
        for param_col, ref_col in zip(param_dict['param_col'], param_dict['ref_col']):
            residual_vs_pressure(
                btl_df[param_col],
                btl_df[ref_col],
                btl_df["CTDPRS"],
                stn=btl_df["STNNBR"],
                xlim=(-param_dict['lim'], param_dict['lim']),
                xlabel=f"{param_name} residual ({param_dict['units']})",
                f_out=f"{outdir}/{ref_col}-{param_col}_vs_p{ext}",
            )
            residual_vs_station(
                btl_df[param_col],
                btl_df[ref_col],
                btl_df["CTDPRS"],
                btl_df["STNNBR"],
                ylim=(-param_dict['lim'], param_dict['lim']),
                ylabel=f"{param_name} residual ({param_dict['units']})",
                f_out=f"{outdir}/{ref_col}-{param_col}_vs_stn{ext}",
            )
            residual_vs_station(
                btl_df[param_col],
                btl_df[ref_col],
                btl_df["CTDPRS"],
                btl_df["STNNBR"],
                ylim=(-param_dict['lim'], param_dict['lim']),
                ylabel=f"{param_name} residual ({param_dict['units']})",
                deep=True,
                f_out=f"{outdir}/{ref_col}-{param_col}_vs_stn_deep{ext}",
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
    plt.xlabel("CTDTMP1-CTDTMP2 residual (T90 C)", fontsize=12)
    plt.ylabel("CTDCOND1-CTDCOND2 residual (mS/cm)", fontsize=12)
    plt.title("CTDCOND1-CTDCOND2 vs. CTDTMP1-CTDTMP2", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{outdir}/c_t_coherence{ext}")
    plt.close()

    # salinity plots
    log.info("Generating salinity residual plots")
    residual_vs_pressure(
        btl_df["CTDSAL"],
        btl_df["SALNTY"],
        btl_df["CTDPRS"],
        stn=btl_df["STNNBR"],
        xlabel="CTDSAL residual (PSU)",
        f_out=f"{outdir}/btlsal-sal_vs_p{ext}",
    )
    residual_vs_station(
        btl_df["CTDSAL"],
        btl_df["SALNTY"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        ylabel="CTDSAL residual (PSU)",
        f_out=f"{outdir}/btlsal-sal_vs_stn{ext}",
    )
    residual_vs_station(
        btl_df["CTDSAL"],
        btl_df["SALNTY"],
        btl_df["CTDPRS"],
        btl_df["STNNBR"],
        ylabel="CTDSAL residual (PSU)",
        deep=True,
        f_out=f"{outdir}/btlsal-sal_vs_stn_deep{ext}",
    )


def pressure_offset(infile, pgroups):

    data = pd.read_csv(infile, dtype={'cast_id': str})
    for i, group in enumerate(pgroups):
        print("Pressure group %s" % i)
        print(
                f"Average deck pressure:\n{data.loc[data['cast_id'].isin(group)].describe().loc[['min', 'max', 'mean']]}\n"
        )
        print(f"Average offset:\n{data.loc[data['cast_id'].isin(group)]['pressure_start'].describe()}\n")


def calculate_residuals(btl_df):

    log.info("Calculating T/C/O residuals\n")

    low_grad_rows = (btl_df["CTDTMP1"] - btl_df["CTDTMP2"]).abs() < 0.002
    deep_rows = btl_df["CTDPRS"] > 2000

    t_flag2 = (btl_df["CTDTMP1_FLAG_W"] == 2) | (btl_df["CTDTMP2_FLAG_W"] == 2)
    c_flag2 = (btl_df["CTDCOND1_FLAG_W"] == 2) | (btl_df["CTDCOND2_FLAG_W"] == 2)
    s_flag2 = (btl_df["CTDSAL_FLAG_W"] == 2) | (btl_df["CTDSAL_FLAG_W"] == 2)
    o_flag2 = btl_df["CTDOXY_FLAG_W"] == 2
    r_flag2 = btl_df["CTDRINKO_FLAG_W"] == 2

    def fmt(x):
        return np.format_float_positional(x, precision=5)

    # Temperature
    low_grad = btl_df[low_grad_rows & t_flag2]
    deep = btl_df[deep_rows & t_flag2]
    print(f"REFT-T1 = {fmt((low_grad['REFTMP'] - low_grad['CTDTMP1']).std() * 2)}")
    print(f"REFT-T2 = {fmt((low_grad['REFTMP'] - low_grad['CTDTMP2']).std() * 2)}")
    print(f"T1-T2 = {fmt((low_grad['CTDTMP1'] - low_grad['CTDTMP2']).std() * 2)}")
    print(f"REFT-T1 (deep) = {fmt((deep['REFTMP'] - deep['CTDTMP1']).std() * 2)}")
    print(f"REFT-T2 (deep) = {fmt((deep['REFTMP'] - deep['CTDTMP2']).std() * 2)}")
    print(f"T1-T2 (deep) = {fmt((deep['CTDTMP1'] - deep['CTDTMP2']).std() * 2)}")
    print("")

    # Conductivity
    low_grad = btl_df[low_grad_rows & c_flag2]
    deep = btl_df[deep_rows & c_flag2]
    print(f"REFC-C1 = {fmt((low_grad['BTLCOND'] - low_grad['CTDCOND1']).std() * 2)}")
    print(f"REFC-C2 = {fmt((low_grad['BTLCOND'] - low_grad['CTDCOND2']).std() * 2)}")
    print(f"C1-C2 = {fmt((low_grad['CTDCOND1'] - low_grad['CTDCOND2']).std() * 2)}")
    print(f"REFC-C1 (deep) = {fmt((deep['BTLCOND'] - deep['CTDCOND1']).std() * 2)}")
    print(f"REFC-C2 (deep) = {fmt((deep['BTLCOND'] - deep['CTDCOND2']).std() * 2)}")
    print(f"C1-C2 (deep) = {fmt((deep['CTDCOND1'] - deep['CTDCOND2']).std() * 2)}")
    print("")

    # Salinity
    low_grad = btl_df[low_grad_rows & s_flag2]
    deep = btl_df[deep_rows & s_flag2]
    print(f"REFS-SAL = {fmt((low_grad['SALNTY'] - low_grad['CTDSAL']).std() * 2)}")
    print(f"REFS-SAL (deep) = {fmt((deep['SALNTY'] - deep['CTDSAL']).std() * 2)}")
    print("")

    # SBE43
    # low_grad = btl_df[low_grad_rows & o_flag2]
    flag2 = btl_df[o_flag2]
    deep = btl_df[deep_rows & o_flag2]
    print(f"OXY-SBE43 = {fmt((flag2['OXYGEN'] - flag2['CTDOXY']).std() * 2)}")
    print(f"OXY-SBE43 (deep) = {fmt((deep['OXYGEN'] - deep['CTDOXY']).std() * 2)}")
    print("")

    # RINKO
    # low_grad = btl_df[low_grad_rows & r_flag2]
    flag2 = btl_df[r_flag2]
    deep = btl_df[deep_rows & r_flag2]
    print(f"OXY-RINKO = {fmt((flag2['OXYGEN'] - flag2['CTDRINKO']).std() * 2)}")
    print(f"OXY-RINKO (deep) = {fmt((deep['OXYGEN'] - deep['CTDRINKO']).std() * 2)}")
    print("")


def cruise_report(report_dir, outdir, pgroups):
    """Do all the things for cruise report"""
    btl_file = Path(report_dir, "report_data.csv")
    btl_df = pd.read_csv(btl_file)

    plot_residuals(btl_df, outdir)
    pressure_offset(Path(report_dir, 'ondeck_pressure.csv'), pgroups)
    calculate_residuals(btl_df)


# if __name__ == "__main__":
#     make_cruise_report_plots()

"""
odf_quickplot

Script for generating pre/postfit figures

* Temperature primary and secondary continuous data, including residual plot
* Conductivity primary and secondary continuous data, including residual plot
* Oxygen (SBE43) and RINKO prefits -> Voltages
* Oxygen (SBE43) and RINKO postfits -> umol/kg

Written 2023 - DMB cruise
"""

import gsw
import pandas as pd
import numpy as np
import logging
from pathlib import Path
from datetime import datetime

from ctdcal import (
    get_ctdcal_config,
    io,
    odf_io,
    process_bottle,
    process_ctd,
    convert,
    ctd_plots,
)

log = logging.getLogger(__name__)
cfg = get_ctdcal_config()


def odf_quickplot(type):
    #   Import the SSSCC list
    if type == "ssscc":
        ssscc_list = process_ctd.get_ssscc_list()
        rosette = "Mixed"
    elif type == "odf-only":
        ssscc_list = process_ctd.get_ssscc_list(fname="data/ssscc/ssscc_odf.csv")
        rosette = "ODF"
    elif type == "gtc-only":
        ssscc_list = process_ctd.get_ssscc_list(fname="data/ssscc/ssscc_gtc.csv")
        rosette = "GTC"
    else:
        raw_files = Path(cfg.dirs["raw"]).glob("*.hex")
        ssscc_list = sorted([f.stem for f in raw_files])

    #   Import fitting coefs for the figure captions

    i = 0
    odf_io.printProgressBar(i, len(ssscc_list), prefix="Progress:", length=50)

    time_data_all = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all  = process_bottle.load_all_btl_files(ssscc_list)
    btl_data_all[cfg.column["refC"]] = convert.CR_to_cond(
            btl_data_all["CRavg"],
            btl_data_all["BathTEMP"],
            btl_data_all[cfg.column["t1"]],
            btl_data_all[cfg.column["p"]],
        )
    
    time_data_fit = pd.read_pickle(cfg.dirs["pressure"] + "ct1.pkl")

    #   Make a directory, where figures can be written out to (a folder called 00101)
    if not Path(cfg.dirs["figs"]).exists(): #   If the parent folder does not exist
        Path(cfg.dirs["figs"]).mkdir()

    #   For loop for each SSSCC
    for ssscc in ssscc_list:
        if not Path(cfg.dirs["figs"] + ssscc).exists():
            Path(cfg.dirs["figs"] + ssscc).mkdir()
            title_lead = f"{ssscc}: {rosette}"
            # pre = time_data_all[time_data_all.SSSCC == ssscc]
            pre_b = btl_data_all[btl_data_all.SSSCC == ssscc]
            #   Import the converted time-series .pkl
            pre = pd.read_pickle(cfg.dirs["time"] + ssscc + "_time.pkl")
            # pre = pd.read_pickle(cfg.dirs["converted"] + ssscc + ".pkl")
            #   Import the postfit ct1 file
            # post = io.load_exchange_ctd(cfg.dirs["pressure"] + ssscc + "_ct1.csv")[1]
            post = time_data_fit[time_data_fit.SSSCC == ssscc]

            #   On all postfit figures, write out text underneath for the fitting equation (if cond, add in the cc term)

            #   Plot the temperature prefit w/ residuals
            ctd_plots.two_element(
                pre.CTDTMP1,
                pre.CTDTMP2,
                pre.CTDPRS,
                ssscc,
                f_out=cfg.dirs["figs"] + ssscc + "/T-before",
            )
            ctd_plots.two_element(
                pre.CTDCOND1,
                pre.CTDCOND2,
                pre.CTDPRS,
                ssscc,
                f_out=cfg.dirs["figs"] + ssscc + "/C-before",
            )
            if type == "gtc-only":
                continue
            else:
                ctd_plots.two_element(
                    pre.CTDOXYVOLTS,
                    pre.U_DEF_poly1,
                    pre.CTDPRS,
                    ssscc,
                    f_out=cfg.dirs["figs"] + ssscc + "/oxygen-before",
                )
            ctd_plots.TCcoherence_plot(pre,outdir=cfg.dirs["figs"] + ssscc + "/TC-coherence-before",ext=".png")

            ctd_plots.conductivity_overlap(ssscc, pre_b, pre, title_lead=title_lead)

            ctd_plots.conductivity_overlap(ssscc, pre_b, pre, time_df2=post, title_lead=title_lead)
            #   Plot the temperature T1 pre-vs-post w/ residuals (if same length)

            #   Do the above using the conductivity sensor

            #   Do that using the rinko and SBE43 (voltages, then umol/kg)

        i += 1
        odf_io.printProgressBar(
            i, len(ssscc_list), prefix="Progress:", suffix=ssscc, length=50
        )

if __name__ == "__main__":

    odf_quickplot()

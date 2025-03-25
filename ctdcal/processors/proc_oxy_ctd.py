"""
Processes dissolved oxygen data from CTD sensors.
"""
import logging
import xml.etree.cElementTree as ET
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from ctdcal import get_ctdcal_config
from ctdcal.fitting.fit_common import get_node, NodeNotFoundError
from ctdcal.flagging.flag_common import nan_values
from ctdcal.common import get_ssscc_list
from ctdcal.processors.functions_oxy import oxy_ml_to_umolkg
from ctdcal.processors.proc_oxy_odf import calculate_bottle_oxygen

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


# TODO: remove underscore in front of name from functions accessed by other modules
def _get_sbe_coef(idx=0):
    """
    Get SBE oxygen coefficients from raw .xmlcon files.
    Defaults to using first station in ssscc.csv file.

    Returns the following tuple of coefficients: Soc, offset, Tau20, Tcor, E
    """
    station = get_ssscc_list()[idx]
    xmlfile = cfg.dirs["raw"] + station + ".XMLCON"

    tree = ET.parse(xmlfile)
    root_eq0 = tree.find(".//CalibrationCoefficients[@equation='0']")  # Owens-Millard
    root_eq1 = tree.find(".//CalibrationCoefficients[@equation='1']")  # SBE equation

    coefs = {c.tag: float(c.text) for c in root_eq1}
    coefs["Tcor"] = float(root_eq0.find("Tcor").text)  # only coef needed from eq0
    keep_keys = ["Soc", "offset", "Tau20", "Tcor", "E"]

    return tuple(coefs[key] for key in keep_keys)

def _PMEL_oxy_eq(coefs, inputs, cc=[1.92634e-4, -4.64803e-2]):
    """
    Modified oxygen equation for SBE 43 used by NOAA/PMEL
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E
    """
    Soc, Voff, Tau20, Tcorr, E = coefs
    oxyvolts, pressure, temp, dvdt, os = inputs
    o2 = (
        Soc
        * (
            oxyvolts
            + Voff
            + Tau20 * np.exp(cc[0] * pressure + cc[1] * (temp - 20)) * dvdt
        )
        * os
        * np.exp(Tcorr * temp)
        * np.exp((E * pressure) / (temp + 273.15))
    )

    return o2

def prepare_oxy(btl_df, time_df, ssscc_list, user_cfg, ref_node):
    """
    Calculate oxygen-related variables needed for calibration:
    sigma, oxygen solubility (OS), and bottle oxygen

    Parameters
    ----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data
    ssscc_list : list of str
        List of stations to process
    user_cfg : Munch object
        Munch dictionary of user-defined parameters
    ref_node : str
        Name of reference parameter

    Returns
    -------

    """
    # Calculate SA and CT
    btl_df["SA"] = gsw.SA_from_SP(
        btl_df[cfg.column["sal"]],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["lon"]],
        btl_df[cfg.column["lat"]],
    )
    btl_df["CT"] = gsw.CT_from_t(
        btl_df["SA"],
        btl_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p"]],
    )
    time_df["SA"] = gsw.SA_from_SP(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )
    time_df["CT"] = gsw.CT_from_t(
        time_df["SA"],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
    )

    # calculate sigma
    btl_df["sigma_btl"] = gsw.sigma0(btl_df["SA"], btl_df["CT"])
    time_df["sigma_btl"] = gsw.sigma0(time_df["SA"], time_df["CT"])

    # Calculate oxygen solubility in µmol/kg
    btl_df["OS"] = gsw.O2sol(
        btl_df["SA"],
        btl_df["CT"],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["lon"]],
        btl_df[cfg.column["lat"]],
    )
    time_df["OS"] = gsw.O2sol(
        time_df["SA"],
        time_df["CT"],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )
    # Convert CTDOXY units
    btl_df["CTDOXY"] = oxy_ml_to_umolkg(btl_df["CTDOXY1"], btl_df["sigma_btl"])
    # Calculate bottle oxygen
    # TODO: this part should be in proc_oxy_btl module
    btl_df[cfg.column["refO"]] = calculate_bottle_oxygen(
        ssscc_list,
        btl_df["SSSCC"],
        btl_df["TITR_VOL"],
        btl_df["TITR_TEMP"],
        btl_df["FLASKNO"],
    )
    btl_df[cfg.column["refO"]] = oxy_ml_to_umolkg(
        btl_df[cfg.column["refO"]], btl_df["sigma_btl"]
    )
    btl_df["OXYGEN_FLAG_W"] = nan_values(btl_df[cfg.column["refO"]])

    # Load manual OXYGEN flags
    flag_file = Path(user_cfg.datadir, "flag", user_cfg.bottleflags_man)
    oxy_flags_manual = None
    if flag_file.exists():
        try:
            oxy_flags_manual = get_node(flag_file, ref_node)
        except NodeNotFoundError:
            log.info(
                "No previously flagged values for %s found in flag file." % ref_node
            )
    else:
        log.info("No pre-existing flag file found.")

    if oxy_flags_manual is not None:
        log.info("Merging previously flagged values for %s." % ref_node)
        oxy_flags_manual_df = pd.DataFrame.from_dict(oxy_flags_manual)
        oxy_flags_manual_df = oxy_flags_manual_df.rename(
            columns={
                "cast_id": "SSSCC",
                "bottle_num": "btl_fire_num",
                "value": "OXYGEN_FLAG_W",
            }
        )
        btl_df.set_index(["SSSCC", "btl_fire_num"], inplace=True)
        btl_df.update(oxy_flags_manual_df.set_index(["SSSCC", "btl_fire_num"]))
        btl_df.reset_index(inplace=True)

    return True

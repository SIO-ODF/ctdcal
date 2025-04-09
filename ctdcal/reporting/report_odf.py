"""
Formats data and generates visuals for ODF-style cruise reports.
"""
import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from ctdcal import get_ctdcal_config
from ctdcal.common import validate_dir
from ctdcal.flagging import flag_common as flagging


cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def export_report_data(df):
    """
    Write out the data used for report generation as a csv.

    Params
    ------
    df : Pandas DataFrame
        Fit bottle data

    """
    df["STNNBR"] = [int(x[0:3]) for x in df["SSSCC"]]
    df["CTDPRS"] = df["CTDPRS"].round(1)
    cruise_report_cols = [
        "STNNBR",
        "CTDPRS",
        "CTDTMP1",
        "CTDTMP1_FLAG_W",
        "CTDTMP2",
        "CTDTMP2_FLAG_W",
        "REFTMP",
        "CTDCOND1",
        "CTDCOND1_FLAG_W",
        "CTDCOND2",
        "CTDCOND2_FLAG_W",
        "BTLCOND",
        "CTDSAL",
        "CTDSAL_FLAG_W",
        "SALNTY",
        "CTDOXY",
        "CTDOXY_FLAG_W",
        "CTDRINKO",
        "CTDRINKO_FLAG_W",
        "OXYGEN",
    ]

    # add in missing flags
    df["CTDTMP1_FLAG_W"] = flagging.by_residual(
        df["CTDTMP1"], df["REFTMP"], df["CTDPRS"]
    )
    df["CTDTMP2_FLAG_W"] = flagging.by_residual(
        df["CTDTMP1"], df["REFTMP"], df["CTDPRS"]
    )
    df["CTDCOND1_FLAG_W"] = flagging.by_residual(
        df["CTDCOND1"], df["BTLCOND"], df["CTDPRS"]
    )
    df["CTDCOND2_FLAG_W"] = flagging.by_residual(
        df["CTDCOND2"], df["BTLCOND"], df["CTDPRS"]
    )
    df["CTDOXY_FLAG_W"] = flagging.by_percent_diff(df["CTDOXY"], df["OXYGEN"])
    df["CTDRINKO_FLAG_W"] = flagging.by_percent_diff(df["CTDRINKO"], df["OXYGEN"])

    df[cruise_report_cols].to_csv("data/report_data.csv", index=False)

    return


def make_depth_log(time_df, cast_id_col='SSSCC', report_dir=cfg.dirs["logs"], threshold=80):
    """
    Create depth log file from maximum depth of each station/cast in time DataFrame.
    If rosette does not get within the threshold distance of the bottom, returns NaN.

    Parameters
    ----------
    time_df : DataFrame
        DataFrame containing continuous CTD data
    threshold : int, optional
        Maximum altimeter reading to consider cast "at the bottom" (defaults to 80)

    """
    report_dir = validate_dir(report_dir, create=True)

    df = time_df[[cast_id_col, "CTDPRS", "GPSLAT", "ALT"]].copy().reset_index()
    df_group = df.groupby(cast_id_col, sort=False)
    idx_p_max = df_group["CTDPRS"].idxmax()
    bottom_df = pd.DataFrame(
        data={
            cast_id_col: df[cast_id_col].unique(),
            "max_p": df.loc[idx_p_max, "CTDPRS"],
            "lat": df.loc[idx_p_max, "GPSLAT"],
            "alt": df.loc[idx_p_max, "ALT"],
        }
    )
    bottom_df.loc[bottom_df["alt"] > threshold, "alt"] = np.nan
    # pandas 1.2.1 ufunc issue workaround with pd.to_numpy()
    bottom_df["DEPTH"] = (
        (
            bottom_df["alt"]
            + np.abs(gsw.z_from_p(bottom_df["max_p"], bottom_df["lat"].to_numpy()))
        )
        .fillna(value=-999)
        .round()
        .astype(int)
    )
    outfile = Path(report_dir, 'depth_log.csv')
    bottom_df[[cast_id_col, "DEPTH"]].to_csv(outfile, index=False)

    return True


def export_bottom_bottle_details(report_dir):
    report_dir = validate_dir(report_dir, create=True)
    report_file = Path(report_dir, 'bottom_bottle_details')

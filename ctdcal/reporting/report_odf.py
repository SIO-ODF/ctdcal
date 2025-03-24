"""
Formats data and generates visuals for ODF-style cruise reports.
"""
from ctdcal.flagging import flag_common as flagging


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

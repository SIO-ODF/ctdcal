import gsw
import numpy as np
import pandas as pd

"""code_pruning: is there anything worth keeping in this module?
only imported in outdated scripts"""


def merge_bottle_trip_dfs(file_ssscc):
    """Merge the bottle trip dataframes, and add flag columns.
    Input:
    file_ssscc - file path
    Output:
    df_all - dataframe, holding all relevant dataframe info
    """
    df_all = pd.DataFrame()

    ### load in all station/casts to be used for ctd processing
    all_casts = []
    with open(file_ssscc, "r") as filename:
        all_casts = [line.strip() for line in filename]

    for ssscc in all_casts:
        station = int(ssscc[0:3])
        cast = int(ssscc[3:5])

        ### rewrite this line to be more portable
        df = pd.read_pickle(f"data/bottle/{ssscc}_btl_mean.pkl")

        ### chop dataframe shorter for ease of use
        df = df[
            ["CTDPRS", "CTDTMP1", "CTDTMP2", "CTDCOND1", "CTDCOND2", "btl_fire_num"]
        ]
        df.rename(index=str, columns={"btl_fire_num": "SAMPNO"}, inplace=True)

        ### add columns
        df["CASTNO"] = cast
        df["STNNBR"] = station
        df["CTDTMP1_FLAG_W"] = 2
        df["CTDTMP2_FLAG_W"] = 2
        df["CTDCOND1_FLAG_W"] = 2
        df["CTDCOND2_FLAG_W"] = 2

        ### finally, join the frames together to have one giant dataframe
        df_all = pd.concat([df_all, df], axis=0)

    return df_all


def prelim_ctd_bottle_df(file_ssscc, file_whp_bottle):
    """Merge bottle trip and bottle ctd info."""
    ### load bottle trip data
    df_bottle_trip = merge_bottle_trip_dfs(file_ssscc)
    ### load whp_bottle file from odf_db
    df_whp_bottle = load_exchange_bottle_file(file_whp_bottle)
    df_whp_bottle = df_whp_bottle[
        [
            "STNNBR",
            "CASTNO",
            "SAMPNO",
            "SALNTY",
            "SALNTY_FLAG_W",
            "CTDSAL",
            "CTDSAL_FLAG_W",
            "REFTMP",
            "REFTMP_FLAG_W",
            "CTDOXY",
            "CTDOXY_FLAG_W",
            "OXYGEN",
            "OXYGEN_FLAG_W",
        ]
    ]

    ### merge both dataframes together
    df_merged = df_bottle_trip.merge(df_whp_bottle, on=["STNNBR", "CASTNO", "SAMPNO"])
    return df_merged


def ctd_residuals_df(df_bottle):
    """Compute residuals and add to dataframe.
    Operate in place.
    """
    # Salinity should really be grabbed straight from the files, but...
    try:
        df_bottle["BTLCOND"] = gsw.C_from_SP(
            df_bottle["SALNTY"], df_bottle["CTDTMP1"], df_bottle["CTDPRS"]
        )
    except ValueError:
        print("WHOOPS SOMETHING WENT WRONG")
        df_bottle["SALNTY"] = df_bottle["SALNTY"].replace(to_replace=-999, value=np.nan)
        df_bottle["BTLCOND"] = gsw.C_from_SP(
            df_bottle["SALNTY"], df_bottle["CTDTMP1"], df_bottle["CTDPRS"]
        )

    df_bottle["BTL_O"] = df_bottle["OXYGEN"] - df_bottle["CTDOXY"]

    df_bottle["BTL_C1"] = df_bottle["BTLCOND"] - df_bottle["CTDCOND1"]
    df_bottle["BTL_C2"] = df_bottle["BTLCOND"] - df_bottle["CTDCOND2"]
    df_bottle["C1_C2"] = df_bottle["CTDCOND1"] - df_bottle["CTDCOND2"]

    df_bottle["BTL_T1"] = df_bottle["REFTMP"] - df_bottle["CTDTMP1"]
    df_bottle["BTL_T2"] = df_bottle["REFTMP"] - df_bottle["CTDTMP2"]
    df_bottle["T1_T2"] = df_bottle["CTDTMP1"] - df_bottle["CTDTMP2"]

    df_bottle["BTL_SAL_UP"] = df_bottle["SALNTY"] - gsw.SP_from_C(
        df_bottle["CTDCOND2"], df_bottle["CTDTMP2"], df_bottle["CTDPRS"]
    )
    df_bottle["BTL_SAL"] = df_bottle["SALNTY"] - df_bottle["CTDSAL"]
    return None


def get_qc_t_from_df(df):
    """Limit data to only good data (WOCE discrete sample code 2)
    Returns a new copy of a dataframe.
    Input:
    df - dataframe of all data
    Output:
    df limited to only good temperature data
    """
    df = (
        df.groupby(["REFTMP_FLAG_W", "CTDTMP1_FLAG_W", "CTDTMP2_FLAG_W"])
        .get_group((2, 2, 2))
        .copy()
    )
    return df


def get_qc_c_from_df(df):
    """Limit data to only good data (WOCE discrete sample code 2)
    Returns a new copy of a dataframe.
    Input:
    df - dataframe of all data
    Output:
    df limited to only good conductivity data
    """
    df = (
        df.groupby(["SALNTY_FLAG_W", "CTDCOND1_FLAG_W", "CTDCOND2_FLAG_W"])
        .get_group((2, 2, 2))
        .copy()
    )
    return df


def get_qc_s_from_df(df):
    """Limit data to only good data (WOCE discrete sample code 2)
    Returns a new copy of a dataframe.
    Input:
    df - dataframe of all data
    Output:
    df limited to only good salinity data
    """
    df = df.groupby(["SALNTY_FLAG_W", "CTDSAL_FLAG_W"]).get_group((2, 2)).copy()
    return df


def get_qc_o_from_df(df):
    """Limit data to only good data (WOCE discrete sample code 2)
    Returns a new copy of a dataframe.
    Input:
    df - dataframe of all data
    Output:
    df limited to only good oxygen data
    """
    df = df.groupby(["OXYGEN_FLAG_W", "CTDOXY_FLAG_W"]).get_group((2, 2)).copy()
    return df


def get_qc_all_from_df(df):
    """Limit data to only good data (WOCE discrete sample code 2)
    Returns a new copy of a dataframe.
    Input:
    df - dataframe of all data
    Output:
    df limited to only good ctdo data
    """
    try:
        df = (
            df.groupby(
                [
                    "REFTMP_FLAG_W",
                    "CTDTMP1_FLAG_W",
                    "CTDTMP2_FLAG_W",
                    "CTDCOND1_FLAG_W",
                    "CTDCOND2_FLAG_W",
                    "SALNTY_FLAG_W",
                    "CTDSAL_FLAG_W",
                    "OXYGEN_FLAG_W",
                    "CTDOXY_FLAG_W",
                ]
            )
            .get_group((2, 2, 2, 2, 2, 2, 2, 2, 2))
            .copy()
        )
    except KeyError:
        df = (
            df.groupby(
                [
                    "REFTMP_FLAG_W",
                    "CTDTMP1_FLAG_W",
                    "CTDTMP2_FLAG_W",
                    "CTDCOND1_FLAG_W",
                    "CTDCOND2_FLAG_W",
                    "SALNTY_FLAG_W",
                    "CTDSAL_FLAG_W",
                    "OXYGEN_FLAG_W",
                    "CTDOXY_FLAG_W",
                ]
            )
            .get_group((2, 2, 2, 2, 2, 2, 2, 2, 1))
            .copy()
        )
    return df


def residual_stddev(
    df,
    param=[
        "BTL_O",
        "BTL_C1",
        "BTL_C2",
        "C1_C2",
        "BTL_T1",
        "BTL_T2",
        "T1_T2",
        "BTL_SAL",
        "BTL_SAL_UP",
    ],
):
    """Calculate standard deviations of parameters or residuals.
    Input:
    df - dataframe containing all data from casts
    param - list of parameter names
    Output:
    output - dictionary, key = is parameter name or parameter deep, value is 2 stddevs
    """
    output = {}
    for x in param:
        output[f"{x}"] = df[f"{x}"].std() * 2
        df_deep = df[df["CTDPRS"] > 2000]
        output[f"{x}_DEEP"] = df_deep[f"{x}"].std() * 2
    return output

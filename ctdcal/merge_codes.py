import os
import sys

import gsw
import numpy as np
import pandas as pd

"""code_pruning: is there anything worth keeping in this module?
only imported in outdated scripts"""


def string_converter(value):
    """To deal with Courtney CTD codes"""
    return value.split(":")[-1].strip()


def int_converter(value):
    """To deal with Courtney CTD codes"""
    return int(float(value.split(":")[-1].strip()))


def float_converter(value):
    """To deal with Courtney CTD codes"""
    return float(value.split(":")[-1])


def merge_courtney_ctd_codes(bottle_df, qc_df):
    """Take in ctd bottle file to be given to database, qc logs, and combine together.
    For use with Courtney style codes (P06W 2017).
    Input:
    bottle_file - dataframe
    qc_file - dataframe

    Output:
    b_test - dataframe with updated codes
    """

    flags = {
        "T": "REFTMP_FLAG_W",
        "T1": "CTDTMP1_FLAG_W",
        "T2": "CTDTMP2_FLAG_W",
        "C": "SALNTY_FLAG_W",
        "C1": "CTDCOND1_FLAG_W",
        "C2": "CTDCOND2_FLAG_W",
    }

    index = ["STNNBR", "CASTNO", "SAMPNO"]

    ### too lazy to fix it all in jupyter, do later
    a = qc_df
    b = bottle_df

    ### setup indicies for writing to
    a["CASTNO"] = a["stacast"] % 100
    a["STNNBR"] = (a["stacast"] - a["CASTNO"]) // 100
    a["new_params"] = None

    ### update new_params to flag_parameter for df.groupby() later
    for param, flag_param in flags.items():
        a.loc[a["parameter"] == param, "new_params"] = flag_param
    a_gb = a.groupby("new_params")

    ### df.update() requires the index to be set before being called
    b_test = b.set_index(index)

    ### loop through all flag_parameter names
    for param, flag_param in flags.items():
        try:
            a_gb_project = a_gb.get_group(flag_param).copy()
            a_gb_project.rename(columns={"flag": flag_param}, inplace=True)

            ### df.update() requires the index to be set before being called
            a_test = a_gb_project[index + [flag_param]].set_index(index)

            ### df.update() returns inplace, so return signature = None
            b_test.update(a_test)
        except KeyError:
            print(f"Parameter not found: {flag_param}")

    ### reset to work with again
    b_test = b_test.reset_index()

    return b_test


def load_courtney_ctd_codes(qc_file):
    """Load Courtney style ctd codes, and pass out dataframe.
    Input:
    qc_file - filepath to a courtney style ctd code
    Output:
    dataframe
    """
    converters = {
        0: int_converter,
        1: int_converter,
        2: string_converter,
        3: float_converter,
        4: int_converter,
        5: float_converter,
        6: float_converter,
        7: float_converter,
    }
    names = [
        "stacast",
        "SAMPNO",
        "parameter",
        "pressure",
        "flag",
        "primary_diff",
        "secondary_diff",
        "p_minus_s",
    ]
    return pd.read_csv(qc_file, header=None, converters=converters, names=names)


def load_to_odf_db_bottle_file(bottle_file):
    """File to be given to odf db maintainer to be loaded
    Input:
    bottle_file - filepath
    Output:
    dataframe
    """
    return pd.read_csv(bottle_file, skiprows=[1], skipfooter=1, engine="python")


def load_exchange_bottle_file(bottle_file):
    """WHP-exchange bottle file from odf db
    Input:
    bottle_file - filepath
    Output:
    dataframe
    """
    return pd.read_csv(
        bottle_file,
        skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11],
        skipfooter=1,
        engine="python",
    )


def write_to_odf_db_bottle_file(df_bottle, output_file):
    """Write WHP_exchange file for odf_db to load.
    Input:
    df_bottle - dataframe with bottlefile to write out
    output_file - filepath to write to
    Output:
    None
    """
    stringy = "CASTNO,CTDBACKSCATTER,CTDBACKSCATTER_FLAG_W,CTDFLUOR,CTDFLUOR_FLAG_W,"
    stringy += (
        "CTDOXY,CTDOXY_FLAG_W,CTDPRS,CTDPRS_FLAG_W,CTDRINKO,CTDRINKO_FLAG_W,CTDSAL,"
    )
    stringy += "CTDSAL_FLAG_W,CTDTMP,CTDTMP_FLAG_W,CTDXMISS,CTDXMISS_FLAG_W,REFTMP,REFTMP_FLAG_W,"
    stringy += "SAMPNO,STNNBR\n,VOLTS,,VOLTS,,UMOL/KG,,DBAR,,VOLTS,,PSS-78,,ITS-90,,VOLTS,,ITS-90,,,\n"
    stringy += df_bottle.to_string(
        header=False,
        index=False,
        formatters=[
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x) + ",",
            lambda x: str(x),
        ],
    )
    stringy += "\nEND_DATA\n"

    with open(output_file, "w") as f_handle:
        f_handle.write(stringy)
    return None


def merge_ctd_codes(log_folder, df_bottle):
    """
    Chopped version of old merge_ctd_codes to deal with dfs
    Input:
    log_folder - file path
    bottle_file - df

    Output:
    btl -  dataframe
    """
    ### Expecting all codes to be in one directory for now. Must change later
    log_files = os.listdir(log_folder)

    ### Merge quality codes, and store a copy in-memory in btl
    for f in log_files:
        x = load_courtney_ctd_codes(f"{log_folder}{f}")
        df_bottle = merge_courtney_ctd_codes(df_bottle, x)

    return df_bottle


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


def main(argv):
    """Example run"""
    qual_codes_filepath = "data/quality_codes/"
    cruise_dir = "data/"
    log_dir = "data/logs"

    file_ssscc = f"{cruise_dir}ssscc.csv"
    bottle_file = f"{qual_codes_filepath}320620180309_hy1.csv"

    ### compile bottle trips into one file, then merge with odf_db bottle file
    df_bottle = prelim_ctd_bottle_df(file_ssscc, bottle_file)
    ### compute residuals for plotting
    ctd_residuals_df(df_bottle)
    load_courtney_ctd_codes()
    df_bottle = merge_ctd_codes(log_dir, df_bottle)
    df_bottle.to_pickle(f"{qual_codes_filepath}merged_ctd_codes.pkl")

    # df_coded = get_qc_all_from_df(df_bottle)
    ### write file out
    return None


if __name__ == "__main__":
    main(sys.argv[1:])

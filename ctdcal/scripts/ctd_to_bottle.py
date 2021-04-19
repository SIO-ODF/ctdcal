import sys

import numpy as np
import pandas as pd


def main(argv):
    """Creates a bottle file with CTD downcast information.

    Needs to be cleaned up from script form to other form, hopefully
    """
    # general fileloading area
    output_file = "data/scratch_folder/ctd_to_bottle.csv"
    file_ssscc = "data/ssscc.csv"
    df_all = pd.DataFrame()

    all_casts = []
    with open(file_ssscc, "r") as filename:
        all_casts = [line.strip() for line in filename]

    for ssscc in all_casts:
        station = int(ssscc[0:3])
        cast = int(ssscc[3:5])
        # bottle handling section
        dir_bottle = "data/bottle/"
        bottle_postfix = "_btl_mean.pkl"

        df_bottle = pd.read_pickle(f"{dir_bottle}{ssscc}{bottle_postfix}")
        # next line not strictly necessaryas we don't move every column, but left just in case
        df_bottle.rename(
            index=str,
            columns={
                "FREE1": "CTDFLUOR",
                "FREE2": "CTDBACKSCATTER",
                "FREE3": "CTDRINKO",
                "btl_fire_num": "SAMPNO",
                "CTDTMP1": "CTDTMP",
                "CTDOXY1": "CTDOXY",
            },
            inplace=True,
        )
        df_bottle = df_bottle.loc[:, ["CTDPRS", "SAMPNO"]]
        df_bottle["bins"] = pd.cut(
            df_bottle["CTDPRS"],
            range(0, int(np.ceil(df_bottle["CTDPRS"].max())) + 5, 2),
            right=False,
            include_lowest=True,
        )

        # ctd handling section
        dir_ctd = "data/pressure/"
        ctd_postfix = "_ct1.csv"
        ctd_skiprows = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]

        df_c = pd.read_csv(
            dir_ctd + ssscc + ctd_postfix,
            skiprows=ctd_skiprows,
            skipfooter=1,
            engine="python",
        )
        df_c["bins"] = pd.cut(
            df_c["CTDPRS"],
            range(0, int(np.floor(df_c["CTDPRS"].max())) + 5, 2),
            right=False,
            include_lowest=True,
        )

        df_merged_btl_ctd_4 = pd.merge(
            df_c, df_bottle[["bins"]], on="bins", how="outer"
        )

        df_merged_btl_ctd_6 = pd.merge(
            df_merged_btl_ctd_4.iloc[:, 1:], df_bottle, on="bins", how="inner"
        )
        df_merged_btl_ctd_6.groupby("SAMPNO").first()

        # reference temperature section
        # there should not be any 9 values, unless the reftemp dies? then we have a PROBLEM
        dir_reft = "data/reft/"
        reft_postfix = "_reft.csv"
        try:
            df_reft = pd.read_csv(
                dir_reft + ssscc + reft_postfix,
                names=["SAMPNO", "REFTMP"],
                skiprows=1,
                usecols=[2, 5],
            )
            df_reft["REFTMP_FLAG_W"] = 2
            df_incomplete_reft = pd.merge(
                df_reft, df_merged_btl_ctd_6, on="SAMPNO", how="outer"
            )
            # this line breaks code; seem redundant? MK
            # df_incomplete_reft.groupby('SAMPNO').fillna(inplace=True, method='pad')
            df_complete_reft = df_incomplete_reft.fillna(
                value={"REFTMP": -999, "REFTMP_FLAG_W": 5}
            )
        # if there are no SBE35 files for the cast, flag as 5
        except FileNotFoundError:
            df_merged_btl_ctd_6["REFTMP"] = -999
            df_merged_btl_ctd_6["REFTMP_FLAG_W"] = 5
            df_complete_reft = df_merged_btl_ctd_6

        df_complete_reft["STNNBR"] = station
        df_complete_reft["CASTNO"] = cast
        df_complete_reft = df_complete_reft.groupby("SAMPNO").first()
        df_complete_reft.reset_index(inplace=True)
        df_complete_reft = df_complete_reft.sort_index(axis=1)
        df_all = df_all.append(df_complete_reft)
        # turn any nans to flag 9 values
        fillna_values = {
            "CTDBACKSCATTER": -999,
            "CTDBACKSCATTER_FLAG_W": 9,
            "CTDFLUOR": -999,
            "CTDFLUOR_FLAG_W": 9,
            "CTDOXY": -999,
            "CTDOXY_FLAG_W": 9,
            "CTDPRS": -999,
            "CTDPRS_FLAG_W": 9,
            "CTDRINKO": -999,
            "CTDRINKO_FLAG_W": 9,
            "CTDSAL": -999,
            "CTDSAL_FLAG_W": 9,
            "CTDTMP": -999,
            "CTDTMP_FLAG_W": 9,
            "CTDXMISS": -999,
            "CTDXMISS_FLAG_W": 9,
            "REFTMP": -999,
            "REFTMP_FLAG_W": 9,
        }
        df_all = df_all.fillna(fillna_values)

    # create the stupid string
    # stringy = "CASTNO,CTDBACKSCATTER,CTDBACKSCATTER_FLAG_W,CTDFLUOR,CTDFLUOR_FLAG_W,"
    # MK 03/27/21: no backscatter on A20/22
    stringy = "CASTNO,CTDFLUOR,CTDFLUOR_FLAG_W,"
    stringy += (
        # "CTDOXY,CTDOXY_FLAG_W,CTDPRS,CTDRINKO,CTDRINKO_FLAG_W,CTDSAL,"
        "CTDOXY,CTDOXY_FLAG_W,CTDPRS,CTDSAL,"
    )
    stringy += "CTDSAL_FLAG_W,CTDTMP,CTDTMP_FLAG_W,CTDXMISS,CTDXMISS_FLAG_W,REFTMP,REFTMP_FLAG_W,"
    # stringy += "SAMPNO,STNNBR\n,VOLTS,,VOLTS,,UMOL/KG,,DBAR,,VOLTS,,PSS-78,,ITS-90,,VOLTS,,ITS-90,,,\n"
    stringy += "SAMPNO,STNNBR\n,VOLTS,,UMOL/KG,,DBAR,,VOLTS,,PSS-78,,ITS-90,,VOLTS,,ITS-90,,,\n"
    # the next line is incredibly stupid. find a better way to create
    stringy += df_all.iloc[:, 0:-1].to_string(
        header=False,
        index=False,
        formatters=[
            # lambda x: str(x) + ",",
            # lambda x: str(x) + ",",
            # lambda x: str(x) + ",",
            # lambda x: str(x) + ",",
            # lambda x: str(x) + ",",
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


if __name__ == "__main__":
    main(sys.argv[1:])

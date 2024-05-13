import pandas as pd

from ctdcal import process_bottle

# def test_add_btlnbr_cols():
#     df = {
#         "btl_fire_num": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
#         "other_id":     [1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15]    #   Like the LADCP is mounted on it
#     }
#     df = pd.DataFrame(df)
#     btl_num_col = "btl_fire_num"

#     #   Call the function
#     result_df = process_bottle.add_btlnbr_cols(df, btl_num_col)

#     #   Check if the columns are added right
#     assert "BTLNBR" in result_df.columns
#     assert "BTLNBR_FLAG_W" in result_df.columns
#     assert all(result_df["BTLNBR_FLAG_W"] == 2)
#     assert result_df["BTLNBR"].equals(df[btl_num_col])

#     #   Cross compare
#     assert result_df["BTLNBR"].iloc[-1] != result_df["other_id"].iloc[-1]

def test_load_hy_file(tmp_path):
    # Create a sample file
    sample_data = """BOTTLE, SomecruiseODF
EXPOCODE,SECT_ID,STNNBR,CASTNO,SAMPNO,BTLNBR,BTLNBR_FLAG_W,DATE,TIME,LATITUDE,LONGITUDE,DEPTH,CTDPRS,CTDTMP,CTDSAL,CTDSAL_FLAG_W,SALNTY,SALNTY_FLAG_W,CTDOXY,CTDOXY_FLAG_W,OXYGEN,OXYGEN_FLAG_W,REFTMP,REFTMP_FLAG_W
,,,,,,,,,,,METERS,DBAR,ITS-90,PSS-78,,PSS-78,,UMOL/KG,,UMOL/KG,,ITS-90,
Somecruise,GOSHIP,1,1,36,36,2,19990722,1345,-32.69,114.97,114,3.5,20.30,35.4,2,-999.0,9,283.8,2,-999.0,9,-999.0,9
Somecruise,GOSHIP,1,1,35,35,2,19990722,1345,-32.69,114.97,114,4.4,20.30,35.4,2,-999.0,9,311.7,2,-999.0,9,20.3,2
Somecruise,GOSHIP,1,1,34,34,2,19990722,1345,-32.69,114.97,114,2.8,20.29,35.4,2,-999.0,9,274.8,2,-999.0,9,-999.0,9

END_DATA
    """
    fname = tmp_path / "sample_hyfile.csv"
    with open(fname, "w") as f:
        f.write(sample_data)

    df = process_bottle.load_hy_file(fname)

    #   Make sure everything is inside that we're expecting
    assert not df.empty
    assert "EXPOCODE" in df.columns
    assert "SECT_ID" in df.columns
    assert "STNNBR" in df.columns
    assert "CASTNO" in df.columns
    assert "DEPTH" in df.columns

    #   Check if the final row with EXPOCODE "END_DATA" is dropped
    assert "END_DATA" not in df["EXPOCODE"].values
    #   Value check
    assert df["BTLNBR"].iloc[-1] == 34

def test_load_btl_data(tmp_path):
    fname = tmp_path / "90909_btl_mean.pkl"
    df = pd.DataFrame(
        {
            "CTDTMP1": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            "CTDTMP2": [0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
            "CTDPRS": [5000, 4000, 3000, 2000, 1000, 0],
            "btl_fire_num": [1, 2, 3, 4, 5, 6],
        }
    )
    df.to_pickle(fname)

    # test read all columns
    btl_data = process_bottle._load_btl_data(fname)
    assert df.equals(btl_data[df.columns])
    assert btl_data.shape == (6, 5)
    assert all(btl_data["SSSCC"] == "90909")

    # test read select columns
    cols = ["CTDTMP1", "CTDPRS", "btl_fire_num"]
    btl_data_select = process_bottle._load_btl_data(fname, cols=cols)
    assert "CTDTMP2" not in btl_data_select.columns
    assert all(btl_data["SSSCC"] == "90909")


def test_load_reft_data(tmp_path):
    fname = tmp_path / "90909_reft.csv"
    df = pd.DataFrame(
        {
            "btl_fire_num": [1, 2, 3, 4, 5, 6],
            "T90": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            "REFTMP_FLAG_W": 2,
        }
    )
    df.to_csv(fname)

    # test read
    reft_data = process_bottle._load_reft_data(fname)
    assert df.equals(reft_data[df.columns])
    assert all(reft_data["SSSCC_TEMP"] == "90909")
    assert reft_data["REFTMP"].equals(reft_data["T90"])


def test_load_salt_data(tmp_path):
    fname = tmp_path / "90909_salts.csv"
    df = pd.DataFrame(
        {
            "SAMPNO": [1, 2, 3, 4, 5, 6],
            "SALNTY": [30.0, 31.0, 32.0, 33.0, 34.0, 35.0],
            "BathTEMP": 24,
            "CRavg": [0.9, 0.91, 0.92, 0.93, 0.94, 0.95],
        }
    )
    df.to_csv(fname)

    # test read
    salt_data = process_bottle._load_salt_data(fname)
    assert all(salt_data["SSSCC_SALT"] == "90909")

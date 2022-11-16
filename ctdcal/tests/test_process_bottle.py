import pandas as pd

from ctdcal import process_bottle


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

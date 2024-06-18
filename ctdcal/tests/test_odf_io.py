import io
import logging
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from ctdcal import odf_io


def check_type(to_check, sub_dtype):
    """
    Type comparisons can fail on Windows due to int32/int64 behavior.
    This function checks subtypes so that tests don't fail on Windows.

    See https://stackoverflow.com/q/64901822/13155619 for more info
    """
    if sub_dtype is int:
        np_type = np.integer
    elif sub_dtype is float:
        np_type = np.floating

    return all([np.issubdtype(dtype, np_type) for dtype in to_check.dtypes])


def make_salt_file(stn=1, cast=1, comment=None, flag=False, to_file=None):
    #  seed RNG for Reading3, comments, and flags
    rng = np.random.default_rng(seed=100)

    # build dummy DataFrame
    constants = {
        "STNNBR": f"{stn:04.0f}",
        "CASTNO": f"{cast:02.0f}",
        "BathTemp": 21,
        "unk": 5907,
    }
    salts = pd.DataFrame(data=constants, index=np.arange(12))
    salts.insert(2, "SAMPNO", [f"{x:02.0f}" for x in salts.index])
    salts.insert(4, "CRavg", np.linspace(1.95, 1.98, 12)[::-1])
    salts.insert(5, "autosalSAMPNO", salts["SAMPNO"].astype(int))
    salts["autosalSAMPNO"] = salts["autosalSAMPNO"].astype(object)  #   Change dtype for integers and strings
    salts.loc[[0, 11], "autosalSAMPNO"] = "worm"
    times = pd.date_range(start="18:32:04", end="19:36:03", periods=24)
    salts["StartTime"] = times[::2].strftime("%H:%M:%S")
    salts["EndTime"] = times[1::2].strftime("%H:%M:%S")

    # assign CR values around mean (give only some samples 3 readings)
    salts["Reading1"] = salts["CRavg"] + 1e-4
    salts["Reading2"] = salts["CRavg"] - 1e-4
    salts["Reading3"] = salts["CRavg"] * rng.choice([1, np.nan], 12)
    attempts = salts[["Reading1", "Reading2", "Reading3"]].count(axis=1)
    salts.insert(9, "Attempts", attempts.map("{:02.0f}".format))

    # add comment marker (#, x)
    if comment is not None:
        salts["STNNBR"] = rng.choice(["", comment], 12, p=[0.8, 0.2]) + salts["STNNBR"]

    # add LabView "quality"(?) flag
    if flag:
        salts["EndTime"] = salts["EndTime"] + rng.choice(["", "*"], 12, p=[0.8, 0.2])

    # final formatting, remove blank Reading3 cells to match Autosal
    header = "12-345 operator: ABC box: S batch: P678 k15: 0.91011 std dial 408"
    string_df = salts.to_string(
        header=False, index=False, float_format="{:.5f}".format, na_rep=""
    )
    text_out = "\n".join([header, string_df.replace("        ", "")])

    if to_file is not None:
        with open(to_file, "w+") as f:
            f.write(text_out)
    else:
        return text_out


def test_salt_loader(caplog, tmp_path):
    # check salt file loads in correctly
    salt_file = make_salt_file()
    saltDF, refDF = odf_io._salt_loader(io.StringIO(salt_file))

    assert saltDF.shape == (10, 14)
    assert all(saltDF[["StartTime", "EndTime"]].dtypes == object)
    assert check_type(saltDF[["CRavg", "Reading1", "Reading2", "Reading3"]], float)
    assert check_type(saltDF[["STNNBR", "CASTNO", "SAMPNO", "autosalSAMPNO"]], int)
    assert check_type(saltDF[["BathTEMP", "Unknown", "Attempts", "IndexTime"]], int)
    assert all(saltDF.index == np.arange(1, 11))
    assert saltDF["Reading3"].isna().sum() == 5

    assert refDF.shape == (2, 2)
    assert all(refDF.dtypes == float)
    assert all(refDF.index == [0, 11])

    # check commented lines are ignored ("bottles" 1, 4, 10)
    salt_file = make_salt_file(comment="#")
    with caplog.at_level(logging.DEBUG):
        saltDF, refDF = odf_io._salt_loader(io.StringIO(salt_file))
        assert "(#, x)" in caplog.messages[0]
        assert "test_odf_io" in caplog.messages[0]
    assert saltDF.shape == (7, 14)
    assert all(saltDF[["StartTime", "EndTime"]].dtypes == object)
    assert check_type(saltDF[["CRavg", "Reading1", "Reading2", "Reading3"]], float)
    assert check_type(saltDF[["STNNBR", "CASTNO", "SAMPNO", "autosalSAMPNO"]], int)
    assert check_type(saltDF[["BathTEMP", "Unknown", "Attempts", "IndexTime"]], int)
    assert all(saltDF.index == [2, 3, 5, 6, 7, 8, 9])
    assert saltDF["Reading3"].isna().sum() == 3

    assert refDF.shape == (2, 2)
    assert all(refDF.dtypes == float)
    assert all(refDF.index == [0, 11])

    # check flagged EndTimes are added to flag_file ("bottles" 1, 4, 10)
    salt_file = make_salt_file(flag=True)
    with caplog.at_level(logging.DEBUG):
        f_out = io.StringIO()
        saltDF, refDF = odf_io._salt_loader(io.StringIO(salt_file), flag_file=f_out)
        flagged = f_out.getvalue().split("\n")
        assert "Found * in test_odf_io" in caplog.messages[1]
        assert "test_odf_io,1,,3,Auto-flagged" in flagged[0]
        assert "test_odf_io,4,,3,Auto-flagged" in flagged[1]
        assert "test_odf_io,10,,3,Auto-flagged" in flagged[2]
        assert flagged[3] == ""
    assert saltDF.shape == (10, 14)
    assert all(saltDF[["StartTime", "EndTime"]].dtypes == object)
    assert check_type(saltDF[["CRavg", "Reading1", "Reading2", "Reading3"]], float)
    assert check_type(saltDF[["STNNBR", "CASTNO", "SAMPNO", "autosalSAMPNO"]], int)
    assert check_type(saltDF[["BathTEMP", "Unknown", "Attempts", "IndexTime"]], int)
    assert all(saltDF.index == np.arange(1, 11))
    assert saltDF["Reading3"].isna().sum() == 5

    assert refDF.shape == (2, 2)
    assert all(refDF.dtypes == float)
    assert all(refDF.index == [0, 11])

    # check behavior with filename (saves time not writing full file to disk 3x)
    d = tmp_path / "salt"
    d.mkdir()
    fake_file = d / "90909"
    fake_file.write_text("\n1 2 3 4 5 6 7 00:01:00 00:02:00 10 11") #   0001 06 13 24 1.99187   13 5427 16:31:39  16:32:16  02 1.99186 1.99188
    saltDF, refDF = odf_io._salt_loader(fake_file)
    assert all(saltDF[["StartTime", "EndTime"]].dtypes == object)
    assert check_type(saltDF[["CRavg", "Reading1"]], float)
    assert check_type(saltDF[["STNNBR", "CASTNO", "SAMPNO", "autosalSAMPNO"]], int)
    assert check_type(saltDF[["BathTEMP", "Unknown", "Attempts", "IndexTime"]], int)
    assert saltDF.shape == (1, 12)
    assert refDF.empty


def test_remove_autosal_drift(caplog):
    index_time = np.linspace(0, 300, 7)
    CRavg = np.array([2.0] * 6 + [1.4])
    saltDF = pd.DataFrame(data={"STNNBR": 1, "CASTNO": 1}, index=np.arange(1, 6))
    saltDF["IndexTime"], saltDF["CRavg"] = (index_time[1:-1], CRavg[1:-1])
    refDF = pd.DataFrame(data={"IndexTime": index_time[::6], "CRavg": CRavg[::6]})

    # check that drift is removed correctly
    removed = odf_io.remove_autosal_drift(saltDF, refDF)
    np.testing.assert_allclose(removed["CRavg"], np.arange(1.5, 2, 0.1)[::-1])
    assert "IndexTime" not in removed.columns

    # check input DataFrames have not been modified
    assert all(saltDF["CRavg"] == CRavg[1:-1])
    assert all(refDF["CRavg"] == CRavg[::6])

    # return unmodified saltDF if refDF is wrong size
    original = odf_io.remove_autosal_drift(saltDF, refDF.iloc[0])
    assert all(original["CRavg"] == saltDF["CRavg"])
    assert "IndexTime" not in original.columns
    assert "start/end reference" in caplog.messages[0]
    assert "00101" in caplog.messages[0]


def test_salt_exporter(tmp_path, caplog):
    saltDF = pd.DataFrame(data={"STNNBR": 1, "CASTNO": 1}, index=np.arange(5))
    saltDF["CRavg"] = np.ones(5)

    # check test data completes round trip
    assert not Path(tmp_path / "00101_salts.csv").exists()
    odf_io._salt_exporter(saltDF, outdir=str(tmp_path))
    assert Path(tmp_path / "00101_salts.csv").exists()
    with open(Path(tmp_path / "00101_salts.csv")) as f:
        assert saltDF.equals(pd.read_csv(f))

    # check "file already exists" message
    with caplog.at_level(logging.INFO):
        odf_io._salt_exporter(saltDF, outdir=str(tmp_path))
        assert "00101_salts.csv already exists" in caplog.messages[0]

    # check file write with multiple stations
    saltDF["STNNBR"] = [1, 1, 2, 2, 2]
    odf_io._salt_exporter(saltDF, outdir=str(tmp_path))
    assert Path(tmp_path / "00201_salts.csv").exists()

    # check file write with multiple casts
    saltDF["CASTNO"] = [1, 2, 1, 2, 3]
    odf_io._salt_exporter(saltDF, outdir=str(tmp_path))
    assert Path(tmp_path / "00102_salts.csv").exists()
    assert Path(tmp_path / "00202_salts.csv").exists()
    assert Path(tmp_path / "00203_salts.csv").exists()

    # check empty (all NaN) Reading# columns are dropped
    saltDF["STNNBR"], saltDF["CASTNO"] = 3, 1
    saltDF["Reading1"] = 1
    saltDF["Reading2"] = np.nan
    assert not Path(tmp_path / "00301_salts.csv").exists()
    odf_io._salt_exporter(saltDF, outdir=str(tmp_path))
    with open(Path(tmp_path / "00301_salts.csv")) as f:
        empty = pd.read_csv(f)
        assert "Reading1" in empty.columns
        assert "Reading2" not in empty.columns


def test_process_salts(tmp_path, caplog):
    # check missing salt files are logged and skipped
    odf_io.process_salts(["90909"], salt_dir=str(tmp_path))
    assert "90909 does not exist" in caplog.messages[0]

    # check salt file is processed through to .csv without errors
    f_path = tmp_path / "90909"
    assert not Path(tmp_path / "90909").exists()
    make_salt_file(stn=909, cast=9, to_file=f_path)  # outfile name depends on stn/cast
    with caplog.at_level(logging.INFO):
        assert not Path(tmp_path / "90909_salts.csv").exists()
        odf_io.process_salts(["90909"], salt_dir=str(tmp_path))
        assert len(caplog.messages) == 1
        assert Path(tmp_path / "90909_salts.csv").exists()

    # check nothing happens if .csv already exists
    with caplog.at_level(logging.INFO):
        odf_io.process_salts(["90909"], salt_dir=str(tmp_path))
        assert "90909_salts.csv already exists" in caplog.messages[1]

    def test_print_progress_bar():
        # Test parameters
        iteration = 3
        total = 10
        prefix = "Progress"
        suffix = "Complete"
        decimals = 1
        length = 20
        fill = "#"
        printEnd = "\r"

        # Expected output
        expected_output = "\rProgress |###-------| 30.0% Complete"

        # Redirect stdout to mock_stdout
        with patch('sys.stdout', io.StringIO()):
            # Call function
            odf_io.print_progress_bar(iteration, total, prefix, suffix, decimals, length, fill, printEnd)

        # Get value from mock_stdout
        actual_output = io.StringIO().getvalue()

        # Check if the actual output matches the expected output
        assert actual_output == expected_output

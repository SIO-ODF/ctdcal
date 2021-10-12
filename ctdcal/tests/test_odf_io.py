import logging
from pathlib import Path

import numpy as np
import pandas as pd

from ctdcal import odf_io


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

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
    assert "start/end reference" in caplog.records[0].message
    assert "00101" in caplog.records[0].message

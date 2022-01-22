import logging

import numpy as np
import requests

from ctdcal import io


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


def test_load_exchange_btl(caplog, tmp_path, monkeypatch):
    # make fake/empty Exchange file
    content = [
        "BOTTLE,20210101ODFSIO\n",
        "#\n#\n#\n#\n",
        "EXPOCODE,SECT_ID,STNNBR,CASTNO,SAMPNO\n",
        ",,,,\n",
        "  123420210101,   ABC,     1,  2,      3\n",
        "  123420210101,   ABC,     1,  2,      4\n",
        "END_DATA",
    ]

    # write fake data and check read from file
    with open(tmp_path / "test_hy1.csv", "w+") as f:
        for line in content:
            f.write(line)

    with caplog.at_level(logging.INFO):
        hy1 = io.load_exchange_btl(tmp_path / "test_hy1.csv")
    assert hy1.shape == (2, 5)
    assert check_type(hy1[["EXPOCODE", "STNNBR", "CASTNO", "SAMPNO"]], int)
    assert all(hy1["SECT_ID"] == "ABC")  # should not have any leading spaces
    assert "from local file" in caplog.messages[0]

    # check read works with str path as well
    with caplog.at_level(logging.INFO):
        assert hy1.equals(io.load_exchange_btl(f"{str(tmp_path)}/test_hy1.csv"))
        assert "from local file" in caplog.messages[1]

    # mock downloading data from CCHDO
    class MockResponse(object):
        def __init__(self):
            self.text = "".join(content)

    def mock_get(*args, **kwargs):
        return MockResponse()

    # "overwrite" requests.get() call in ctdcal.io with mock_get()
    monkeypatch.setattr(requests, "get", mock_get)

    # check read works from URL
    with caplog.at_level(logging.INFO):
        assert hy1.equals(io.load_exchange_btl("http://fakeurl"))
        assert "from http link" in caplog.messages[2]


def test_write_pressure_details(tmp_path):
    f_path = tmp_path / "prs_log.csv"

    # check file is created if it doesn't exist
    assert not (f_path).exists()
    io.write_pressure_details("00101", f_path, "00:00:01", "00:04:01")
    assert (f_path).exists()

    # check only new data are appended (not another header)
    io.write_pressure_details("00201", f_path, "00:05:01", "00:09:01")
    with open(f_path, "rb") as f:
        contents = f.readlines()
        assert len(contents) == 3
        assert b"SSSCC" in contents[0]
        assert b"00101" in contents[1]
        assert b"00201" in contents[2]


def test_write_cast_details(tmp_path):
    f_path = tmp_path / "cast.csv"

    # check file is created if it doesn't exist
    assert not f_path.exists()
    io.write_cast_details("00101", f_path, 1.0, 2.0, 3.0, 0.0, 100.0, 5.0, -70.0, 170.0)
    assert f_path.exists()

    # check only new data are appended (not another header)
    io.write_cast_details("00201", f_path, 2.0, 3.0, 4.0, 0.0, 99.0, 6.0, -70.0, 170.0)
    with open(f_path, "rb") as f:
        contents = f.readlines()
        assert len(contents) == 3
        assert b"SSSCC" in contents[0]
        assert b"00101" in contents[1]
        assert b"00201" in contents[2]

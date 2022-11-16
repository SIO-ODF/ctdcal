import logging
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo
from zipimport import ZipImportError

import numpy as np
import pytest
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


def test_load_cnv(tmp_path):
    # make fake/empty .cnv file
    content = [
        "* header rows\n* header rows\n* header rows\n* header rows\n",
        "# nquan = 3\n# nvalues = 4\n# units = specified\n",
        "# name 0 = prDM\n# name 1 = depSM\n# name 2 = t090C\n",
        "# bad_flag = -9.990e-29\n",
        "# comment rows\n# comment rows\n# comment rows\n",
        "*END*\n",
        "      4.000      3.962     8.9121\n",
        "      6.000      5.943     8.9125\n",
        "      8.000      7.924     8.9123\n",
        "     10.000      9.905     8.9126\n",
    ]

    # write fake data and check read from file
    with open(tmp_path / "test_1.cnv", "w+") as f:
        for line in content:
            f.write(line)

    cnv = io.load_cnv(tmp_path / "test_1.cnv")
    assert cnv.shape == (4, 3)
    assert check_type(cnv, float)
    assert cnv.columns.tolist() == ["prDM", "depSM", "t090C"]


def test_load_exchange_btl(caplog, tmp_path, monkeypatch):
    # make fake/empty Exchange file
    content = [
        "BOTTLE,20210101ODFSIO\n",
        "# comment about EXPOCODE\n",
        "#\n#\n#\n#\n",
        "EXPOCODE,SECT_ID,STNNBR,CASTNO,SAMPNO\n",
        ",,,,\n",
        "  123420210101,   ABC,     1,  2,      3\n",
        "  123420210101,   ABC,     1,  2,      4\n",
        "END_DATA",
    ]

    # write fake data and check read from file
    fname = "test_hy1.csv"
    with open(tmp_path / fname, "w+") as f:
        for line in content:
            f.write(line)

    with caplog.at_level(logging.INFO):
        hy1 = io.load_exchange_btl(tmp_path / fname)
    assert hy1.shape == (2, 5)
    assert check_type(hy1[["EXPOCODE", "STNNBR", "CASTNO", "SAMPNO"]], int)
    assert all(hy1["SECT_ID"] == "ABC")  # should not have any leading spaces
    assert f"{fname} from local file" in caplog.messages[0]

    # check read works with str path as well
    with caplog.at_level(logging.INFO):
        assert hy1.equals(io.load_exchange_btl(f"{str(tmp_path)}/{fname}"))
        assert f"{fname} from local file" in caplog.messages[1]

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
        assert hy1.equals(io.load_exchange_btl(f"http://fakeurl/{fname}"))
        assert f"{fname} from http link" in caplog.messages[2]


def test_load_exchange_ctd(caplog, tmp_path, monkeypatch):
    # make fake/empty Exchange file
    content = [
        "CTD,20210101ODFSIO\n",
        "# comment about NUMBER_HEADERS\n",
        "# comment about CTDPRS\n",
        "#\n#\n#\n#\n",
        "NUMBER_HEADERS = 11\n",
        "EXPOCODE = 012345678910\n",
        "SECT_ID = ABC\n",
        "STNNBR = 1\n",
        "CASTNO = 1\n",
        "DATE = 20220101\n",
        "TIME = 1234\n",
        "LATITUDE = 45.6789\n",
        "LONGITUDE = -50.1234\n",
        "DEPTH = 96\n",
        "INSTRUMENT_ID = 987\n",
        "CTDPRS,CTDTMP,CTDTMP_FLAG_W,CTDSAL,CTDSAL_FLAG_W,CTDOXY,CTDOXY_FLAG_W\n",
        "DBAR,ITS-90,,PSS-78,,UMOL/KG,\n",
        "      0.0,   2.8000,1,  32.5320,2,    330.0,3\n",
        "      2.0,   2.8070,2,  33.0000,3,    200.0,4\n",
        "END_DATA",
    ]

    # write fake data and check read from file
    fname = "test_ct1.csv"
    with open(tmp_path / fname, "w+") as f:
        for line in content:
            f.write(line)

    with caplog.at_level(logging.INFO):
        header, ct1 = io.load_exchange_ctd(tmp_path / fname)
    assert ct1.shape == (2, 7)
    assert len(header) == 11
    assert check_type(ct1[["CTDPRS", "CTDTMP", "CTDSAL", "CTDOXY"]], float)
    assert check_type(ct1[["CTDTMP_FLAG_W", "CTDSAL_FLAG_W", "CTDOXY_FLAG_W"]], int)
    assert f"{fname} from local file" in caplog.messages[0]

    # check read works with str path as well
    caplog.clear()
    with caplog.at_level(logging.INFO):
        header_str, ct1_str = io.load_exchange_ctd(f"{str(tmp_path)}/{fname}")
        assert ct1.equals(ct1_str)
        assert header_str == header
        assert f"{fname} from local file" in caplog.messages[0]

    # mock downloading data from CCHDO (single cast file; .zip tested separately)
    class MockResponse(object):
        def __init__(self):
            self.content = "".join(content).encode("utf8")

    def mock_get(*args, **kwargs):
        return MockResponse()

    # "overwrite" requests.get() call in ctdcal.io with mock_get()
    monkeypatch.setattr(requests, "get", mock_get)

    # check read works from URL
    caplog.clear()
    with caplog.at_level(logging.INFO):
        header_url, ct1_url = io.load_exchange_ctd(f"http://fakeurl/{fname}")
        assert ct1.equals(ct1_url)
        assert header_url == header
        assert f"{fname} from http link" in caplog.messages[0]

    # write fake data (x3 files) to .zip
    zname = "test_ctd.zip"
    with ZipFile(tmp_path / zname, "w") as zf:
        for n in [1, 2, 3]:
            zf.writestr(
                ZipInfo(filename=f"CTD_{n}_ct1.csv"),
                "".join(content).encode("utf8"),
                compress_type=ZIP_DEFLATED,
            )

    # check read works from .zip
    caplog.clear()
    with caplog.at_level(logging.INFO):
        header_zip, ct1_zip = io.load_exchange_ctd(tmp_path / zname)
    assert "from local file" in caplog.messages[0]
    assert "from .zip" in caplog.messages[1]
    assert all(["open file object" in caplog.messages[n] for n in [2, 3, 4]])
    assert len(ct1_zip) == 3
    assert all([ct1_zip[n].equals(ct1) for n in [0, 1, 2]])
    assert all([header_zip[n] == header for n in [0, 1, 2]])

    # check reading specific # of files from .zip works
    caplog.clear()
    with caplog.at_level(logging.INFO):
        header_zip, zipped = io.load_exchange_ctd(tmp_path / zname, n_files=2)
    assert "from local file" in caplog.messages[0]
    assert "from .zip" in caplog.messages[1]
    assert all(["open file object" in caplog.messages[n] for n in [2, 3]])
    assert len(zipped) == 2
    assert all([zipped[n].equals(ct1) for n in [0, 1]])
    assert all([header_zip[n] == header for n in [0, 1]])

    # write nested zip (including "real" data is not important)
    with ZipFile(tmp_path / "level1.zip", "w") as zf1:
        zf1.writestr(ZipInfo(filename="test.csv"), " ".encode("utf8"))
    with ZipFile(tmp_path / "level0.zip", "w") as zf0:
        zf0.write(tmp_path / "level1.zip")

    # check error on recursive .zip
    with pytest.raises(ZipImportError, match="Recursive .zip files"):
        io.load_exchange_ctd(tmp_path / "level0.zip")


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

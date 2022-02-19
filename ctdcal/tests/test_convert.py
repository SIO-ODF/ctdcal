import logging
import random
from importlib import resources
from shutil import copy

import numpy as np
import pandas as pd
import pytest

from ctdcal import convert


def make_hex_file(path, hex_chars=88, n_lines=50, seed=100):
    # generate large integer and string format into hex
    random.seed(a=seed)
    raw_hex = [f"{random.randrange(16 ** hex_chars):X}"]
    if not path.parent.exists():
        path.parent.mkdir()
    with open(path, "w") as f:
        f.write("\n".join(raw_hex * n_lines))


def make_converted_pkl(path, rows=4001, n_deck=10, n_soak=200, p_outliers=0):
    # make fake data
    downcast = np.concatenate(
        (
            [0] * n_deck,  # on deck
            np.linspace(0, 20, n_soak),
            [20] * n_soak * 3,  # fake soak
            np.linspace(20, 0, n_soak),
            np.linspace(0, 500, rows - n_soak * 5),
        )
    )
    upcast = np.linspace(downcast[-1], 0, n_deck + rows - 1)
    cond = np.concatenate(([0] * n_deck, np.linspace(30, 35, rows)))
    temp = np.concatenate(([10] * n_deck, np.linspace(5, 0, rows)))
    pump = np.concatenate(([0] * n_deck, [1] * rows))

    # construct and export DataFrame
    df = pd.DataFrame()
    df["CTDPRS"] = np.concatenate((downcast, upcast))
    df["CTDCOND1"] = np.concatenate((cond, cond[:-1][::-1]))
    df["CTDCOND2"] = np.concatenate((cond, cond[:-1][::-1]))
    df["CTDSAL"] = np.concatenate((cond, cond[:-1][::-1]))
    df["CTDTMP1"] = np.concatenate((temp, temp[:-1][::-1]))
    df["CTDTMP2"] = np.concatenate((temp, temp[:-1][::-1]))
    n_data = len(df)

    # add pressure outlier
    if p_outliers:
        bad_rows = np.random.randint(0, n_data, 2)
        df.loc[bad_rows, "CTDPRS"] = 8000

    # extra weirdness to make remove_on_deck not fail
    # (TODO: remove need for this junk!!!)
    df["pump_on"] = np.concatenate((pump, pump[:-1][::-1]))
    df["scan_datetime"] = 1645231099 + np.linspace(0, 1e3, n_data)
    df["GPSLAT"] = np.array([45] * n_data)
    df["GPSLON"] = 1645231099 + np.linspace(0, 1e3, n_data)
    df["ALT"] = 590 - df["CTDPRS"]  # altimeter goes 100m->10m and back
    df.loc[df["ALT"] > 100, "ALT"] = 100

    df.to_pickle(path)


def test_hex_to_ctd(caplog, tmp_path):
    # make fake file/path (but don't yet copy over test .xmlcon)
    raw_dir = tmp_path / "raw"
    cnv_dir = tmp_path / "converted"
    make_hex_file(tmp_path / "raw" / "90909.hex")
    cnv_dir.mkdir()

    # check warning if .xmlcon is missing
    convert.hex_to_ctd("90909", raw_dir=raw_dir, cnv_dir=cnv_dir)
    assert "90909.XMLCON does not exist" in caplog.messages[0]
    caplog.clear()

    # check file is converted correctly (with correct columns) if .xmlcon exists
    with resources.path("ctdcal.tests.data", "90909.XMLCON") as test_xml:
        copy(test_xml, raw_dir)
    convert.hex_to_ctd("90909", raw_dir=raw_dir, cnv_dir=cnv_dir)
    cnv = pd.read_pickle(cnv_dir / "90909.pkl")
    assert cnv.shape == (50, 24)
    assert cnv.columns.str.contains("CTD").sum() == 11
    assert cnv.columns.str.contains("U_DEF").sum() == 2
    assert cnv.columns.str.contains("FREE").sum() == 2
    assert cnv.columns.str.contains("GPS").sum() == 2
    misc_cols = [
        "ALT",
        "REF_PAR",
        "new_fix",
        "pressure_temp_int",
        "pump_on",
        "btl_fire",
        "scan_datetime",
    ]
    for col in misc_cols:
        assert col in cnv.columns

    # check behavior if file already exists
    with caplog.at_level(logging.INFO):
        convert.hex_to_ctd("90909", raw_dir=raw_dir, cnv_dir=cnv_dir)
    assert "Converting 1 .hex file" in caplog.messages[0]
    assert "90909.pkl already exists" in caplog.messages[1]
    caplog.clear()

    # error if .hex file does not exist
    with caplog.at_level(logging.INFO):
        # check that function works on str or list[str]
        convert.hex_to_ctd("00101", raw_dir=raw_dir, cnv_dir=cnv_dir)
        convert.hex_to_ctd(["00101"], raw_dir=raw_dir, cnv_dir=cnv_dir)
    assert "Converting 1 .hex file" in caplog.messages[0]
    assert "00101.hex does not exist" in caplog.messages[1]
    assert caplog.messages[:2] == caplog.messages[2:]


# ignore NumPy warnings caused by fake data
@pytest.mark.filterwarnings("ignore: Mean of empty slice")
@pytest.mark.filterwarnings("ignore: invalid value encountered in double_scalars")
def test_make_time_files(caplog, tmp_path):
    cnv_dir = tmp_path / "converted"
    time_dir = tmp_path / "time"
    logs_dir = tmp_path / "logs"
    for d in [cnv_dir, time_dir, logs_dir]:
        d.mkdir()

    # make fake DataFrame/.pkl and push through make_time_files
    make_converted_pkl(cnv_dir / "90909.pkl", p_outliers=2)
    with caplog.at_level(logging.DEBUG):
        convert.make_time_files(
            "90909",
            filter_cols=[],
            cnv_dir=cnv_dir,
            time_dir=time_dir,
            logs_dir=logs_dir,
        )
    assert "Generating 1 time.pkl file" in caplog.messages[0]
    assert "90909: 2 bad pressure points removed" in caplog.messages[1]
    assert "90909: Only 0.4 seconds of start pressure" in caplog.messages[2]
    assert "90909: Only 0.4 seconds of end pressure" in caplog.messages[3]
    assert (time_dir / "90909_time.pkl").exists()
    assert (logs_dir / "cast_details.csv").exists()
    assert (logs_dir / "ondeck_pressure.csv").exists()

    # check that function works on str or list[str]
    caplog.clear()
    with caplog.at_level(logging.INFO):
        convert.make_time_files("00101", cnv_dir=cnv_dir, time_dir=time_dir)
        convert.make_time_files(["00101"], cnv_dir=cnv_dir, time_dir=time_dir)

    # error if .hex file does not exist
    assert "Generating 1 time.pkl file" in caplog.messages[0]
    assert "00101.pkl does not exist" in caplog.messages[1]
    assert caplog.messages[:2] == caplog.messages[2:]


def test_make_btl_files(caplog, tmp_path):
    cnv_dir = tmp_path / "converted"
    btl_dir = tmp_path / "bottle"

    # check that function works on str or list[str]
    with caplog.at_level(logging.INFO):
        convert.make_btl_mean("00101", cnv_dir=cnv_dir, btl_dir=btl_dir)
        convert.make_btl_mean(["00101"], cnv_dir=cnv_dir, btl_dir=btl_dir)

    # error if .hex file does not exist
    assert "Generating 1 btl_mean.pkl file" in caplog.messages[0]
    assert "00101.pkl does not exist" in caplog.messages[1]
    assert caplog.messages[:2] == caplog.messages[2:]

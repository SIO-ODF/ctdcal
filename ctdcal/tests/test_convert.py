import logging
import random
from importlib import resources
from shutil import copy

import pandas as pd

from ctdcal import convert


def make_hex_file(path, hex_chars=88, n_lines=50, seed=100):
    # generate large integer and string format into hex
    random.seed(a=seed)
    raw_hex = [f"{random.randrange(16 ** hex_chars):X}"]
    if not path.parent.exists():
        path.parent.mkdir()
    with open(path, "w") as f:
        f.write("\n".join(raw_hex * n_lines))


def test_hex_to_ctd(caplog, tmp_path):
    raw_dir = tmp_path / "raw"
    cnv_dir = tmp_path / "converted"

    # make fake file/path and copy over test .xmlcon
    make_hex_file(tmp_path / "raw" / "90909.hex")
    cnv_dir.mkdir()

    # check warning if .xmlcon  is missing
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


def test_make_time_files(caplog, tmp_path):
    cnv_dir = tmp_path / "converted"
    time_dir = tmp_path / "time"

    # check that function works on str or list[str]
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

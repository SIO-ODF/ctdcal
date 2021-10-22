import logging

from ctdcal import convert


def test_hex_to_ctd(caplog, tmp_path):
    raw_dir = tmp_path / "raw"
    cnv_dir = tmp_path / "converted"

    # check that function works on str or list[str]
    with caplog.at_level(logging.INFO):
        convert.hex_to_ctd("00101", raw_dir=raw_dir, cnv_dir=cnv_dir)
        convert.hex_to_ctd(["00101"], raw_dir=raw_dir, cnv_dir=cnv_dir)

    # error if .hex file does not exist
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

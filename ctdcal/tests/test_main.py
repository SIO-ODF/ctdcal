import logging
from pathlib import Path, PosixPath, WindowsPath

from click.testing import CliRunner

import ctdcal.__main__ as main


def test_cli():
    runner = CliRunner()
    result = runner.invoke(main.cli)
    assert "cruise-report" in result.output
    assert "import" in result.output
    assert "init" in result.output
    assert "process" in result.output
    assert "qc" in result.output


def test_debug(tmp_path):
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        # default options
        result = runner.invoke(main.cli, "process")
        assert "Debug mode off" in result.output

        # no debug
        result_no_debug = runner.invoke(main.cli, ["--no-debug", "process"])
        assert "Debug mode off" in result_no_debug.output

        # debug (run last so Logger("ctdcal").level=NOTSET for other tests)
        result_debug = runner.invoke(main.cli, ["--debug", "process"])
        assert "Debug mode on" in result_debug.output


def test_init(tmp_path, caplog):
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        # check first init call creates dirs successfully
        with caplog.at_level(logging.INFO):
            result = runner.invoke(main.init)
        assert "Building default /data/ directories" in caplog.messages[0]
        assert result.exit_code == 0
        assert not result.exception

        # check "folders already exist" error
        caplog.clear()
        with caplog.at_level(logging.INFO):
            result_rerun = runner.invoke(main.init)
        assert "Building default /data/ directories" in caplog.messages[0]
        assert result_rerun.exit_code == 1
        assert isinstance(result_rerun.exception, FileExistsError)
        if isinstance(Path.cwd(), PosixPath):
            assert "data/ssscc" in result_rerun.exception.filename
        elif isinstance(Path.cwd(), WindowsPath):
            assert "data\\ssscc" in result_rerun.exception.filename


def test_import_data():
    runner = CliRunner()
    result = runner.invoke(main.import_data)
    assert result.exit_code == 1
    assert isinstance(result.exception, NotImplementedError)


def test_process(tmp_path, caplog):
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        # default options
        with caplog.at_level(logging.INFO):
            result = runner.invoke(main.process)
        assert result.exit_code == 1
        # assert isinstance(result.exception, FileNotFoundError)
        # assert "ssscc.csv" in result.exception.filename
        # assert "generating from .hex file list" in caplog.messages[-1]

        # PMEL option
        with caplog.at_level(logging.INFO):
            result_PMEL = runner.invoke(main.process, ["-g", "PMEL"])
        assert result_PMEL.exit_code == 1
        assert isinstance(result_PMEL.exception, NotImplementedError)


def test_cruise_report(tmp_path, caplog):
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        with caplog.at_level(logging.INFO):
            result = runner.invoke(main.cruise_report)
        assert isinstance(result.exception, FileNotFoundError)
        assert "report_data.csv" in result.exception.filename


def test_quick_convert(tmp_path, caplog):
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        with caplog.at_level(logging.INFO):
            result = runner.invoke(main.quick_convert)
        assert result.exit_code == 0
        assert "Found 0 .cnv files" in caplog.messages[0]

        caplog.clear()
        with caplog.at_level(logging.INFO):
            result_ctd = runner.invoke(main.quick_convert, ["-f", "ct1"])
        assert result_ctd.exit_code == 0
        assert "Found 0 .cnv files" in caplog.messages[0]

        caplog.clear()
        with caplog.at_level(logging.INFO):
            result_btl = runner.invoke(main.quick_convert, ["-f", "hy1"])
        assert result_btl.exit_code == 1
        assert isinstance(result_btl.exception, NotImplementedError)

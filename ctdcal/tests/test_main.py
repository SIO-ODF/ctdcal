import logging

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
        assert "data/ssscc" in result_rerun.exception.filename


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
        assert isinstance(result.exception, FileNotFoundError)
        assert "ssscc.csv" in result.exception.filename
        assert "generating from .hex file list" in caplog.messages[0]

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

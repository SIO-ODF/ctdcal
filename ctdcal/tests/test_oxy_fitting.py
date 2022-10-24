import logging
from pathlib import Path

from ctdcal import oxy_fitting


def test_gather_oxy_params(caplog, tmp_path):
    # make fake file path
    f_path = tmp_path / "90909"
    # breakpoint()
    assert not Path(tmp_path / "90909").exists()
    with caplog.at_level(logging.INFO):
        oxy_params = oxy_fitting.gather_oxy_params(f_path)
        assert "Failed to load" in caplog.messages[0]
        assert oxy_params.isnull().values.all()

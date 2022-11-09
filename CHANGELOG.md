# Change Log
This project adheres to [Semantic Versioning](https://semver.org/).

## v0.1.4b (unreleased)

### Added
* Logging messages of all levels are now logged to `ctdcal.log` file
* `ctdcal` CLI now has a `--debug` flag to display all logging levels (instead of WARNING and above)
* `ctdcal.io` has two new functions for reading Exchange data files from local filesystem or URL, `load_exchange_btl` and `load_exchange_ctd` (see [CCHDO File Formats](https://cchdo.ucsd.edu/formats))
* Tutorial notebooks on how to use `ctdcal` can be found in the [docs](https://ctdcal.readthedocs.io/en/latest/)
* `ctdcal quick-convert` CLI command to create Exchange ct1 files from Sea-Bird .cnv files, without performing any calibration
* `ctdcal qc` CLI command to launch interactive data flagging tool
* Data flagging tool now has an exit button to halt the local server process
* Added Wetlabs Fluorometer to `convert.py`
* Start testing `process_bottle` and `oxy_fitting` modules

### Changed
* By default, only logging levels WARNING and above will be displayed in terminal (see `--debug` addition above)
* Fixed bug when loading single station `ssscc.csv` files into `fit_ctd`

## v0.1.3b (2021-10-21)

### Added
* Created `ctdcal.io` module to hold all non-ODF-specific reading/writing functions
* `isort` pre-commit hook
* Set up [Codecov](https://app.codecov.io/gh/cchdo/ctdcal/)
* `equations_sbe` conversion functions will now return a list of missing coefficients (instead of first missing in equation)
* `equations_sbe` now NaNs out zero-values in all input frequency arrays (applies to sbe3, sbe4, and sbe9 conversions)

### Changed
* Renamed `master` branch to `main`
* Fix sphinx/RTD version bug by importing ctdcal._version instead of importlib.metadata.version("ctdcal")
* Fix GitHub Action `run-tests` bug which triggered twice when pushing tagged commits
* `ctd_plots` functions now return axis handle if filename is not given
* `_intermediate_residual_plot()` is now wrapped around `residual_vs_pressure()` instead of being a duplicate function
* `report_ctd` merged into new `ctdcal.io` module

### Removed
* `cmocean` is no longer a package dependency
* Remove hardcoded cruise report plot code from `ctd_plots` module in favor of `ctdcal/scripts/cruise_report.py`
* Outdated code in top-level `old` folder has been removed in favor of `SBEReader` class and `equations_sbe` module
* Outdated `merge_codes` module has been removed

## v0.1.2b (2021-09-29)

### Added
* Initialize testing framework using `pytest`
* Add testing for `fit_ctd` and `flagging` modules
* Run tests on multiple OS and Python version using GitHub Actions
* More status badges

### Removed
* Deleted outdated `unitTest.sh` testing file

## v0.1.1b (2021-09-23)

### Added
* Basic guide to making package releases
* GitHub Action which publishes tagged versions to [TestPyPI](https://test.pypi.org/project/ctdcal/)

### Changed
* Fixed configuration loader bug by including `user_settings.yaml` in packaged files
* Update installation steps to use `pip` instead of `git clone`
* Added `pip install ctdcal` step to `publish-to-pypi` workflow

## v0.1.0b (2021-09-23)

### Added
* Github Action which publishes released versions to [PyPI](https://pypi.org/project/ctdcal/)

### Changed
* Switch from using `versioneer` to `setuptools_scm`

### Removed
* Deleted `.flake8` settings file in lieu of `setup.cfg`

## v0.1.0a (2021-09-22)

### Added
* Example Jupyter notebook functionality in documentation
* Generate stub pages for all functions/classes/etc.

### Changed
* Improved ReadTheDocs documentation
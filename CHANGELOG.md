# Change Log

## v0.1.3b (unreleased)

## v0.1.2b (09/29/2021)

### Added
* Initialize testing framework using `pytest`
* Add testing for `fit_ctd` and `flagging` modules
* Run tests on multiple OS and Python version using GitHub Actions
* More status badges

### Removed
* Deleted outdated `unitTest.sh` testing file

## v0.1.1b (09/23/2021)

### Added
* Basic guide to making package releases
* GitHub Action which publishes tagged versions to `TestPyPI <https://test.pypi.org/project/ctdcal/>`_

### Changed
* Fixed configuration loader bug by including `user_settings.yaml` in packaged files
* Update installation steps to use `pip` instead of `git clone`
* Added `pip install ctdcal` step to `publish-to-pypi` workflow

## v0.1.0b (09/23/2021)

### Added
* Github Action which publishes released versions to `PyPI <https://pypi.org/project/ctdcal/>

### Changed
* Switch from using `versioneer` to `setuptools_scm`

### Removed
* Deleted `.flake8` settings file in lieu of `setup.cfg`

## v0.1.0a (09/22/2021)

### Added
* Example Jupyter notebook functionality in documentation
* Generate stub pages for all functions/classes/etc.

### Changed
* Improved ReadTheDocs documentation
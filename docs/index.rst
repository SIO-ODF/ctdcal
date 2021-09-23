.. ctdcal documentation master file, created by
   sphinx-quickstart on Wed Dec 20 15:09:16 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:

   Jupyter tutorials <tutorials>
   API reference <api>
   Contribute on GitHub <https://github.com/cchdo/ctdcal>

************************************************
ctdcal: Process and calibrate CTD data in Python
************************************************

This package is under development by the `Oceanographic Data Facility (ODF) <https://scripps.ucsd.edu/ships/shipboard-technical-support/odf>`_ group at `Scripps Institution of Oceanography (SIO) <https://scripps.ucsd.edu/>`_. ODF uses data processing software and techniques that were developed internally and provide exceptional fitting of CTD and hydrographic data. The software is regularly updated to take full advantage of the latest coding standards and improvements.

The software is all-encompassing with the ability to process raw data, perform quality control with automated data flagging, easily visualize the finalized data, and export to `netCDF <https://www.unidata.ucar.edu/software/netcdf/>`_ using standardized conventions. Users have the ability to get the highest quality data out of their measurements.

Installation
============
``ctdcal`` can be installed using pip::
   
   $ pip install ctdcal


CLI usage
=========
Initialize data folders
-----------------------
Initialize default `/data/` folders by running::

   $ ctdcal init

(Future versions of ctdcal are planned have more robust init options/flags/etc.)

Import and process data
-----------------------
To process data, copy over raw `.hex` and `.xmlcon` files into `/data/raw/` and reference data to their appropriate folder (`oxygen`, `reft`, `salt`).

Users can process their data with individual ctdcal functions or try::

   $ ctdcal process [--group ODF]

to process using ODF procedures.

Package usage
=============
Explore user settings
---------------------
Most ctdcal functions get settings from `user_settings.yaml` and subsequently `config.py`. Call the configuration loader to explore default settings::

   from ctdcal import get_ctdcal_config
   cfg = get_ctdcal_config()

   # directories for I/O purposes
   print(cfg.dirs)
   print(cfg.fig_dirs)

   # experiment-specific settings (e.g., expocode, CTD serial number) from user_settings.yaml
   print(cfg.settings)

   # dictionary mapping of short/long column names
   print(cfg.columns)

As ctdcal continues to be developed, more robust `tutorials <https://ctdcal.readthedocs.io/en/latest/tutorials.html>`_ will be added to `our documentation <https://ctdcal.readthedocs.io/en/latest/>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

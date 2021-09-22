.. ctdcal documentation master file, created by
   sphinx-quickstart on Wed Dec 20 15:09:16 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:

   Jupyter tutorials <tutorials>
   API reference <api>
   Contribute on GitHub <https://github.com/cchdo/ctdcal>

ctdcal: Process and calibrate CTD data in Python
================================================

This package is under development by the `Oceanographic Data Facility (ODF) <https://scripps.ucsd.edu/ships/shipboard-technical-support/odf>`_ group at `Scripps Institution of Oceanography (SIO) <https://scripps.ucsd.edu/>`_. ODF uses data processing software and techniques that were developed internally and provide exceptional fitting of CTD and hydrographic data. The software is regularly updated to take full advantage of the latest coding standards and improvements.

The software is all-encompassing with the ability to process raw data, perform quality control with automated data flagging, easily visualize the finalized data, and export to `netCDF <https://www.unidata.ucar.edu/software/netcdf/>`_ using standardized conventions. Users have the ability to get the highest quality data out of their measurements.

Installation
============
``ctdcal`` can be installed using Git and pip as follows:

Pull down the latest version of ctdcal:
   >>> git clone https://github.com/cchdo/ctdcal.git

Change directories to the top-level ctdcal:
   >>> cd ctdcal

Create a new virtual environment with your preferred environment manager and install dependencies with pip:
   >>> pip install .

.. note::
   Note: there is an occasional (conda?) bug where CLI tools are not immediately accessible after install – this can usually be remedied by deactivating and reactiving the virtual environment.

Initialize default `/data/` folders by running:
   >>> ctdcal init

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

# ctdcal

This is an initial draft of the ctdcal installation process in its current form.

### Clone repository
Pull down the latest version of ctdcal:
```
git clone https://github.com/cchdo/ctdcal.git
```

### Install ctdcal and dependencies
Change directories to the top-level ctdcal:
```
cd ctdcal
```

Create a new virtual environment with your preferred environment manager and install with pip:
```
pip install .
```

Note: there is an occasional (conda?) bug where CLI tools are not immediately accessible after install – this can usually be remedied by deactivating and reactiving the virtual environment.

Initialize default `/data/` folders by running:
```
ctdcal init
```

(Future versions of ctdcal are planned have more robust init options/flags/etc.)

### Import and process data
To process data, copy over raw `.hex` and `.xmlcon` files into `/data/raw/` and reference data to their appropriate folder (`oxygen`, `reft`, `salt`).

Users can process their data with individual ctdcal functions or try:
```
ctdcal process [--group ODF]
```
to process using ODF procedures.
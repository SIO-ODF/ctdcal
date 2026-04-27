"""
Extracts parameters and cal coeffs from instrument xml config file.

2025-04-08
Currently only used for getting SBE43 coeffs for oxy calibration stage.
"""
import xml.etree.cElementTree as ET
from pathlib import Path

from ctdcal.common import list_to_file, validate_dir


def parse_coeffs(casts, rawdir, caldir):
    validate_dir(caldir, create=True)
    get_sbe43_coef(casts, rawdir, caldir)


def get_sbe43_coef(casts, rawdir, outdir):
    """
    Get SBE oxygen coefficients from raw .xmlcon files.
    Defaults to using first station in ssscc.csv file.

    Returns the following tuple of coefficients: Soc, offset, Tau20, Tcor, E
    """
    for cast in casts:
        infile = Path(rawdir, '%s.XMLCON' % cast)

        tree = ET.parse(infile)
        root_eq0 = tree.find(".//CalibrationCoefficients[@equation='0']")  # Owens-Millard
        root_eq1 = tree.find(".//CalibrationCoefficients[@equation='1']")  # SBE equation

        coefs = {c.tag: float(c.text) for c in root_eq1}
        coefs["Tcor"] = float(root_eq0.find("Tcor").text)  # only coef needed from eq0
        keep_keys = ["Soc", "offset", "Tau20", "Tcor", "E"]

        outfile = Path('%s_sbe43_coeffs.csv' % cast)
        list_to_file(outfile, outdir, tuple(str(coefs[key]) for key in keep_keys))

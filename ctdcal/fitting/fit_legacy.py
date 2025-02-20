import logging
from pathlib import Path

import yaml

from ctdcal import get_ctdcal_config

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def load_fit_yaml(fname=f"{cfg.dirs['logs']}fit_coefs.yaml", to_object=False):
    """Load polynomial fit order information from .yaml file."""

    if not Path(fname).exists():
        log.warning("Warning: Coefficients fit order YAML does not exist. Generating from scratch...")
        generate_yaml()

    with open(fname, "r") as f:
        ymlfile = yaml.safe_load(f)

    if to_object:
        return type("ymlfile", (object,), ymlfile)
    else:
        return ymlfile


def generate_yaml(fname="fit_coefs.yaml", outdir=cfg.dirs['logs']):
    """
    Create a default coeff. yaml file.
    """
    data = {
        't1': {
            'ssscc_t1': {
                'P_order': 1,
                'T_order': 0,
                'zRange': "1000:6000"
            }
        },
        'c1': {
            'ssscc_c1': {
                'P_order': 1,
                'T_order': 0,
                'C_order': 0,
                'zRange': "1000:6000"
            }
        },
        't2': {
            'ssscc_t1': {
                'P_order': 1,
                'T_order': 0,
                'zRange': "1000:6000"
            }
        },
        'c2': {
            'ssscc_c1': {
                'P_order': 1,
                'T_order': 0,
                'C_order': 0,
                'zRange': "1000:6000"
            }
        }
    }

    # Write the data to a YAML file
    with open(outdir + fname, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)


def write_fit_yaml():
    """For future use with automated fitting routine(s).
    i.e., iterate to find best fit parameters, save to file"""
    pass

"""
Minimal setup.py for building ctd processing package.
"""

import sys
import glob

from setuptools import setup

import versioneer


# Check Python version.
if sys.version_info < (3, 5):
    pip_message = ('This may be due to an out of date pip. '
                   'Make sure you have pip >= 9.0.1.')
    try:
        import pip
        pip_version = tuple([int(x) for x in pip.__version__.split('.')[:3]])
        if pip_version < (9, 0, 1):
            pip_message = ('Your pip version is out of date, '
                           'please install pip >= 9.0.1. '
                           'pip {} detected.').format(pip.__version__)
        else:
            # pip is new enough - it must be something else.
            pip_message = ''
    except Exception:
        pass

    error = """
Latest ctd does not support Python < 3.5.
Python {py} detected.
{pip}
""".format(py=sys.version_info, pip=pip_message)

    print(error, file=sys.stderr)
    sys.exit(1)

scripts = glob.glob('scripts/*.py') + glob.glob('scripts/institutions/*/*.py')

config = dict(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    download_url='https://pypi.python.org/pypi/gsw/',
    scripts=scripts,
)

setup(**config)

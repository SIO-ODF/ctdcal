"""
Minimal setup.py for building ctd processing package.
"""

import glob

from setuptools import setup

import versioneer

# scripts = glob.glob("scripts/*.py") + glob.glob("scripts/institutions/*/*.py")

config = dict(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    download_url="https://pypi.python.org/pypi/gsw/",
    # scripts=scripts,
)

setup(**config)

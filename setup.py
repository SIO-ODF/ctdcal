"""
Minimal setup.py for building ctd processing package.
"""

import os
import sys
import glob

import pkg_resources
from setuptools import Extension, setup

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


rootpath = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    return open(os.path.join(rootpath, *parts), 'r').read()


LICENSE = read('LICENSE')
long_description = read('README.md')


scripts = glob.glob('scripts/*.py')

config = dict(
    name='ctdcal',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=['ctdcal'],
    author=['Joseph Gum', 'Andrew Barna'],
    author_email='jgum@ucsd.edu',
    description='CTD and bottle data processing package from UCSD ODF',
    long_description=long_description,
    license=LICENSE,
    url='https://github.com/somts/odf-ctd-proc',
    download_url='https://pypi.python.org/pypi/gsw/',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
    ],
    python_requires='>=3.5',
    platforms='any',
    keywords=['oceanography', 'seawater', 'TEOS-10', 'ctd', 'calibration'],
    install_requires=[
                      'numpy',
                      'matplotlib',
                      'scipy',
                      'pandas==0.20.3',
                      'gsw>=3.2',
                      ],
    scripts=scripts,
)

setup(**config)

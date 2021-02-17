"""
This script will initialize the data directory for ctdcal according to ODF usage.
Joseph Gum, January 22, 2019
"""
import pathlib
import sys


def build_dirs():
    data_directory = "./data/"
    directories_to_be_made = [
        "raw",
        "salt",
        "bottle",
        "converted",
        "logs",
        "oxygen",
        "pressure",
        "quality_codes",
        "reft",
        "time",
    ]

    # will make parents directory if not missing, will not throw an error if directory already exists
    for x in directories_to_be_made:
        pathlib.Path(data_directory + x).mkdir(parents=True, exist_ok=True)


def main(argv):
    build_dirs()


if __name__ == "__main__":
    main(sys.argv[1:])

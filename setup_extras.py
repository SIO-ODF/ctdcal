from pathlib import Path

DIR_NAMES = [
        "raw",
        "reft",
        "salt",
        "oxygen",
        ]


def make_data_dirs():
    """Build folders for raw and reference data"""
    base_name = Path("./data/")
    for sub_dir in DIR_NAMES:
        path = base_name / sub_dir
        path.mkdir(parents=True)

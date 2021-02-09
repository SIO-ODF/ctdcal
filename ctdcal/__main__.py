import click
from pathlib import Path


@click.group()
def cli():
    """The ctdcal command creates and manipulates data directories

    Documentation: tbd
    """
    pass


@cli.command()
def init():
    """Setup data folder with appropriate subfolders"""

    DEFAULT_DIRS = [
        "raw",
        "converted",
        "time",
        "bottle",
        "reft",
        "salt",
        "oxygen",
        "logs",
        "ssscc",
        "pressure",
    ]

    print(f"Building default data directories: \n {*DEFAULT_DIRS,}")

    base_name = Path("./data/")
    for sub_dir in DEFAULT_DIRS:
        path = base_name / sub_dir
        path.mkdir(parents=True)


@cli.command("import")  # click workaround to get a command named 'import'
def import_data():
    """Import data from given folder into ctdcal for processing"""
    # something like this?
    # ctdcal import _path_

    # TODO: smart imports based on file ext? .hex, .xmlcon, .cap
    # NOTE: ODF file types vs. others (oxygen, salt)

    pass


@cli.command()
def process():
    """Process data using certain group's methodology"""
    # something like this?
    # ctdcal process -[odf / pmel / ...] -[bottle / ctd / all]
    pass


if __name__ == "__main__":
    cli()

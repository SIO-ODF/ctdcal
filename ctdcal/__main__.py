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


if __name__ == "__main__":
    cli()

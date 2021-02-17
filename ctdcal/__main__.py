import click
from pathlib import Path
from .scripts.odf_process_all import odf_process_all


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
@click.option("-g", "--group", type=click.Choice(["ODF", "PMEL"], case_sensitive=False))
@click.option(
    "-t", "--type", type=click.Choice(["bottle", "ctd", "all"], case_sensitive=False)
)
def process(group, type):
    """Process data using a particular group's methodology"""

    if group is None:
        print("No group specificed, default to processing as ODF")
        group = "ODF"

    if group == "ODF":
        odf_process_all()
    elif group == "PMEL":
        # pmel_process()
        pass
    else:
        print("Group not specified")
        pass


if __name__ == "__main__":
    cli()

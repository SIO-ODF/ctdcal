from pathlib import Path

import click
import logging

# Rich handling
# from rich.logging import RichHandler
# from rich.console import Console
# handler = logging.StreamHandler()
# handler.addFilter(logging.Filter("ctdcal"))  # filter out msgs from other modules
# FORMAT = "%(funcName)s: %(message)s"
# logging.basicConfig(
#     level="INFO",
#     format=FORMAT,
#     datefmt="[%X]",
#     handlers=[handler, RichHandler(console=Console(stderr=True))],
# )

# log = logging.getLogger(__name__)

handler = logging.StreamHandler()
handler.addFilter(logging.Filter("ctdcal"))  # filter out msgs from other modules
FORMAT = "%(funcName)s: %(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[handler],
    # handlers=[RichHandler(console=Console(stderr=True))],
)

log = logging.getLogger(__name__)


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

    log.info(f"Building default data directories: \n {*DEFAULT_DIRS,}")

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
@click.option(
    "-g",
    "--group",
    type=click.Choice(["ODF", "PMEL"], case_sensitive=False),
    default="ODF",
)
@click.option(
    "-t", "--type", type=click.Choice(["bottle", "ctd", "all"], case_sensitive=False)
)
def process(group, type):
    """Process data using a particular group's methodology"""

    if group == "ODF":
        from .scripts.odf_process_all import odf_process_all

        odf_process_all()
    elif group == "PMEL":
        # pmel_process()
        pass


@cli.command()
def cruise_report():
    """Generate bottle residual figures for cruise report"""

    from .scripts.cruise_report import cruise_report_residuals
    cruise_report_residuals()


if __name__ == "__main__":
    cli()

from pathlib import Path

import click
import logging
import time
from . import get_ctdcal_config

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
    level="NOTSET",
    format=FORMAT,
    datefmt="[%X]",
    handlers=[handler],
    # handlers=[RichHandler(console=Console(stderr=True))],
)

log = logging.getLogger(__name__)
cfg = get_ctdcal_config()


@click.group()
def cli():
    """The ctdcal command creates and manipulates data directories

    Documentation: tbd
    """
    pass


@cli.command()
def init():
    """Setup data folder with appropriate subfolders"""

    log.info(f"Building default /data/ directories: \n {*cfg.dirs.keys(),}")

    for sub_dir in cfg.dirs.values():
        Path(sub_dir).mkdir(parents=True)
    import shutil
    shutil.copy2("fit_coefs.yaml", "./data/logs/fit_coefs.yaml")

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

    t = time.time()
    if group == "ODF":
        from .scripts.odf_process_all import odf_process_all

        odf_process_all()
    elif group == "PMEL":
        # pmel_process()
        pass
    elapsed = time.time() - t
    log.info("Processing complete: " + str(round(elapsed/60)) + " minutes and " + str(round(elapsed)) + " seconds.")

@cli.command()
def cruise_report():
    """Generate bottle residual figures for cruise report"""

    from .scripts.cruise_report import cruise_report_residuals

    cruise_report_residuals()

@cli.command()  #   MK
def qc():
    """Launch interactive data flagging web app for QA/QC"""
    from bokeh.application import Application
    from bokeh.application.handlers.script import ScriptHandler
    from bokeh.server.server import Server
    from tornado.ioloop import IOLoop
    io_loop = IOLoop.current()
    bokeh_app = Application(ScriptHandler(filename="ctdcal/tools/data_qc.py"))
    server = Server(bokeh_app, io_loop=io_loop)
    server.start()
    server.show("/")    #   Blank tab, no output
    io_loop.start()


if __name__ == "__main__":
    cli()

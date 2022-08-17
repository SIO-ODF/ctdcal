import logging
from importlib import resources
from pathlib import Path

import click
import time
from bokeh.application import Application
from bokeh.application.handlers.script import ScriptHandler
from bokeh.server.server import Server
from tornado.ioloop import IOLoop

from . import get_ctdcal_config

## logging settings
# terminal output
stream = logging.StreamHandler()
stream.setLevel(logging.WARNING)
stream.addFilter(logging.Filter("ctdcal"))  # filter out msgs from other modules

# ctdcal.log output
logfile_FORMAT = "%(asctime)s | %(funcName)s |  %(levelname)s: %(message)s"
logfile = logging.FileHandler("ctdcal.log")
logfile.setLevel(logging.NOTSET)
logfile.addFilter(logging.Filter("ctdcal"))  # filter out msgs from other modules
logfile.setFormatter(logging.Formatter(logfile_FORMAT))

# global configs
FORMAT = "%(funcName)s: %(levelname)s: %(message)s"
logging.basicConfig(
    level=logging.NOTSET,
    format=FORMAT,
    datefmt="[%X]",
    handlers=[stream, logfile],
)

log = logging.getLogger(__name__)
cfg = get_ctdcal_config()


@click.group()
@click.option("--debug/--no-debug", default=False)
def cli(debug):
    """The ctdcal command creates and manipulates data directories

    Documentation: tbd
    """
    if debug:
        click.echo("Debug mode on (logging all levels)")
        logging.getLogger("ctdcal").setLevel(logging.NOTSET)
    else:
        click.echo("Debug mode off (logging 'WARNING' and higher levels)")
        logging.getLogger("ctdcal").setLevel(logging.WARNING)


@cli.command()
@click.option(
    "-r",
    "--rebuild",
    default=0,
    help="Option to save the current /data folder and rebuild from scratch. 0 (default) = No duplication. 1 = Fresh data folder. 2 = Copy raw only. 3 = Copy all non-processed files.",
)
def init(rebuild=0):
    """Setup data folder with appropriate subfolders"""

    import shutil

    if rebuild != 0:
        import os

        source = "data/"
        target = "".join(map(str, time.localtime()[0:6])) + "/"
        log.info(f"Storing data within {target}")
        os.rename(source, target)

    log.info(f"Building default /data/ directories: \n {*cfg.dirs.keys(),}")
    for sub_dir in cfg.dirs.values():
        if rebuild == 2:
            if sub_dir == "data/raw/":
                log.info(f"Duplicating raw data only")
                shutil.copytree(target + "/raw", source + "raw")
            else:
                Path(sub_dir).mkdir(parents=True)
        elif rebuild == 3:
            ignore_csv = lambda d, files: [
                f
                for f in files
                if os.path.isfile(os.path.join(d, f)) and f[-4:] == ".csv"
            ]
            if sub_dir == "data/ssscc/":
                log.info(f"Duplicating all data for a fresh run")
                shutil.copytree(target + "/ssscc", source + "ssscc")
                shutil.copy2(
                    target + "ssscc.csv", source + "ssscc.csv"
                )  #   Also copy ssscc file in parent
            elif sub_dir == "data/raw/":
                log.info("Getting raw folder...")
                shutil.copytree(target + "/raw", source + "raw")
            elif sub_dir == "data/salt/":
                log.info("Getting unprocessed salt files...")
                shutil.copytree(target + "/salt", source + "salt", ignore=ignore_csv)
            elif sub_dir == "data/reft/":
                log.info("Getting unprocessed reference temperature files...")
                shutil.copytree(target + "/reft", source + "reft", ignore=ignore_csv)
            elif sub_dir == "data/oxygen/":
                log.info("Getting oxygen files...")
                shutil.copytree(target + "/oxygen", source + "oxygen")
            elif sub_dir == "data/logs/":
                log.info("Getting core logs...")
                Path(sub_dir).mkdir(parents=True)
                shutil.copy2(
                    target + "/logs/manual_depth_log.csv",
                    source + "logs/manual_depth_log.csv",
                )
            else:
                Path(sub_dir).mkdir(parents=True)
        else:
            Path(sub_dir).mkdir(parents=True)

    shutil.copy2("fit_coefs.yaml", "./data/logs/fit_coefs.yaml")
    shutil.copy2("o2flasks.vol", "./data/oxygen/o2flasks.vol")
    log.info("All data directories successfully created")


@cli.command("import")  # click workaround to get a command named 'import'
def import_data():
    """Import data from given folder into ctdcal for processing"""
    # something like this?
    # ctdcal import _path_

    # TODO: smart imports based on file ext? .hex, .xmlcon, .cap
    # NOTE: ODF file types vs. others (oxygen, salt)

    raise NotImplementedError


@cli.command()
@click.option(
    "-g",
    "--group",
    type=click.Choice(["ODF", "PMEL", "WHOI"], case_sensitive=False),
    default="ODF",
)
# @click.option(
#     "-t",
#     "--type",
#     type=click.Choice(["bottle", "ctd", "all"], case_sensitive=False),
#     default="all",
# )
def process(group):
    """Process data using a particular group's methodology"""

    if group == "ODF":
        from .scripts.odf_process_all import odf_process_all

        odf_process_all()
    elif group == "PMEL":
        # pmel_process()
        raise NotImplementedError

    elif group == "WHOI":
        from .scripts.whoi_process_all import whoi_process_all

        whoi_process_all()


@cli.command()
def cruise_report():
    """Generate bottle residual figures for cruise report"""

    #   2022 OSNAP: Need additional figures for report

    from .scripts.cruise_report import cruise_report_residuals

    cruise_report_residuals()


@cli.command()
def qc():  # pragma: no cover
    """Launch interactive data flagging web app for QA/QC"""
    io_loop = IOLoop.current()
    with resources.path("ctdcal.tools", "data_qc.py") as fname:
        bokeh_app = Application(ScriptHandler(filename=fname))
    server = Server(bokeh_app, io_loop=io_loop)
    server.start()
    server.show("/")
    io_loop.start()


@click.option(
    "-f",
    "--file",
    type=click.Choice(["hy1", "ct1", "all"], case_sensitive=False),
    default="all",
)
@cli.command()
def quick_convert(file):
    """Convert Sea-Bird .cnv files to Exchange CTD (ct1) and bottle (hy1) files."""
    from .scripts.quick_convert import cnv_to_ct1

    if file == "ct1":
        cnv_to_ct1()
    elif file == "hy1":
        # cnv_to_hy1()
        raise NotImplementedError
    elif file == "all":
        cnv_to_ct1()
        # cnv_to_hy1()

    # should we have non-Exchange output options?


if __name__ == "__main__":
    cli()

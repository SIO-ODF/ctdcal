"""
Process all CTD and bottle data using ODF routines.
"""
import logging
from pathlib import Path

from ctdcal.common import load_user_config, load_fit_groups, make_cast_id_list
from ctdcal.fitting.fit_ctd import calibrate_temp, calibrate_cond, load_time_all, calibrate_pressure
from ctdcal.fitting.fit_oxy import calibrate_oxy
from ctdcal.fitting.fit_oxy_rinko import calibrate_oxy as calibrate_rinko
from ctdcal.formats.exchange import export_exchange
from ctdcal.parsers.parse_ctd_xmlcon import parse_coeffs
from ctdcal.processors.cast_tools import make_time_files
from ctdcal.processors.convert_legacy import hex_to_ctd
from ctdcal.processors.proc_bottle import make_btl_files, load_btl_all
from ctdcal.processors.proc_oxy_ctd import prepare_oxy
from ctdcal.processors.proc_oxy_odf import proc_oxy
from ctdcal.processors.proc_reft import proc_reft
from ctdcal.processors.proc_reft_rbr import proc_reft as proc_rbr
from ctdcal.processors.proc_salt_odf import proc_salt
from ctdcal.reporting.report_odf import make_depth_log, export_report_data, make_odf_report_table
from ctdcal.scripts.cruise_report_odf import cruise_report

# setup logging if not using ctdcal CLI
logger = logging.getLogger(__name__)
stream = logging.StreamHandler()
stream.setLevel(logging.WARNING)
stream.addFilter(logging.Filter("ctdcal"))  # filter out msgs from other modules
logfile_FORMAT = "%(asctime)s | %(funcName)s |  %(levelname)s: %(message)s"
logfile = logging.FileHandler("ctdcal.log")
logfile.setLevel(logging.NOTSET)
logfile.addFilter(logging.Filter("ctdcal"))  # filter out msgs from other modules
logfile.setFormatter(logging.Formatter(logfile_FORMAT))
FORMAT = "%(funcName)s: %(levelname)s: %(message)s"
logging.basicConfig(
    level=logging.DEBUG,
    format=FORMAT,
    encoding='utf-8',
    datefmt="[%X]",
    handlers=[stream, logfile],
)

USERCONFIG = '/Users/als026/data/i09n_2025/cfg_i09n.yaml'
EXCHANGECONFIG = '/Users/als026/data/i09n_2025/i09n_exchange.yaml'
cfg = load_user_config(USERCONFIG)

# Runtime flags:
# if this flag is set to True, the calibration routines will be bypassed
skip_calibrate = True
# if this and the above flags are set to True, the export routines will be bypassed
skip_export = True
# if this flag is set to True, only the cruise report processing will execute
process_cruise_report = True


def odf_process_all():

    #####
    # Step 0: Load and define necessary variables
    #####

    datadir = cfg.datadir
    fit_coeffs = cfg.fit_coeffs

    # single CTD setup
    inst = 'ctd'
    logdir = Path(datadir, 'log')
    ssscc_dir = Path(datadir, 'ssscc', inst)
    rawdir = Path(datadir, 'raw', inst)
    caldir = Path(datadir, 'cal', inst)
    cnvdir = Path(datadir, 'cnv', inst)
    timedir = Path(datadir, 'time', inst)
    btldir = Path(datadir, 'btl', inst)
    reportdir = Path(datadir, 'report', inst)
    plotdir = Path(datadir, 'plots', inst)
    outdir = Path(datadir, 'save')
    flagfile = Path(datadir, 'flag', cfg.bottleflags_man)
    salt_rawdir = Path(datadir, 'raw', 'salt')
    salt_cnvdir = Path(datadir, 'cnv', 'salt')
    salt_figdir = Path(plotdir, 'cond')
    reft_rawdir = Path(datadir, 'raw', 'reft')
    reft_parsedir = Path(datadir, 'parsed', 'reft')
    reft_cnvdir = Path(datadir, 'cnv', 'reft')
    reft_figdir = Path(plotdir, 'temp')
    oxy_rawdir = Path(datadir, 'raw', 'oxy')
    oxy_cnvdir = Path(datadir, 'cnv', 'oxy')
    oxy_figdir = Path(plotdir, 'oxy')
    rinko_figdir = Path(plotdir, 'rinko')

    rbr_rawdir = Path('/Users/als026/data/i09n_2025/rbr/raw')
    rbr_parseddir = Path(datadir, 'parsed', 'rbr')
    rbr_cnvdir = Path(datadir, 'cnv', 'rbr')

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file or make it from all raw data files
    ssscc_list = make_cast_id_list(rawdir)

    # convert raw .hex files
    hex_to_ctd(ssscc_list, rawdir, cnvdir)
    parse_coeffs(ssscc_list, rawdir, caldir)

    # process time files
    make_time_files(ssscc_list, timedir, cnvdir, reportdir, cfg.filter_params, cfg.soak_det_params, cfg.inst_params.freq)

    # process bottle file
    make_btl_files(ssscc_list, rawdir, btldir, cnvdir)

    # generate salt .csv files
    proc_salt(ssscc_list, salt_rawdir, salt_cnvdir, btldir, flagfile)

    # generate reftemp .csv files
    proc_reft(ssscc_list, reft_rawdir, reft_parsedir, reft_cnvdir)

    # parse rbr reft
    # proc_rbr(ssscc_list, rbr_rawdir, cnvdir, btldir, rbr_cnvdir, sync_times=True)
    # proc_rbr(ssscc_list, rbr_rawdir, cnvdir, btldir, rbr_cnvdir, export_parsed=True, parsed_dir=rbr_parseddir)

    # generate oxygen .csv files
    proc_oxy(ssscc_list, oxy_rawdir, oxy_cnvdir)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####
    # load fit groups by parameter
    fit_groups = load_fit_groups(ssscc_list, cfg.fit_groups)
    if not skip_calibrate:

        # load in all bottle and time data into DataFrame
        time_data_all = load_time_all(ssscc_list, timedir)
        btl_data_all = load_btl_all(ssscc_list, btldir, reft_cnvdir, salt_cnvdir, oxy_cnvdir)

        # process pressure offset
        # TODO: these functions return an updated dataframe, which we aren't
        #   assigning or reassigning to anything. Instead we trust that the
        #   updates which happen in the other module are visible by this one
        #   too (they  indeed seem to be). Is this a safe assumption?
        calibrate_pressure(btl_data_all, time_data_all, fit_groups.pressure, reportdir)

        # create cast depth log file
        make_depth_log(time_data_all, cast_id_col='cast_id', report_dir=reportdir)

        # calibrate temperature against reference
        calibrate_temp(
                btl_data_all, time_data_all, fit_groups.temperature, fit_coeffs.temperature,
                reft_figdir, reportdir, cast_id='cast_id'
        )

        # calibrate conductivity against reference
        btl_data_all, time_data_all = calibrate_cond(
                btl_data_all, time_data_all, fit_groups.conductivity, fit_coeffs.conductivity,
                salt_figdir, reportdir, flagfile, cast_id='cast_id'
        )

        # calculate params needs for oxy/rinko calibration
        oxy_cast_list = make_cast_id_list(oxy_rawdir, pattern='?????')
        prepare_oxy(btl_data_all, time_data_all, ssscc_list, flagfile, oxy_rawdir, 'oxygen', 'cast_id')

        # calibrate oxygen against reference
        calibrate_oxy(btl_data_all, time_data_all, oxy_figdir, reportdir, caldir, oxy_cast_list, cast_id_col='cast_id')
        calibrate_rinko(
                btl_data_all, time_data_all,
                rinko_figdir, reportdir,
                oxy_cast_list,
                fit_groups.rinko,
                cast_id_col='cast_id',
                full_cast_list=ssscc_list,
        )

    #####
    # Step 3: export data
    #####
    if process_cruise_report is True:
        # export files and figs for the cruise report
        cruise_report_dir = Path(outdir, 'cruise_report')
        cruise_report(reportdir, cruise_report_dir, fit_groups.pressure)
        return

    if skip_export is False and skip_calibrate is False:
        make_odf_report_table(btl_data_all, reportdir)

        # export to Exchange format
        exchange_settings = load_user_config(EXCHANGECONFIG)
        export_exchange(time_data_all, btl_data_all, ssscc_list, exchange_settings, outdir, reportdir, 'cast_id')


if __name__ == "__main__":
    odf_process_all()

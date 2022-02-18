import logging
from pathlib import Path

from ctdcal import flagging, get_ctdcal_config, io

log = logging.getLogger(__name__)

cfg = get_ctdcal_config()


def cnv_to_ct1():
    """
    A script for converting Sea-Bird .cnv files to uncalibrated Exchange CTD files
    for plotting and QA/QC purposes (in ODV, etc.).
    """
    cnv_files = list(Path(cfg.dirs["converted"]).glob("*.cnv"))
    log.info(f"Found {len(cnv_files)} .cnv files, attempting to convert")

    # make CTD files
    sbe_to_woce = {  # map col names from SBE to WOCE
        "prDM": "CTDPRS",
        "t090C": "CTDTMP1",
        "t190C": "CTDTMP2",
        "c0mS/cm": "CTDCOND1",
        "c1mS/cm": "CTDCOND2",
        "sal00": "CTDSAL1",
        "sal11": "CTDSAL2",
        "sbeox0V": "CTDOXYVOLTS",
        "sbeox0ML/L": "CTDOXY",
        "flECO-AFL": "CTDFLUOR",
        "CStarTr0": "CTDXMISS",
    }
    for f in cnv_files:
        df = io.load_cnv(f).rename(mapper=sbe_to_woce, axis=1)
        df = df[sbe_to_woce.values()]

        # give everything WOCE-named uncalibrated flags
        for idx, col in enumerate(df.columns):
            flags = flagging.nan_values(df[col], flag_good=1, flag_nan=9)
            df.insert(idx * 2 + 1, col + "_FLAG_W", flags)

        # export to pressure folder
        # TODO: replace this with a to_exchange() function (should live in ctdcal.io)
        df.to_csv(f"{cfg.dirs['pressure']}{f.stem}_ct1.csv", na_rep="-999")


def cnv_to_hy1():
    """
    A script for converting Sea-Bird .cnv files to an uncalibrated Exchange bottle file
    for plotting and QA/QC purposes (in ODV, etc.).
    """
    # TODO: make bottle file
    # this will probably have significant overlap with cnv_to_ct1()
    # maybe there should be three functions in this file?:
    # _load_cnv_files (load/map names, add flags)
    # cnv_to_ct1 (bin to 2db and export file)
    # cnv_to_hy1 (match to data from .bl and export file)

    # something like this to start:

    # bl_files = list(Path(cfg.dirs["converted"]).glob("*.bl"))
    # if not bl_files:
    #     log.warning("Not .bl files found - cannot make Exchange bottle files")

    # find data to average using .bl file
    pass

"""
Attempt to write a cleaner processing script from scratch.
"""

# import necessary packages
import os
import sys
import subprocess
import pandas as pd
import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import ctdcal.oxy_fitting as oxy_fitting
import ctdcal.rinko as rinko
import ctdcal.odf_io as odf_io
import scripts.odf_sbe_metadata as odf_sbe_metadata


def process_all():

    # TODO: document which functions/scripts produce which files

    #####
    # Step 0: Load and define necessary variables
    #####

    import config as cfg

    #####
    # Step 1: Generate intermediate file formats (.pkl, _salts.csv, _reft.csv)
    #####

    # load station/cast list from file
    ssscc_list = []
    with open(cfg.directory["ssscc_file"], "r") as filename:
        ssscc_list = [line.strip() for line in filename]

    # convert hex to ctd (TODO: convert this to function form)
    cnv_dir_list = os.listdir("data/converted/")
    for ssscc in ssscc_list:
        if "{}.pkl".format(ssscc) not in cnv_dir_list:
            subprocess.run(
                [
                    "odf_convert_sbe.py",
                    "data/raw/" + ssscc + ".hex",
                    "data/raw/" + ssscc + ".XMLCON",
                    "-o",
                    "data/converted",
                ],
                stdout=subprocess.PIPE,
            )
            print("odf_convert_sbe.py SSSCC: " + ssscc + " done")

    # first half of CTD data processing
    # TODO: clean up and document better
    # TODO: export to odf_io instead of standalone script?
    # this produces #_ct1.csv (preliminary), #_time.pkl, and "ondeck_pressure.csv"
    time_dir_list = os.listdir("data/time/")
    for ssscc in ssscc_list:
        if "{}_time.pkl".format(ssscc) not in time_dir_list:
            odf_sbe_metadata.main("data/converted/" + ssscc + ".pkl")
            print("odf_sbe_metadata.py SSSCC: " + ssscc + " done")

    # process bottle file (TODO: convert this to function form)
    btl_dir_list = os.listdir("data/bottle/")
    for ssscc in ssscc_list:
        if "{}_btl_mean.pkl".format(ssscc) not in btl_dir_list:
            subprocess.run(
                [
                    "odf_process_bottle.py",
                    "data/converted/" + ssscc + ".pkl",
                    "-o",
                    "data/bottle/",
                ],
                stdout=subprocess.PIPE,
            )
            print("odf_process_bottle.py SSSCC: " + ssscc + " done")

    # generate salt .csv files
    odf_io.process_salts(ssscc_list)

    # generate reftemp .csv files
    process_ctd.process_reft(ssscc_list)

    #####
    # Step 2: calibrate pressure, temperature, conductivity, and oxygen
    #####

    # load in all bottle and time data into DataFrame
    # TODO: clean up process_ctd.load_all_ctd_files
    btl_data_all = process_ctd.load_all_ctd_files(ssscc_list, "bottle", cfg.btl_cols)
    time_data_all = process_ctd.load_all_ctd_files(ssscc_list, "time", None)

    # process pressure offset
    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)

    # create cast depth log file
    process_ctd.make_depth_log(time_data_all)

    # calibrate temperature against reference
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)

    # calibrate temperature against reference
    fit_ctd.calibrate_cond(btl_data_all, time_data_all)

    # calculate params needs for oxy/rinko calibration
    oxy_fitting.prepare_oxy(btl_data_all, time_data_all, ssscc_list)

    # calibrate oxygen against reference
    oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)

    # TODO: calibrate rinko against reference, similar to oxy_fitting.calibrate_oxy()
    # rinko.calibrate_oxy()  # or something

    #####
    # Step 3: export data
    #####

    # export time data to _ct1.csv format
    # TODO: clean this up more
    process_ctd.export_ct1(time_data_all, ssscc_list)

    # run: ctd_to_bottle.py

    breakpoint()

    # TODO: write/abstract code to a function (in process_ctd?)
    # TODO: look at Argo .nc files for inspiration on structuring
    #
    # experiment with xarray
    # import xarray  # not needed apparently? at least as currently coded

    # da_out = df.to_xarray()  # xarray calls them DataArrays instead of DataFrames

    # # set attributes
    # da_out["CTDTMP"].attrs["long_name"] = "Temperature"
    # da_out["CTDTMP"].attrs["units"] = "ITS-90"
    # da_out["CTDTMP"].attrs["description"] = "Continuous temperature from CTD downcast"

    # # can set attrs from dict
    # prs_attrs = dict(
    #     long_name="Pressure", units="dbar", description="Continuous pressure"
    # )
    # tmp_attrs = dict(
    #     long_name="Temperature", units="ITS-90", description="Continuous temperature"
    # )
    # sal_attrs = dict(
    #     long_name="Salinity", units="PSS-78", description="Continuous salinity",
    # )

    # # can do one at at time
    # da_out["CTDSAL"].attrs = sal_attrs

    # # maybe do nested dicts in config.py and loop? e.g.:
    # ctd_attrs = dict(
    #     CTDPRS=prs_attrs,
    #     CTDTMP=tmp_attrs,
    #     CTDSAL=sal_attrs,
    #     # CTDOXY=oxy_attrs,
    #     # CTDRINKO=rinko_attrs,
    #     # CTDXMISS=xmiss_attrs,
    #     # CTDFLUOR=fluor_attrs,
    # )

    # for var in da_out.keys():
    #     if var == "SSSCC":
    #         continue
    #     if not var.endswith("_FLAG_W"):
    #         da_out[var].attrs = ctd_attrs[var]

    # # output files
    # # don't actually run bc this is 24Hz data...
    # # da_out.to_netcdf('example.nc')


def main(argv):
    """Run everything.
    """
    process_all()


if __name__ == "__main__":
    main(sys.argv[1:])

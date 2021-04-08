### configuration file for odf_process_all_MK.py
# this is an example to fill with your own values
#
# TODO: organize these by editable/fixed variables

# Cruise specifics
# not really needed for processing, move elsewhere
cruise = dict(
    cruise_title="GO-SHIP TN389",
    cruise_name="TN389",
    cruisedb="TN389",
    vessel="R/V Thomas G Thompson",
    ship_code="3250",
    expocode="325020210316",
    chief_sci="Ryan Woosley",
    sectionid="A20",
    start_date="2021-03- 14:00:00",
    start_port="Woods Hole, Massachussetts, USA",
    start_latlon="42.8821 S 147.3272 W",
    end_date="2018-05-14 10:00:00",
    end_port="St. Thomas, USVI",
    end_latlon="53.1638 S 70.9171 W",
)

ctd_serial = 914

# CTD variables/flags/units
# move elsewhere when xarray is implemented
ctd_outputs = dict(
    press=["CTDPRS", "CTDPRS_FLAG_W", "DBAR", ""],
    temp=["CTDTMP", "CTDTMP_FLAG_W", "ITS-90", ""],
    salt=["CTDSAL", "CTDSAL_FLAG_W", "PSS-78", ""],
    doxy=["CTDOXY", "CTDOXY_FLAG_W", "UMOL/KG", ""],
    rinko=["CTDRINKO", "CTDRINKO_FLAG_W", "UMOL/KG", ""],
    # rinko=["CTDOXY", "CTDOXY_FLAG_W", "UMOL/KG", ""],  # reporting Rinko as primary oxy
    xmiss=["CTDXMISS", "CTDXMISS_FLAG_W", "0-5VDC", ""],
    fluor=["CTDFLUOR", "CTDFLUOR_FLAG_W", "0-5VDC", ""],
    # backscatter=["CTDBACKSCATTER", "CTDBACKSCATTER_FLAG_W", "0-5VDC", ""],
    # bbp = ['CTDBBP700RAW', 'CTDBBP700RAW_FLAG_W', '0-5VDC', ''],
)

ctd_col_names = []
ctd_col_units = []
for i in range(len(ctd_outputs)):
    param = list(ctd_outputs.keys())[i]
    param_list = ctd_outputs[param]
    ctd_col_names.append(param_list[0])
    ctd_col_names.append(param_list[1])
    ctd_col_units.append(param_list[2])
    ctd_col_units.append(param_list[3])

ctd_time_output = dict(col_names=ctd_col_names, col_units=ctd_col_units)

# T/C fitting parameters
# temperature params are (P_order, T_order, zRange)
# conductivity params are (P_order, T_order, C_order, zRange)
fit_orders1 = {
    "ssscc_t1": (1, 0, "750:6000"),
    "ssscc_t2": (1, 0, "1000:6000"),
    "ssscc_t3": (1, 0, "1500:6000"),
    "ssscc_c1": (1, 0, 0, "1000:6000"),
    "ssscc_c2": (2, 0, 1, "900:6000"),
    "ssscc_c3": (2, 0, 1, "900:6000"),
}
fit_orders2 = {
    "ssscc_t1": (1, 0, "750:6000"),
    "ssscc_t2": (1, 0, "1000:6000"),
    "ssscc_t3": (1, 0, "1500:6000"),
    "ssscc_c1": (1, 0, 0, "1000:6000"),
    "ssscc_c2": (1, 1, 0, "1000:6000"),
    "ssscc_c3": (2, 1, 1, "1200:6000"),
}

# List of directories for I/O purposes
directory = {
    "ssscc_file": "data/ssscc.csv",
    "ssscc": "data/ssscc/",
    "raw": "data/raw/",
    "converted": "data/converted/",
    "time": "data/time/",
    "pressure": "data/pressure/",
    "bottle": "data/bottle/",
    "salt": "data/salt/",
    "reft": "data/reft/",
    "oxy": "data/oxygen/",
    "logs": "data/logs/",
    "qual_temp_primary": "quality_code/temp_primary/",
    "qual_temp_secondary": "quality_code/temp_secondary/",
    "qual_cond_primary": "quality_code/cond_primary/",
    "qual_cond_secondary": "quality_code/cond_secondary/",
    "t1_fit_figs": "data/logs/fitting_figs/temp_primary/",
    "t2_fit_figs": "data/logs/fitting_figs/temp_secondary/",
    "c1_fit_figs": "data/logs/fitting_figs/cond_primary/",
    "c2_fit_figs": "data/logs/fitting_figs/cond_secondary/",
    "ox_fit_figs": "data/logs/fitting_figs/oxy_primary/",
    "rinko_fit_figs": "data/logs/fitting_figs/oxy_rinko/",
}

# remnant of old system, will be pushed into xarray metadata/attrs
# Labels for CTD columns
column = {
    "p": "CTDPRS",
    "p_btl": "CTDPRS",
    "t1": "CTDTMP1",
    "t2": "CTDTMP2",
    "t1_btl": "CTDTMP1",
    "t2_btl": "CTDTMP2",
    "c1": "CTDCOND1",
    "c2": "CTDCOND2",
    "c1_btl": "CTDCOND1",
    "c2_btl": "CTDCOND2",
    "reft": "T90",
    "refc": "BTLCOND",
    "sal": "CTDSAL",
    "sal_btl": "SALNTY",
    "rinko_oxy": "U_DEF_poly1",  # CHECK THIS!
    "oxy_btl": "OXYGEN",
    "oxyvolts": "CTDOXYVOLTS",
    "lat": "GPSLAT",
    "lon": "GPSLON",
    "lat_btl": "GPSLAT",
    "lon_btl": "GPSLON",
}

# List of columns to filter
filter_cols = []
for x in ["p", "t1", "t2", "c1", "c2", "sal", "rinko_oxy", "oxyvolts", "lat", "lon"]:
    filter_cols.append(column[x])

# List of bottle columns to be read during calibration
btl_cols = [
    "CTDTMP1",
    "CTDTMP2",
    "CTDPRS",
    "CTDCOND1",
    "CTDCOND2",
    "CTDSAL",
    "CTDOXY1",
    "CTDOXYVOLTS",
    "U_DEF_poly1",
    "CTDXMISS",
    "CTDFLUOR",
    "ALT",
    "REF_PAR",
    "GPSLAT",
    "GPSLON",
    "new_fix",
    "pressure_temp_int",
    "pump_on",
    "btl_fire",
    "scan_datetime",
    "btl_fire_num",
]

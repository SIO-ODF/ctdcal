"""
P02 Process the BIO casts, which run before the full cast, on
full cast coeffs.

Bio casts:
* Typically run about 1000 m
    * The altimeter is almost never used
* Have no Winkler/Autosal procedures done
* Tend to have most bottles fired at the surface in quick bursts
    * Therefore the refT can be sporadic in later bottles (<15 sec)

Order of operations:
* ctdcal process -> processes the full cast and generates fit coeffs
* ctdcal bio     -> processes bio casts using the previous ctdcal run
"""

# from sqlalchemy import true
from time import time
from ctdcal import (
    convert,
    fit_ctd,
    get_ctdcal_config,
    odf_io,
    oxy_fitting,
    process_bottle,
    process_ctd,
    rinko,
    sbe_reader as sbe_rd,
    flagging,
)

from pathlib import Path
import logging
import pandas as pd
import numpy as np
import gsw
from datetime import datetime, timezone

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)

def odf_process_bio():
    """
    Rebuilding based on procedure outlined in odf_process_all.py
    """
    #   Step 0: Initialize

    #   Step 1: Import the fit coeffs for T, C, oxy, and Rinko
    #   Note: The refT for the bio casts exist and can be refit.
    # fname_t1    = "data/logs/fit_coef_t1.csv"
    # t1_coefs    = pd.read_csv(fname_t1)
    # fname_t2    = "data/logs/fit_coef_t2.csv"
    # t2_coefs    = pd.read_csv(fname_t2)
    fname_c1    = "data/logs/fit_coef_c1.csv"
    c1_coefs    = pd.read_csv(fname_c1)
    fname_c2    = "data/logs/fit_coef_c2.csv"
    c2_coefs    = pd.read_csv(fname_c2)
    fname_sbe4  = "data/logs/sbe43_coefs.csv"
    sbe43_coefs = pd.read_csv(fname_sbe4)
    fname_rinko = "data/logs/rinko_coefs.csv"
    rinko_coefs = pd.read_csv(fname_rinko)
    log.info("All coefs loaded.")

    #   Step 2: Make the bio_ssscc
    #   Should keep this separate from the full ssscc_list
    #   to prevent anything from trying to be fit without having
    #   data that can be used to constrain the fit
    fname_ssscc = "data/ssscc_bio.csv"
    ssscc_list = []
    with open(fname_ssscc, "r") as filename:
        ssscc_list = [line.strip() for line in filename]

    #   Step 3: Convert the hex into bottle averages and time .pkl
    log.info("Converting .hex files")
    for ssscc in ssscc_list:
        #   Stash .pkl files in the normal folders
        if not Path(cfg.dirs["converted"] + ssscc + ".pkl").exists():
            hexFile = "data/raw_bio/" + ssscc + ".hex"  #   Keep the raw files separated for safety
            xmlconFile = "data/raw_bio/" + ssscc + ".XMLCON"
            sbeReader = sbe_rd.SBEReader.from_paths(hexFile, xmlconFile)
            converted_df = convert.convertFromSBEReader(sbeReader, ssscc)
            converted_df.to_pickle(cfg.dirs["converted"] + ssscc + ".pkl")

    convert.make_time_files(ssscc_list)
    convert.make_btl_mean(ssscc_list)
    log.info('All .pkl files generated.')

    #   Might as well process the refT files for the bio casts, saves trouble of calibrate temp
    process_bottle.process_reft(ssscc_list)

    #   Step 4: Get the time and bottle data into the workspace
    time_data_all   = process_ctd.load_all_ctd_files(ssscc_list)
    btl_data_all    = process_bottle.load_all_btl_files(ssscc_list) #   Extra cols should be full of NaNs
    if '05301' in ssscc_list:
        bottles_05301 = [float(bnum) for bnum in [1, 2, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21,
        22, 23, 25, 27, 28, 29, 31, 33, 35]]    #   Custom bottle positions
        btl_data_all["btl_fire_num"].loc[btl_data_all["SSSCC"] == "05301"] = bottles_05301

    log.info('Time and Bottle data loaded.')
    # print(time_data_all.head())
    # print(btl_data_all.head())

    #   Add SSSCC2 column, incrementing up from bio num to full num
    btl_data_all["SSSCC2"] = btl_data_all["SSSCC"].str.replace(r"01", "02") #   Sorry Mike
    time_data_all["SSSCC2"] = time_data_all["SSSCC"].str.replace(r"01", "02")
    #   Station 032 is the only one where the bio SSSCC ends in a 2
    btl_data_all["SSSCC2"].loc[btl_data_all["SSSCC"]=='03202'] = '03203'
    time_data_all["SSSCC2"].loc[time_data_all["SSSCC"]=='03202'] = '03203'

    #   Step 5: Apply the pressure offset
    process_ctd.apply_pressure_offset(btl_data_all)
    process_ctd.apply_pressure_offset(time_data_all)
    print("Pressure offset successfully applied.")

    #   Step 6: Make the depth log
    #   Note: Bio is really shallow, so all the depths will be -999
    #   To fix: Look for int(SSSCC)+1 in the depth log list and use that
    # process_ctd.make_depth_log(time_data_all) #   (Overwrites existing log)
    depth_log_bio(time_data_all)

    #   Step 7: Calibrate the temperature
    #   Note: The refT for the bio casts exist and can be refit.
    # calibrate_temp_bio(btl_data_all, time_data_all, t1_coefs, t2_coefs)
    fit_ctd.calibrate_temp(btl_data_all, time_data_all)
    print("Temperature calibration successful...?")

    #   Step 8: Calibrate the cond w/ pre-existing coeffs
    btl_data_all, time_data_all = calibrate_cond_bio(
        btl_data_all, time_data_all, c1_coefs, c2_coefs)
    print("Conductivity has been 'calibrated'.")

    #   Step 9: Calibrate the SBE43 w/ pre-existing coeffs
    btl_data_all["OXYGEN_FLAG_W"] = 9
    time_data_all["OXYGEN_FLAG_W"] = 9
    prepare_oxy_bio(btl_data_all, time_data_all, ssscc_list)

    ### Coming back to this.
    #   oxy_fitting.calibrate_oxy(btl_data_all, time_data_all, ssscc_list)
    #   No plotting with NaNs
    # calibrate_sbe43_bio(btl_data_all, time_data_all,
    #     sbe43_coefs,
    # )
    print("Oxygen columns are up.")

    #   Step 10: Calibrate the Rinko w/ pre-existing coeffs
    prepare_rinko_bio(btl_data_all, time_data_all, rinko_coefs)
    print('RINKO columns are up.')


    #   Step 11: Export the ct1 files
    export_ct1_bio(time_data_all)

    #   Step 12: Export/merge the hy1 file
    #   Barna would prefer the hy1 files to be one big one that he doesn't need to merge himself
    export_hy1_bio(btl_data_all)

def calibrate_cond_bio(btl_df, time_df, c1_coefs, c2_coefs):
    # import gsw
    log.info("Calibrating conductivity without AutoSal for bio casts...")
    btl_df["SALNTY_FLAG_W"] = 9 
    # ssscc_subsets = sorted(Path(cfg.dirs["ssscc"]).glob("ssscc_c*.csv"))
    # for cN, tN in zip(["c1", "c2"], ["t1", "t2"]):
    #     for f in ssscc_subsets:
    #         ssscc_sublist = pd.read_csv(f, header=None, dtype="str").squeeze(axis=1).to_list()
    #         #   Since no salts are processed, the ssscc should be pulled from the full cast
    #         btl_rows = btl_df["SSSCC2"].isin(ssscc_sublist).values
    #         time_rows = time_df["SSSCC2"].isin(ssscc_sublist).values
            
    #         #   Prepare data for fitting: I think I need the df_good, bad
    #         #   Can't do that - no refC therefore no df_good
    #         if cN == 'c1':
    #             P_coefs = c1_coefs.cp1
    #             T_coefs = c1_coefs.ct1
    #             C_coefs = c1_coefs.c0
    #         elif cN == 'c2':
    #             P_coefs = c2_coefs.cp1
    #             T_coefs = c2_coefs.ct1
    #             C_coefs = c2_coefs.c0
    #         btl_df.loc[btl_rows, cfg.column[cN]] = fit_ctd.apply_polyfit(
    #             btl_df.loc[btl_rows, cfg.column[cN]],
    #             C_coefs,
    #             (btl_df.loc[btl_rows, cfg.column["p"]], P_coefs),
    #             (btl_df.loc[btl_rows, cfg.column[tN]], T_coefs),
    #         )
    #         time_df.loc[time_rows, cfg.column[cN]] = fit_ctd.apply_polyfit(
    #             time_df.loc[time_rows, cfg.column[cN]],
    #             C_coefs,
    #             (time_df.loc[time_rows, cfg.column["p"]], P_coefs),
    #             (time_df.loc[time_rows, cfg.column[tN]], T_coefs),
    #         )
    time_df[cfg.column["sal"]] = gsw.SP_from_C(
        time_df[cfg.column["c1"]],
        time_df[cfg.column["t1"]],
        time_df[cfg.column["p"]],
    )
    btl_df[cfg.column["sal"]] = gsw.SP_from_C(
        btl_df[cfg.column["c1"]],
        btl_df[cfg.column["t1"]],
        btl_df[cfg.column["p"]],
    )
    time_df["SALNTY"] = np.nan
    btl_df["SALNTY"] = np.nan
    time_df[cfg.column["sal"] + "_FLAG_W"] = 2 #   Questionable, since fitting with other coeffs
    btl_df[cfg.column["sal"] + "_FLAG_W"] = 2

    return btl_df, time_df

def prepare_oxy_bio(btl_df, time_df, ssscc_list):
    btl_df["SA"] = gsw.SA_from_SP(
        btl_df[cfg.column["sal"]],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["lon"]],
        btl_df[cfg.column["lat"]],
    )
    btl_df["CT"] = gsw.CT_from_t(
        btl_df["SA"],
        btl_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p"]],
    )
    time_df["SA"] = gsw.SA_from_SP(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )
    time_df["CT"] = gsw.CT_from_t(
        time_df["SA"],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
    )

    # calculate sigma
    btl_df["sigma_btl"] = gsw.sigma0(btl_df["SA"], btl_df["CT"])
    time_df["sigma_btl"] = gsw.sigma0(time_df["SA"], time_df["CT"])

    # Calculate oxygen solubility in µmol/kg
    btl_df["OS"] = gsw.O2sol(
        btl_df["SA"],
        btl_df["CT"],
        btl_df[cfg.column["p"]],
        btl_df[cfg.column["lon"]],
        btl_df[cfg.column["lat"]],
    )
    time_df["OS"] = gsw.O2sol(
        time_df["SA"],
        time_df["CT"],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )
    time_df["dv_dt"] = oxy_fitting.calculate_dV_dt(time_df["CTDOXYVOLTS"], time_df["scan_datetime"])
    # Convert CTDOXY units
    btl_df["CTDOXY"] = oxy_fitting.oxy_ml_to_umolkg(btl_df["CTDOXY1"], btl_df["sigma_btl"])
    btl_df[cfg.column["refO"]] = np.nan
    time_df["CTDOXY"] = oxy_fitting.oxy_ml_to_umolkg(time_df["CTDOXY1"], time_df["sigma_btl"])
    btl_df["CTDOXY_FLAG_W"] = 2
    time_df["CTDOXY_FLAG_W"] = 2

    #   Not sure where NAN rows are coming from on the tail end of 053.
    if any(btl_df["SSSCC"].isnull()):
        nan_rows = btl_df[btl_df["SSSCC"].isnull()]
        print("There are rows of NaNs around station 053?")

    return True

def calibrate_sbe43_bio(btl_df, time_df, sbe43_coefs):
    log.info("Calibrating oxygen ")
    all_sbe43_merged = pd.DataFrame()
    sbe43_dict = {}
    all_sbe43_fit = pd.DataFrame()
    ssscc_list = btl_df["SSSCC"].dropna().unique()
    ssscc_list2 = btl_df["SSSCC2"].dropna().unique() #   For grabbing other coefs

    #   No reference oxy, I think most of this will be skipped
    #   No density matching
    # for ssscc in ssscc_list:
    #     time_data = time_df[time_df["SSSCC"] == ssscc].copy()
    #     btl_data = btl_df[btl_df["SSSCC"] == ssscc].copy()
    #     # can't calibrate without bottle oxygen ("OXYGEN")
    #     if (btl_data["OXYGEN_FLAG_W"] == 9).all():
    #         sbe43_dict[ssscc] = np.full(5, np.nan)
    #         log.warning(ssscc + " skipped, all oxy data is NaN")
    #         continue

    # sbe_coef0 = sbe43_coefs

    ### Come back to this ###
    #   We need to get sbe_df, passing in just btl_df isn't enough

    #   I can't quite follow what's going on with sbe43_merged coming from sigma matching.

    # for ssscc in ssscc_list:
    #     sbe_coef, sbe_df = oxy_fitting.sbe43_oxy_fit(
    #         all_sbe43_merged.loc[all_sbe43_merged["SSSCC"] == ssscc].copy(),
    #         sbe_coef0=sbe_coef0,
    #         f_suffix=f"_{ssscc}",
    #     )
    #     # build coef dictionary
    #     if ssscc not in sbe43_dict.keys():  # don't overwrite NaN'd stations
    #         sbe43_dict[ssscc] = sbe_coef
    #     # all non-NaN oxygen data with flags
    #     all_sbe43_fit = pd.concat([all_sbe43_fit, sbe_df])

def prepare_rinko_bio(btl_df, time_df, rinko_coefs):
    #   This may be easier that the SBE43, since Uchida_DO_eq allows for us to feed coefs in
    #   Get SSSCC coefs in, get the voltage for each SSSCC, pass to Uchida
    ssscc_list2 = btl_df["SSSCC2"].dropna().unique()
    varnames = ["c0", "c1", "c2", "d0", "d1", "d2", "cp"]
    btl_df["CTDRINKO"] = np.nan
    time_df["CTDRINKO"] = np.nan
    for ssscc in ssscc_list2:
        ssscc_coefs = rinko_coefs[rinko_coefs["Unnamed: 0"].str.contains(ssscc)]
        ssscc_coefs = ssscc_coefs.iloc[0]
        ssscc_coefs = tuple(ssscc_coefs[1:])

        btl_df["CTDRINKO"].loc[btl_df["SSSCC2"] == ssscc] = rinko._Uchida_DO_eq(
            ssscc_coefs, 
            (
                btl_df["U_DEF_poly1"].loc[btl_df["SSSCC2"] == ssscc],
                btl_df["CTDPRS"].loc[btl_df["SSSCC2"] == ssscc],
                btl_df["CTDTMP1"].loc[btl_df["SSSCC2"] == ssscc],
                btl_df["CTDSAL"].loc[btl_df["SSSCC2"] == ssscc],
                btl_df["OS"].loc[btl_df["SSSCC2"] == ssscc],
            ))
        time_df["CTDRINKO"].loc[time_df["SSSCC2"] == ssscc] = rinko._Uchida_DO_eq(
            ssscc_coefs, 
            (
                time_df["U_DEF_poly1"].loc[time_df["SSSCC2"] == ssscc],
                time_df["CTDPRS"].loc[time_df["SSSCC2"] == ssscc],
                time_df["CTDTMP1"].loc[time_df["SSSCC2"] == ssscc],
                time_df["CTDSAL"].loc[time_df["SSSCC2"] == ssscc],
                time_df["OS"].loc[time_df["SSSCC2"] == ssscc],
            ))
    btl_df["CTDRINKO_FLAG_W"] = 2
    time_df["CTDRINKO_FLAG_W"] = 2

    return True

def depth_log_bio(time_data_all):
    df = time_data_all[["SSSCC", "CTDPRS", "GPSLAT", "ALT"]].copy().reset_index()
    df_group = df.groupby("SSSCC", sort=False)
    idx_p_max = df_group["CTDPRS"].idxmax()
    bottom_df = pd.DataFrame(
        data={
            "SSSCC": df["SSSCC"].unique(),
            "max_p": df.loc[idx_p_max, "CTDPRS"],
            "lat": df.loc[idx_p_max, "GPSLAT"],
            "alt": df.loc[idx_p_max, "ALT"],
        }
    )
    #   While bio casts aren't too deep, they *could* be shallow and see bottom in SD
    #   Keep this in.
    bottom_df.loc[bottom_df["alt"] > 80, "alt"] = np.nan
    #   Read in the depth log to get actual depths. Fix as per normal later on in hy1 writeout.
    auto_depth_df   = pd.read_csv(cfg.dirs["logs"]+"depth_log.csv")
    #   Bio casts always preceed a full cast. Look for ssscc + 1 in the auto depth log for the depth to use
    bottom_df["SSSCC2"] = bottom_df["SSSCC"].astype(int)+1
    pull = auto_depth_df["DEPTH"].loc[auto_depth_df["SSSCC"].isin(bottom_df["SSSCC2"])]
    pull.index = bottom_df.index    #   Reindex to make replacement work
    bottom_df["DEPTH"] = pull
    #   The manual depth log is pulled when the ct1 and hy1 files are made.
    #   Will need to tweak those functions to read from the bio_depth_log
    #   I don't want to overwrite the existing log.
    if any(bottom_df["DEPTH"] == -999):
        # manual_depth_df = pd.read_csv(
        #     cfg.dirs["logs"] + "manual_depth_log.csv", dtype={"SSSCC": str}
        # )
        # pull = manual_depth_df.loc[manual_depth_df.SSSCC.astype(float).isin(bottom_df.SSSCC2)]
        #   Why are pandas indices so awful
        bottom_df["DEPTH"].iloc[2] = 6221.0
        bottom_df["DEPTH"].iloc[3] = 6185.0

    bottom_df[["SSSCC", "DEPTH"]].to_csv(
        cfg.dirs["logs"] + "depth_log_bio.csv", index=False
    )

def export_ct1_bio(time_df):
    time_df["CTDFLUOR_FLAG_W"] = 1
    time_df["CTDXMISS_FLAG_W"] = 1
    # rename outputs as defined in user_settings.yaml
    for param, attrs in cfg.ctd_outputs.items():
        if param not in time_df.columns:
            print(param, 'being renamed from', attrs["sensor"])
            time_df.rename(columns={attrs["sensor"]: param}, inplace=True)

    # check that all columns are there
    time_df[cfg.ctd_col_names]

    time_df["SSSCC"] = time_df["SSSCC"].astype(str).copy()
    cast_details = pd.read_csv(
        # cfg.dirs["logs"] + "cast_details.csv", dtype={"SSSCC": str}
        cfg.dirs["logs"] + "bottom_bottle_details.csv",
        dtype={"SSSCC": str},
    )
    depth_df = pd.read_csv(
        cfg.dirs["logs"] + "depth_log_bio.csv", dtype={"SSSCC": str}, na_values=-999
    ).dropna()
    #   depth_df actually looks good
    # manual_depth_df = pd.read_csv(
    #     cfg.dirs["logs"] + "manual_depth_log.csv", dtype={"SSSCC": str}
    # )
    # full_depth_df = pd.concat([depth_df, manual_depth_df])
    # full_depth_df.drop_duplicates(subset="SSSCC", keep="first", inplace=True)

    #   Iteratively writing the ct1 files out
    for ssscc in time_df["SSSCC"].unique():   
        time_data = time_df[time_df["SSSCC"] == ssscc].copy()
        time_data = process_ctd.pressure_sequence(time_data)
        
        bad_list = ["01001", "01101", "01201"]  #   P02
        if ssscc in bad_list:
            print(f"Using SBE43 as CTDOXY for {ssscc}")
        else:   # switch oxygen primary sensor to rinko
            print(f"Using Rinko as CTDOXY for {ssscc}")
            time_data.loc[:, "CTDOXY"] = time_data["CTDRINKO"]
            time_data.loc[:, "CTDOXY_FLAG_W"] = time_data["CTDRINKO_FLAG_W"]
        
        time_data = time_data[cfg.ctd_col_names]
        # time_data = time_data.round(4)
        time_data = time_data.where(~time_data.isnull(), -999)  # replace NaNs with -999

        # force flags back to int (TODO: make flags categorical)
        for col in time_data.columns:
            if col.endswith("FLAG_W"):
                time_data[col] = time_data[col].astype(int)

        try:
            depth = depth_df.loc[depth_df["SSSCC"] == ssscc, "DEPTH"].iloc[0]
        except IndexError:
            log.warning(f"No depth logged for {ssscc}, setting to -999")
            depth = -999

        # get cast_details for current SSSCC
        cast_dict = cast_details[cast_details["SSSCC"] == ssscc].to_dict("records")[0]
        b_datetime = (
            datetime.fromtimestamp(cast_dict["bottom_time"], tz=timezone.utc)
            .strftime("%Y%m%d %H%M")
            .split(" ")
        )
        
        btm_lat = cast_dict["latitude"]
        btm_lon = cast_dict["longitude"]

        now = datetime.now(timezone.utc)
        file_datetime = now.strftime("%Y%m%d")  # %H:%M")
        file_datetime = file_datetime + "ODFSIO"

        #   Write out to a place separate from the full casts so I can see them
        with open(f"{cfg.dirs['bio']}{ssscc}_ct1.csv", "w+") as f:
            ctd_header = (  # this is ugly but prevents tabs before label
                f"CTD,{file_datetime}\n"
                f"NUMBER_HEADERS = 11\n"
                f"EXPOCODE = {cfg.expocode}\n"
                f"SECT_ID = {cfg.section_id}\n"
                f"STNNBR = {ssscc[:3]}\n"  # STNNBR = SSS
                f"CASTNO = {ssscc[3:]}\n"  # CASTNO = CC
                f"DATE = {b_datetime[0]}\n"
                f"TIME = {b_datetime[1]}\n"
                f"LATITUDE = {btm_lat:.4f}\n"
                f"LONGITUDE = {btm_lon:.4f}\n"
                f"INSTRUMENT_ID = {cfg.ctd_serial}\n"
                f"DEPTH = {depth:.0f}\n"
            )
            f.write(ctd_header)
            np.asarray(cfg.ctd_col_names).tofile(f, sep=",", format="%s")
            f.write("\n")
            np.asarray(cfg.ctd_col_units).tofile(f, sep=",", format="%s")
            f.write("\n")
            time_data.to_csv(f, header=False, index=False)
            f.write("END_DATA")

def export_hy1_bio(df, out_dir=cfg.dirs["bio"], org="ODF"):
    log.info("Exporting bottle file")
    btl_data = df.copy()
    now = datetime.now()
    file_datetime = now.strftime("%Y%m%d")

    # TODO: move to config; integrate Barna's "params" package instead?
    btl_columns = {
        "EXPOCODE": "",
        "SECT_ID": "",
        "STNNBR": "",
        "CASTNO": "",
        "SAMPNO": "",
        "BTLNBR": "",
        "BTLNBR_FLAG_W": "",
        "DATE": "",
        "TIME": "",
        "LATITUDE": "",
        "LONGITUDE": "",
        "DEPTH": "METERS",
        "CTDPRS": "DBAR",
        "CTDTMP": "ITS-90",
        "CTDSAL": "PSS-78",
        "CTDSAL_FLAG_W": "",
        "SALNTY": "PSS-78",
        "SALNTY_FLAG_W": "",
        "CTDOXY": "UMOL/KG",
        "CTDOXY_FLAG_W": "",
        # "OXYGEN": "UMOL/KG",  # P02 Barna won't use Oxy
        # "OXYGEN_FLAG_W": "",
        "REFTMP": "ITS-90",
        "REFTMP_FLAG_W": "",
        "CTDFLUOR": "0-5VDC",   # P02 Barna wants FLUOR + XMISS
        "CTDFLUOR_FLAG_W": "",
        "CTDXMISS": "0-5VDC",
        "CTDXMISS_FLAG_W": "",
    }

    # rename outputs as defined in user_settings.yaml
    for param, attrs in cfg.ctd_outputs.items():
        if param not in btl_data.columns:
            print(param, 'renamed to', attrs["sensor"])
            btl_data.rename(columns={attrs["sensor"]: param}, inplace=True)

    if any(btl_data["SSSCC"].isnull()):
        print("Still need to look at the NaNs in the SSSCC..")
        #   Drop the rows where the SSSCC is NaN (which shouldn't be happening)
        btl_data = btl_data.dropna(subset=['SSSCC'])

    btl_data["EXPOCODE"] = cfg.expocode
    btl_data["SECT_ID"] = cfg.section_id
    btl_data["STNNBR"] = [int(x[0:3]) for x in btl_data["SSSCC"]]
    btl_data["CASTNO"] = [int(x[3:]) for x in btl_data["SSSCC"]]
    btl_data["SAMPNO"] = btl_data["btl_fire_num"].astype(int)
    btl_data = process_bottle.add_btlnbr_cols(btl_data, btl_num_col="btl_fire_num")

    # sort by decreasing sample number (increasing pressure) and reindex
    btl_data = btl_data.sort_values(
        by=["STNNBR", "SAMPNO"], ascending=[True, False], ignore_index=True
    )

    #   P02 add in FLUOR, XMISS flags per Barna's instructions
    btl_data["CTDFLUOR_FLAG_W"] = 1
    btl_data["CTDXMISS_FLAG_W"] = 1
    # round data
    # for col in ["CTDTMP", "CTDSAL", "SALNTY", "REFTMP"]:
    #     btl_data[col] = btl_data[col].round(4)
    # for col in ["CTDPRS", "CTDOXY", "OXYGEN"]:
    #     btl_data[col] = btl_data[col].round(1)

    # add depth
    depth_df = pd.read_csv(
        cfg.dirs["logs"] + "depth_log_bio.csv", dtype={"SSSCC": str}, na_values=-999
    ).dropna()
    btl_data["DEPTH"] = -999
    for index, row in depth_df.iterrows():
        btl_data.loc[btl_data["SSSCC"] == row["SSSCC"], "DEPTH"] = int(row["DEPTH"])

    # deal with nans
    # TODO: missing REFTMP not obvious til loading data - where to put this?
    # _reft_loader() is not the right place
    # maybe during loading step flag missing OXYGEN, REFTMP, BTLCOND?
    btl_data["REFTMP_FLAG_W"] = flagging.nan_values(
        btl_data["REFTMP_FLAG_W"], old_flags=btl_data["REFTMP_FLAG_W"]
    )
    btl_data = btl_data.where(~btl_data.isnull(), -999)

    # check columns
    try:
        btl_data[btl_columns.keys()]
    except KeyError as err:
        log.info("Column names not configured properly... attempting to correct")
        bad_cols = err.args[0].split("'")[1::2]  # every other str is a column name
        for col in bad_cols:
            if col.endswith("FLAG_W"):
                log.warning(col + " missing, flagging with 9s")
                btl_data[col] = 9
            else:
                log.warning(col + " missing, filling with -999s")
                btl_data[col] = -999

    btl_data = btl_data[btl_columns.keys()]
    time_stamp = file_datetime + org
    with open(out_dir + cfg.expocode + "_bio_hy1.csv", mode="w+") as f:
        f.write("BOTTLE, %s\n" % (time_stamp))
        f.write(",".join(btl_columns.keys()) + "\n")
        f.write(",".join(btl_columns.values()) + "\n")
        btl_data.to_csv(f, header=False, index=False)
        f.write("\n" + "END_DATA")

    return


if __name__ == "__main__":
    odf_process_bio()
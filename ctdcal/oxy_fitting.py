#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:04:19 2018

@author: k3jackson
"""

import csv
from pathlib import Path

import config as cfg
import gsw
import numpy as np
import pandas as pd
import scipy

import ctdcal.ctd_plots as ctd_plots
import ctdcal.process_ctd as process_ctd
import ctdcal.sbe_equations_dict as sbe_eq
import ctdcal.sbe_reader as sbe_rd


def oxy_loader(oxyfile):

    f = open(oxyfile,newline='')
    oxyF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
    oxy_Array = []

    for row in oxyF:
        if len(row) > 9:
            row = row[:9]
        oxy_Array.append(row)

    params = oxy_Array[0]

    del oxy_Array[0]

    header = ['STNNO_OXY','CASTNO_OXY','BOTTLENO_OXY','FLASKNO','TITR_VOL','TITR_TEMP','DRAW_TEMP','TITR_TIME','END_VOLTS']

    O2 = np.array(oxy_Array)
    df = pd.DataFrame(O2,columns=header)

    f.close()

    df = df.apply(pd.to_numeric, errors='ignore')

    # Remove "Dummy Data"
    #df['BOTTLENO_OXY'] = df['BOTTLENO_OXY'].astype(int,errors='ignore')
    df = df[df['BOTTLENO_OXY']!=99]

    # Remove "ABORTED DATA"

    df = df[df['TITR_VOL'] > 0]

    # Sort and reindex
    df = df.sort_values('BOTTLENO_OXY')
    #df = df.reset_index()

    # Get necessary columns for output
    df['STNNO_OXY'] = df['STNNO_OXY'].astype(str)
    df['CASTNO_OXY'] = df['CASTNO_OXY'].astype(str)
    df['FLASKNO']= df['FLASKNO'].astype(str)

    df['SSSCC_OXY'] = df['STNNO_OXY']+'0'+df['CASTNO_OXY']

    return df ,params #DF

def flask_load(flaskfile=cfg.directory["oxy"] + 'o2flasks.vol', skip_rows=12):
    """Load information from flask.vol file

    """
## Convert to Dataframe
    with open(flaskfile, 'r') as f:
        flasks = {}
        for l in f:
            if 'Volume' not in l:
                if l.strip().startswith("#"):
                    continue
                row = l.strip().split()
                flasks[str(row[0])] = float(row[1])

    flasks = pd.DataFrame.from_dict(flasks,orient='index')
    flasks = flasks.rename(columns={0:'FLASK_VOL'})
    flasks['FLASKNO'] = flasks.index.values
    flasks = flasks[['FLASKNO','FLASK_VOL']]
    return flasks

#    f = open(flaskfile,newline='')
#    flaskF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
#
#    flask_array = []
#    for row in flaskF:
#        flask_array.append(row)
#
#    f.close
#    del flask_array[0:skip_rows]
#
#    df = pd.DataFrame(flask_array)
#    df['FLASKNO'] = df[0]
#    df['FLASK_VOL'] = df[1]
#    df = df[['FLASKNO','FLASK_VOL']]
#    df['FLASKNO'] = df['FLASKNO'].astype(float,errors='ignore')

#    return df


# code formerly in fit_ctd, need to sort out this vs match_flask
def get_flask_vol(flask, o2flasks, t=20, glass="borosilicate"):
    _coef = {
            "borosilicate": 0.00001,
            "soft": 0.000025,
            }
    _t = 20
    flasks = flask_load()
    fv = flasks[flask]
    coef = _coef[glass]
    return fv * (1.0 + coef * (t - _t));


def match_flask(flask_nums, flask_list, flask_vols, t=20, glass="borosilicate"):
    """
    flask_nums = flasks used in winkler titrations
    flasks_list = entire list of flask numbers from .vol file
    flasks_vols = entire list of flask volumes from .vol file

    """


    ### NOTE: LOOK AT ANDREW DICKSON'S SOP 13: GRAVIMETRIC CALIBRATION OF A VOLUME CONTAINED USING WATER!!!
    ### These formulas/values may be out of date
    _coef = {
            "borosilicate": 0.00001,
            "soft": 0.000025,
            }
    _t = 20

    # Make flask_nums into DataFrame for merging

    #flask_num_df = pd.DataFrame({'BOTTLENO_OXY':sample_nums, 'FLASKNO':flask_nums})
    flask_num_df = pd.DataFrame({'FLASKNO':flask_nums})
    flask_list_vols = pd.DataFrame({'FLASKNO':flask_list,'FLASK_VOL':flask_vols})



    coef = _coef[glass]
    merge_df = flask_num_df.merge(flask_list_vols, how='left')

    # code formerly in fit_ctd, need to sort this out w/ above NOTE
    # flasks = flask_load()
    # fv = flasks[flask]
    # coef = _coef[glass]
    # return fv * (1.0 + coef * (t - _t));

    volumes = merge_df['FLASK_VOL'].values
    #volumes = merge_df['FLASK_VOL']

    #bottle_num = merge_df['BOTTLENO_OXY'].values
#    merge_df = merge_df[merge_df[ssscc_col].notnull()]

#    merge_df = merge_df.sort_values(by='master_index')
#    merge_df = merge_df.set_index(df.index.values)

    return volumes #bottle_num

def gather_oxy_params(oxyfile):

    """
    Collects winkler oxygen measurement parameters into a DataFrame (First line
    of output from Labview program)

    Parameters
    ----------

    oxyfile :string
             Path to labview output file for winkler titration

    Returns
    -------

    df :dataframe
        DataFrame containing oxygen measurement parameters: TITR, BLANK, KIO3_N
        KIO3_V, KIO3_T, and THIO_T
    """

    f = open(oxyfile,newline='')
    oxyF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
    oxy_Array = []

    for row in oxyF:
        oxy_Array.append(row)

    f.close()

    params = oxy_Array[0]
    df = pd.DataFrame(params)
    df = df.transpose()

    df = df[[0,1,2,3,4,5]]

# If column names are needed later on for the parameters DF
#    cols = ['TITR','BLANK','KIO3_N','KIO3_V','KIO3_T','THIO_T']
#
#    df.rename(columns={0:cols[0],1:cols[1],2:cols[2],3:cols[3],4:cols[4],5:cols[5]})

    return df

def gather_all_oxy_params(ssscc,oxyfile_prefix=cfg.directory["oxy"],oxyfile_postfix=''):
    """
    Collects all oxygen parameters for a given SSSCC (multiple stations)

    Parameters
    ----------

    ssscc :list
           List containing the station/casts to have parameters collected from.

    oxyfile_prefix :string
                    Path string prefix to file
    oxyfile_postfix :string
                     Postfix for oxygen_file (any extensions)

    Returns
    -------

    df_data_all :Pandas DataFrame
                 DataFrame containing all of the oxygen parameters for each
                 station.
    """

    df_data_all = pd.DataFrame()

    for x in ssscc:
        oxyfile = oxyfile_prefix + x + oxyfile_postfix
        params = gather_oxy_params(oxyfile)
        params['SSSCC'] = x

        df_data_all = pd.concat([df_data_all,params])

    df_data_all = df_data_all.reset_index()
    df_data_all = df_data_all.set_index('SSSCC')
    df_data_all.drop('index', axis=1, inplace=True)
    return df_data_all

def thio_n_calc(params,ssscc):#(titr, blank, kio3_n, kio3_v, kio3_t, thio_t):

    """
    Calculates normality of thiosulfate used in Winkler oxygen titrations

    Parameters
    ----------

    params :DataFrame
            Table listing calibration parameters used for each station winkler
            titrations:
                Titr, Blank, KIO3_Norm, KIO3_Vol, KIO3_Temp, and Thio temp
            ex:
                         0         1          2        3      4      5
                SSSCC
                14402  0.53388  -0.00012  0.0123631  10.0019  22.69  22.64
                14501  0.53388  -0.00012  0.0123631  10.0019  22.69  22.64
                14601  0.53378  -0.00011  0.0123631  10.0019  21.56  22.08
                14701  0.53378  -0.00011  0.0123631  10.0019  21.56  22.08

    ssscc :array-like
           Column of station/cast values FOR EACH measurement taken:
           (The length of ssscc should be the same as the amount of oxygen
           titrations taken)

    Returns
    -------

    thio_n :array-like
            Thiosulfate normality

    """

    params = params.loc[ssscc]
    params = params.apply(pd.to_numeric)

    titr   = params[0]
    blank  = params[1]
    kio3_n = params[2]
    kio3_v = params[3]
    kio3_t = params[4]
    thio_t = params[5]

    rho_stp  = gsw.rho_t_exact(0, 20, 0)
    rho_kio  = gsw.rho_t_exact(0, kio3_t, 0)
    rho_thio = gsw.rho_t_exact(0, thio_t, 0)

    kio3_v_20c = kio3_v * (rho_kio/rho_stp)
    thio_v_20c = titr * (rho_thio/rho_stp) - blank

    thio_n = kio3_v_20c * kio3_n / thio_v_20c

    return thio_n.values

def titr_20_calc(titr, titr_temp):


    rho_titr = gsw.rho_t_exact(0, titr_temp, 0)
    rho_stp = gsw.rho_t_exact(0, 20, 0)

    #convert the titr ammount to 20c equivalent
    titr_20c = titr * rho_titr/rho_stp
    return titr_20c

def calculate_bottle_oxygen(ssscc_list, ssscc_col, titr_vol, titr_temp, flask_nums, d_temp='DRAW_TEMP'):
    """
    Calculates oxygen values from Winkler titrations

    Parameters
    ----------

    ssscc_list :array-like
                list containing all of the stations to be calculated


    """



    params = gather_all_oxy_params(ssscc_list)

    thio_n = thio_n_calc(params, ssscc_col)

    titr_20C = titr_20_calc(titr_vol, titr_temp)

    flask_df = flask_load()
    flask_list = flask_df['FLASKNO']
    flask_vols = flask_df['FLASK_VOL']

    volumes = match_flask(flask_nums, flask_list, flask_vols)

    params = params.loc[ssscc_col]
    params = params.apply(pd.to_numeric)

    blank  = params[1]

    blank = blank.reset_index()
    blank = blank[1].apply(pd.to_numeric)

    oxy_mlL = oxygen_eq(titr_20C, blank.values, thio_n, volumes)

    return oxy_mlL

def flag_winkler_oxygen(oxygen):
    flag = pd.Series(oxygen).copy()
    flag.loc[flag.notnull()] = 2
    flag.loc[flag.isnull()] = 9
    flag = flag.astype(int)
    return flag

def oxygen_eq(titr,blank,thio_n,flask_vol):

    oxy_mlL = (((titr - blank) * thio_n * 5.598 - 0.0017) / ((flask_vol - 2.0) * 0.001))

    return oxy_mlL

def oxy_from_ctd_eq(coef, oxyvolts, pressure, temp, dvdt, os, cc=[1.92634e-4,-4.64803e-2]):
    """
    Modified oxygen equation for SBE 43 used by NOAA/PMEL

    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E
    """


   # ox=c(1)*(doxvo+c(2)+c(3)*exp(cc(1)*dpres+cc(2)*dtepr).*dsdvdt).*os.*exp(c(4)*dtepr).*exp(c(5)*dpres./(273.15+dtepr));
    oxy_mlL = coef[0] * (oxyvolts + coef[1] + coef[2] * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt) * os \
            * np.exp(coef[3] * temp) \
            * np.exp((coef[4] * pressure) \
            / (temp + 273.15))


    return oxy_mlL

def SB_oxy_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc=[1.92634e-4,-4.64803e-2]):
    """
    Oxygen (ml/l) = Soc * (V + Voffset) * (1.0 + A * T + B * T^2 + C * T^3 ) * OxSol(T,S) * exp(E * P / K)
    Original equation from calib sheet dated 2014:

    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = A
    coef[3] = B
    coef[4] = C
    coef[5] = E
    coef[6] = Tau20

    """

    oxygen = (coef[0]) * ((oxyvolts + coef[1] + (coef[6] * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt))
                 * (1.0 + coef[2] * temp + coef[3] * np.power(temp,2) + coef[4] * np.power(temp,3) )
                 * os
                 * np.exp((coef[5] * pressure) / (temp + 273.15)))

    return oxygen

def get_SB_coef(hexfile,xmlfile):

    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1,
                        'ranking': 5, 'column': 'CTDOXYVOLTS',
                        'sensor_info': rawConfig['Sensors'][i]}


    Soc = oxy_meta['sensor_info']['Soc'] # Oxygen slope
    Voff = oxy_meta['sensor_info']['offset'] # Sensor output offset voltage
    Tau20 = oxy_meta['sensor_info']['Tau20'] #Sensor time constant at 20 deg C and 1 Atm
    Tcorr = oxy_meta['sensor_info']['Tcor'] # Temperature correction
    E  = oxy_meta['sensor_info']['E'] #Compensation Coef for pressure effect on membrane permeability (Atkinson et al)
    A = oxy_meta['sensor_info']['A'] # Compensation Coef for temp effect on mem perm
    B = oxy_meta['sensor_info']['B'] # Compensation Coef for temp effect on mem perm
    C = oxy_meta['sensor_info']['C'] # Compensation Coef for temp effect on mem perm
    D = [oxy_meta['sensor_info']['D0'],oxy_meta['sensor_info']['D1'],oxy_meta['sensor_info']['D2']] # Press effect on time constant Usually not fitted for.

    #coef = {'Soc':Soc,'Voff':Voff,'Tau20':Tau20,'Tcorr':Tcorr,'E':E,'A':A,'B':B,'C':C,'D':D}

    # Make into an array in order to do keast squares fitting easier
#        coef0s:
#    coef[0] = Soc
#    coef[1] = Voffset
#    coef[2] = Tau20
#    coef[3] = Tcorr
#    coef[4] = E
#
#    cc[0] = D1
#    cc[1] = D2
    coef = [Soc,Voff,A,B,C,E,Tau20]

    return coef


def apply_oxy_coef(oxyvolts, press, temp, sal, time, lon, lat, coef, cc=[1.92634e-4,-4.64803e-2]):

    sigma0 = sigma_from_CTD(sal, temp, press, lon, lat)

    os = sbe_eq.OxSol(temp, sal)
    dv_dt = calculate_dVdT(oxyvolts, time)

    new_OXY = coef[0] * (oxyvolts + coef[1] + coef[2] \
            * np.exp(cc[0] * press + cc[0] \
            * temp) \
            * dv_dt) * os \
            * np.exp(coef[3] * temp) \
            * np.exp((coef[4] * press) \
            / (temp + 273.15))

    oxy_uMolkg = oxy_ml_to_umolkg(new_OXY,sigma0)

    return oxy_uMolkg

def sigma_from_CTD(sal, temp, press, lon, lat, ref=0):
    """
    Calculates potential density from CTD parameters at various reference pressures

    Parameters
    ----------

    sal :array-like
         Salinity in PSU (PSS-78)

    temp :array_like
          In-situ temperature in deg C

    press :array_like
           Pressure in dbar

    lon :array_like
         longitude in decimal degrees

    lat :array_like
         latitute in decimal degrees

    ref :int
         reference pressure point for caluclateing sigma0 (in ref * 1000 dBar ref[0-4])

    Returns
    -------

    simga :array-like
           Potential density calculated at a reference pressure of ref * 1000 dBar

    """

    CT = gsw.CT_from_t(sal, temp, press)

    SA = gsw.SA_from_SP(sal, press, lon, lat)

    # Reference pressure in ref*1000 dBars
    if ref == 0:
        sigma = gsw.sigma0(SA, CT)
    elif ref == 1:
        sigma = gsw.sigma1(SA, CT)
    elif ref == 2:
        sigma = gsw.sigma2(SA, CT)
    elif ref == 3:
        sigma = gsw.sigma3(SA, CT)
    elif ref == 4:
        sigma = gsw.sigma4(SA, CT)

    return sigma

def oxy_ml_to_umolkg(oxy_mL_L, sigma0):
    """Convert dissolved oxygen from units of mL/L to micromol/kg.

    Parameters
    ----------
    oxy_mL_L : array-like
        Dissolved oxygen in units of [mL/L]
    sigma0 : array-like
        Potential density anomaly (i.e. sigma - 1000) referenced to 0 dbar [kg/m^3]

    Returns
    -------
    oxy_umol_kg : array-like
        Dissolved oxygen in units of [umol/kg]

    Notes
    -----
    Conversion value 44660 is exact for oxygen gas and derived from the ideal gas law.
    (c.f. Sea-Bird Application Note 64, pg. 6)
    """

    oxy_umol_kg = oxy_mL_L * 44660 / (sigma0 + 1000)

    return oxy_umol_kg

def calculate_dVdT(oxyvolts, time):
#
     doxyv = np.diff(oxyvolts)


     dt = np.diff(time)

     # Replace 0 values with median value
     m = np.median(dt[dt > 0])
     dt[dt == 0] = m


     dv_dt = doxyv/dt
     dv_dt = np.insert(dv_dt,0,0)
     dv_dt_mean = np.mean(dv_dt[(dv_dt != np.inf) & (dv_dt !=-np.inf)])
     dv_dt[(dv_dt == np.inf) | (dv_dt == -np.inf)] = np.NaN
     #dv_dt = dv_dt.replace([np.inf,-np.inf],np.NaN)
     #dv_dt = dv_dt.fillna(dv_dt_mean)

     np.nan_to_num(dv_dt_mean)

     #dv_dt_conv = np.convolve(dv_dt,[0.5,0.5])

#     a = 1
#     windowsize = 5
#     b = (1/windowsize)*np.ones(windowsize)
#     #filtered_dvdt = scipy.signal.filtfilt(b,a,dv_dt_conv)
#     filtered_dvdt = scipy.signal.filtfilt(b,a,dv_dt)


     return dv_dt#filtered_dvdt


def _get_sbe_coef(idx=0):
    """
    Get SBE oxygen coefficients from raw data files.
    Defaults to using first station in ssscc.csv file.
    """
    ssscc_list = process_ctd.get_ssscc_list()
    station = ssscc_list[idx]
    hexfile = cfg.directory["raw"] + station + ".hex"
    xmlfile = cfg.directory["raw"] + station + ".XMLCON"

    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1,
                        'ranking': 5, 'column': 'CTDOXYVOLTS',
                        'sensor_info': rawConfig['Sensors'][i]}


    Soc = oxy_meta['sensor_info']['Soc'] # Oxygen slope
    Voff = oxy_meta['sensor_info']['offset'] # Sensor output offset voltage
    Tau20 = oxy_meta['sensor_info']['Tau20'] #Sensor time constant at 20 deg C and 1 Atm
    Tcorr = oxy_meta['sensor_info']['Tcor'] # Temperature correction
    E  = oxy_meta['sensor_info']['E'] #Compensation Coef for pressure effect on membrane permeability (Atkinson et al)
    A = oxy_meta['sensor_info']['A'] # Compensation Coef for temp effect on mem perm
    B = oxy_meta['sensor_info']['B'] # Compensation Coef for temp effect on mem perm
    C = oxy_meta['sensor_info']['C'] # Compensation Coef for temp effect on mem perm
    D = [oxy_meta['sensor_info']['D0'],oxy_meta['sensor_info']['D1'],oxy_meta['sensor_info']['D2']] # Press effect on time constant Usually not fitted for.

    coef = [Soc,Voff,Tau20,Tcorr,E]#,A,B,C,D]

    return coef

def calculate_weights(pressure):

    eps=1e-5
    wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
    wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]#[20,20,25,25,50,50,100,100,200,200,500,500]
    wgt = scipy.interpolate.interp1d(wrow1,wrow2)

    weights = wgt(pressure)

    return weights


def oxy_equation(X, Soc, Voffset, A, B, C, E, Tau20):
    # eq (3) in Uchida CTD manual
    cc=[1.92634e-4,-4.64803e-2]
    oxyvolts, pressure, temp, dvdt, os = X

    oxygen = (Soc * ((oxyvolts + Voffset + (Tau20 * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt))
                 * (1.0 + A * temp + B * np.power(temp,2) + C * np.power(temp,3) )
                 * os
                 * np.exp((E * pressure) / (temp + 273.15))))

    return oxygen


def _PMEL_oxy_eq(coefs,inputs,cc=[1.92634e-4,-4.64803e-2]):
    """
    Modified oxygen equation for SBE 43 used by NOAA/PMEL
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E
    """
    Soc, Voff, Tau20, Tcorr, E = coefs
    oxyvolts, pressure, temp, dvdt, os = inputs
    o2 = Soc * (oxyvolts + Voff + Tau20 * np.exp(cc[0] * pressure + cc[1] * (temp - 20)) * dvdt) * os \
            * np.exp(Tcorr * temp) \
            * np.exp((E * pressure) / (temp + 273.15))

    return o2

# TODO: optionally include other residual types
# (abstracted from PMEL code oxygen_cal_ml.m)
# unweighted L2: sum((ref - oxy)^2)  # if weighted fails
# unweighted L4: sum((ref - oxy)^4)  # unsure of use case
# unweighted L1: sum(abs(ref - oxy))  # very far from ideal
# anything else? genericize with integer "norm" function input?
def PMEL_oxy_weighted_residual(coefs,weights,inputs,refoxy):
    return np.sum((weights*(refoxy-_PMEL_oxy_eq(coefs, inputs))**2))/np.sum(weights**2)

def match_sigmas(btl_prs, btl_oxy, btl_tmp, btl_SA, ctd_os, ctd_prs, ctd_tmp, ctd_SA, ctd_oxyvolts, ctd_time):

    # Construct Dataframe from bottle and ctd values for merging
    btl_data = pd.DataFrame(data={
        "CTDPRS": btl_prs,
        "REFOXY": btl_oxy,
        "CTDTMP": btl_tmp,
        "SA": btl_SA,
    })
    time_data = pd.DataFrame(data={
        "CTDPRS": ctd_prs,
        "OS": ctd_os,
        "CTDTMP": ctd_tmp,
        "SA": ctd_SA,
        "CTDOXYVOLTS": ctd_oxyvolts,
        "CTDTIME": ctd_time,
    })
    time_data["dv_dt"] = calculate_dVdT(time_data["CTDOXYVOLTS"], time_data["CTDTIME"])

    # Merge DF
    merged_df = pd.DataFrame(
        columns=["CTDPRS", "CTDOXYVOLTS", "CTDTMP", "dv_dt", "OS"], dtype=float
    )
    merged_df["REFOXY"] = btl_data["REFOXY"].copy()

    # calculate sigma referenced to multiple depths
    for idx, p_ref in enumerate([0, 1000, 2000, 3000, 4000, 5000, 6000]):
        btl_data[f"sigma{idx}"] = (
            gsw.pot_rho_t_exact(
                btl_data["SA"],
                btl_data["CTDTMP"],
                btl_data["CTDPRS"],
                p_ref,
            )
            - 1000  # subtract 1000 to get potential density *anomaly*
        ) + 1e-8*np.random.standard_normal(btl_data["SA"].size)
        time_data[f"sigma{idx}"] = (
            gsw.pot_rho_t_exact(
                time_data["SA"],
                time_data["CTDTMP"],
                time_data["CTDPRS"],
                p_ref,
            )
            - 1000  # subtract 1000 to get potential density *anomaly*
        ) + 1e-8*np.random.standard_normal(time_data["SA"].size)
        rows = (btl_data["CTDPRS"] > (p_ref - 500)) & (btl_data["CTDPRS"] < (p_ref + 500))
        time_sigma_sorted = time_data[f"sigma{idx}"].sort_values().to_numpy()
        sigma_min = np.min([np.min(btl_data.loc[rows, f"sigma{idx}"]), np.min(time_sigma_sorted)])
        sigma_max = np.max([np.max(btl_data.loc[rows, f"sigma{idx}"]), np.max(time_sigma_sorted)])
        time_sigma_sorted = np.insert(time_sigma_sorted, 0, sigma_min - 1e-4)
        time_sigma_sorted = np.append(time_sigma_sorted, sigma_max + 1e-4)
        # TODO: can this be vectorized?
        cols = ["CTDPRS", "CTDOXYVOLTS", "CTDTMP", "dv_dt", "OS"]
        inds = np.concatenate(([0], np.arange(0, len(time_data)), [len(time_data) - 1]))
        for col in cols:
            merged_df.loc[rows, col] = np.interp(
                btl_data.loc[rows, f"sigma{idx}"],
                time_sigma_sorted,
                time_data[col].iloc[inds],
            )

    # Apply coef and calculate CTDOXY
    sbe_coef0 = _get_sbe_coef() # initial coefficient guess
    merged_df['CTDOXY'] = _PMEL_oxy_eq(sbe_coef0, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']))

    return merged_df


def sbe43_oxy_fit(merged_df, sbe_coef0=None, f_suffix=None):

    # Plot data to be fit together
    f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual{f_suffix}_prefit.pdf"
    ctd_plots._intermediate_residual_plot(
        merged_df['REFOXY'] - merged_df['CTDOXY'],
        merged_df["CTDPRS"],
        merged_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10,10)
    )

    # Create DF for good and questionable values
    bad_df = pd.DataFrame()
    good_df = pd.DataFrame()

    if sbe_coef0 is None:
        sbe_coef0 = _get_sbe_coef()  # load initial coefficient guess

    p0 = sbe_coef0[0], sbe_coef0[1], sbe_coef0[2], sbe_coef0[3], sbe_coef0[4]
    
    # Curve fit (weighted)
    weights = calculate_weights(merged_df['CTDPRS'])
    res = scipy.optimize.minimize(PMEL_oxy_weighted_residual,x0=p0,args=(weights, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']), merged_df['REFOXY']), bounds=[(None,None),(None,None),(0,None),(None,None),(None,None)])
    cfw_coefs = res.x
    merged_df['CTDOXY'] = _PMEL_oxy_eq(cfw_coefs, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']))        

    merged_df['residual'] = merged_df['REFOXY'] - merged_df['CTDOXY']
    stdres = np.std(merged_df['residual'])
    cutoff = stdres * 2.8

    thrown_values = merged_df[np.abs(merged_df['residual']) > cutoff]
    bad_df = pd.concat([bad_df, thrown_values])
    merged_df = merged_df[np.abs(merged_df['residual']) <= cutoff]

    while not thrown_values.empty: # runs as long as there are thrown_values

        p0 = cfw_coefs[0], cfw_coefs[1], cfw_coefs[2], cfw_coefs[3], cfw_coefs[4]
        weights = calculate_weights(merged_df['CTDPRS'])
        res = scipy.optimize.minimize(PMEL_oxy_weighted_residual,x0=p0,args=(weights, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']), merged_df['REFOXY']), bounds=[(None,None),(None,None),(0,None),(None,None),(None,None)])
        cfw_coefs = res.x
        merged_df['CTDOXY'] = _PMEL_oxy_eq(cfw_coefs, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS'], merged_df['CTDTMP'], merged_df['dv_dt'], merged_df['OS']))

        merged_df['residual'] = merged_df['REFOXY'] - merged_df['CTDOXY']
        stdres = np.std(merged_df['residual'])
        cutoff = stdres * 2.8
        thrown_values = merged_df[np.abs(merged_df['residual']) > cutoff]
        print(len(thrown_values))
        print(p0)
        bad_df = pd.concat([bad_df, thrown_values])
        merged_df = merged_df[np.abs(merged_df['residual']) <= cutoff]

    # intermediate plots to diagnose data chunks goodness
    # TODO: implement into bokeh/flask dashboard
    if f_suffix is not None:
        f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual{f_suffix}.pdf"
        ctd_plots._intermediate_residual_plot(
            merged_df["residual"],
            merged_df["CTDPRS"],
            merged_df["SSSCC"],
            xlabel="CTDOXY Residual (umol/kg)",
            f_out=f_out,
            xlim=(-10,10)
        )

    # good_df = pd.concat([good_df, merged_df])
    merged_df['CTDOXY_FLAG_W'] = 2
    # bad_df = pd.concat([bad_df, bad_values])
    bad_df['CTDOXY_FLAG_W'] = 3
    df = pd.concat([merged_df,bad_df])

    # df['SSSCC_int'] = df['SSSCC_sbe43'].astype(int)
    # df.sort_values(by=['SSSCC_int','btl_fire_num'],ascending=[True,True],inplace=True)
    # df.drop('SSSCC_int', axis=1, inplace=True)

    return cfw_coefs, df

def prepare_oxy(btl_df, time_df, ssscc_list):
    """
    Calculate oxygen-related variables needed for calibration:
    sigma, oxygen solubility (OS), and bottle oxygen

    Parameters
    ----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data
    ssscc_list : list of str
        List of stations to process

    Returns
    -------

    """
    # Calculate SA and CT
    btl_df["SA"] = gsw.SA_from_SP(
        btl_df[cfg.column["sal"]],
        btl_df[cfg.column["p_btl"]],
        btl_df[cfg.column["lon_btl"]],
        btl_df[cfg.column["lat_btl"]],
    )
    btl_df["CT"] = gsw.CT_from_t(
        btl_df["SA"],
        btl_df[cfg.column["t1_btl"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p_btl"]],
    )
    time_df["SA"] = gsw.SA_from_SP(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon_btl"]],
        time_df[cfg.column["lat_btl"]],
    )
    time_df["CT"] = gsw.CT_from_t(
        time_df["SA"],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
    )
    # calculate sigma
    btl_df["sigma_btl"] = sigma_from_CTD(
        btl_df[cfg.column["sal"]],
        btl_df[cfg.column["t1_btl"]],  # oxygen sensor is on primary line (ie t1)
        btl_df[cfg.column["p_btl"]],
        btl_df[cfg.column["lon_btl"]],
        btl_df[cfg.column["lat_btl"]],
    )
    time_df["sigma_ctd"] = sigma_from_CTD(
        time_df[cfg.column["sal"]],
        time_df[cfg.column["t1"]],  # oxygen sensor is on primary line (ie t1)
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon_btl"]],
        time_df[cfg.column["lat_btl"]],
    )
    # Calculate oxygen solubility in µmol/kg
    btl_df["OS"] = gsw.O2sol(
        btl_df["SA"],
        btl_df["CT"],
        btl_df[cfg.column["p_btl"]],
        btl_df[cfg.column["lon_btl"]],
        btl_df[cfg.column["lat_btl"]],
    )
    time_df["OS"] = gsw.O2sol(
        time_df["SA"],
        time_df["CT"],
        time_df[cfg.column["p"]],
        time_df[cfg.column["lon"]],
        time_df[cfg.column["lat"]],
    )
    # Convert CTDOXY units
    btl_df["CTDOXY"] = oxy_ml_to_umolkg(btl_df["CTDOXY1"], btl_df["sigma_btl"])
    # Calculate bottle oxygen
    btl_df[cfg.column["oxy_btl"]] = calculate_bottle_oxygen(
        ssscc_list,
        btl_df["SSSCC"],
        btl_df["TITR_VOL"],
        btl_df["TITR_TEMP"],
        btl_df["FLASKNO"],
    )
    btl_df[cfg.column["oxy_btl"]] = oxy_ml_to_umolkg(
        btl_df[cfg.column["oxy_btl"]], btl_df["sigma_btl"]
    )
    btl_df["OXYGEN_FLAG_W"] = flag_winkler_oxygen(btl_df[cfg.column["oxy_btl"]])
    # Load manual OXYGEN flags
    if Path("data/oxygen/manual_oxy_flags.csv").exists():
        manual_flags = pd.read_csv(
            "data/oxygen/manual_oxy_flags.csv", dtype={"SSSCC": str}
        )
        for _, flags in manual_flags.iterrows():
            df_row = (btl_df["SSSCC"] == flags["SSSCC"]) & (
                btl_df["btl_fire_num"] == flags["Bottle"]
            )
            btl_df.loc[df_row, "OXYGEN_FLAG_W"] = flags["Flag"]

    return True


def calibrate_oxy(btl_df, time_df, ssscc_list):
    """
    Non-linear least squares fit chemical sensor oxygen against bottle oxygen.

    Parameters
    ----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data
    ssscc_list : list of str
        List of stations to process

    Returns
    -------

    """
    # Plot all pre fit data
    f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual_all_prefit.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df['OXYGEN'] - btl_df['CTDOXY'],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10,10)
    )
    # Prep vars, dfs, etc.
    all_sbe43_merged = pd.DataFrame()
    sbe43_dict = {}
    all_sbe43_fit = pd.DataFrame()

    btl_df["dv_dt"] = np.nan  # initialize column
    # Density match time/btl oxy dataframes
    for ssscc in ssscc_list:
        time_data = time_df[time_df["SSSCC"] == ssscc].copy()
        btl_data = btl_df[btl_df["SSSCC"] == ssscc].copy()
        # can't calibrate without bottle oxygen ("OXYGEN")
        if (btl_data["OXYGEN_FLAG_W"] == 9).all():
            sbe43_dict[ssscc] = np.full(5, np.nan)
            print(ssscc + " skipped, all oxy data is NaN")
            continue
        sbe43_merged = match_sigmas(
            btl_data[cfg.column["p_btl"]],
            btl_data[cfg.column["oxy_btl"]],
            btl_data["CTDTMP1"],
            btl_data["SA"],
            time_data["OS"],
            time_data[cfg.column["p"]],
            time_data[cfg.column["t1"]],
            time_data["SA"],
            time_data[cfg.column["oxyvolts"]],
            time_data["scan_datetime"],
        )
        sbe43_merged = sbe43_merged.reindex(btl_data.index)  # add nan rows back in
        btl_df.loc[btl_df["SSSCC"] == ssscc, ["CTDOXYVOLTS","dv_dt","OS"]] = sbe43_merged[["CTDOXYVOLTS","dv_dt","OS"]]
        sbe43_merged["SSSCC"] = ssscc
        all_sbe43_merged = pd.concat([all_sbe43_merged, sbe43_merged])
        print(ssscc + " density matching done")

    # Only fit using OXYGEN flagged good (2)
    all_sbe43_merged = all_sbe43_merged[btl_df["OXYGEN_FLAG_W"] == 2].copy()

    # Fit ALL oxygen stations together to get initial coefficient guess
    (sbe_coef0, _) = sbe43_oxy_fit(all_sbe43_merged, f_suffix="_ox0")

    # Fit each cast individually
    for ssscc in ssscc_list:
        sbe_coef, sbe_df = sbe43_oxy_fit(
            all_sbe43_merged.loc[all_sbe43_merged["SSSCC"] == ssscc],
            sbe_coef0=sbe_coef0,
            f_suffix=f"_{ssscc}",
        )
        # build coef dictionary
        if ssscc not in sbe43_dict.keys():  # don't overwrite NaN'd stations
            sbe43_dict[ssscc] = sbe_coef
        # all non-NaN oxygen data with flags
        all_sbe43_fit = pd.concat([all_sbe43_fit, sbe_df])

    # TODO: save outlier data from fits?
    # TODO: secondary oxygen flagging step (instead of just taking outliers from fit routine)

    # apply coefs
    time_df["CTDOXY"] = np.nan
    for ssscc in ssscc_list:
        if np.isnan(sbe43_dict[ssscc]).all():
            print(ssscc + " missing oxy data, leaving nan values and flagging as 9")
            time_df.loc[time_df["SSSCC"] == ssscc, "CTDOXY_FLAG_W"] = 9
            time_df.loc[time_df["SSSCC"] == ssscc, "RINKO_FLAG_W"] = 9
            continue
        btl_rows = (btl_df["SSSCC"] == ssscc).values
        time_rows = (time_df["SSSCC"] == ssscc).values
        btl_df.loc[btl_rows, "CTDOXY"] = _PMEL_oxy_eq(
            sbe43_dict[ssscc],
            (
                btl_df.loc[btl_rows, cfg.column["oxyvolts"]],
                btl_df.loc[btl_rows, cfg.column["p_btl"]],
                btl_df.loc[btl_rows, cfg.column["t1_btl"]],
                btl_df.loc[btl_rows, "dv_dt"],
                btl_df.loc[btl_rows, "OS"],
            ),
        )
        print(ssscc + " btl data fitting done")
        time_df.loc[time_rows, "CTDOXY"] = _PMEL_oxy_eq(
            sbe43_dict[ssscc],
            (
                time_df.loc[time_rows, cfg.column["oxyvolts"]],
                time_df.loc[time_rows, cfg.column["p"]],
                time_df.loc[time_rows, cfg.column["t1"]],
                time_df.loc[time_rows, "dv_dt"],
                time_df.loc[time_rows, "OS"],
            ),
        )
        print(ssscc + " time data fitting done")

    # TODO: flag oxy data here? compare w/ T/C routines

    # Plot all post fit data
    f_out = f"{cfg.directory['ox_fit_figs']}sbe43_residual_all_postfit.pdf"
    ctd_plots._intermediate_residual_plot(
        btl_df['OXYGEN'] - btl_df['CTDOXY'],
        btl_df["CTDPRS"],
        btl_df["SSSCC"],
        xlabel="CTDOXY Residual (umol/kg)",
        f_out=f_out,
        xlim=(-10,10)
    )

    # export fitting coefs
    sbe43_coefs = pd.DataFrame.from_dict(
        sbe43_dict, orient="index", columns=["Soc", "Voffset", "Tau20", "Tcorr", "E"]
    )
    sbe43_coefs.to_csv(cfg.directory["logs"] + "sbe43_coefs.csv")
    
    return True

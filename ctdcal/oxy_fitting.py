#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:04:19 2018

@author: k3jackson
"""



import scipy
import numpy as np
#import sys
#sys.path.append('ctdcal/')
#import ctdcal.process_ctd as process_ctd
import ctdcal.fit_ctd as fit_ctd
import ctdcal.sbe_reader as sbe_rd
import ctdcal.sbe_equations_dict as sbe_eq
import gsw
import pandas as pd
import csv

#Line 342 module isopycnals

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

    header = ['STNNO_OXY','CASTNO_OXY','BOTTLENO_OXY','FLASKNO','TITR_VOL','TITR_TEMP','DRAW_TEMP','TITER_TIME','END_VOLTS']

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

def flask_load(flaskfile='data/oxygen/o2flasks.vol',skip_rows=12):
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

def gather_all_oxy_params(ssscc,oxyfile_prefix='data/oxygen/',oxyfile_postfix=''):
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

    rho_stp  = fit_ctd.rho_t(20)
    rho_kio  = fit_ctd.rho_t(kio3_t)
    rho_thio = fit_ctd.rho_t(thio_t)

    kio3_v_20c = kio3_v * (rho_kio/rho_stp)
    thio_v_20c = titr * (rho_thio/rho_stp) - blank

    thio_n = kio3_v_20c * kio3_n / thio_v_20c
    thio_n = thio_n.values

    return thio_n

def titr_20_calc(titr, titr_temp):


    rho_titr = fit_ctd.rho_t(titr_temp)
    rho_stp = fit_ctd.rho_t(20)

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

def os_umol_kg( sal, PT):
    """
    Calculates oxygen solubility in umol/kg as found in Gordon and Garcia 1992
    (Taken from GSW, use this until gsw releases a version in their toolbox)

    From GSW:
    Note that this algorithm has not been approved by IOC and is not work
    from SCOR/IAPSO Working Group 127. It is included in the GSW
    Oceanographic Toolbox as it seems to be oceanographic best practice.

    Paramteters
    -----------

    sal :array-like
         Practical Salinity (PSS-78)

    PT :array-like
        Potential Temperature

    Returns
    -------

    O2sol_umol :array_like
                Oxygen Solubility in mirco-moles per kilogram (µmol/kg)

    """



    a0 = 5.80871
    a1 = 3.20291
    a2 = 4.17887
    a3 = 5.10006
    a4 = -9.86643e-2
    a5 = 3.80369
    b0 = -7.01577e-3
    b1 = -7.70028e-3
    b2 = -1.13864e-2
    b3 = -9.51519e-3
    c0 = -2.75915e-7


    # Calculate potential temp

    PT68 = PT * 1.00024 # the 1968 International Practical Temperature Scale IPTS-68.

    y = np.log((298.15 - PT68) / (273.15 + PT))

    O2sol_umol = np.exp(a0 + y * (a1 + y * (a2 + y * (a3 + y * (a4 + a5 * y)))) + sal * (b0 + y * (b1 + y * (b2 + b3 * y)) + c0 * sal))

    return O2sol_umol

def oxy_ml_to_umolkg(oxy_mlL,sigma0):

    oxy_uMolkg = oxy_mlL * 44660 / (sigma0 + 1000)

    return oxy_uMolkg

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

def merge_parameters(btl_df,time_df,l_param='sigma0_btl',r_param='sigma0_ctd'):

    merge_df = pd.merge_asof(btl_df,time_df,left_on=l_param,right_on=r_param,direction='nearest')

    # Clean up merged dataframe
    # These column names can be taken as a variable from the ini file
    keep_columns = ['CTDTMP1_y','CTDTMP2_y', 'CTDPRS_y', 'CTDCOND1_y', 'CTDCOND2_y', 'CTDSAL_y',
       'CTDOXY1_y', 'CTDOXYVOLTS_y', 'FREE1_y', 'FREE2_y', 'FREE3_y',
       'FREE4_y', 'FLUOR_y', 'CTDBACKSCATTER_y', 'CTDXMISS_y', 'ALT_y',
       'REF_PAR_y', 'GPSLAT_y', 'GPSLON_y', 'new_fix_y', 'pressure_temp_int_y',
       'pump_on_y', 'btl_fire_y', 'scan_datetime_y', 'SSSCC_y',
       'master_index_y','OS_y','dv_dt_y','CTDOXY',r_param]

    merge_df = merge_df[keep_columns]

    col_names = ['CTDTMP1','CTDTMP2', 'CTDPRS', 'CTDCOND1', 'CTDCOND2','CTDSAL',
       'CTDOXY1', 'CTDOXYVOLTS', 'FREE1', 'FREE2', 'FREE3',
       'FREE4', 'FLUOR', 'CTDBACKSCATTER', 'CTDXMISS', 'ALT',
       'REF_PAR', 'GPSLAT', 'GPSLON', 'new_fix', 'pressure_temp_int',
       'pump_on', 'btl_fire', 'scan_datetime', 'SSSCC',
       'master_index','OS','dv_dt','CTDOXY', r_param]

    merge_df.columns = [col_names]


    return merge_df

def get_sbe_coef(hexfile, xmlfile):


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
    coef = [Soc,Voff,Tau20,Tcorr,E]#,A,B,C,D]

    return coef

def calculate_weights(pressure):

    eps=1e-5
    wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
    wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]#[20,20,25,25,50,50,100,100,200,200,500,500]
    wgt = scipy.interpolate.interp1d(wrow1,wrow2)

    weights = wgt(pressure)

    return weights

def interpolate_sigma(ctd_sigma,btl_sigma):

#    all_sigma0 = time_data_clean['sigma0_time'].append(btl_data_clean['sigma0_btl'])
#    fsg = pd.Series(np.min(all_sigma0)-1e-4)
#    lsg = pd.Series(np.max(all_sigma0)+1e-4)
#    new_sigma0 = fsg.append(time_data_clean['sigma0_time'])
#    new_sigma0 = new_sigma0.append(lsg)

    all_sigma = ctd_sigma.append(btl_sigma)
    fsg = pd.Series(np.min(all_sigma)-1e-4)
    lsg = pd.Series(np.max(all_sigma)+1e-4)
    new_sigma = fsg.append(ctd_sigma)
    new_sigma = new_sigma.append(lsg)
    new_sigma = new_sigma.reset_index()
    new_sigma = new_sigma[0]


    #x = np.arange(np.size(ctd_sigma))
#    x_inter = np.arange(np.size(new_sigma0))
#    inter_sigma2 = scipy.interpolate.interp1d(x_inter,new_sigma0)
#    new_x = np.linspace(0,np.max(x_inter),np.size(ctd_sigma))
#    new_sigma0=inter_sigma2(new_x)
#    new_sigma0 = pd.Series(new_sigma0)

    return new_sigma

def interpolate_pressure(ctd_pres,btl_pres):

    all_press = ctd_pres.append(btl_pres)
    fsg = pd.Series(np.min(all_press) - 1)
    lsg = pd.Series(np.max(all_press) + 1)
    new_pres = fsg.append(ctd_pres)
    new_pres = new_pres.append(lsg)
    new_pres = new_pres.reset_index()
    new_pres = new_pres[0]

    return new_pres

def interpolate_param(param):
    """
    First extends
    Generates function using x_inter
    calculates values over a new range "new_x"
    x_inter is the x-range of the extended parameter
    new_x is the x-range of the "interpolated-down" parameter

    """
    ex_param = pd.Series(param.iloc[0])
    ex_param = ex_param.append(param)#,index=[])
    ex_param = ex_param.append(pd.Series(param.iloc[-1]))
    ex_param = ex_param.reset_index()
    ex_param = ex_param[0]

#    x_inter = np.arange(len(ex_param))
#
#    inter_param = scipy.interpolate.interp1d(x_inter,ex_param)
#
#    new_x = np.arange(len(param))
#    dparam = inter_param(new_x)

    return ex_param

def oxy_equation(X, Soc, Voffset, A, B, C, E, Tau20):

    cc=[1.92634e-4,-4.64803e-2]
    oxyvolts, pressure, temp, dvdt, os = X

    oxygen = (Soc * ((oxyvolts + Voffset + (Tau20 * np.exp(cc[0] * pressure + cc[1] * temp) * dvdt))
                 * (1.0 + A * temp + B * np.power(temp,2) + C * np.power(temp,3) )
                 * os
                 * np.exp((E * pressure) / (temp + 273.15))))

    return oxygen


def sbe43_oxy_fit(btl_prs, btl_oxy, btl_sigma, ctd_sigma, ctd_os, ctd_prs, ctd_tmp, ctd_oxyvolts, ctd_time, coef0,btl_ssscc=None):


    # Create DF for good and questionable values

    bad_df = pd.DataFrame()
    good_df = pd.DataFrame()

    # Construct Dataframe from bottle and ctd values for merging
    if 'btl_ssscc' in  locals():
        btl_dict = {'CTDPRS_sbe43_btl':btl_prs, 'REFOXY_sbe43':btl_oxy, 'sigma_sbe43_btl':btl_sigma, 'SSSCC_sbe43':btl_ssscc}
    else:
        btl_dict = {'CTDPRS_sbe43_btl':btl_prs, 'REFOXY_sbe43':btl_oxy, 'sigma_sbe43_btl':btl_sigma}
    btl_data = pd.DataFrame(btl_dict)
    time_dict = {'CTDPRS_sbe43_ctd':ctd_prs, 'sigma_sbe43_ctd':ctd_sigma, 'OS_sbe43_ctd':ctd_os, 'CTDTMP_sbe43_ctd':ctd_tmp, 'CTDOXYVOLTS':ctd_oxyvolts, 'CTDTIME':ctd_time}
    time_data = pd.DataFrame(time_dict)

    # Sort DataFrames by sigma0
    time_data.sort_values('sigma_sbe43_ctd', inplace=True)
    btl_data.sort_values('sigma_sbe43_btl', inplace=True)
    btl_data.dropna(subset=['REFOXY_sbe43'], inplace=True)
    # Merge DF
    merged_df = pd.merge_asof(btl_data, time_data, left_on='sigma_sbe43_btl', right_on='sigma_sbe43_ctd', direction='nearest', suffixes=['_btl','_ctd'])

    #Calculate dv_dt
    merged_df['dv_dt'] = calculate_dVdT(merged_df['CTDOXYVOLTS'], merged_df['CTDTIME'])

    # Apply coef and calculate CTDOXY
    merged_df['CTDOXY'] = SB_oxy_eq(coef0, merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd'])
    #merged_df.dropna(subset=['REFOXY_sbe43'], inplace=True)


    # Curve fit (weighted)
    p0 = coef0[0], coef0[1], coef0[2], coef0[3], coef0[4], coef0[5], coef0[6]
    weights = 1/(np.sqrt(merged_df['CTDPRS_sbe43_ctd']))

    try:
        cfw_coef , cov = scipy.optimize.curve_fit(oxy_equation, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd']), merged_df['REFOXY_sbe43'], p0, sigma=weights, absolute_sigma=False)
        merged_df['CTDOXY'] = SB_oxy_eq(cfw_coef, merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd'])

        merged_df['res_sbe43'] = merged_df['REFOXY_sbe43'] - merged_df['CTDOXY']
        stdres = np.std(merged_df['res_sbe43'])
        cutoff = stdres * 2.8

        thrown_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
        bad_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
        bad_df = pd.concat([bad_df, bad_values])
        merged_df = merged_df[np.abs(merged_df['res_sbe43']) <= cutoff]

        while not thrown_values.empty:

            p0 = cfw_coef[0], cfw_coef[1], cfw_coef[2], cfw_coef[3], cfw_coef[4], cfw_coef[5], cfw_coef[6]
            weights = 1/((merged_df['CTDPRS_sbe43_ctd']))
            cfw_coef , cov = scipy.optimize.curve_fit(oxy_equation, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd']), merged_df['REFOXY_sbe43'], p0, sigma=weights, absolute_sigma=False)
            merged_df['CTDOXY'] = SB_oxy_eq(cfw_coef, merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd'])
            merged_df['res_sbe43'] = merged_df['REFOXY_sbe43'] - merged_df['CTDOXY']
            stdres = np.std(merged_df['res_sbe43'])
            cutoff = stdres * 2.8
            thrown_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
            bad_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
            merged_df = merged_df[np.abs(merged_df['res_sbe43']) <= cutoff]

    except RuntimeError:

        try:#Nested try/except could be better
            print('Weighted Curve fitting failed...using Unweighted Fitting')
            cfw_coef , cov = scipy.optimize.curve_fit(oxy_equation, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd']), merged_df['REFOXY_sbe43'], p0)
            merged_df['CTDOXY'] = SB_oxy_eq(cfw_coef, merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd'])

            merged_df['res_sbe43'] = merged_df['REFOXY_sbe43'] - merged_df['CTDOXY']
            stdres = np.std(merged_df['res_sbe43'])
            cutoff = stdres * 2.8

            thrown_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
            bad_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
            bad_df = pd.concat([bad_df, bad_values])
            merged_df = merged_df[np.abs(merged_df['res_sbe43']) <= cutoff]

            while not thrown_values.empty:

                p0 = cfw_coef[0], cfw_coef[1], cfw_coef[2], cfw_coef[3], cfw_coef[4], cfw_coef[5], cfw_coef[6]
                cfw_coef , cov = scipy.optimize.curve_fit(oxy_equation, (merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd']), merged_df['REFOXY_sbe43'], p0)
                merged_df['CTDOXY'] = SB_oxy_eq(cfw_coef, merged_df['CTDOXYVOLTS'], merged_df['CTDPRS_sbe43_ctd'], merged_df['CTDTMP_sbe43_ctd'], merged_df['dv_dt'], merged_df['OS_sbe43_ctd'])
                merged_df['res_sbe43'] = merged_df['REFOXY_sbe43'] - merged_df['CTDOXY']
                stdres = np.std(merged_df['res_sbe43'])
                cutoff = stdres * 2.8
                thrown_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
                bad_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
                merged_df = merged_df[np.abs(merged_df['res_sbe43']) <= cutoff]

        except:
            print('Curve fitting failed...using SBE coef')
            cfw_coef = coef0
            merged_df['res_sbe43'] = merged_df['REFOXY_sbe43'] - merged_df['CTDOXY']
            stdres = np.std(merged_df['res_sbe43'])
            cutoff = stdres * 2.8
            thrown_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
            bad_values = merged_df[np.abs(merged_df['res_sbe43']) > cutoff]
            merged_df = merged_df[np.abs(merged_df['res_sbe43']) <= cutoff]


    good_df = pd.concat([good_df, merged_df])
    good_df['CTDOXY_FLAG_W'] = 2
    bad_df = pd.concat([bad_df, bad_values])
    bad_df['CTDOXY_FLAG_W'] = 3
    df = pd.concat([good_df,bad_df])
    df.sort_values(by='CTDPRS_sbe43_btl',ascending=False,inplace=True)

    return cfw_coef, df

def apply_oxygen_coef_ctd(df, coef_df, ssscc, ssscc_col='SSSCC',oxyvo_col='CTDOXYVOLTS',
                          time_col='scan_datetime',prs_col='CTDPRS',tmp_col='CTDTMP1',
                          os_col='OS_ctd'):

    df['CTDOXY'] = -999
    for station in ssscc:
        coef = coef_df.loc[station[0:3]].values
        mask = (df[ssscc_col] == station)
        time_mask = df[mask].copy()
        time_mask['dv_dt'] = calculate_dVdT(time_mask[oxyvo_col],time_mask[time_col])
        df.loc[mask, 'CTDOXY'] = oxy_equation((time_mask[oxyvo_col],time_mask[prs_col],time_mask[tmp_col],
                                              time_mask['dv_dt'],time_mask[os_col]),coef[0],coef[1],coef[2],
                                                coef[3],coef[4],coef[5],coef[6])

    df['CTDOXY_FLAG_W'] = 2

    return df


def merge_oxy_df(df,oxy_df,btl_stn_col='SSSCC',btl_prs_col='CTDPRS',
                 oxy_stn_col='SSSCC',oxy_prs_col='CTDPRS_btl'):

    df = df.merge(oxy_df,left_on=[btl_stn_col, btl_prs_col], right_on=[oxy_stn_col, oxy_prs_col],how='outer')

    mask = (df['CTDOXY_y'].notna())
    df.loc[mask,'CTDOXY_x']= df['CTDOXY_y'][mask]

    mask = (df['CTDOXY_FLAG_W_y'].notna())
    df.loc[mask,'CTDOXY_FLAG_W_x']= df['CTDOXY_FLAG_W_y'][mask]

    #Rename Columns

    df.rename(columns={'CTDPRS_x':'CTDPRS','OXYGEN_x':'OXYGEN','CTDOXY_x':'CTDOXY','CTDOXY_FLAG_W_x':'CTDOXY_FLAG_W'},inplace=True)

    drop_list = []
    for key in df.keys():
        if '_x'in key:
            drop_list.append(key)
        elif '_y' in key:
            drop_list.append(key)


    df.drop(columns=drop_list, inplace=True)

    return df

def clean_oxygen_df(df):

    df.rename(columns={'CTDPRS_x':'CTDPRS','OXYGEN_x':'OXYGEN','CTDOXYVOLTS_y':'CTDOXYVOLTS'},inplace=True)

    drop_list = []
    for key in df.keys():
        if '_x'in key:
            drop_list.append(key)
        elif '_y' in key:
                drop_list.append(key)

    df.drop(columns=drop_list, inplace=True)

    return df

def create_coef_df(coef_dict):

    coef_df = pd.DataFrame(coef_dict)
    coef_df = coef_df.transpose()
    coef_df.rename(columns={0:'Soc',1:'Voffset',2:'A',3:'B',4:'C',5:'E',6:'Tau20'})

    return coef_df

def flag_oxy_data(df,oxy_col='OXYGEN',ctd_oxy_col='CTDOXY',p_col ='CTDPRS',flag_col='CTDOXY_FLAG_W'):
    df['res_sbe43'] = df[oxy_col] - df[ctd_oxy_col]
    df_deep = df.loc[(df[p_col] >= 2000) & (df[flag_col]==2)]
    std = df_deep['res_sbe43'].std()
    df.loc[(df[p_col] >= 2000) & (df[flag_col]==2) & (np.abs(df['res_sbe43']) >= (std * 2.8)),flag_col] = 3

    df_shallow = df.loc[(df[p_col] < 2000) & (df[p_col] >= 500) & (df[flag_col]==2)]
    std_shal = df_shallow['res_sbe43'].std()
    df.loc[(df[p_col] < 2000) & (df[p_col] >= 500) & (df[flag_col]==2) & (np.abs(df['res_sbe43']) >= (std_shal * 2.8)),flag_col] = 3

    return df


def least_squares_resid(coef0,oxyvolts,pressure,temp,dvdt,os,ref_oxy,switch,cc=[1.92634e-4,-4.64803e-2]):

    coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef0,
                                     args=(oxyvolts,pressure,temp,dvdt,os,ref_oxy,switch))
    return coef

def oxygen_cal_ml(coef,oxyvolts,pressure,temp,dvdt,os,ref_oxy,switch,cc=[1.92634e-4,-4.64803e-2]):#temp,dvdt,os,

    """"

    NOAA's oxygen fitting routine using the equation:
    OXY(ml/L) = SOC * (doxy_volts+Voffset+Tau20*exp(cc1*PRES+cc2*TEMP)*dvdt) \
                *os*exp(Tcorr*TEMP)*exp(E*PRESS/TEMP_K)

    coef0s:
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E

    cc[0] = D1
    cc[1] = D2

    """
#
#    ctd_oxy_mlL = apply_oxy_coef(time_data,coef,cc=[1.92634e-4,-4.64803e-2],oxyvo_col='CTDOXYVOLTS',p_col='CTDPRS',t_col='CTDTMP1',sal_col='CTDSAL')
    ctd_oxy_mlL = oxy_from_ctd_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc)
#    ctd_oxy_mlL = SB_oxy_eq(coef,oxyvolts,pressure,temp,dvdt,os,cc)
    #Weight Determination
    if switch == 1:

        #weights = calculate_weights(pressure)
        weights = 1/(np.power(pressure,1))
        #resid = ((weights * (ref_oxy - ctd_oxy_mlL))**2) / (np.sum(weights)**2) #Original way (np.sum(weights)**2)
        resid = ((weights * (ref_oxy - ctd_oxy_mlL))**2) #/ #(np.sum(weights)**2) #Original way (np.sum(weights)**2)
    elif switch == 2:
        #L2 Norm
        resid = (ref_oxy - ctd_oxy_mlL)**2

    elif switch == 3:
        #ODF residuals

        resid = np.sqrt(((ref_oxy - ctd_oxy_mlL)**2) / (np.std(ctd_oxy_mlL)**2))

    elif switch == 4:
        # Weighted ODF residuals

        #weights = calculate_weights(pressure)
        weights = 1/(np.power(pressure,1))
        resid = np.sqrt(weights * ((ref_oxy - ctd_oxy_mlL)**2) / (np.sum(weights)**2))#(np.std(ctd_oxy_mlL)**2))

    elif switch == 5:

        #weights = calculate_weights(pressure)
        weights = 1/(np.power(pressure,1))
        resid = ((weights * (ref_oxy - ctd_oxy_mlL))**2) / (np.sum(weights**2))

    return resid

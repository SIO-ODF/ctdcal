#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:04:19 2018

@author: k3jackson
"""


#import matplotlib
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import scipy
import numpy as np
import ctdcal.process_ctd as process_ctd
import ctdcal.sbe_reader as sbe_rd
import ctdcal.sbe_equations_dict as sbe_eq
import gsw
import pandas as pd
#import ctdcal.oxy_fitting as oxy_fitting
import configparser
import os
#import report_ctd
import datetime


#Line 342 module isopycnals

def oxy_fit(time_data, btl_data, ssscc, hexfile, xmlfile, method = 1, 
            sal_col='CTDSAL', t_col='CTDTMP1', p_col='CTDPRS', 
            lon_col='GPSLON', lat_col='GPSLAT', dov_col='CTDOXYVOLTS', 
            cond_col='CTDCOND1', sal_btl_col='CTDSAL', t_btl_col='CTDTMP1',
            p_btl_col='CTDPRS', lat_btl_col='GPSLAT', lon_btl_col='GPSLON',
            dov_btl_col='CTDOXYVOLTS', oxy_btl_col='CTDOXY1',
            cond_btl_col='CTDCOND1'):
    """
    Fits CTD oxygen data to bottle oxygen data and adds to output dataframe.
    
    Default arguments:
        
        time_data -- CTD downcast data (Pandas DataFrame)
        
        btl_data -- upcast bottle data (Pandas DataFrame)
        
        hexfile -- path to the station/cast hex file (string) 
        
        xmlfile -- path to the station/cast xml file (string)
        
    Keyword arguments:
        
        method -- least squares residual calculation method (int):
            1. Weighted L2 Norm. (Default)
            2. Vanilla L2 Norm.
            3. Normalized L2 Norm.
            
        sal_col -- DF column name for CTD downcast salinity in PSU
            (default 'CTDSAL')
        
        t_col -- DF column name for CTD downcast temperature in celcius 
            (default 'CTDTMP1')
            
        p_col -- DF column name for CTD downcast pressure in decibar
            (default 'CTDPRS')
            
        lon_col -- DF column name for CTD downcast longitude (default 'GPSLON')
        
        lat_col -- DF column name for CTD downcast latitude (default 'GPSLAT')
        
        dov_col -- DF column name for CTD downcast dissolved oxygen voltage
            (default 'CTDOXYVOLTS')
        
        cond_col -- DF column name for CTD downcast conductivity in mS/cm
            (default 'CTDCOND1')
        
        sal_btl_col -- DF column name for bottle salinity in PSU
            (default 'CTDSAL')
        
        t_btl_col -- DF column name for bottle temperature in celcius 
            (default 'CTDTMP1')
            
        p_btl_col -- DF column name for bottle pressure in decibar
            (default 'CTDPRS')
            
        lon_btl_col -- DF column name for bottle longitude (default 'GPSLON')
        
        lat_btl_col -- DF column name for bottle latitude (default 'GPSLAT')
        
        dov_btl_col -- DF column name for bottle dissolved oxygen voltage
            (default 'CTDOXYVOLTS')
        
        cond_btl_col -- DF column name for bottle conductivity in mS/cm
            (default 'CTDCOND1')
        
            
    
    """
    
    
    ##NORMALLY INT#####
    stn_nbr = int(ssscc[0:3])
    cst_nbr = int(ssscc[-2:])
    
    #Calculate Practical salinity from conductivity and put into Dataframe
    time_data[sal_col] = gsw.SP_from_C(time_data[cond_col],time_data[t_col],
                                       time_data[p_col])
    
    btl_data[sal_btl_col] = gsw.SP_from_C(btl_data[cond_col],
                                          btl_data[t_btl_col],
                                          btl_data[p_btl_col])
    
    
#    Calculate Conserv Temp, Absolute Sal, and Sigma0 for CTD data and clean up
#    the dataframe
    
    time_data['CT_time'] = gsw.CT_from_t(time_data[sal_col],time_data[t_col],
                                         time_data[p_col])
    
    time_data['SA_time'] = gsw.SA_from_SP(time_data[sal_col],time_data[p_col],
                                          time_data[lon_col],time_data[lat_col])
    
    time_data['sigma0_time'] = gsw.sigma0(time_data['SA_time'],
                                         time_data['CT_time'])
    
    time_data_clean = process_ctd.binning_df(time_data)
    time_data_clean = time_data_clean.dropna(how='all')
    time_data_clean.index = pd.RangeIndex(len(time_data_clean.index))#reindex
    
    #Preserve dataframe for writing out
    time_data_write = time_data_clean
    
#   Rename Dataframe columns to avoid confusion later on
    time_data_clean['CTDOXYVOLTS_time'] = time_data_clean[dov_col] 
    time_data_clean['TMP_KELVIN_CTD'] = time_data_clean[t_col] + 273.15
    time_data_clean['PRESSURE_CTD'] = time_data_clean[p_col]
    time_data_clean['TEMPERATURE_CTD'] = time_data_clean[t_col]
    time_data_clean['SALINITY_CTD'] = time_data_clean[sal_col]
    
#   Reorder by sigma0    
    time_data_clean = time_data_clean.sort_values(by='sigma0_time')
    
#    Calculate Conserv Temp, Absolute Sal, and Sigma0 for BTL data and clean up
#    the dataframe    
    btl_data_clean = btl_data.dropna(how='all')
    btl_data_clean['BTLNBR'] = btl_data['btl_fire_num'] 
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index)) #reindex
    
    btl_data_clean['CT_btl'] = gsw.CT_from_t(btl_data_clean[sal_btl_col],
                                             btl_data_clean[t_btl_col],
                                             btl_data_clean[p_btl_col])
    
    btl_data_clean['SA_btl'] = gsw.SA_from_SP(btl_data_clean[sal_btl_col],
                                              btl_data_clean[p_btl_col],
                                              btl_data_clean[lon_btl_col],
                                              btl_data_clean[lat_btl_col])
    
    btl_data_clean['sigma0_btl'] = gsw.sigma0(btl_data_clean['SA_btl'],
                                              btl_data_clean['CT_btl'])
    
#   Sort bottles by sigma0
    btl_data_clean = btl_data_clean.sort_values(by='sigma0_btl') 
    
#   Calculate dV/dT for bottle data and convolve
    doxyv_btl = np.diff(btl_data_clean['CTDOXYVOLTS'])
    dt_btl = np.diff(btl_data_clean['scan_datetime'])
    dv_dt_btl = doxyv_btl/dt_btl
    dv_dt_conv_btl = np.convolve(dv_dt_btl,[0.5,0.5])
    btl_data_clean['dv_dt_conv_btl'] = dv_dt_conv_btl
    
#   Calculate dV/dT for CTD data, convolve, and filter
    doxyv = np.diff(time_data_clean['CTDOXYVOLTS_time'])
    dt = np.diff(time_data_clean['scan_datetime'])
    
    dv_dt = doxyv/dt
    dv_dt_conv = np.convolve(dv_dt,[0.5,0.5])
    
    a = 1
    windowsize = 5
    b = (1/windowsize)*np.ones(windowsize)
    
    #filt = scipy.signal.lfilter(b,a,dv_dt_conv)
    filtfilt = scipy.signal.filtfilt(b,a,dv_dt_conv)
    
    time_data_clean['dv_dt_time'] = filtfilt
    
#   Create Array containing all BTL and CTD sigma0 values and add increase the
#   range of the CTD sigma0 values in both directions by 1e-4 from the largest
#   and smallest simga0 values
    
    all_sigma0 = time_data_clean['sigma0_time'].append(btl_data_clean['sigma0_btl'])
    fsg = pd.Series(np.min(all_sigma0)-1e-4)
    lsg = pd.Series(np.max(all_sigma0)+1e-4)
    new_sigma0 = fsg.append(time_data_clean['sigma0_time'])
    new_sigma0 = new_sigma0.append(lsg)

#   Increase the size of the CTD press., sal.,temp., and oxyvo., to match the
#   length of new_sigma0
    
    new_pressure = pd.Series(time_data_clean['CTDPRS'].iloc[0])
    new_pressure = new_pressure.append(time_data_clean['CTDPRS'])
    new_pressure = new_pressure.append(pd.Series(time_data_clean['CTDPRS'].iloc[-1]))
    #len(new_pressure)

    new_sal = pd.Series(time_data_clean['CTDSAL'].iloc[0])
    new_sal = new_sal.append(time_data_clean['CTDSAL'])
    new_sal = new_sal.append(pd.Series(time_data_clean['CTDSAL'].iloc[-1]))
    #np.size(new_sal)

    new_temp = pd.Series(time_data_clean[t_col].iloc[0])
    new_temp = new_temp.append(time_data_clean[t_col])
    new_temp = new_temp.append(pd.Series(time_data_clean[t_col].iloc[-1]))
    #np.size(new_temp)

    new_oxyvo = pd.Series(time_data_clean['CTDOXYVOLTS_time'].iloc[0])
    new_oxyvo = new_oxyvo.append(time_data_clean['CTDOXYVOLTS_time'])
    new_oxyvo = new_oxyvo.append(pd.Series(time_data_clean['CTDOXYVOLTS_time'].iloc[-1]))
    #np.size(new_oxyvo)
    
#   Interpolate sigma0 for CTD Data 
    x = np.arange(np.size(time_data_clean['sigma0_time']))
    x_inter = np.arange(np.size(new_sigma0))
    inter_sigma2 = scipy.interpolate.interp1d(x_inter,new_sigma0)
    new_x = np.linspace(0,np.max(x_inter),np.size(time_data_clean['sigma0_time']))
    new_sigma0=inter_sigma2(new_x)
    
#   Interpolate press, sal, temp, and oxyvo values

    inter_pres2 = scipy.interpolate.interp1d(x_inter,new_pressure)
    inter_sal2 = scipy.interpolate.interp1d(x_inter,new_sal)
    inter_temp2 = scipy.interpolate.interp1d(x_inter,new_temp)
    inter_oxyvo2 = scipy.interpolate.interp1d(x_inter,new_oxyvo)
       
    dpres = inter_pres2(new_x)
    dsal = inter_sal2(new_x)
    dtemp = inter_temp2(new_x)
    doxyvo = inter_oxyvo2(new_x)
    
#   Add interpolated values to dataframe

    time_data_clean['sigma0_time'] = new_sigma0
    time_data_clean['PRESSURE_CTD'] = dpres
    time_data_clean['SALINITY_CTD'] = dsal
    time_data_clean['TEMPERATURE_CTD'] = dtemp
    time_data_clean['CTDOXYVOLTS_time'] = doxyvo

#   Interpolate dV/dT down to length of the bottle size
    btl_len_x = np.arange(np.size(btl_data_clean['dv_dt_conv_btl']))
    dv_dt_inter_x = np.arange(np.size(time_data_clean['dv_dt_time']))
    dv_dt_inter = np.interp(btl_len_x,dv_dt_inter_x,time_data_clean['dv_dt_time'])
    
#   Recalculate Simga0
    btl_data_clean = btl_data.dropna(how='all')
    btl_data_clean['BTLNBR'] = btl_data['btl_fire_num']
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index))
   
    btl_data_clean['CT_btl'] = gsw.CT_from_t(btl_data_clean[sal_btl_col],
                                             btl_data_clean[t_btl_col],
                                             btl_data_clean[p_btl_col])
    
    btl_data_clean['SA_btl'] = gsw.SA_from_SP(btl_data_clean[sal_btl_col],
                                              btl_data_clean[p_btl_col],
                                              btl_data_clean[lon_btl_col],
                                              btl_data_clean[lat_btl_col])
    
    btl_data_clean['sigma0_btl'] = gsw.sigma0(btl_data_clean['SA_btl'],
                                              btl_data_clean['CT_btl'])

#   Sort bottles for matching
    btl_data_clean = btl_data_clean.sort_values(by='sigma0_btl') #reindex
    
#   Merge BTL and TIME dataframes on sigma0 and reindex
    matched_df = pd.merge_asof(btl_data_clean,time_data_clean,
                               left_on='sigma0_btl',right_on='sigma0_time',
                               direction='nearest')
    
    btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index))

    time_data_matched = matched_df
    
#   Match CTD dV/dT direction to that of BTL data
    time_data_matched['dv_dt_time'] = dv_dt_inter*-1
    
#   Calculate Oxygen Solubility for both BTL and CTD data
    btl_data_clean['OS'] = sbe_eq.OxSol(btl_data_clean[t_btl_col],
                                        btl_data_clean['CTDSAL'])
    
    time_data_matched['OS'] = sbe_eq.OxSol(time_data_matched['TEMPERATURE_CTD'],
                                           time_data_matched['SALINITY_CTD'])

#   Determine Initial Guess from hex files    
    coef0 = None
    

    sbeReader = sbe_rd.SBEReader.from_paths(hexfile, xmlfile)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
        sensor_id = rawConfig['Sensors'][i]['SensorID']
        if str(sensor_id) == '38':
            oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1, 
                        'ranking': 5, 'column': 'CTDOXYVOLTS', 
                        'sensor_info': rawConfig['Sensors'][i]}

    if coef0 is None:
        coef0 = [oxy_meta['sensor_info']['Soc'], 
                 oxy_meta['sensor_info']['offset'], oxy_meta['sensor_info']['A'], 
                 oxy_meta['sensor_info']['B'], oxy_meta['sensor_info']['C'], 
                 oxy_meta['sensor_info']['E']]
        print(coef0)
        
    Tau20=oxy_meta['sensor_info']['Tau20']
    Tcorr = oxy_meta['sensor_info']['Tcor']
    
    coef0 =[oxy_meta['sensor_info']['Soc'],oxy_meta['sensor_info']['offset'],
            oxy_meta['sensor_info']['Tau20'],oxy_meta['sensor_info']['Tcor'],
            oxy_meta['sensor_info']['E']]
    
    cc=[1.92634e-4,-4.64803e-2]
    
#   Calculate oxygen in ml/L using NOAA's equation
    
    btl_data_clean['dv_dt_conv_btl'] = dv_dt_conv_btl
    
    btl_data_clean['NOAA_oxy_mlL'] = coef0[0]*(btl_data_clean['CTDOXYVOLTS'] \
        + coef0[1] + Tau20 * np.exp(cc[0] * btl_data_clean['CTDPRS'] + cc[0] 
        * btl_data_clean[t_btl_col]) \
        * btl_data_clean['dv_dt_conv_btl']) * btl_data_clean['OS'] \
        * np.exp(Tcorr * btl_data_clean[t_btl_col]) \
        * np.exp((coef0[4]*btl_data_clean['CTDPRS']) / (btl_data_clean[t_btl_col] \
        + 273.15))

    time_data_matched['NOAA_oxy_mlL'] = coef0[0] \
        * (time_data_matched['CTDOXYVOLTS_time'] \
        + coef0[1]+Tau20 * np.exp(cc[0] \
        * time_data_matched['CTDPRS_y'] + cc[0] \
        * time_data_matched['TEMPERATURE_CTD']) \
        * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
        * np.exp(Tcorr * time_data_matched['TEMPERATURE_CTD']) \
        * np.exp((coef0[4] * time_data_matched['CTDPRS_y']) \
        / (time_data_matched['TEMPERATURE_CTD'] + 273.15)) 
        
#   Calculate weights for fitting    
    eps=1e-5
    wrow1 = [0,100,100+eps,300,300+eps,500,500+eps,1200,1200+eps,2000,2000+eps,7000]
    wrow2 = [20,20,25,25,50,50,100,100,200,200,500,500]
    wgt = scipy.interpolate.interp1d(wrow1,wrow2)
    
    time_data_matched['weights'] = wgt(time_data_matched['CTDPRS_y'])
    
#   Least Squares fitting to determine new coefficients   
    coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef0,
                                     args=(time_data_matched,btl_data_clean,
                                           method))

#   Apply new coeficients to recalculate Oxygen in ml/L 
    new_OXY = coef[0] * (time_data_matched['CTDOXYVOLTS_y'] + coef[1] + coef[2] \
        * np.exp(cc[0] * time_data_matched['CTDPRS_y'] + cc[0] \
        * time_data_matched['TEMPERATURE_CTD']) \
        * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
        * np.exp(coef[3] * time_data_matched['TEMPERATURE_CTD']) \
        * np.exp((coef[4] * time_data_matched['CTDPRS_y']) \
        / (time_data_matched['TEMPERATURE_CTD']+273.15))
    
#   Convert Oxygen from ml/L to uMol/kg and put CTD data into BTL dataframe
    time_data_matched['Fitted_OXY_uMOLKG_time'] = new_OXY * 44660 \
        / (time_data_matched['sigma0_time'] + 1000)
    
    btl_data_clean['OXY_uMOLKG_btl'] = btl_data_clean['NOAA_oxy_mlL'] * 44660 \
        / (btl_data_clean['sigma0_btl'] + 1000)
    
    btl_data_clean['CTDOXY'] = time_data_matched['Fitted_OXY_uMOLKG_time']
    btl_data_output = btl_data_clean
    
#   Calculate Residual and determine cutoff based upon standard deviation 
#   outlier values stored in bad_values dataframe
    btl_data_output['residual'] = np.abs(btl_data_clean['OXY_uMOLKG_btl'] 
                            - time_data_matched['Fitted_OXY_uMOLKG_time'])
 
    cutoff = np.std(btl_data_output['residual']) * 2.8
    btl_data_output['index'][btl_data_output['residual']<=cutoff]
    btl_data_clean=btl_data_output[btl_data_output['residual']<=cutoff].copy()
    thrown_values = btl_data_output[btl_data_output['residual']>cutoff]
    bad_values = btl_data_output[btl_data_output['residual']>cutoff]
    btl_data_clean
#    print(thrown_values['BTLNBR'])
    
#   Loop until there are no outliers
    
    while not thrown_values.empty:
        
#       Re-match new BTL dataframe with original time dataframe and reindex
        
        matched_df = pd.merge_asof(btl_data_clean,time_data_clean,
                                   left_on='sigma0_btl',right_on='sigma0_time',
                                   direction='nearest')
        
        btl_data_clean.index = pd.RangeIndex(len(btl_data_clean.index)) #reindex DF
        time_data_matched = matched_df.copy()
        
#       Re-Interpolate dV/dT
        btl_len_x = np.arange(np.size(btl_data_clean['dv_dt_conv_btl']))
        dv_dt_inter_x = np.arange(np.size(time_data_clean['dv_dt_time']))
        
        dv_dt_inter = np.interp(btl_len_x,dv_dt_inter_x,
                                time_data_clean['dv_dt_time'])
        
        time_data_matched['dv_dt_time'] = dv_dt_inter*-1
        
#       Recalculate Oxygen Solubility
        
        btl_data_clean.loc[:,'OS'] = sbe_eq.OxSol(btl_data_clean[t_btl_col].values,
                                            btl_data_clean['CTDSAL'].values)
        
        time_data_matched['OS'] = sbe_eq.OxSol(time_data_matched['TEMPERATURE_CTD'],
                                               time_data_matched['SALINITY_CTD'])
        
#       Recalculate Oxygen (ml/L)       
        time_data_matched['NOAA_oxy_mlL'] = coef[0] \
                * (time_data_matched['CTDOXYVOLTS_time'] \
                + coef[1]+Tau20 * np.exp(cc[0] \
                * time_data_matched['CTDPRS_y'] + cc[0] \
                * time_data_matched['TEMPERATURE_CTD']) \
                * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
                * np.exp(Tcorr * time_data_matched['TEMPERATURE_CTD']) \
                * np.exp((coef[4] * time_data_matched['CTDPRS_y']) \
                / (time_data_matched['TEMPERATURE_CTD'] + 273.15))   
                
#       Recalculate coeficients 
        coef,flag=scipy.optimize.leastsq(oxygen_cal_ml,coef,
                                         args=(time_data_matched,btl_data_clean,
                                               method))
        
#       Re-apply coeficients
        new_OXY = coef[0] * (time_data_matched['CTDOXYVOLTS_y'] + coef[1] + coef[2] \
                  * np.exp(cc[0] * time_data_matched['CTDPRS_y'] + cc[0] \
                  * time_data_matched['TEMPERATURE_CTD']) \
                  * time_data_matched['dv_dt_time']) * time_data_matched['OS'] \
                  * np.exp(coef[3] * time_data_matched['TEMPERATURE_CTD']) \
                  * np.exp((coef[4] * time_data_matched['CTDPRS_y']) \
                  / (time_data_matched['TEMPERATURE_CTD']+273.15))
                  
#       Convert oxygen ml/L to uMol/kg
        time_data_matched['Fitted_OXY_uMOLKG_time'] = new_OXY * 44660 \
                  / (time_data_matched['sigma0_time'] + 1000)
                  
        btl_data_clean['OXY_uMOLKG_btl'] = btl_data_clean['NOAA_oxy_mlL'] \
                  * 44660 / (btl_data_clean['sigma0_btl'] + 1000)
                  
#       Rename Columns
        btl_data_clean['OXYGEN'] = btl_data_clean['OXY_uMOLKG_btl']
        btl_data_clean['CTDOXY'] = time_data_matched['Fitted_OXY_uMOLKG_time']
        
#       Recalculate residuals and residual standard deviation
        btl_data_output = btl_data_clean
        btl_data_clean['residual'] = np.abs(btl_data_output['OXYGEN']-btl_data_output['CTDOXY'])
        std_res = np.std(btl_data_output['residual'])
        cutoff = np.std(btl_data_output['residual'])*2.8
                  
#       Use residuals to determine cutoff and collect outliers
        btl_data_clean=btl_data_clean[btl_data_output['residual']<=cutoff].copy()
        bad_values2 = btl_data_output[btl_data_output['residual']>cutoff]
        thrown_values = btl_data_clean[btl_data_output['residual']>cutoff]
#        print(thrown_values['BTLNBR'])
#        print(btl_data_clean['BTLNBR'])
        bad_values = pd.concat([bad_values,bad_values2])
#        print('NEW STANDARD DEVIATION IS:',std_res)
        
#      Add Station and Cast Numbers to DataFrame
    time_data_matched['STNNBR']=int(stn_nbr)
    time_data_matched['CASTNO']=int(cst_nbr)

    btl_data_clean['STNNBR']=int(stn_nbr)
    btl_data_clean['CASTNO']=int(cst_nbr)
        
##   Sanity PLOT
#    plt.plot(btl_data_clean['residual'],btl_data_clean['CTDPRS']*-1,'bx')
#    plt.xlim(xmax=10)
#    plt.show()
    #dataframe_concat = pd.concat([dataframe_concat,btl_data_clean])
        
#   Flag data    
    bad_values['CTDOXY_FLAG_W'] = 4
    bad_values['STNNBR'] = int(stn_nbr)
    bad_values['CASTNO'] = int(cst_nbr)
    btl_data_clean['CTDOXY_FLAG_W'] = 2
    btl_data_clean = pd.concat([btl_data_clean,bad_values])
    btl_data_clean=btl_data_clean.sort_values(by='BTLNBR')
    
#   Compose writeable dataframe for exportation 
    btl_data_clean = btl_data_clean.reset_index(drop=True)
    btl_data_write = pd.DataFrame()
    btl_data_write['STNNBR'] = btl_data_clean['STNNBR'].values
    btl_data_write['CASTNO'] = btl_data_clean['CASTNO'].values
    btl_data_write['SAMPNO'] = btl_data_clean['BTLNBR'].values
    btl_data_write['CTDOXY'] = btl_data_clean['CTDOXY'].values
    btl_data_write['CTDOXY_FLAG_W'] = btl_data_clean['CTDOXY_FLAG_W'].values
    
    print('Final coefficients: ',coef)
#    btl_data_write.to_csv(filestring,index=False)
#    btl_write_concat = pd.concat([btl_write_concat,btl_data_write])
    
    return btl_data_write, btl_data_clean, coef

def oxygen_cal_ml(coef0,time_data,btl_data,switch):

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

    cc = [1.92634e-4,-4.64803e-2]


    #MODIFIED CODE
    time_data['NOAA_oxy_mlL'] = coef0[0] * (time_data['CTDOXYVOLTS_time'] 
            + coef0[1] + coef0[2] * np.exp(cc[0] * time_data['PRESSURE_CTD'] \
            + cc[0] * time_data['TEMPERATURE_CTD']) * time_data['dv_dt_time']) \
            * time_data['OS'] * np.exp(coef0[3] * time_data['TEMPERATURE_CTD']) \
            * np.exp((coef0[4] * time_data['PRESSURE_CTD']) \
            / (time_data['TEMPERATURE_CTD'] + 273.15))

    #Weight Determination
    if switch == 1:
        
        eps = 1e-5
        
        wrow1 = [0, 100, 100 + eps, 300, 300 + eps, 500, 500 + eps, 1200, 1200 
                 + eps, 2000, 2000 + eps, 7000]
        
        wrow2 = [20, 20, 25, 25, 50, 50, 100, 100, 200, 200, 500, 500]
        
        wgt = scipy.interpolate.interp1d(wrow1,wrow2)

        #Modified CODE
        time_data['weights'] = wgt(time_data['PRESSURE_CTD'])

        resid = (time_data['weights'] * (btl_data['NOAA_oxy_mlL'] 
            - time_data['NOAA_oxy_mlL'])) / (np.sum(time_data['weights'])**2) 
        
    elif switch == 2:
       
        resid =np.sqrt(((btl_data['NOAA_oxy_mlL'] - time_data['NOAA_oxy_mlL'])**2))

                
    
    elif switch == 3:
              
        
        resid = np.sqrt(((btl_data['NOAA_oxy_mlL'] 
                - time_data['NOAA_oxy_mlL'])**2) 
                / (np.std(time_data['NOAA_oxy_mlL'])**2))

    return resid    
    
def apply_oxy_coef(df,coef,oxyvo_col = 'CTDOXYVOLTS',p_col = 'CTDPRS',
                   t_col = 'CTDTMP1',dvdt_col = 'dv_dt_time',sal_col = 'CTDSAL',
                   date_time = 'scan_datetime',lon_col='GPSLON',lat_col='GPSLAT'):
   """ Apply determined coefficients to time data """     
   ###COEF[2] is TAU20
   
   Tau20 = coef[2]
   Tcorr = coef[3]
   cc=[1.92634e-4,-4.64803e-2]
   
   #calculate sigma0
   df['CT'] = gsw.CT_from_t(df[sal_col],df[t_col],
                                         df[p_col])
    
   df['SA'] = gsw.SA_from_SP(df[sal_col],df[p_col],
                                 df[lon_col],df[lat_col])
    
   df['sigma0'] = gsw.sigma0(df['SA'], df['CT'])
   
   
   
   #calculate OS
   df['OS'] = sbe_eq.OxSol(df[t_col], df[sal_col])
   
   #Calculate dv_dt and filter
   doxyv = pd.Series(np.diff(df[oxyvo_col]))
   dt = pd.Series(np.diff(df[date_time]))
   #Remove 0s and replace with means
   dt = dt.replace(0,dt.mean())
   
    
   dv_dt = doxyv/dt
   #df[dvdt_col] = dv_dt
   dv_dt_conv = np.convolve(dv_dt,[0.5,0.5])
    
   a = 1
   windowsize = 5
   b = (1/windowsize)*np.ones(windowsize)
   filtfilt = scipy.signal.filtfilt(b,a,dv_dt_conv)
    
   df[dvdt_col] = filtfilt
   
   #calculate oxy in ml/L using new coef
   df['NOAA_oxy_fitted'] = coef[0] * (df[oxyvo_col] \
                + coef[1]+Tau20 * np.exp(cc[0] \
                * df[p_col] + cc[0] \
                * df[t_col]) \
                * df[dvdt_col]) * df['OS'] \
                * np.exp(Tcorr * df[t_col]) \
                * np.exp((coef[4] * df[p_col]) \
                / (df[t_col] + 273.15))  
    
    #convert to umol/kg
   df['CTDOXY'] = df['NOAA_oxy_fitted'] * 44660 \
                 / (df['sigma0'] + 1000)
  
   return df

def write_oxy_coef(coef,sta_cast):
    """
    
    coefs:
    coef[0] = Soc
    coef[1] = Voffset
    coef[2] = Tau20
    coef[3] = Tcorr
    coef[4] = E

    cc[0] = D1
    cc[1] = D2
    
    
    """
    
    oxy_coef = [{'SSSCC': int(sta_cast), 'Soc': coef[0], 'Voffset': coef[1], 'Tau20': coef[2], 'Tcorr': coef[3], 'E': coef[4]}]
    df = pd.DataFrame(oxy_coef)
    df = df[['SSSCC', 'Soc', 'Voffset', 'Tau20','Tcorr','E']]
    df.set_index('SSSCC',inplace=True)
    
    return df

def get_oxy_coef(ssscc,log_file ='data/logs/oxy_fit_coefs.csv',ind_col = 'SSSCC'):
    
    df = pd.read_csv(log_file,index_col='SSSCC')
    ssscc = int(ssscc)
    coef = np.zeros(5)
    coef[0] = df.loc[ssscc]['Soc']
    coef[1] = df.loc[ssscc]['Voffset']
    coef[2] = df.loc[ssscc]['Tau20']
    coef[3] = df.loc[ssscc]['Tcorr']
    coef[4] = df.loc[ssscc]['E']
    
    
    return coef

def load_ct1_file(ssscc, dir_ctd = 'data/pressure/',ctd_postfix = '_ct1.csv',
                  ctd_skiprows = [0,1,2,3,4,5,6,7,8,9,10,11,13]):
    """ Loads in ct1 pressure File.
    
    This function loads in the SSSCC_ct1.csv file from the pressure directory 
    and converts it to a Pandas DataFrame
    
    Parameters
    ----------
    
    ssscc : str
        The name of the station/cast pressure file to be loaded in. 
        Ex. '00101'
        
    dir_ctd : str
        The name of the directory that the ssscc file is located
        
    ctd_postfix : str
        The the postfix for the ssscc file. Default '_ct1.csv'
        
    ctd_skiprows : array
        Header rows to skip in .csv file
        
    Returns
    -------
    DataFrame
        DataFrame containing the CTD pressure data found in the `ssscc` 
        csv file.
    
    Example
    -------
    >>> df = oxy_fitting.load_ct1_file('00101')
            index  CTDTMP1  CTDTMP2       CTDPRS  CTDCOND1  CTDCOND2     CTDSAL  \
    0         0     17.1067  17.1071     2.933453   45.5290   45.5339  35.456234   
    1         1     17.1068  17.1071     2.933497   45.5291   45.5339  35.456233   
    2         2     17.1068  17.1071     2.935130   45.5291   45.5339  35.456233   
    3         3     17.1069  17.1072     2.938334   45.5292   45.5339  35.456231   
    4         4     17.1070  17.1072     2.942899   45.5292   45.5339  35.456141   
    5         5     17.1070  17.1073     2.948284   45.5293   45.5340  35.456226      
    
    """
    
    df = pd.read_csv(dir_ctd + ssscc + ctd_postfix, skiprows=ctd_skiprows, skipfooter=1, engine='python')
    
    
    
    return df

def load_time_data(time_file):
    """ Loads in CTD time file.
    
    This function returns a pandas DataFrame of the CTD time data (in .pkl format)
    using the path string `time_file`.
    
    Parameters
    ----------
    time_file : str
        String containing the full path including filename with extension of
        the CTD time pickle file you wish to load and convert to a pandas DataFrame
    
    Returns
    -------
    DataFrame
            DataFrame containing the CTD time data found in the `time_file` 
            pickle file.
        
    
    Example
    -------
    
    >>> df = oxy_fitting.load_time_data('data/time/00101_time.pkl')
            index  CTDTMP1  CTDTMP2       CTDPRS  CTDCOND1  CTDCOND2     CTDSAL  \
    0         0     17.1067  17.1071     2.933453   45.5290   45.5339  35.456234   
    1         1     17.1068  17.1071     2.933497   45.5291   45.5339  35.456233   
    2         2     17.1068  17.1071     2.935130   45.5291   45.5339  35.456233   
    3         3     17.1069  17.1072     2.938334   45.5292   45.5339  35.456231   
    4         4     17.1070  17.1072     2.942899   45.5292   45.5339  35.456141   
    5         5     17.1070  17.1073     2.948284   45.5293   45.5340  35.456226  

    
    
    """
    
    time_data = process_ctd.dataToNDarray(time_file,float,True,',',1)
    time_data = pd.DataFrame.from_records(time_data)  
    #time_data = process_ctd.pressure_sequence(time_data)
    
    return time_data

def apply_time_data(ssscc,p_col = 'CTDPRS',t_col = 'CTDTMP1',dvdt_col = 'dv_dt_time',
                    sal_col = 'CTDSAL',date_time = 'scan_datetime',lon_col='GPSLON',
                    lat_col='GPSLAT',oxyvo_col = 'CTDOXYVOLTS'):
    
    time_file = 'data/time/'+ssscc+'_time.pkl'
    time_data = load_time_data(time_file)
    coef = get_oxy_coef(ssscc)
    
    time_data['CT'] = gsw.CT_from_t(time_data[sal_col],time_data[t_col],time_data[p_col])
    time_data['SA'] = gsw.SA_from_SP(time_data[sal_col],time_data[p_col],
                              time_data[lon_col],time_data[lat_col])
    time_data['sigma0'] = gsw.sigma0(time_data['SA'], time_data['CT'])
    time_data['OS'] = sbe_eq.OxSol(time_data[t_col], time_data[sal_col])
    
    time_data = pressure_sequence(time_data)
    
    time_data = apply_oxy_coef(time_data,coef)
    # Fill NA for dv_dt at beginning of the cast with mean
    time_data['dv_dt_time'] = time_data['dv_dt_time'].replace([np.inf, -np.inf], np.nan)
    time_data['dv_dt_time'] = time_data['dv_dt_time'].fillna(time_data['dv_dt_time'].mean())
    
    
    Tau20 = coef[2]
    Tcorr = coef[3]
    cc=[1.92634e-4,-4.64803e-2]
    
    time_data['NOAA_oxy_fitted'] = coef[0] * (time_data[oxyvo_col] \
             + coef[1]+Tau20 * np.exp(cc[0] \
             * time_data[p_col] + cc[0] \
             * time_data[t_col]) \
             * time_data[dvdt_col]) * time_data['OS'] \
             * np.exp(Tcorr * time_data[t_col]) \
             * np.exp((coef[4] * time_data[p_col]) \
             / (time_data[t_col] + 273.15))

    time_data['CTDOXY'] = time_data['NOAA_oxy_fitted'] * 44660 \
              / (time_data['sigma0'] + 1000)
    
    
    
    return time_data['CTDOXY']
#def oxyfit_to_ctd(ssscc,df):
    
def write_ct1_file_oxy(ssscc, ctd_oxy, dir_ctd = 'data/pressure/',ctd_postfix = '_ct1.csv'):
    
    
    
    ### Gather Cruise/Station information
    FILE_EXT = 'csv'
    filename_base = ssscc
    
    depth = -999
    
    iniFile = 'data/ini-files/configuration.ini'
    config = configparser.RawConfigParser()
    config.read(iniFile)
   
    #Initialise Configuration Parameters
    expocode = config['cruise']['expocode']
    sectionID = config['cruise']['sectionid']
    pressure_directory = config['ctd_processing']['pressure_data_directory']
    log_directory = config['ctd_processing']['log_directory']

    dopkg_col = config['analytical_inputs']['dopkg']
    ctd = config['ctd_processing']['ctd_serial']
    p_column_one = list(config['pressure_series_output']['q1_columns'].split(','))
    p_column_data = config['pressure_series_output']['data'].split(',')
    p_column_names = config['pressure_series_output']['column_name'].split(',')
    p_column_units = config['pressure_series_output']['column_units'].split(',')
    p_column_qual = config['pressure_series_output']['qual_columns'].split(',')
    
        # Collect Cast Details from Log
    logfileName = str('cast_details' + '.' + FILE_EXT)
    logfilePath = os.path.join(log_directory, logfileName)

    cast_details = process_ctd.dataToNDarray(logfilePath,str,None,',',0)
    for line in cast_details:
        if filename_base in line[0]:
            for val in line:
                if 'at_depth' in val: btime = float(str.split(val, ':')[1])
                if 'latitude' in val: btm_lat = float(str.split(val, ':')[1])
                if 'longitude' in val: btm_lon = float(str.split(val, ':')[1])
                if 'altimeter_bottom' in val: btm_alt = float(str.split(val, ':')[1])
            break
   

    ct1_df = load_ct1_file(ssscc)
    
    qual_pseq_data = ct1_df['CTDPRS_FLAG_W']
    
    ct1_df['CTDOXY'] = ctd_oxy
    ct1_df['CTDOXY_FLAG_W'] = 2
       
    
    report_pressure_series_data_oxy(filename_base, expocode, sectionID, btime, btm_lat, btm_lon, depth, btm_alt, ctd, pressure_directory, p_column_names, p_column_units, p_column_data, qual_pseq_data, ct1_df, ct1_df)

    return

def report_pressure_series_data_oxy(stacast, expocode, section_id, btime=-999, btm_lat=-999, btm_lon=-999, depth=-999, btm_alt=-999, ctd=-999, p_dir=None, p_column_names=None, p_column_units=None, p_column_data=None, qualMat=None, doArr=None, df=None):
    """report_pressure_series_data function

    Function takes full NUMPY ndarray with predefined dtype array
    and writes time series data with quality codes to csv file.

    Args:
        param1 (str): stacast, actually file base name from hex data.
        param2 (str): expocode, vessel id and startd date cruise identification code
        param3 (str): section_id, US hydrographic section id
        param4 (str): btime, UTC date time stamp at bottom of cast
        param5 (str): btm_lat, latitude value at bottom of cast.
        param6 (str): btm_lon, longitude value at bottom of cast.
        param7 (str): depth, ctd depth value + altimeter value.
        param8 (str): btm_alt, altimeter value at bottom of cast.
        param9 (str): ctd, serial number of ctd inst.
        param10 (str): p_dir, directory to write file.
        param11 (str): p_column_names, header column names.
        param12 (str): p_column_units, header column units.
        param13 (dataframe): qualMat, input dataframe.
        param13 (ndarray): doArr, input do data.
        param14 (ndarray): inMat, input data ndarray.

    Prints formatted csv file.

    Returns:
        No return
    """

    if df is None:
        print("In report_pressure_series_data: No data array.")
        return
    else:
        out_col = []
        h_num = 11
        now = datetime.datetime.now()
        file_datetime = now.strftime("%Y%m%d %H:%M")
        bdt = datetime.datetime.fromtimestamp(btime).strftime('%Y%m%d %H%M').split(" ")
        b_date = bdt[0]
        b_time = bdt[1]

        s_num = stacast[-5:-2]
        c_num = stacast[-2:]

        outfile = open(p_dir+stacast+'_ct1.csv', "w+")
        outfile.write("CTD, %s\nNUMBER_HEADERS = %s \nEXPOCODE = %s \nSECT_ID = %s\nSTNNBR = %s\nCASTNO = %s\n DATE = %s\nTIME = %s\nLATITUDE = %f\nLONGITUDE = %f\nDEPTH = %s\nINSTRUMENT_ID = %s\n" % (file_datetime, h_num, expocode, section_id, s_num, c_num, b_date, b_time, btm_lat, btm_lon, depth, ctd))
        cn = np.asarray(p_column_names)
        cn.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')
        cu = np.asarray(p_column_units)
        cu.tofile(outfile,sep=',', format='%s')
        outfile.write('\n')

        # Need to rewrite to remove hardcoded column data.
        # No ideal method to print formatted output to csv in python ATM
        for i in range(0,len(df)):
            ### Alter configuration.ini in pressure_series in order to change output
            ### Will need to be altered in order to put out RINKO in UMOL/KG
            ### CTDOXY is slipped in and piggybacks on CTDXMISS qual code
            outfile.write("%8.1f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d,%10.4f,%d\n" % (df['CTDPRS'][i], df['CTDPRS_FLAG_W'][i], df['CTDTMP'][i], df['CTDTMP_FLAG_W'][i], df['CTDSAL'][i], df['CTDSAL_FLAG_W'][i], df['CTDOXY'][i], df['CTDOXY_FLAG_W'][i], df['CTDXMISS'][i], df['CTDXMISS_FLAG_W'][i], df['CTDFLUOR'][i], df['CTDFLUOR_FLAG_W'][i],df['CTDBACKSCATTER'][i], df['CTDBACKSCATTER_FLAG_W'][i], df['CTDRINKO'][i], df['CTDRINKO_FLAG_W'][i]))
        outfile.write('END_DATA')
    return
    

###### TAKEN FROM PROCESS_CTD ##############
def _roll_filter(df, pressure_column="CTDPRS", direction="down"):
    #fix/remove try/except once serialization is fixed
    try:
        if direction == 'down':
            monotonic_sequence = df[pressure_column].expanding().max()
        elif direction == 'up':
            monotonic_sequence = df[pressure_column].expanding().min()
        else:
            raise ValueError("direction must be one of (up, down)")
    except KeyError:
        pressure_column = 'CTDPRS'
        if direction == 'down':
            monotonic_sequence = df[pressure_column].expanding().max()
        elif direction == 'up':
            monotonic_sequence = df[pressure_column].expanding().min()
        else:
            raise ValueError("direction must be one of (up, down)")

    return df[df[pressure_column] == monotonic_sequence]


def roll_filter(df, p_col='CTDPRS', up='down', frames_per_sec=24, search_time=15, start=0):
    """roll_filter function

    Function takes full NUMPY ndarray with predefined dtype array
    and subsample arguments to return a roll filtered ndarray.

    
    Jackson Notes: calculate sampling frequency 
        Two types of filters: differential filter, alias filt
        Get index of second diff function and do discrete analysis

    Args:
        param1 (str): stacast, station cast info
        param2 (ndarray): inMat, numpy ndarray with dtype array
        param3 (str): up, direction to filter cast (up vs down)
        param4 (int): frames_per_sec, subsample selection rate
        param5 (int): seach_time, search time past pressure inversion

    Returns:
        Narray: The return value ndarray of data with ship roll removed
    """
    #When the "pressure sequence" code is fixed, uncomment and use this instead
    #start = kwargs.get("start", 0)
    if df is None:
        print("Roll filter function: No input data.")
        return
    
    end = df[p_col].idxmax()
    full_matrix = df
    tmp_df = df[start:end]
    tmp_df = _roll_filter(tmp_df)

    return tmp_df

#    remove = []
#    frequency = 24 # Hz of package
#
#    if (frames_per_sec > 0) & (frames_per_sec <= 24):
#        sample = int(frequency/frames_per_sec) # establish subsample rate to time ratio
#    else: sample = frequency
#
#    # Adjusted search time with subsample rate
#    search_time = int(sample*frequency*int(search_time))


#    else:
#        P = df[p_col]
#        dP = P.diff()
#        
#        if up is 'down':
#            #index_to_remove = np.where(dP < 0)[0] # Differential filter Use DIff command
#            #subMat = np.delete(df, index_to_remove, axis=0)# use dataframe boolean
#            
#            P = P[(dP>0) | (dP.isna()==True)]#Remove pressure value increases and recover first element with or
#            P = P.reset_index(drop=True)
#            dP2 = P.diff()
#            P2 = P[(dP2<0)]# Questionable data points
#            indicies = P2.index
#
#            tmp = np.array([])
#            for i in range(0,len(P)-1):#Lose If Statement
#               if P[i] > P[i+1]:# Use another diff command to find indicies

#                   deltaP = P[i+1] + abs(P[i] - P[i+1])
#                   # Remove aliasing
#                   k = np.where(P == min(P[i+1:i+search_time], key=lambda x:abs(x-deltaP)))[0]
#                   tmp = np.arange(i+1,k[0]+1,1)
#               remove = np.append(remove,tmp)
#               deltaP = 0
#        elif up is 'up':
#            index_to_remove = np.where(dP > 0)[0] # Differential filter
#            subMat = np.delete(inMat, index_to_remove, axis=0)
#
#            P = subMat[p_col]
#            tmp = np.array([])
#            for i in range(0,len(P)-1):
#               if P[i] < P[i+1]:
#                   deltaP = P[i+1] - abs(P[i] - P[i+1])
#                   # Remove aliasing
#                   k = np.where(P == min(P[i+1:i+search_time], key=lambda x:abs(x-deltaP)))[0]
#                   tmp = np.arange(i+1,k[0]+1,1)
#               remove = np.append(remove,tmp)
#               deltaP = 0
#
#        subMat = np.delete(subMat,remove,axis=0)
#
#    return subMat

def pressure_sequence(df, p_col='CTDPRS', intP=2.0, startT=-1.0, startP=0.0, up='down', sample_rate=12, search_time=15):
    """pressure_sequence function

    Function takes a dataframe and several arguments to return a pressure 
    sequenced data ndarray.

    Pressure sequencing includes rollfilter.

    Necessary inputs are input Matrix (inMat) and pressure interval (intP).
    The other inputs have default settings. The program will figure out
    specifics for those settings if left blank.
    Start time (startT), start pressure (startP) and up are mutually exclusive.
    If sensors are not not fully functional when ctd starts down cast
    analyst can select a later start time or start pressure but not both.
    There is no interpolation to the surface for other sensor values.
    'up' indicates direction for pressure sequence. If up is set startT and startP
    are void.

    Args:
        param1 (Dataframe: Dataframe containing measurement data
        param2 (str): p_col, pressure column name
        param3 (float): starting pressure interval
        param5 (float): start time (startT) for pressure sequence
        param6 (float): start pressure (startP) for pressure sequence
        param7 (str): pressure sequence direction (down/up)
        param8 (int): sample_rate, sub sample rate for roll_filter. Cleans & speeds processing.
        param9 (int): search_time, truncate search index for the aliasing part of ship roll.
        param10 (ndarray): inMat, input data ndarray

    Returns:
        Narray: The return value is a matrix of pressure sequenced data

    todo: deep data bin interpolation to manage empty slices
    """
    # change to take dataframe with the following properties
    # * in water data only (no need to find cast start/end)
    # * The full down and up time series (not already split since this method will do it)
    # New "algorithm" (TODO spell this right)
    # * if direction is "down", use the input as is
    # * if direction is "up", invert the row order of the input dataframe
    # Use the "roll filter" method to get only the rows to be binned
    # * the roll filter will treat the "up" part of the cast as a giant roll to be filtered out
    # * the reversed dataframe will ensure we get the "up" or "down" part of the cast
    # * there is no need to reverse the dataframe again as the pressure binning process will remove any "order" information (it doesn't care about the order)
    # That's basically all I (barna) have so far TODO Binning, etc...
    # pandas.cut() to do binning

    #lenP, prvPrs not used
    # Passed Time-Series, Create Pressure Series

    start = 0

    # Roll Filter
    roll_filter_matrix = roll_filter(df, p_col, up, sample_rate, search_time, start=start)

    df_roll_surface = fill_surface_data(roll_filter_matrix, bin_size=2)
    #bin_size should be moved into config
    binned_df = binning_df(df_roll_surface, bin_size=2)
    binned_df = binned_df.reset_index(drop=True)
    return binned_df
### Once serialization has been fixed, fix try/except to compact code

def binning_df(df, **kwargs):
    '''Bins records according to bin_size, then finds the mean of each bin and returns a df.
    '''
    bin_size = kwargs.get("bin_size", 2)
    try:
        labels_in = [x for x in range(0,int(np.ceil(df['CTDPRS_DBAR'].max())),2)]
        df['bins'] = pd.cut(df['CTDPRS_DBAR'], range(0,int(np.ceil(df['CTDPRS_DBAR'].max()))+bin_size,bin_size), right=False, include_lowest=True, labels=labels_in)
        df['CTDPRS_DBAR'] = df['bins'].astype('float64')
        df_out = df.groupby('bins').mean()
        return df_out
    except KeyError:
        labels_in = [x for x in range(0,int(np.ceil(df['CTDPRS'].max())),2)]
        df['bins'] = pd.cut(df['CTDPRS'], range(0,int(np.ceil(df['CTDPRS'].max()))+bin_size,bin_size), right=False, include_lowest=True, labels=labels_in)
        df['CTDPRS'] = df['bins'].astype('float64')
        df_out = df.groupby('bins').mean()
        return df_out

def fill_surface_data(df, **kwargs):
    '''Copy first scan from top of cast, and propgate up to surface
    '''
    surface_values = []
    bin_size = kwargs.get("bin_size", 2)
    try:
        for x in range(1, int(np.floor(df.iloc[0]['CTDPRS_DBAR'])), bin_size):
            surface_values.append(x)
        df_surface = pd.DataFrame({'CTDPRS_DBAR': surface_values})
        df_merged = pd.merge(df_surface, df, on='CTDPRS_DBAR', how='outer')
    except KeyError:
        for x in range(1, int(np.floor(df.iloc[0]['CTDPRS'])), bin_size):
            surface_values.append(x)
        df_surface = pd.DataFrame({'CTDPRS': surface_values})
        df_merged = pd.merge(df_surface, df, on='CTDPRS', how='outer')

    return df_merged.fillna(method='bfill')

#Load time_data to df
    
    #####           GRAVEYARD             ##########
        
    #NOAA Calucaltion (unmodified):
        #ORIGINAL CODE
        #time_data['NOAA_oxy_mlL'] =coef0[0]*(time_data['CTDOXYVOLTS_y']+\
        #coef0[1]+coef0[2]*np.exp(cc[0]*time_data['CTDPRS_y']+cc[0]*time_data['TEMPERATURE_CTD'])\
        #*time_data['dv_dt_time'])*time_data['OS']*np.exp(coef0[3]*time_data['TEMPERATURE_CTD'])\
        #*np.exp((coef0[4]*time_data['CTDPRS_y'])/(time_data['TEMPERATURE_CTD']+273.15))
    
    #SWITCH 1 code:
    
        #Original Code
        #time_data['weights'] = wgt(time_data['CTDPRS_y'])
        #resid = np.sum((time_data['weights']*(btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL']))/(np.sum(time_data['weights'])**2))
        #resid = (time_data['weights']*(btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL']))/(np.sum(time_data['weights'])**2) Working
        
    #SWITCH 3 code:
    
        #Original Code
        #resid = np.sqrt(((btl_data['CTDOXY1']-time_data['NOAA_oxy_mlL'])**2)/(np.std(time_data['NOAA_oxy_mlL'])**2))
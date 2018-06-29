#!/usr/bin/env python
import math
import scipy
import numpy as np
import pandas as pd
#import ctdcal.process_ctd as process_ctd
import ctdcal.sbe_reader as sbe_rd
import ctdcal.sbe_equations_dict as sbe_eq
from scipy.optimize import leastsq
import gsw
import csv
#import requests
import os
import json

M = 31.9988   #Molecular weight of O2
R = 831.432    #/* Gas constant, X 10 J Kmole-1 K-1 */
D = 1.42905481 #  /* density O2 g/l @ 0 C */


def offset(offset, inArr):
    """offset column of data

    Input:
        - inMat, 1d numpy array with return True np.isnans()
    Output:
        - Mat with offset applied to column of data
    Example:
        >>> outArray = offset(offset, col, inMat)
    """
    #for i in range(0, len(inArr)):
    #    inArr[i] = float(inArr[i]) + offset

    #return inArr
    return inArr + offset


def IESRho(df,sal_col='CTDSAL',t_col='CTDTMP1',p_col='CTDPRS'):
    """
        Calculate the IES density
    
    Parameters
    ----------
    df : Dataframe
         dataframe containing salinity, temperature and pressure data

    sal_col : string
             dataframe column name of salinity data [PSU]
    
    t_col : string
            dataframe colum,n name of temperature data [C (ITS-90)]
    
    p_col : string
            dataframe column name of pressure data [dBar]

    Returns
    -------
    df : Dataframe
         input dataframe but with an addition 'rho' column

    
    
    """
    
    
    bars  = df[p_col] * 0.1
#/* pressure in bars */
#/*
#**  Calculate Rho(s,t,0.0)
#*/
    
    rhow = (999.842594 + df[t_col] *( 0.06793952 + df[t_col] * (-0.00909529 
            + df[t_col] * (1.001685e-4 + df[t_col] * (-1.120083e-6 + df[t_col]
            * 6.536332e-9)))))
    
    #/* rhow = density of pure water kg/m**3 */
    
    kw = (((df[t_col] * (-0.0040899+ df[t_col] *( 7.6438e-5+ df[t_col]
         * (-8.2467e-7 + df[t_col] *  5.3875e-9)))) + 0.824493) * df[sal_col])
    #/* pure water secant bulk modulus */
    
    termc = df[sal_col] * np.sqrt(df[sal_col]);
    kst0  = ((-0.00572466 + df[t_col] * ( 1.0227e-4  + df[t_col] 
             * (-1.6546e-6))) * termc)
    #/* k(s,t,0) */
    rho   = rhow + kw + kst0 + 4.8314e-4 * df[sal_col] * df[sal_col]
#      /* rho(s,t,0.0)  kg/m**3 */
#      /*
#      **  Calculate pressure effects.
#      */
#  	/*
#  	**                 Rho(s,t,0)
#  	**  Rho(s,t,p) = -------------
#  	**                        p
#  	**               1.0 - -------
#  	**                     k(s,t,p)
#  	*/
    kw    = (df[t_col] * (148.4206 + df[t_col] * ( -2.327105+ df[t_col] 
            * (0.01360477 + df[t_col] * ( -5.155288e-5)))) + 1.965221e4)
    
    kst0  =  ((54.6746 + df[t_col] * (-0.603459 + df[t_col] * (0.0109987 
               + df[t_col] * (-6.167e-5)))) * df[sal_col] + kw + (0.07944 
              + df[t_col] * (0.016483 + df[t_col] * (-5.3009e-4))) * termc)
#  	/*
#  	**  Calculate pressure terms.
#  	*/
    terma = (3.239908 + df[t_col] * (0.00143713 + df[t_col] * (1.16092e-4
             + df[t_col] * (-5.77905e-7))) + (0.0022838 + df[t_col]
            * (-1.0981e-5 + df[t_col] * (-1.6078e-6))) * df[sal_col]
            + 1.91075e-4 * termc)
    
    termb = (8.50935e-5  + df[t_col] * (-6.12293e-6  + df[t_col] *  5.2787e-8)
             + (-9.9348e-7 + df[t_col] * (2.0816e-8 + df[t_col] * 9.1697e-10))
             * df[sal_col])
    
    kstp  = kst0 + bars * (terma + bars * termb)
    df['rho']   = rho/(1.0-bars/kstp)
    return df

def get_flasks(o2flasks):
    with open(o2flasks, 'r') as f:
        flasks = {}
        for l in f:
            if 'Volume' not in l:
                if l.strip().startswith("#"):
                    continue
                row = l.strip().split()
                flasks[int(row[0])] = float(row[1])
    return flasks

def get_flask_vol(flask, o2flasks, t=20, glass="borosilicate"):
    _coef = {
            "borosilicate": 0.00001,
            "soft": 0.000025,
            }
    _t = 20
    flasks = get_flasks(o2flasks)
    fv = flasks[flask]
    coef = _coef[glass]
    return fv * (1.0 + coef * (t - _t));

#def stp_rho(rho_func=IESRho):
#    return rho_func(0, 20, 0)
#
def rho_t(t):
    z0       =  9.9983952e2
    z1       =  1.6945176e1
    z2       = -7.9870401e-3
    z3       = -4.6170461e-5
    z4       =  1.0556302e-7
    z5       = -2.8054253e-10
    z6       =  1.6879850e-2
    Rho20Deg =     998.2041
    return ((z0+(t)*(z1+(t)*((z2+(t)*(z3+(t)*(z4+(t)*z5))))))/(1.0+z6*(t)))

def thio_n_calc(titr, blank, kio3_n, kio3_v, kio3_t, thio_t):

    rho_stp  = rho_t(20)
    rho_kio  = rho_t(kio3_t)
    rho_thio = rho_t(thio_t)

    kio3_v_20c = kio3_v * (rho_kio/rho_stp)
    thio_v_20c = titr * (rho_thio/rho_stp) - blank

    thio_n = kio3_v_20c * kio3_n / thio_v_20c
    return thio_n

def titr_20_calc(titr, titr_temp,):
    rho_titr = rho_t(titr_temp)
    rho_stp = rho_t(20)

    #convert the titr ammount to 20c equivalent
    titr_20c = titr * rho_titr/rho_stp
    return titr_20c

def mll_to_umolkg(o2ml, s, t, rho_func=IESRho):
    o2kg = o2ml / ((M/D * 0.001) * rho_func(s, t, 0)/1000)
    return o2kg

# Find nearest value to argument in array
# Return the index of that value
def find_isopycnals(p_btl_col, t_btl_col, sal_btl_col, dov_btl_col, lat_btl_col, lon_btl_col, btl_data, p_col, t_col, sal_col, dov_col, lat_col, lon_col, time_data):
    """find_iscopycnals

    p_btl_col:   Pressure column for bottle data
    t_btl_col:   Temperature column for bottle data
    sal_btl_col: Salinity column for bottle data
    dov_btl_col: Oxygen voltage column for bottle data
    lat_btl_col: Latitude bottle column for bottle data
    lon_btl_col: Longitude bottle column for bottle data
    btl_data:    Bottle data ndarray
    p_col:       Pressure column for bottle data
    t_col:       Temperature column for bottle data
    sal_col:     Salinity column for bottle data
    dov_col:     Oxygen voltage column for bottle data
    lat_col:     Latitude column for bottle data
    lon_col:     Longitude column for bottle data
    time_data:   Time data ndarray

    """

    time_sigma = []
    CT = gsw.CT_from_t(time_data[sal_col],time_data[t_col],time_data[p_col])
    SA = gsw.SA_from_SP(time_data[sal_col],time_data[p_col],time_data[lon_col],time_data[lat_col])
    time_sigma = gsw.sigma0(SA,CT)

    # Pressure-bounded isopycnal search.
    # Based on maximum decent rate and package size.
    CT = gsw.CT_from_t(btl_data[sal_btl_col],btl_data[t_btl_col],btl_data[p_btl_col])
    SA = gsw.SA_from_SP(btl_data[sal_btl_col],btl_data[p_btl_col],btl_data[lon_btl_col],btl_data[lat_btl_col])
    btl_sigma = gsw.sigma0(SA,CT)
    for i in range(0,len(btl_data[p_btl_col])):
        #CT = gsw.CT_from_t(btl_data[sal_btl_col][i],btl_data[t_btl_col][i],btl_data[p_btl_col][i])
        #SA = gsw.SA_from_SP(btl_data[sal_btl_col][i],btl_data[p_btl_col][i],btl_data[lon_btl_col][i],btl_data[lat_btl_col][i])
        #btl_sigma = gsw.sigma0(SA,CT)
        p_indx = find_nearest(time_data[p_col], btl_data[p_btl_col][i])
        indx = find_nearest(time_sigma[p_indx-int(24*1.5):p_indx+int(24*1.5)], btl_sigma[i])

        #print('Bottle:')
        #print('Sigma: '+str(btl_sigma))
        #print('Pres: '+str(btl_data[p_btl_col][i]))
        #print('Temp: '+str(btl_data[t_btl_col][i]))
        #print('Salt: '+str(btl_data[sal_btl_col][i]))
        #print('Pressure: '+str(p_indx)+' '+str(indx+p_indx))
        #print('Sigma: '+str(time_sigma[indx+p_indx]))
        #print('Pres: '+str(time_data[p_col][indx+p_indx]))
        #print('Temp: '+str(time_data[t_col][indx+p_indx]))
        #print('Salt: '+str(time_data[sal_col][indx+p_indx]))

        if indx+p_indx > len(time_sigma):
            btl_data[t_btl_col][i] = time_data[t_col][len(time_data)-1]
            btl_data[sal_btl_col][i] = time_data[sal_col][len(time_data)-1]
            btl_data[dov_btl_col][i] = time_data[dov_col][len(time_data)-1]
        else:
            btl_data[t_btl_col][i] = time_data[t_col][indx+p_indx]
            btl_data[sal_btl_col][i] = time_data[sal_col][indx+p_indx]
            btl_data[dov_btl_col][i] = time_data[dov_col][indx+p_indx]

    return btl_data


# Find nearest value to argument in array
# Return the index of that value
def find_nearest(yarr, val):
    """find_nearest

    """
    indx = (np.abs(yarr-val)).argmin()
    return indx


# Residual calculation
def calibration(independent_arr, dependent_diff_arr, order):
    """calibration

    """
    return np.polyfit(independent_arr, dependent_diff_arr, order)


# Oxygen Coef
def find_oxy_coef(o2pl, p, t, salt, dov, hexfilePath, xmlfilePath):
    """fit_oxy fits CTD dissolved oxygen

    """
    kelvin = []
    for i in range(0,len(t)):
        kelvin.append(t[i] + 273.15)

    # Retrieve Config data
    sbeReader = sbe_rd.SBEReader.from_paths(hexfilePath, xmlfilePath)
    rawConfig = sbeReader.parsed_config()
    for i, x in enumerate(rawConfig['Sensors']):
       sensor_id = rawConfig['Sensors'][i]['SensorID']
       if str(sensor_id) == '38':
           oxy_meta = {'sensor_id': '38', 'list_id': 0, 'channel_pos': 1, 'ranking': 5, 'column': 'CTDOXYVOLTS', 'sensor_info': rawConfig['Sensors'][i]}


    coef0 = [oxy_meta['sensor_info']['Soc'], oxy_meta['sensor_info']['offset'], oxy_meta['sensor_info']['A'], oxy_meta['sensor_info']['B'], oxy_meta['sensor_info']['C'], oxy_meta['sensor_info']['E']]
    oxy_data = oxy_dict(coef0, p, kelvin, t, salt, dov)
    # Non Linear fit routine
    coefs, flag = leastsq(residualO2, coef0, args=(o2pl.astype(float),p,kelvin,t,salt,dov))

    return coefs

def oxy_dict(calib, P, K, T, S, V):
    """SBE equation for converting engineering units to oxygen (ml/l).
    SensorID: 38

    calib is a dict holding Soc, Voffset, Tau20, A, B, C, E
    The following are single or list/tuple:
    P is pressure in decibars
    K is temperature in Kelvin
    T is temperature in Celcius
    S is Practical Salinity Units
    V is Voltage from instrument

    Original equation from calib sheet dated 2014:
    Oxygen (ml/l) = Soc * (V + Voffset) * (1.0 + A * T + B * T + C * T ) * OxSol(T,S) * exp(E * P / K)

    """
    #array mode
    try:
       oxygen = []
       for P_x, K_x, T_x, S_x, V_x in zip(P, K, T, S, V):
           #print(T_x)
           temp = (calib[0] * (V_x + calib[1])
                   * (1.0 + calib[2] * T_x + calib[3] * math.pow(T_x,2) + calib[4] * math.pow(T_x,3) )
                   * sbe_eq.OxSol(T_x,S_x)
                   * math.exp(calib[5] * P_x / K_x)) #foo
           temp = round(temp,4)
           oxygen.append(temp)
    #Single mode.
    except:
       oxygen = (calib[0] * (V + calib[1])
                 * (1.0 + calib[2] * T + calib[3] * math.pow(T,2) + calib[4] * math.pow(T,3) )
                 * sbe_eq.OxSol(T,S)
                 * math.exp(calib[5] * P / K))
    return oxygen
    # return (calib[0] * (V + calib[1])
    #           * (1.0 + calib[2] * T + calib[3] * math.pow(T,2) + calib[4] * math.pow(T,3) )
    #           * sbe_eq.OxSol(T,S)
    #           * np.exp(calib[5] * P / K))


# Residual calculation
def residualO2(calib, o2pl, P, K, T, S, V):
    """residual weighted difference of dissolved oxygen bottle data
       vs dissolved oxygen CTD data.

    This conversion is included for least squares fitting routine.

    calib is a dict holding Soc, Voffset, Tau20, A, B, C, E
    The following are single or list/tuple:
    calib is a list of oxy_dict coefficients to be optimized
    o2pl is dissolved oxygen winkler titrated data
    P is pressure in decibars
    K is temperature in Kelvin
    T is temperature in Celcius
    S is Practical Salinity Units
    V is Voltage from instrument
    """
    weight = []
    ctd_o2pl = oxy_dict(calib, P, K, T, S, V)
    sig = np.std(ctd_o2pl)

    # Least sq residual
    for i in range(0, len(o2pl)):
        if o2pl[i] > 0:
            weight.append(scipy.sqrt((o2pl[i] - oxy_dict(calib, P[i], K[i], T[i], S[i], V[i]))**2/sig**2))
    return weight


#def conductivity_polyfit(C, P, T, cond):
#    """Polynomial used to fit conductivity data with pressure effect.
#    The following are single or list/tuple:
#    C is starting estimate for coefficients
#    P is pressure in decibars
#    T is temperature in Celcius
#    cond is conductivity in mS/cm
#
#    Original equation from ...
#    Conductivity mS/cm = cond + C0 * P^2 + C1 * P + C2 * T^2 + C3 * T + C4 * cond^2 + C5 * cond + C6
#
#    Another time based fit must be run at the end of cruise to account for time dependent drift.
#
#    """
#    try:
#        c_arr = []
#        for P_x, T_x, cond_x in zip(P, T, cond):
#            tmp = cond_x + C[0] * np.power(P_x,2) + C[1] * P_x + C[2] * np.power(T_x,2) + C[3] * T_x + C[4] * np.power(cond_x,2) + C[5] * cond_x + C[6]
#            c_arr.append(round(tmp,4))
#    #Single mode.
#    except:
#        tmp = cond + C[0] * np.power(P,2) + C[1] * P + C[2] * np.power(T,2) + C[3] * T + C[4] * np.power(cond_x,2) + C[5] * cond_x + C[6]
#        c_arr = round(tmp,4)
#    #tmp = cond + C[0] * math.pow(P,2) + C[1] * P + C[2] * math.pow(T,2) + C[3] * T + C[4] * math.pow(cond_x,2) + C[5] * cond_x + C[6]
#    #c_arr = tmp.round(decimals=4)
#    #c_arr = round(tmp,4)
#
#    return c_arr

def conductivity_polyfit(df,coef,t_col,cond_col,p_col='CTDPRS'):
#    
     fitted_cond = df[cond_col] + (coef[0] * (df[p_col]**2) + coef[1] * df[p_col] + coef[2] * (df[t_col]**2) \
                    + coef[3] * df[t_col] + coef[4] * (df[cond_col]**2) + coef[5] * df[cond_col] + coef[6])
     fitted_cond = fitted_cond.round(4)
     fitted_sal =  gsw.SP_from_C(fitted_cond,df[t_col],df[p_col])
     return fitted_cond, fitted_sal
     
     
#def temperature_polyfit(C, P, T):
#    """Polynomial used to fit data with pressure effect.
#
#    The following are single or list/tuple:
#    fit declares fit type
#    C is starting estimate for coefficients
#    P is pressure in decibars
#    T is temperature in Celcius
#
#    Original equation from ...
#    Temperature degC ITS-90 = T + C0 * P^2 + C1 * P + C2 * T^2 + C3 * T + C4
#
#    Another time based fit must be run at the end of cruise to account for time dependent drift.
#
#    """
#    try:
#        t_arr = []
#        for P_x, T_x in zip(P, T):
#            tmp = T_x + C[0] * np.power(P_x,2) + C[1] * P_x + C[2] * np.power(T_x,2) + C[3] * T_x + C[4]
#            t_arr.append(round(tmp,4))
#       #Single mode.
#    except:
#        tmp = T + C[0] * np.power(P,2) + C[1] * P + C[2] * np.power(T,2) + C[3] * T + C[4]
#        t_arr = round(tmp,4)
#
#    return t_arr

def temperature_polyfit(df,coef,t_col,p_col='CTDPRS'):
    
    fitted_temp = df[t_col] + coef[0] * (df[p_col]**2) + coef[1] * df[p_col] + coef[2] * (df[t_col]**2) + coef[3] * df[t_col] + coef[4]
    fitted_temp = fitted_temp.round(4)
    
    return fitted_temp

#def load_qual(path):
#    comment_dict = {}
#    with open(path) as f:
#        reader = csv.reader(f)
#        # ignore first line
#        next(reader)
#        for line in reader:
#            sta = int(line[0])
#            cast = int(line[1])
#            bottle = int(line[2])
#            param = line[3]
#            flag = int(line[4])
#            comment = line[5]
#
#            comment_dict[(sta, cast, bottle, param)] = [flag, comment]
#
#    return comment_dict
#
#
#
#salts = requests.get("http://go-ship.rrevelle.sio.ucsd.edu/api/salt").json()
#def o2_calc(path, o2_payload, thio_ns):

def o2_calc(o2flasks, o2path, btl_num): #, salt
#    qual = load_qual("/Volumes/public/O2Backup/o2_codes_001-083.csv")

    btl_num.astype(int)
    o2ml = np.zeros(shape=(len(btl_num),), dtype=[('BTLNUM', np.int),('OXYGEN',np.float)])
    o2kg = np.zeros(shape=(len(btl_num),), dtype=[('BTLNUM', np.int),('OXYGEN',np.float)])

    with open(o2path, 'r') as f:
#        rho = IESRho
        params = next(f).strip().split()

        titr   = float(params[0])
        blank  = float(params[1])
        kio3_n = float(params[2])
        kio3_v = float(params[3])
        kio3_t = float(params[4])
        thio_t = float(params[5])

        thio_n = thio_n_calc(titr, blank, kio3_n, kio3_v, kio3_t, thio_t)
#        rho_stp = rho_t(20)

        btl_counter = 0
        #try:
        for l in f:
            row = l.split()
            if "ABORT" in row[-1]:
                continue
            station = int(row[0])
            cast = int(row[1])
            bottle = int(row[2])
            if (bottle == 99) or (bottle > 36):
                continue
            flask = int(row[3])
            titr = float(row[4])
            titr_temp = float(row[5])
            draw_temp = float(row[6])
            flask_vol = get_flask_vol(flask, o2flasks, draw_temp)
            titr_20c = titr_20_calc(titr, titr_temp)

            if bottle in btl_num:
                #print(btl_num[bottle])
                btl_counter += 1
                o2ml['BTLNUM'][bottle-1] = int(bottle)
                o2ml['OXYGEN'][bottle-1] = (((titr_20c - blank) * thio_n * 5.598 - 0.0017)/((flask_vol - 2.0) * 0.001))
                o2kg['BTLNUM'][bottle-1] = int(bottle)
                #o2kg['OXYGEN'][bottle-1] = mll_to_umolkg(o2ml['OXYGEN'][bottle-1], salt[bottle-1], draw_temp,rho)
            else:
                btl_counter += 1
                o2ml['BTLNUM'][bottle-1] = btl_counter
                o2ml['OXYGEN'][bottle-1] = '-999'
                o2kg['OXYGEN'][bottle-1] = '-999'
                o2kg['BTLNUM'][bottle-1] = btl_counter
#           row_dict = {"station": str(station), "cast": str(cast),"bottle": str(bottle), "o2": o2kg}
        #except ValueError:
            #print('File probably malformed. Check datafile for problems.')
    return o2ml #o2kg,


def salt_calc(saltpath, btl_num_col, btl_tmp_col, btl_p_col, btl_data):
    
    f = open(saltpath, newline='')
    saltF = csv.reader(f,delimiter=' ', quoting=csv.QUOTE_NONE, skipinitialspace='True')
    
    saltArray = []
    for row in saltF:
        saltArray.append(row)
    del saltArray[0]
         
    header = ['STNNBR','CASTNO','SAMPNO','BathTEMP','CRavg','autosalSAMPNO',\
              'Unknown','StartTime','EndTime','Attempts','Reading1','Reading2',\
              'Reading3', 'Reading4', 'Reading5','Reading6','Reading7','Reading8',\
              'Reading9', 'Reading10','Reading11','Reading12']
    f.close()
    # make all rows of Salt files the same length as header   
    for row in saltArray:
        if len(row) < len(header):
            row.extend([np.NaN]*(len(header)-len(row)))
            
    saltArray = np.array(saltArray) # change to np array
    
    saltDF = pd.DataFrame(saltArray,columns=header) # change to DataFrame
    saltDF = saltDF.apply(pd.to_numeric, errors='ignore')
    
    cond = saltDF[['autosalSAMPNO','SAMPNO','CRavg','BathTEMP']]
    # Remove standard measurements
    cond = cond[(cond['autosalSAMPNO']!='worm')]
    #if all(cond['autosalSAMPNO'].values != cond['SAMPNO'].values):
    #    raise ValueError('Problem with sample numbers in salt file (check file)')
    cond = cond.drop('autosalSAMPNO',axis=1)
    cond = cond.apply(pd.to_numeric) # For some reason doesn't completely work the first time
    #cond = cond[(cond['SAMPNO']!=0) & (cond['SAMPNO']!=99)]
    # Filter unmeansured bottle data from btl_data
    data = btl_data[btl_data[btl_num_col].isin(cond['SAMPNO'].tolist())]
    
    salinity = SP_salinometer((cond['CRavg']/2.0),cond['BathTEMP'])
    try:
        cond['BTLCOND'] = gsw.C_from_SP(salinity,data[btl_tmp_col],data[btl_p_col])
        #cond = cond.drop('SAMPNO',1)
    except ValueError:
        print('Possible mis-entered information in salt file (Check salt file)')
        raise ValueError
        #data = btl_data[btl_data[btl_num_col].isin(cond['SAMPNO'].tolist())]
        #salinity = SP_salinometer((cond['CRavg']/2.0),cond['BathTEMP'])
        #cond['BTLCOND'] = gsw.C_from_SP(salinity,data[btl_tmp_col],data[btl_p_col])
        # Drop and rename autosalSAMPNO
        #cond = cond.drop('autosalSAMPNO',1)
        #cond = cond.rename(columns={'SAMPNO':'autosalSAMPNO'})
    # Create 36-place DF
    cond
    DF = pd.DataFrame(data=np.arange(1,37),columns=['SAMPNO'],index=range(1,37))
    # Merge
    DF = DF.merge(cond,on="SAMPNO",how='outer')
    DF = DF.set_index(np.arange(1,37))
    
    
    return DF
    
def CR_to_cond(cond_ratio,bath_temp,btl_temp,btl_press):

    ### Clean up to avoid runtimewarning ###
    cond_df = cond_ratio.copy()
    bath_df = bath_temp.copy()
    nans = np.isnan(cond_df)
    bnans = np.isnan(bath_df)
    cond_df.loc[nans] = 0
    bath_df.loc[bnans] = 0
    
    
    salinity = SP_salinometer((cond_df / 2.0),bath_df)
    cond = gsw.C_from_SP(salinity,btl_temp,btl_press)  
    
    cond[cond<=1] = np.nan
    
    #cond = pd.series(cond)
    
    
    return cond
    
##    qual = load_qual("/Volumes/public/O2Backup/o2_codes_001-083.csv")
#
#    salt_file_name = os.path.basename(saltpath)
#    salt_sta = int(salt_file_name[0:3])
#    salt_cst = int(salt_file_name[3:5])
#
#    btl_num = btl_data[btl_num_col].astype(int)
#
#    psu = np.zeros(shape=(len(btl_num),), dtype=[(btl_num_col, np.int),('SALNTY',np.float)])
#    mspcm = np.zeros(shape=(len(btl_num),), dtype=[(btl_num_col, np.int),('BTLCOND',np.float)])
#    tmp_tmp = np.zeros(len(btl_num), dtype=np.float)
#    tmp_p = np.zeros(len(btl_num), dtype=np.float)
#    bath_tmp = np.zeros(len(btl_num), dtype=np.float) #patched in - JIG 2017-07-15
#
#    with open(saltpath, 'r') as f:
#        params = next(f).strip().split()
#        #std   = float(params[8])
#
#        params = next(f).strip().split()
#        worm1   = float(params[4])
#
#        #btl_counter = 0
#        for l in f:
#            row = l.split()
#            station = int(row[0])
#            cast = int(row[1])
#            if (station == salt_sta) and (cast == salt_cst):
#                if (row[5] == 'worm'):
#                    worm2 = float(row[4])
#                else:
#                   bottle = int(row[5])
#                   cond = float(row[4])
#                   bath = int(row[3]) #patched in - JIG 2017-07-15
#                   bath_tmp = np.full_like(bath_tmp, bath)
#
#                if bottle in btl_num:
#                    i = int(np.where(btl_num == bottle)[0][0])
#                    j = int(np.where(btl_data[btl_num_col] == bottle)[0][0])
#                    tmp_tmp[j] = btl_data[btl_tmp_col][j]
#                    tmp_p[j] = btl_data[btl_p_col][j]
#                    psu[btl_num_col][i] = int(bottle)
#                    mspcm[btl_num_col][i] = int(bottle)
#                    psu['SALNTY'][i] = cond / 2.0
#                    #import pdb; pdb.set_trace()
#
#
#        psu['SALNTY'] = SP_salinometer(psu['SALNTY'], bath_tmp)
#        mspcm['BTLCOND'] = gsw.C_from_SP(psu['SALNTY'], tmp_tmp, tmp_p)
#        
#        mspcm = pd.DataFrame(mspcm)
#        psu = pd.DataFrame(psu)
##           row_dict = {"station": str(station), "cast": str(cast),"bottle": str(bottle), "o2": o2kg}
    
    
    
#    return mspcm, psu

def write_calib_coef(ssscc,coef,param):
    """ Write coef to csv
    
    
    """
    df = pd.DataFrame()
    df['SSSCC'] = ssscc  
    
    if param == 'T':

        df['coef_0'] = coef[0]
        df['coef_1'] = coef[1]
        df['coef_2'] = coef[2]
        df['coef_3'] = coef[3]
        df['coef_4'] = coef[4]
        df['coef_5'] = coef[5]
        df['coef_6'] = coef[6]
        
    if param == 'C':
       
        df['coef_0'] = coef[0]
        df['coef_1'] = coef[1]
        df['coef_2'] = coef[2]
        df['coef_3'] = coef[3]
        df['coef_4'] = coef[4]
        

    return df

#
def apply_fit_coef(df,ssscc,coef_frame,param,sensor,t_col = 'CTDTMP',p_col = 'CTDPRS',
                   cond_col = 'CTDCOND'):
    """ Applies Coef to time and bottle Data
    
    
    """
    
    coef = coef_frame.loc[coef_frame['SSSCC']== ssscc]
    
    if sensor == 1:
        t_col = t_col+'1'
        cond_col = cond_col+'1'
        
    elif sensor == 2:
        t_col = t_col+'2'
        cond_col = cond_col+'2'
           
    
    if param == 'T':
        
        df[t_col] = temperature_polyfit(coef,df[p_col],df[t_col])
        
    
    elif param == 'C':
    
        df[cond_col] = conductivity_polyfit(coef,df[p_col],df[t_col],
                                              df[cond_col])
    
    
    
    
    return df

def apply_pressure_offset(df,p_off,p_col='CTDPRS'):
    
    df[p_col] = offset(p_off, df[p_col])
#    df[p_col] = np.around(df[p_col],decimals=3)
    return df
  
##
#            key = (station, cast, bottle, "o2")
#            if key in qual:
#                flag, comment = qual[key]
#                row_dict["o2_qual"] = flag
#                row_dict["comment"] = comment
#
#            o2_payload.append(row_dict)
#
#
#            # stoichiometric relation between mols of thio and mols of o2
#            #print(station, cast, bottle, o2ml, o2kg)
#            print("{:>2}, {:8.1f}".format(bottle, o2kg))
#
#        thio_ns.append([station, thio_n])
#
#o2_payload = []
#thio_ns = []
#for root, dirs, files in os.walk("/Volumes/public/O2Backup/O2/"):
#    for file in files:
#        if file.startswith("9") or file.startswith("8") or file.startswith("."):
#            continue
#        path = os.path.join(root, file)
#        o2_calc(path, o2_payload, thio_ns)
#    break
#
#with open("o2kg.json", "w") as f:
#    json.dump(o2_payload,f, indent=2)
#with open("thios.json", "w") as f:
#    json.dump(thio_ns, f)

### Taken from python-gsw, aka 48-term version. Needs to be rewritten and added
### to GSW-python or this module must be adjusted to deal with generic forms.
def SP_salinometer(Rt, t):
    """
    Calculates Practical Salinity SP from a salinometer, primarily using
    the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical
    Salinity is only valid in the range 2 < SP < 42.  If the PSS-78 algorithm
    produces a Practical Salinity that is less than 2 then the Practical
    Salinity is recalculated with a modified form of the Hill et al. (1986)
    formula. The modification of the Hill et al. (1986) expression is to
    ensure that it is exactly consistent with PSS-78 at SP = 2.
    A laboratory salinometer has the ratio of conductivities, Rt, as an output,
    and the present function uses this conductivity ratio and the temperature t
    of the salinometer bath as the two input variables.
    Parameters
    ----------
    Rt : array
         C(SP,t_68,0)/C(SP=35,t_68,0) [unitless]
         conductivity ratio
         :math:`R = \frac{C(S, t_68, 0)}{C(35, 15(IPTS-68),0)} [unitless]
    t : array
        Temperature of the bath of the salinometer [:math:`^\circ` C (ITS-90)]
    Returns
    -------
    SP : array
         Practical Salinity [psu (PSS-78), unitless]
    Examples
    --------
    TODO
    References
    -----------
    .. [1] Fofonoff, P. and R.C. Millard Jr. 1983: Algorithms for computation
       of fundamental properties of seawater.  Unesco Tech. Pap. in Mar. Sci.,
       44, 53 pp.
    .. [2] Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the
       Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng.,
       11, 109 - 112.
    .. [3] IOC, SCOR and IAPSO, 2010: The international thermodynamic equation
       of seawater - 2010: Calculation and use of thermodynamic properties.
       Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
       UNESCO (English), 196 pp. See appendix E of this TEOS-10 Manual, and in
       particular, Eqns. (E.2.1) and (E.2.6).
    """

    a = (0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081)
    b = (0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144)
    c = (0.6766097, 2.00564e-2, 1.104259e-4, -6.9698e-7, 1.0031e-9)
    d = (3.426e-2, 4.464e-4, 4.215e-1, -3.107e-3)
    e = (2.070e-5, -6.370e-10, 3.989e-15)
    P = (4.577801212923119e-3, 1.924049429136640e-1, 2.183871685127932e-5,
         -7.292156330457999e-3, 1.568129536470258e-4, -1.478995271680869e-6,
         9.086442524716395e-4, -1.949560839540487e-5, -3.223058111118377e-6,
         1.175871639741131e-7, -7.522895856600089e-5, -2.254458513439107e-6,
         6.179992190192848e-7, 1.005054226996868e-8, -1.923745566122602e-9,
         2.259550611212616e-6, 1.631749165091437e-7, -5.931857989915256e-9,
         -4.693392029005252e-9, 2.571854839274148e-10, 4.198786822861038e-12)
    q = (5.540896868127855e-5, 2.015419291097848e-1, -1.445310045430192e-5,
         -1.567047628411722e-2, 2.464756294660119e-4, -2.575458304732166e-7,
         5.071449842454419e-3, -9.081985795339206e-5, -3.635420818812898e-6,
         2.249490528450555e-8, -1.143810377431888e-3, 2.066112484281530e-5,
         7.482907137737503e-7, 4.019321577844724e-8, -5.755568141370501e-10,
         1.120748754429459e-4, -2.420274029674485e-6, -4.774829347564670e-8,
         -4.279037686797859e-9, -2.045829202713288e-10, 5.025109163112005e-12)
    r = (3.432285006604888e-3, 1.672940491817403e-1, 2.640304401023995e-5,
         1.082267090441036e-1, -6.296778883666940e-5, -4.542775152303671e-7,
         -1.859711038699727e-1, 7.659006320303959e-4, -4.794661268817618e-7,
         8.093368602891911e-9, 1.001140606840692e-1, -1.038712945546608e-3,
         -6.227915160991074e-6, 2.798564479737090e-8, -1.343623657549961e-10,
         1.024345179842964e-2, 4.981135430579384e-4, 4.466087528793912e-6,
         1.960872795577774e-8, -2.723159418888634e-10, 1.122200786423241e-12)
    u = (5.180529787390576e-3, 1.052097167201052e-3, 3.666193708310848e-5,
         7.112223828976632, -3.631366777096209e-4, -7.336295318742821e-7,
         -1.576886793288888e+2, -1.840239113483083e-3, 8.624279120240952e-6,
         1.233529799729501e-8, 1.826482800939545e+3, 1.633903983457674e-1,
         -9.201096427222349e-5, -9.187900959754842e-8, -1.442010369809705e-10,
         -8.542357182595853e+3, -1.408635241899082, 1.660164829963661e-4,
         6.797409608973845e-7, 3.345074990451475e-10, 8.285687652694768e-13)
    k = 0.0162

    a, b, c, d, e, P, q, r, u, k = (np.asarray(x)
                                    for x in (a, b, c, d, e, P, q, r, u, k))

    Rt, t = np.broadcast_arrays(Rt, t, subok=True)

    t68 = t * 1.00024
    ft68 = (t68 - 15) / (1 + k * (t68 - 15))

    Rt[Rt < 0] = np.ma.masked
    Rtx = np.sqrt(Rt)

    SP = (a[0] + (a[1] + (a[2] + (a[3] + (a[4] + a[5] * Rtx) * Rtx) * Rtx) *
                  Rtx) * Rtx + ft68 *
          (b[0] + (b[1] + (b[2] + (b[3] + (b[4] + b[5] * Rtx) * Rtx) * Rtx) *
                   Rtx) * Rtx))

    """The following section of the code is designed for SP < 2 based on the
    Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
    exactly equal to the PSS-78 algorithm at SP = 2."""

    I2 = SP < 2
    if I2.any():
        Hill_ratio = Hill_ratio_at_SP2(t[I2])
        x = 400 * Rt[I2]
        sqrty = 10 * Rtx[I2]
        part1 = 1 + x * (1.5 + x)
        part2 = 1 + sqrty * (1 + sqrty * (1 + sqrty))
        SP_Hill_raw = SP[I2] - a[0] / part1 - b[0] * ft68[I2] / part2
        SP[I2] = Hill_ratio * SP_Hill_raw
    # Ensure that SP is non-negative.
    SP = np.maximum(SP, 0)
    return SP

def Hill_ratio_at_SP2(t):
    """
    USAGE:
     Hill_ratio = Hill_ratio_at_SP2(t)
    DESCRIPTION:
     Calculates the Hill ratio, which is the adjustment needed to apply for
     Practical Salinities smaller than 2.  This ratio is defined at a
     Practical Salinity = 2 and in-situ temperature, t using PSS-78. The Hill
     ratio is the ratio of 2 to the output of the Hill et al. (1986) formula
     for Practical Salinity at the conductivity ratio, Rt, at which Practical
     Salinity on the PSS-78 scale is exactly 2.
    INPUT:
     t  =  in-situ temperature (ITS-90)                  [ deg C ]
    OUTPUT:
     Hill_ratio  =  Hill ratio at SP of 2                [ unitless ]
    AUTHOR:
     Trevor McDougall and Paul Barker
    VERSION NUMBER: 3.0 (26th March, 2011)
    """

    SP2 = 2 * np.ones_like(t)

    a0 = 0.0080
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 = 2.7081
    b0 = 0.0005
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 = 0.0636
    b5 = -0.0144
    g0 = 2.641463563366498e-1
    g1 = 2.007883247811176e-4
    g2 = -4.107694432853053e-6
    g3 = 8.401670882091225e-8
    g4 = -1.711392021989210e-9
    g5 = 3.374193893377380e-11
    g6 = -5.923731174730784e-13
    g7 = 8.057771569962299e-15
    g8 = -7.054313817447962e-17
    g9 = 2.859992717347235e-19
    k = 0.0162
    t68 = t * 1.00024
    ft68 = (t68 - 15) / (1 + k * (t68 - 15))
    # -------------------------------------------------------------------------
    # Find the initial estimates of Rtx (Rtx0) and of the derivative dSP_dRtx
    # at SP = 2.
    # -------------------------------------------------------------------------
    Rtx0 = g0 + t68 * (g1 + t68 * (g2 + t68 * (g3 + t68 * (g4 + t68 * (g5
              + t68 * (g6 + t68 * (g7 + t68 * (g8 + t68 * g9))))))))
    dSP_dRtx = (a1 + (2 * a2 + (3 * a3 + (4 * a4 + 5 * a5 * Rtx0) * Rtx0) *
                Rtx0) * Rtx0 + ft68 * (b1 + (2 * b2 + (3 * b3 + (4 * b4 + 5 *
                b5 * Rtx0) * Rtx0) * Rtx0) * Rtx0))
    # -------------------------------------------------------------------------
    # Begin a single modified Newton-Raphson iteration to find Rt at SP = 2.
    # -------------------------------------------------------------------------
    SP_est = (a0 + (a1 + (a2 + (a3 + (a4 + a5 * Rtx0) * Rtx0) * Rtx0) * Rtx0) *
              Rtx0 + ft68 * (b0 + (b1 + (b2 + (b3 + (b4 + b5 * Rtx0) * Rtx0) *
              Rtx0) * Rtx0) * Rtx0))
    Rtx = Rtx0 - (SP_est - SP2) / dSP_dRtx
    Rtxm = 0.5 * (Rtx + Rtx0)
    dSP_dRtx = (a1 + (2 * a2 + (3 * a3 + (4 * a4 + 5 * a5 * Rtxm) * Rtxm) *
                Rtxm) * Rtxm + ft68 * (b1 + (2 * b2 + (3 * b3 + (4 * b4 + 5 *
                b5 * Rtxm) * Rtxm) * Rtxm) * Rtxm))
    Rtx = Rtx0 - (SP_est - SP2) / dSP_dRtx
    # This is the end of one full iteration of the modified Newton-Raphson
    # iterative equation solver.  The error in Rtx at this point is equivalent
    # to an error in SP of 9e-16 psu.
    x = 400 * Rtx * Rtx
    sqrty = 10 * Rtx
    part1 = 1 + x * (1.5 + x)
    part2 = 1 + sqrty * (1 + sqrty * (1 + sqrty))
    SP_Hill_raw_at_SP2 = SP2 - a0 / part1 - b0 * ft68 / part2
    return 2. / SP_Hill_raw_at_SP2

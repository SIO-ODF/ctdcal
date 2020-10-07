"""A module for SBE conversion equations and related helper equations.

Eventual goal is to convert all outputs to numpy arrays to make compatible with
gsw libraries, and remove written wrappers.

SBE internally references each sensor with an assigned SensorID for indexing. The
given sensor numbers are determined empirically from .XMLCON files in 2016-2017 and
not via an official document, and may change according to SBE wishes.

"""

'''code_pruning: this module really needs a new name, and all functions inside it as well.'''
'''code_pruning: Also all functions should be checked to work with vectorized versions and possibly remove single datapoint version.'''
'''code_pruning: calls to the math library/module should be removed in favor of faster numpy (does that work with xarray/dask, or does code need to be rewritten for that?)'''

import math

import numpy as np


def temp_its90_dict(calib, freq, verbose = 0):
    """SBE equation for converting engineering units to Celcius according to ITS-90.

    TO BE DEPRECIATED
    SensorID: 55

    Parameters
    ----------
    calib : dict
        calib is a dict holding G, H, I, J, F0
        G, H, I, J, F0: coefficients given in calibration.
    freq : float
        freq = frequency sampled by sensor, either as a single value or a list or tuple.

    Returns
    -------
    ITS90: arrylike or float
        A float depicting ITS90 temperature.

    Original form from calib sheet dated 2012:
    Temperature ITS-90 = 1/{g+h[ln(f0/f )]+i[ln2(f0/f)]+j[ln3(f0/f)]} - 273.15 (°C)

    """
    f = freq
    #array mode
    try:
        ITS90 = []
        for i, f_x in enumerate(f):
            #Hack for dealing with 0 frequencies, needs additional reporting later
            if f_x == 0:
                f_x = 1
                if verbose > 0:
                    print("Zero (0) frequency temperature record being processed as 1. Record: ", i)
            temp = 1/(calib['G']
                      + calib['H'] * (math.log(calib['F0']/f_x))
                      + calib['I'] * math.pow((math.log(calib['F0']/f_x)),2)
                      + calib['J'] * math.pow((math.log(calib['F0']/f_x)),3)
                     ) - 273.15
            temp = round(temp,4)
            ITS90.append(temp)
    #single mode
    except:
        if f == 0:
            f = 1
            if verbose > 0:
                print("Zero (0) frequency temperature record [singleton] being processed.")
        ITS90 = 1/(calib['G']
                   + calib['H'] * (math.log(calib['F0']/f))
                   + calib['I'] * math.pow((math.log(calib['F0']/f)),2)
                   + calib['J'] * math.pow((math.log(calib['F0']/f)),3)
                  ) - 273.15
        ITS90 = round(ITS90,4)
    return ITS90

def temp_its90(calib, f):
    """SBE equation for converting engineering units to Celcius according to ITS-90.

    SensorID: 55

    Parameters
    ----------
    calib : dict
        calib is a dict holding G, H, I, J, F0
        G, H, I, J, F0: coefficients given in calibration.
    f : float
        f = frequency sampled by sensor, either as a single value or a list or tuple.

    Returns
    -------
    ITS90 : arraylike or float
        A float depicting ITS90 temperature.

    Original form from calib sheet dated 2012:
    Temperature ITS-90 = 1/{g+h[ln(f0/f )]+i[ln2(f0/f)]+j[ln3(f0/f)]} - 273.15 (°C)

    """
    # The array might come in as type=object, which throws AttributeError. Maybe this should be in try/except?
    f = f.astype(float)

    ITS90 = np.around(
            (1/(calib['G']
              + calib['H'] * (np.log(calib['F0']/f))
              + calib['I'] * np.power((np.log(calib['F0']/f)),2)
              + calib['J'] * np.power((np.log(calib['F0']/f)),3)
            ) - 273.15
            ),4)
#     Was used when using raw python and not numpy. Need to figure out how to handle 0s now, it defaults
#     to not producing a value where 0 is used, ending in a -273.15 value being output
    return ITS90

def sbe9(calib, f, t):
    """NOT TESTED

    SBE/STS(?) equation for converting SBE9 pressure frequency to pressure.
    SensorID: 45

    Parameters
    ----------
    calib : dict
        T1: coefficient
        T2: coefficient
        T3: coefficient
        T4: coefficient
        T5: not used
        C1: coefficient
        C2: coefficient
        C3: coefficient
        D1: coefficient
        D2: coefficient
        AD590M: used in digiquartz temeperature correction
        AD590B: used in digiquartz temperature correction

    f : array-like
        f is sensor frequency (usually between 30kHz and 42kHz)
    t : array-like
        t is sensor integer from the digiquartz temperature probe

    Returns
    -------
    p : array-like
        p is pressure in decibar
    """
    #Equation expecting pressure period in microseconds, so divide f by 1,000,000.
    uf = f / 1000000
    t = (calib['AD590M'] * int(t)) + calib['AD590B']
    T0 = calib['T1'] + calib['T2']*t + calib['T3']*np.power(t,2) + calib['T4']*np.power(t,3)
    w = 1-T0*T0*f*f
    pressure = (0.6894759*((calib['C1']+calib['C2']*t+calib['C3']*t*t)*w*(1-(calib['D1']+calib['D2']*t)*w)-14.7))
    return np.around(pressure,4)

def sbe4c(calib, F, t, p):
    """NOT TESTED

    SBE equation for converting SBE4C frequency to conductivity. Calculates mS/cm

    SensorID: 3

    Parameters
    ----------
    calib : dict
        G : coefficient
        H : coefficient
        I : coefficient
        J : coefficient
        CPcor : coefficient (nominal)
        CTcor : coefficient (nominal)
    F : array-like
        F is instrument frequency
    t : array-like
        t is temperature (ITS-90 degrees C)
    p : array-like
        p is pressure (decibars)

    Returns
    -------
    c : array-like
        c is conductivity in mS/cm
    """
    f = F/1000
    c = ((calib['G'] + calib['H'] * np.power(f,2)
                     + calib['I'] * np.power(f,3)
                     + calib['J'] * np.power(f,4))
                    / (1 + calib['CTcor'] * t + calib['CPcor'] * p))
    #S/m to mS/cm
    c = c * 0.1
    c = round(c,5)
    return c

def sbe43(calib, P, K, T, S, V):
    """NOT TESTED
    IS THIS IN ITS-90? DOES IT MATTER?
    SBE equation for converting SBE43 engineering units to oxygen (ml/l).

    SensorID: 38

    Parameters
    ----------

    calib : dictionary
        calib is a dict holding Soc, Voffset, Tau20, A, B, C, E
    P : array-like
        P is pressure in decibars
    K : array-like
        K is temperature in Kelvin
    T : array-like
        T is temperature in Celcius
    S : array-like
        S is Practical Salinity Units, PSS-78
    V : array-like
        V is voltage from instrument

    Returns
    -------
    oxygen : array-like
        oxygen is the oxygen measurement in ml/l

    Original equation from calib sheet dated 2014:
    Oxygen (ml/l) = Soc * (V + Voffset) * (1.0 + A * T + B * T + C * T ) * OxSol(T,S) * exp(E * P / K)
    """

    oxygen = np.around(calib['Soc'] * (V + calib['offset'])
              * (1.0 + calib['A'] * T + calib['B'] * np.power(T,2) + calib['C'] * np.power(T,3) )
              * OxSol(T,S)
              * np.exp(calib['E'] * P / K),4)
    return oxygen

def OxSol_New(T,S):
    """Should work with sbe43(), should be depreciated in favor of gsw in the future

    Eq. 8 from Garcia and Gordon, 1992.
    Harded coded to do ml/l for now. If requeste, add in new mode for ug/l.

    Parameters
    ----------
    T : array-like
        T is the ITS-90 temperature in Celcius
    S : array-like
        S is the Practical Salinity, PSS-78
    """

    x = S
    y = np.log((298.15 - T)/(273.15 + T))

    """umol/kg coefficients
    a0 =  5.80871
    a1 =  3.20291
    a2 =  4.17887
    a3 =  5.10006
    a4 = -9.86643e-2
    a5 =  3.80369
    b0 = -7.01577e-3
    b1 = -7.70028e-3
    b2 = -1.13864e-2
    b3 = -9.51519e-3
    c0 = -2.75915e-7
    """

    """ml/l coefficients"""
    a0 = 2.00907
    a1 = 3.22014
    a2 = 4.05010
    a3 = 4.94457
    a4 = -2.56847e-1
    a5 = 3.88767
    b0 = -6.24523e-3
    b1 = -7.37614e-3
    b2 = -1.03410e-2
    b3 = -8.17083e-3
    c0 = -4.88682e-7

    O2sol = np.exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
    return O2sol


def OxSol(T,S):
    """Eq. 8 from Garcia and Gordon, 1992.
    Harded coded to do ml/l for now. If requeste, add in new mode for ug/l.

    Parameters
    ----------
    T : Temperature
    S : Practical Salinity
    """

    x = S
    y = np.log((298.15 - T)/(273.15 + T))

    """umol/kg coefficients
    a0 =  5.80871
    a1 =  3.20291
    a2 =  4.17887
    a3 =  5.10006
    a4 = -9.86643e-2
    a5 =  3.80369
    b0 = -7.01577e-3
    b1 = -7.70028e-3
    b2 = -1.13864e-2
    b3 = -9.51519e-3
    c0 = -2.75915e-7
    """

    """ml/l coefficients"""
    a0 = 2.00907
    a1 = 3.22014
    a2 = 4.05010
    a3 = 4.94457
    a4 = -2.56847e-1
    a5 = 3.88767
    b0 = -6.24523e-3
    b1 = -7.37614e-3
    b2 = -1.03410e-2
    b3 = -8.17083e-3
    c0 = -4.88682e-7

    O2sol = np.exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
    return O2sol

def OxSol_umolkg(T,S):
    
    """Eq. 8 from Garcia and Gordon, 1992.
    Harded coded to do ml/l for now. If requeste, add in new mode for ug/l.

    Inputs:
    T = Potential Temperau
    S = Practical Salinity
    """

    x = S
    y = np.log((298.15 - T)/(273.15 + T))

    """umol/kg coefficients"""
    a0 =  5.80871
    a1 =  3.20291
    a2 =  4.17887
    a3 =  5.10006
    a4 = -9.86643e-2
    a5 =  3.80369
    b0 = -7.01577e-3
    b1 = -7.70028e-3
    b2 = -1.13864e-2
    b3 = -9.51519e-3
    c0 = -2.75915e-7
    

    """ml/l coefficients"""
#    a0 = 2.00907
#    a1 = 3.22014
#    a2 = 4.05010
#    a3 = 4.94457
#    a4 = -2.56847e-1
#    a5 = 3.88767
#    b0 = -6.24523e-3
#    b1 = -7.37614e-3
#    b2 = -1.03410e-2
#    b3 = -8.17083e-3
#    c0 = -4.88682e-7

    O2sol = np.exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
    return O2sol


def oxy_hysteresis_voltage(calib, voltage, scan_window=48):
    '''SBE equation for computing hysteresis from raw voltage.
    Must be run before oxy_dict.
    Because of looking backwards, must skip i = 0.

    Input:
    calib: a dict holding H1, H2, H3, Voffset
    voltage: a sequence of engineering voltages
    scan_window: an int for scans to skip between, default 48 scans at 24Hz OR 2 seconds
    '''
    output = []

    for i, x in enumerate(voltage):
        if i == 0:
            continue

        D = 1 + calib['H1']*(exp(P(i)/calib['H2']) - 1)
        C = exp(-1 * ())

    return output

def sbe43_hysteresis_voltage(calib, voltage):
    '''NOT TESTED NOT FINISHED

    SBE equation for computing hysteresis from raw voltage in relation to SBE43.
    Must be run before oxy_dict.
    Because of looking backwards, must skip i = 0.

    Parameters
    ----------
    calib : dict
        calib is a dict holding H1, H2, H3, Voffset
    voltage : array-like
        voltage is a sequence of engineering voltages

    Returns
    -------
    output : array-like
        output is the
    '''
    output = []

    for i, x in enumerate(voltage):
        if i == 0:
            continue

        D = 1 + calib['H1']*(exp(P(i)/calib['H2']) - 1)
        C = exp(-1 * ())

    return output


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
            temp = (calib['Soc'] * (V_x + calib['offset'])
                    * (1.0 + calib['A'] * T_x + calib['B'] * math.pow(T_x,2) + calib['C'] * math.pow(T_x,3) )
                    * OxSol(T_x,S_x)
                    * math.exp(calib['E'] * P_x / K_x)) #foo
            temp = round(temp,4)
            oxygen.append(temp)
    #Single mode.
    except:
        oxygen = (calib['Soc'] * (V + calib['offset'])
                  * (1.0 + calib['A'] * T + calib['B'] * math.pow(T,2) + calib['C'] * math.pow(T,3) )
                  * OxSol(T,S)
                  * math.exp(calib['E'] * P / K))
    return oxygen


def cond_dict(calib, F, t, p, units='mS'):
    """SBE equation for converting frequency to conductivity. Calculates mS/cm
    SensorID: 3

    Inputs:
    calib:
        G: coefficient
        H: coefficient
        I: coefficient
        J: coefficient
        CPcor: coefficient (nominal)
        CTcor: coefficient (nominal)

    F: instrument frequency
    t: temperature (ITS-90 degrees C)
    p: pressure (decibars)

    Output:
    sequence or float of mS/cm
    """
    try:
        Conductivity = []
        f = [x/1000 for x in F]
        for F_0, t_0, p_0 in zip(f, t, p):
            temp = ((calib['G'] + calib['H'] * math.pow(F_0,2)
                     + calib['I'] * math.pow(F_0,3)
                     + calib['J'] * math.pow(F_0,4))
                    / (1 + calib['CTcor'] * t_0 + calib['CPcor'] * p_0))
            #S/m to mS/cm
            if units == 'mS':
                temp = temp * 1
            elif units == 'S':
                temp = temp * 0.1
            temp = round(temp, 5)
            Conductivity.append(temp)
    #single mode
    except:
        f = F/1000
        Conductivity = ((calib['G'] + calib['H'] * math.pow(f,2)
                         + calib['I'] * math.pow(f,3)
                         + calib['J'] * math.pow(f,4))
                        / (1 + calib['CTcor'] * t + calib['CPcor'] * p))
        #S/m to mS/cm
        Conductivity = Conductivity * 0.1
        Conductivity = round(Conductivity,5)
    return Conductivity

def pressure_dict(calib, f, t):
    """SBE/STS(?) equation for converting pressure frequency to temperature.
    SensorID: 45

    Inputs:
    calib is a dictionary of coefficients
    T1: coefficient
    T2: coefficient
    T3: coefficient
    T4: coefficient
    T5: not used
    C1: coefficient
    C2: coefficient
    C3: coefficient
    D1: coefficient
    D2: coefficient
    AD590M: used in digiquartz temeperature correction
    AD590B: used in digiquartz temperature correction

    f: sensor frequency (usually between 30kHz and 42kHz)
    t: sensor integer from the digiquartz temperature probe

    """
    #array mode
    try:
        t_converted = []
        for x in t:
            t_converted.append((calib['AD590M'] * int(x)) + calib['AD590B'])
        pressure = []
        """Equation expecting pressure period in microseconds, so divide f by 1,000,000. """
        uf = [x/1000000 for x in f]
        for f_x, t_x in zip(uf, t_converted):
            T0 = calib['T1'] + calib['T2']*t_x + calib['T3']*math.pow(t_x,2) + calib['T4']*math.pow(t_x,3)
            w = 1-T0*T0*f_x*f_x
            temp = (0.6894759*((calib['C1']+calib['C2']*t_x+calib['C3']*t_x*t_x)*w*(1-(calib['D1']+calib['D2']*t_x)*w)-14.7))
            pressure.append(round(temp,4))
    #single mode
    except:
        t = (calib['AD590M'] * int(t)) + calib['AD590B']
        T0 = calib['T1'] + calib['T2']*t + calib['T3']*math.pow(t,2) + calib['T4']*math.pow(t,3)
        w = 1-T0*T0*f*f
        pressure = (0.6894759*((calib['C1']+calib['C2']*t+calib['C3']*t*t)*w*(1-(calib['D1']+calib['D2']*t)*w)-14.7))
    return pressure


def wetlabs_flrtd_chl_dict(calib, counts):
    """Wetlabs

    UNFINISHED

    """
    chl = calib['scale_factor'] * (output - calib['darkcounts'])
    return chl


def wetlabs_transmissometer_cstar_dict(calib, signal):
    """Wetlabs CStar Transmissiometer.
    Equation from calib sheet for S/N#: CST-479-DR, Date: October 31, 2014
    SensorID: 71

    Inputs:
    calib is a dictionary of constants/coefficients
        calib['dark'] = voltage when beam is blocked. V_d on calib sheet
        calib['air'] = voltage with clear beam path in air. V_air on calib sheet
        calib['reference'] = voltage with beam path in clean water. V_ref on calib sheet
    signal: dict/single of signal voltage

    Relationship of transmittance (Tr) to beam attenuation coefficient (c), and pathlength (x, in meters): Tr = e^-ex
    beam attenuation coefficient is determined as: c = -1/x * ln (Tr)

    """

    #array mode
    try:
        tx = []
        for signal_x in signal:
            temp = (signal_x - calib['dark'])/(calib['reference'] - calib['dark'])
            tx.append(temp)
    #single mode
    except:
        tx = (signal - calib['dark'])/(calib['reference'] - calib['dark'])
    return tx


def benthos_psa916_dict(calib, signal):
    """Equation for determining altitude from a Benthos PSA-916 altimeter.
    Equation provided by SBE as AN95, or here: http://www.seabird.com/document/an95-setting-teledyne-benthos-altimeter-sea-bird-profiling-ctd
    Equation stated as: altimeter height = [300 * voltage / scale factor] + offset
    SensorID: 0

    Inputs:
    calib is a dictionary of coefficients
        calib['ScaleFactor']: scaling factor to be applied
        calib['Offset']: offset to be applied

    signal: signal voltage
    """

    #array mode
    try:
        altitude = []
        for signal_x in signal:
            temp = (300 * signal_x / calib['ScaleFactor']) + calib['Offset']
            altitude.append(temp)
    #single mode
    except:
        altitude = (300 * signal / calib['ScaleFactor']) + calib['Offset']
    return altitude


def fluoro_seapoint_dict(calib, signal):
    """
    Raw voltage supplied from fluorometer right now, after looking at xmlcon.
    The method will do nothing but spit out the exact values that came in.
    SensorID: 11

    Inputs:
    calib is a dictionary of coefficients(?)
        GainSetting: the gain applied. according to xmlcon,
            "<!-- The following is an array index, not the actual gain setting. -->"
        Offset: offset applied

    signal: signal voltage
    """
    try:
        fluoro = []
        for signal_x in signal:
            temp = signal_x
            fluoro.append(round(temp,6))
    #single mode
    except:
        fluoro = round(signal,6)
    return fluoro

def altimeter_voltage(calib, volts):
    """
    SBE Equation for converting voltages from an altimeter to meters.

    While the SBE documentation refers to a Teledyne Benthos altimeter, the equation
    works for all altimeters typically found in the wild.

    Sensor ID: 0

    Parameters
    ----------
    calib : dict
        calib is a dict holding SerialNumber, CalibrationDate, ScaleFactor, and Offset
    volts : array-like
        volts is the voltage from the altimeter

    Returns
    -------
    bottom_distance : array-like
        bottom_distance is the distance from the altimeter to an object below it, in meters.

    Equation provided by SBE as AN95, or here:
    http://www.seabird.com/document/an95-setting-teledyne-benthos-altimeter-sea-bird-profiling-ctd
    Equation stated as: altimeter height = [300 * voltage / scale factor] + offset
    """

    # The array might come in as type=object, which throws AttributeError. Maybe this should be in try/except?
    volts = volts.astype(float)

    bottom_distance = np.around((
                                (300 * volts / calib['ScaleFactor'])
                                + calib['Offset']
                                ),1)
    return bottom_distance

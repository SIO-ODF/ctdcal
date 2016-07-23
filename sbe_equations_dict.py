"""Needs truncation implemented."""

"""Following equations taken from SBE3+ calib sheet"""

def temp_its90_dict(calib, freq):
    """SBE used equation for converting engineering units to Celcius
    according to ITS-90.
    calib is a dict holding g,h,i,j
    g, h, i, j, f0: coefficients given in calibration.
    f: frequency sampled by sensor, either as a single value or a list or tuple.

    Original form from calib sheet dated 2012:
    ITS-90 = 1/{g+h[ln(f0/f )]+i[ln2(f0/f)]+j[ln3(f0/f)]} - 273.15 (°C)

    """
    f = freq
    #array mode
    try:
        ITS90 = []
        for f_x in f:
            temp = 1/(calib['g']
                      + calib['h'] * (math.log(calib['f0']/f_x))
                      + calib['i'] * math.pow((math.log(calib['f0']/f_x)),2)
                      + calib['j'] * math.pow((math.log(calib['f0']/f_x)),3)
                     ) - 273.15
            temp = round(temp,4)
            ITS90.append(temp)
    #single mode
    except:
        ITS90 = 1/(calib['g']
                   + calib['h'] * (math.log(calib['f0']/f))
                   + calib['i'] * math.pow((math.log(calib['f0']/f)),2)
                   + calib['j'] * math.pow((math.log(calib['f0']/f)),3)
                  ) - 273.15
        ITS90 = round(ITS90,4)
    return ITS90

def temp_its68_dict(calib, freq):
    """THIS FUNCTION IS NOT FULLY IMPLEMENTED/CORRECT!

    SBE used equation for converting engineering units to Celcius
    according to IPTS-68.
    calib is a dict holding a, b, c, d, f0
    a, b, c, d, f0: coefficients given in calibration.
    f: frequency sampled by sensor, either as a single value or a list or tuple.

    Original form from calib sheet dated 2012:
    IPTS-68 = 1/{a+b[ln(f0/f )]+c[ln2(f0/f)]+d[ln3(f0/f)]} - 273.15 (°C)

    """
    f = freq
    #array mode
    try:
        ITS68 = []
        for f_x in f:
            temp = 1/(calib['a']
                      + calib['b'] * (math.log(calib['f0']/f_x))
                      + calib['c'] * math.pow((math.log(calib['f0']/f_x)),2)
                      + calib['d'] * math.pow((math.log(calib['f0']/f_x)),3)
                     ) - 273.15
            temp = round(temp,4)
            ITS68.append(temp)
    #single mode
    except:
        ITS68 = 1/(calib['a']
                   + calib['b'] * (math.log(calib['f0']/f))
                   + calib['c'] * math.pow((math.log(calib['f0']/f)),2)
                   + calib['d'] * math.pow((math.log(calib['f0']/f)),3)
                  ) - 273.15
        ITS68 = round(ITS68,4)
    return ITS68

"""Following equations taken from SBE 43 calib sheet."""

def oxy_dict(calib, P, K, T, S, V):
    """SBE equation for converting engineering units to oxygen (ml/l).

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

    #Assumes all are arrays, or none are arrays. Need way to test for them.
    #Array mode
    try:
        oxygen = []
        for P_x, K_x, T_x, S_x, V_x in zip(P, K, T, S, V):
            temp = (calib['Soc'] * (V_x + calib['Voffset'])
                    * (1.0 + calib['A'] * T_x
                       + calib['B'] * math.pow(T_x,2)
                       + calib['C'] * math.pow(T_x,3))
                    * OxSol(T_x,S_x)
                    * math.exp(calib['E'] * P_x / K_x)) #foo
            temp = round(temp,4)
            oxygen.append(temp)
    #Single mode.
    except:
        oxygen = (calib['Soc'] * (V + calib['Voffset'])
                  * (1.0 + calib['A'] * T
                     + calib['B'] * math.pow(T,2)
                     + calib['C'] * math.pow(T,3))
                  * OxSol(T,S)
                  * math.exp(calib['E'] * P / K))
    return round(oxygen,4)

def OxSol(T,S):
    """Eq. 8 from Garcia and Gordon, 1992.

    Inputs:
    T = ITS-90 Temperature
    S = Practical Salinity
    """

    x = S
    y = math.log((298.15 - T)/(273.15 + T))

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

    O2sol = math.exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y))))
                     + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
    return O2sol

def cond_dict(calib, F, t, p):
    """SBE equation for converting frequency to conductivity.

    Comes out an order of magnitude too large???? Investigate.

    Inputs:

    G: coefficient
    H: coefficient
    I: coefficient
    J: coefficient
    CPcor: coefficient (nominal)
    CTcor: coefficient (nominal)

    F: instrument frequency, single or list/tuple
    t: temperature (ITS-90 degrees C), single or list/tuple
    p: pressure (decibars), single or list/tuple
    """

    try:
        Conductivity = []
        f = [x/1000 for x in F]
        for F_x, T_x, P_x in zip(f, t, p):
            temp = ((calib['G'] + calib['H'] * math.pow(F_x,2)
                     + calib['I'] * math.pow(F_x,3)
                     + calib['J'] * math.pow(F_x,4))
                    / (1 + calib['CTcor'] * T_x + calib['CPcor'] * P_x))
            temp = round(temp, 5)
            Conductivity.append(temp)
    #single mode
    except:
        f = F/1000
        Conductivity = ((calib['G'] + calib['H'] * math.pow(f,2)
                         + calib['I'] * math.pow(f,3)
                         + calib['J'] * math.pow(f,4))
                        / (1 + calib['CTcor'] * t + calib['CPcor'] * p))
        Conductivity = round(Conductivity,5)
    return Conductivity

def pressure_dict(calib, f, t):
    """SBE/STS(?) equation for converting pressure frequency to temperature.

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
    AD590M: not used
    AD590B: not used

    f: sensor frequency (usually between 30kHz and 42kHz)
    t: sensor temperature in C

    """
    #array mode
    try:
        pressure = []
        """Equation expecting pressure period in microseconds, so divide f by 1,000,000. """
        uf = [x/1000000 for x in f]
        for f_x, t_x in zip(uf, t):
            T0 = calib['T1'] + calib['T2']*t_x + calib['T3']*math.pow(t_x,2) + calib['T4']*math.pow(t_x,3)
            w = 1-T0*T0*f_x*f_x
            temp = (0.6894759*((calib['C1']+calib['C2']*t_x+calib['C3']*t_x*t_x)*w*(1-(calib['D1']+calib['D2']*t_x)*w)-14.7))
            pressure.append(round(temp,2))
    #single mode
    except:
        T0 = calib['T1'] + calib['T2']*t + calib['T3']*math.pow(t,2) + calib['T4']*math.pow(t,3)
        w = 1-T0*T0*f*f
        pressure = (0.6894759*((calib['C1']+calib['C2']*t+calib['C3']*t*t)*w*(1-(calib['D1']+calib['D2']*t)*w)-14.7))
    return pressure

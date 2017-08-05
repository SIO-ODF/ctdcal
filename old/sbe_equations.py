"""Needs truncation implemented."""

"""Following equations taken from SBE3+ calib sheet"""

def temp_its90(g, h, i, j, f0, f):
    """SBE equation for converting engineering units to Celcius according to ITS-90.

    g, h, i, j, f0: coefficients given in calibration.
    f: frequency sampled by sensor.

    Original form from calib sheet dated 2012:
    Temperature ITS-90 = 1/{g+h[ln(f0/f )]+i[ln2(f0/f)]+j[ln3(f0/f)]} - 273.15 (°C)

    """
    ITS90 = 1/(g +
        h * (math.log(f0/f)) +
        i * math.pow((math.log(f0/f)),2) +
        j * math.pow((math.log(f0/f)),3)
        ) - 273.15
    return ITS90

def temp_its68(a, b, c, d, f0, f):
    """SBE equation for converting engineering units to Celcius according to IPTS-68.

    a, b, c, d, f0: coefficients given in calibration.
    f: frequency sampled by sensor.

    Original form from calib sheet dated 2012:
    Temperature IPTS-68 = 1/{a+b[ln(f0/f )]+c[ln2(f0/f)]+d[ln3(f0/f)]} - 273.15 (°C)

    """
    ITS68 = 1/(a +
        b * (math.log(f0/f)) +
        c * math.pow((math.log(f0/f)),2) +
        d * math.pow((math.log(f0/f)),3)
        ) - 273.15
    return ITS68

"""Following equations taken from SBE 43 calib sheet."""

def oxy(Soc, Voffset, Tau20, A, B, C, E, P, K, T, S, V):
    """SBE equation for converting engineering units to oxygen (ml/l).

    asdasdas

    Original equation from calib sheet dated 2014:
    Oxygen (ml/l) = Soc * (V + Voffset) * (1.0 + A * T + B * T + C * T ) * OxSol(T,S) * exp(E * P / K)

    """

    oxygen = Soc * (V + Voffset) * (1.0 + A * T + B * T + C * T ) * OxSol(T,S) * math.exp(E * P / K)
    return oxygen

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

    O2sol = math.exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
    return O2sol

def Cond(g, h, i, j, CPcor, CTcor, F, t, p):
    """SBE equation for converting frequency to conductivity.

    Inputs:
    g: coefficient
    h: coefficient
    i: coefficient
    j: coefficient
    CPcor: coefficient (nominal)
    CTcor: coefficient (nominal)

    F: instrument frequency
    t: temperature (ITS-90 degrees C)
    p: pressure (decibars)
    """
    f = F/1000
    Conductivity = (g + h * round(math.pow(f,2),3) + i * round(math.pow(f,3),3) + j * round(math.pow(f,4),3) ) / (1 + CTcor * t + CPcor * p)
    return round(Conductivity,5)

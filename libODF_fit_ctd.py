#!/usr/bin/env python
from math import sqrt
import numpy as np
import conversions as convert
import density_enthalpy_48 as density
import libODF_process_ctd as process_ctd
import gsw
import csv
import requests
import os
import json

M = 31.9988   #Molecular weight of O2
R = 831.432    #/* Gas constant, X 10 J Kmole-1 K-1 */
D = 1.42905481 #  /* density O2 g/l @ 0 C */

def offset(offset, col, inMat):
    """offset column of data 

    Input:
        - inMat, 1d numpy array with return True np.isnans()
    Output:
        - Mat with offset applied to column of data 
    Example:
        >>> outArray = offset(offset, col, inMat) 
    """
    for i in range(0, len(inMat)):
        inMat[col][i] = inMat[col][i] + offset

    return inMat

def IESRho(s, t, p):
    bars  = p * 0.1
#/* pressure in bars */
#/*
#**  Calculate Rho(s,t,0.0)
#*/
    rhow  =  (999.842594 +t*
       ( 0.06793952 +t*
       (-0.00909529 +t*
       ( 1.001685e-4+t*
       (-1.120083e-6+t*
       6.536332e-9)))))
    #/* rhow = density of pure water kg/m**3 */
    kw    = (((t*(-0.0040899+
       t*( 7.6438e-5+
       t*(-8.2467e-7+
       t*  5.3875e-9))))+
       0.824493)*s)
    #/* pure water secant bulk modulus */
    termc = s * sqrt(s);
    kst0  = ((-0.00572466 +
       t*( 1.0227e-4  +
       t*(-1.6546e-6))) * termc)
    #/* k(s,t,0) */
    rho   = rhow + kw + kst0 + 4.8314e-4 *s*s
#      /* rho(s,t,0.0)  kg/m**3 */
#      /*
#      **  Calculate pressure effects.
#      */
    if (bars > 0.0):
#  	/*
#  	**                 Rho(s,t,0)
#  	**  Rho(s,t,p) = -------------
#  	**                        p
#  	**               1.0 - -------
#  	**                     k(s,t,p)
#  	*/
        kw    = (t*(148.4206          +
           t*( -2.327105        +
           t*(  0.01360477      +
           t*( -5.155288e-5)))) +
           1.965221e4)
        kst0  =  ( (54.6746    +
           t*(-0.603459  +
           t*( 0.0109987 +
           t*(-6.167e-5))))*s + kw +
           ( 0.07944   +
           t*( 0.016483  +
           t*(-5.3009e-4)))*termc)
#  	/*
#  	**  Calculate pressure terms.
#  	*/
        terma = (    3.239908     +
           t*( 0.00143713   +
           t*( 1.16092e-4   +
           t*(-5.77905e-7)))+
           ( 0.0022838    +
           t*(-1.0981e-5    +
           t*(-1.6078e-6)))*s +
           1.91075e-4*termc)
        termb = (8.50935e-5  +
           t*(-6.12293e-6  +
           t*  5.2787e-8)  +
           (-9.9348e-7   +
           t*( 2.0816e-8   +
           t*  9.1697e-10))*s)
        kstp  = kst0 + bars*(terma + bars*termb)
        rho   = rho/(1.0-bars/kstp)
    return rho

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

#
#def H2odVdT(v1, t1):
#    return (v1*rho_t(t1)/rho_t(20))
#
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
def o2_calc(o2flasks, o2path, sal_col, btlpath):
#    qual = load_qual("/Volumes/public/O2Backup/o2_codes_001-083.csv")
    btl_data = process_ctd.dataToNDarray(btlpath,None,True,',',0)
    btls = btl_data['btl_fire_num'][1:]
    btls = btls.astype(float)
    sal = (btl_data[sal_col][1:]) 
    sal = sal.astype(float)
    print(len(sal))

    o2ml = np.ndarray(shape=len(sal), dtype=np.float)
    o2kg = np.ndarray(shape=len(sal), dtype=np.float)

    with open(o2path, 'r') as f:
        rho = IESRho
        params = next(f).strip().split()

        titr   = float(params[0])
        blank  = float(params[1])
        kio3_n = float(params[2])
        kio3_v = float(params[3])
        kio3_t = float(params[4])
        thio_t = float(params[5])

        thio_n = thio_n_calc(titr, blank, kio3_n, kio3_v, kio3_t, thio_t)
        rho_stp = rho_t(20)

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
            if float(bottle) in btls:
                btl_index = np.where(btls == bottle)[0]
                print(btl_index)
                print(sal[btl_index])
                o2ml[btl_index] = (((titr_20c - blank) * thio_n * 5.598 - 0.0017)/((flask_vol - 2.0) * 0.001))
                o2kg[btl_index] = mll_to_umolkg(o2ml[btl_index], float(sal[btl_index]), draw_temp,rho)
            else: o2kg[btl_index] = '-999'
#           row_dict = {"station": str(station), "cast": str(cast),"bottle": str(bottle), "o2": o2kg}
    return o2kg, o2ml
#
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

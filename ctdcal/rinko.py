from collections import namedtuple
import numpy as np
import oxy_fitting
import pandas as pd
import scipy

RinkoO2Cal = namedtuple("RinkoO2Cal", [*"ABCDEFGH"])
RinkoTMPCal = namedtuple("RinkoTMPCal", [*"ABCD"])

def rinko_o2_cal_parameters(**kwargs):
    A = kwargs.get("A", -4.524084e1)
    B = kwargs.get("B",  1.449377e2)
    C = kwargs.get("C", -3.051590e-1)
    D = kwargs.get("D",  1.065300e-2)
    E = kwargs.get("E",  4.000000e-3)
    F = kwargs.get("F",  6.250000e-5)
    G = kwargs.get("G",  0.000000e0)
    H = kwargs.get("H",  1.000000e0)
    return RinkoO2Cal(A,B,C,D,E,F,G,H)

def rinko_temperature_cal_parameters(**kwargs):
    A = kwargs.get("A", -5.305905e0)
    B = kwargs.get("B",  1.666857e1)
    C = kwargs.get("C", -2.142681e0)
    D = kwargs.get("D",  4.582805e-1)
    return RinkoTMPCal(A,B,C,D)

def rinko_temperature(v, tmp_cal:RinkoTMPCal):
    if type(tmp_cal) is not RinkoTMPCal:
        raise ValueError("tmp_cal must be of type RinkoTMPCal")

    A, B, C, D = tmp_cal
    return A + B*v + C*v**2 + D*v**3

def rinko_pprime_aro_cav(v, t, o2_cal:RinkoO2Cal):
    """
    Calculates Rinko P' of the equation P = G + H * P'
    where P is DO physical value IN PERCENT [%]
    """
    A, B, C, D, E, F, G, H = o2_cal

    term_1_denominator = 1 + D*(t-25) + F*(t-25)**2
    term_1 = A/term_1_denominator

    term_2_denominator = v * (1 + D*(t-25) + F*(t-25)**2) + C
    term_2 = B/term_2_denominator

    return term_1 + term_2

def rinko_saturation(pprime, o2_cal:RinkoO2Cal):
    """
        Calculates Rinko P of the equation P = G + H * P'
        where P is DO physical value IN PERCENT [%]
    """

    A, B, C, D, E, F, G, H = o2_cal

    return G + H * pprime

def rinko_correct_for_pressure(p, d, o2_cal:RinkoO2Cal):
    """Note that the pressure term, d, must be in MPa

    1 decibar = 0.01 Mpa
    """
    A, B, C, D, E, F, G, H = o2_cal

    return p*(1 + E*d)

def rinko_saturation(df, film="B", model="ARO-CAV", **kwargs):
    pass

def rinko_oxy_eq( press, temp, oxyvo, os, o2_cal:RinkoO2Cal):

    #Calculate pprime

    pprime = rinko_pprime_aro_cav(oxyvo,temp,o2_cal)

    # Calculate P (DO physical value in %)

    #p = rinko_saturation(pprime, o2_cal)

    # Correct for pressure * d is pressure in Mpa *

    d = press * 0.01

    p_corr = rinko_correct_for_pressure(pprime,d,o2_cal)

    # Divide by 100 to get percents in the form of 0.xx

    p_corr = p_corr / 100

    # Multiply by OS to get DO (os can be in either ml/l or umol/kg)

    DO = p_corr * os

    return DO

def rinko_curve_fit_eq(X, a, b, c, d, e, f, g, h):
    """
    Same as rinko_oxy_eq, but in a form that is more suitible for scipy's curve fit routine
    X contains pressure, temperature, voltage, and OS (the normal arguments for rinko_oxy_eq)
    """

    press, temp, oxyvo, os = X
    o2_cal = RinkoO2Cal(a, b, c, d, e, f, g, h)

    #Calculate pprime

    pprime = rinko_pprime_aro_cav(oxyvo,temp,o2_cal)

    # Calculate P (DO physical value in %)

    #p = rinko_saturation(pprime, o2_cal)

    # Correct for pressure * d is pressure in Mpa *

    d = press * 0.01

    p_corr = rinko_correct_for_pressure(pprime,d,o2_cal)

    # Divide by 100 to get percents in the form of 0.xx

    p_corr = p_corr / 100

    # Multiply by OS to get DO (os can be in either ml/l or umol/kg)

    DO = p_corr * os

    return DO

def rinko_oxygen_cal(o2_cal,pressure,temp,oxyvolts,os,ref_oxy,switch):

    """"

    Rinko oxygen fitting routine using the equation:
        Calculates Rinko P of the equation P = G + H * P'
        where P is DO physical value IN PERCENT [%]

    Fits 7 Coef:
    coef[0] = A
    coef[1] = B
    coef[2] = C
    coef[3] = D
    coef[4] = E
    coef[5] = F
    coef[6] = G


    """

    rinko_oxy = rinko_oxy_eq(pressure, temp, oxyvolts, os, o2_cal)

    #Weight Determination
    if switch == 1:

        weights = oxy_fitting.calculate_weights(pressure)

        resid = ((weights * (ref_oxy - rinko_oxy))**2) / (np.sum(weights)**2) #Original way (np.sum(weights)**2)

    elif switch == 2:
        #L2 Norm
        resid = (ref_oxy - rinko_oxy)**2

    elif switch == 3:
        #ODF residuals

        resid = np.sqrt(((ref_oxy - rinko_oxy)**2) / (np.std(rinko_oxy)**2))

    elif switch == 4:
        # Weighted ODF residuals

        weights = oxy_fitting.calculate_weights(pressure)
        resid = np.sqrt(weights * ((ref_oxy - rinko_oxy)**2) / (np.sum(weights)**2))#(np.std(ctd_oxy_mlL)**2))

    elif switch == 5:

        weights = oxy_fitting.calculate_weights(pressure)

        resid = ((weights * (ref_oxy - rinko_oxy))**2) / (np.sum(weights**2))


    return resid

def rinko_oxygen_fit(btl_prs,btl_oxy, btl_sigma, ctd_sigma, ctd_os, ctd_prs, ctd_tmp, rinkovolts, coef0, btl_ssscc=None):

    # Create DF for good and questionable values

    bad_df = pd.DataFrame()
    good_df = pd.DataFrame()

    # Construct Dataframe from bottle and ctd values for merging
    if 'btl_ssscc' in locals():
        btl_dict = {'CTDPRS_rinko_btl':btl_prs, 'REFOXY_rinko':btl_oxy, 'sigma_rinko_btl':btl_sigma, 'SSSCC_rinko':btl_ssscc}
    else:
        btl_dict = {'CTDPRS_rinko_btl':btl_prs, 'REFOXY_rinko':btl_oxy, 'sigma_rinko_btl':btl_sigma}
    btl_data = pd.DataFrame(btl_dict)
    time_dict = {'CTDPRS_ctd_rinko':ctd_prs, 'sigma_rinko_ctd':ctd_sigma, 'OS_rinko_ctd':ctd_os, 'CTDTMP_rinko':ctd_tmp, 'CTDRINKOVOLTS':rinkovolts}
    time_data = pd.DataFrame(time_dict)

    # Sort DataFrames by sigma0
    time_data.sort_values('sigma_rinko_ctd', inplace=True)
    btl_data.sort_values('sigma_rinko_btl', inplace=True)
    btl_data.dropna(subset=['REFOXY_rinko'], inplace=True)

    # Merge DF
    merged_df = pd.merge_asof(btl_data, time_data, left_on='sigma_rinko_btl', right_on='sigma_rinko_ctd', direction='nearest', suffixes=['_btl','_ctd'])
    merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],coef0)

    # Curve fit (weighted)
    p0 = coef0
    weights = 1/(np.sqrt(merged_df['CTDPRS_ctd_rinko']))

    try:
        cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq,
                                                  (merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']),
                                                  merged_df['REFOXY_rinko'], coef0, sigma=weights, absolute_sigma=False)

        merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)
        merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
        stdres = np.std(merged_df['res_rinko'])
        cutoff = stdres * 2.8

        thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
        bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
        bad_df = pd.concat([bad_df, bad_values])
        merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]


        while not thrown_values.empty:

            p0 = cfw_coef
            weights = 1/(np.power(merged_df['CTDPRS_ctd_rinko'],1/2))
            cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq, (merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']), merged_df['REFOXY_rinko'], p0, sigma=weights, absolute_sigma=False)
            merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)
            merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
            stdres = np.std(merged_df['res_rinko'])
            cutoff = stdres * 2.8
            thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
            bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
            bad_df = pd.concat([bad_df, bad_values])
            merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

    except RuntimeError:

            try:#Nested try/except could be better
                print('Weighted Curve fitting failed...using Unweighted Fitting')
                cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq, (merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']), merged_df['REFOXY_rinko'], coef0)
                merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)

                merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
                stdres = np.std(merged_df['res_rinko'])
                cutoff = stdres * 2.8

                thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
                bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
                bad_df = pd.concat([bad_df, bad_values])
                merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

                while not thrown_values.empty:

                    p0 = cfw_coef
                    cfw_coef , cov = scipy.optimize.curve_fit(rinko_curve_fit_eq, (merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd']), merged_df['REFOXY_rinko'], coef0)
                    merged_df['CTDRINKO'] = rinko_oxy_eq(merged_df['CTDPRS_ctd_rinko'],merged_df['CTDTMP_rinko'],merged_df['CTDRINKOVOLTS'],merged_df['OS_rinko_ctd'],cfw_coef)

                    merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
                    stdres = np.std(merged_df['res_rinko'])
                    cutoff = stdres * 2.8

                    thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
                    bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
                    bad_df = pd.concat([bad_df, bad_values])
                    merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]

            except:
                print('Logging RINKO coef...')
                cfw_coef = coef0
                merged_df['res_rinko'] = merged_df['REFOXY_rinko'] - merged_df['CTDRINKO']
                stdres = np.std(merged_df['res_rinko'])
                cutoff = stdres * 2.8
                thrown_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
                bad_values = merged_df[np.abs(merged_df['res_rinko']) > cutoff]
                merged_df = merged_df[np.abs(merged_df['res_rinko']) <= cutoff]


    #coef_dict[station] = cfw_coef
    good_df = pd.concat([good_df, merged_df])
    good_df['CTDRINKO_FLAG_W'] = 2
    bad_df = pd.concat([bad_df, bad_values])
    bad_df['CTDRINKO_FLAG_W'] = 3
    df = pd.concat([good_df,bad_df])
    df.sort_values(by='CTDPRS_rinko_btl',ascending=False,inplace=True)
    oxy_df = df.copy()

    return cfw_coef, df



if __name__ == "__main__":
    print(rinko_temperature(1.565323565, rinko_temperature_cal_parameters()))
    print(rinko_correct_for_pressure(
    rinko_pprime_aro_cav(1.5421, 0.442, rinko_o2_cal_parameters()),
    2/100,
    rinko_o2_cal_parameters()
    ))

import numpy as np
import pytest

from ctdcal import oxy_fitting, rinko

#   Coefs from CTDCAL run
coefs = (2.4837e+0, 1.7548e-1, -1.0020e-3, 3.1672e-2, -1.8098e-1, 2.685e-1, 9.7101e-2)
#   Fake pressure points
pressure = np.array([50, 150, 250, 400, 800, 1500, 2500, 5000])
refoxy = np.array([400, 160, 170, 200, 280, 335, 360, 400])
#   Fake inputs tuple
#   (raw voltage, pressure, temperature, salinity, oxygen solubility)
inputs = (1.5, 500, 20, 35, 100)


def test_rinko_DO():
    #   Assign test cases
    p_prime = 0.5
    G = 10
    H = 5

    assert rinko.rinko_DO(p_prime, G, H) == 10 + 5 * 0.5

def test_rinko_p_prime():
    #   Shared raw and temperature arrays for test cases
    #   Assuming "raw" are voltages
    N = np.array([1.7, 1.9, 2.0])
    t = np.array([20, 25, 22])

    A, B, C, D, E, F, G, H = 2, 3, 4, 5, 6, 7, 8, 9
    p_solved = (A / (1 + D * (t - 25))) + (B / ((N - F) * (1 + D * (t - 25)) + C + F))
    p_prime = rinko.rinko_p_prime(N, t, A, B, C, D, E, F, G, H)
    assert np.allclose(p_prime, p_solved)

    #   Now with 0297 coeffs
    A, B, C, D, E, F, G, H = -4.367428e+01, 1.376636e+02, -3.647983e-01, 1.044300e-02, 4.300000e-03, 6.810000e-05, 0, 1
    p_solved = (A / (1 + D * (t - 25))) + (B / ((N - F) * (1 + D * (t - 25)) + C + F))
    p_prime = rinko.rinko_p_prime(N, t, A, B, C, D, E, F, G, H)
    assert np.allclose(p_prime, p_solved)

def test_correct_pressure():
    #   Assign test cases
    P = np.array([80, 90, 85])  #   Corrected % oxygen
    d = np.array([0.1, 0.2, 0.15])  #   Assuming 1 MPa = 10 dbar
    E = 0.0043  #   Grabbed from SN 0297
    p_solved = P * (1 + E * d)

    p_corr = rinko.correct_pressure(P, d, E)

    assert np.allclose(p_corr, p_solved)

def test_salinity_correction():
    #   Assign test cases
    DO_c = np.array([400, 200, 180])
    T = np.array([20.5, 29.5, 33.0])
    S = np.array([33.0, 34.5, 35.2])
    B0 = -6.24523e-3
    B1 = -7.37614e-3
    B2 = -1.03410e-2
    B3 = -8.17083e-3
    C0 = -4.88682e-7
    T_scaled = np.log((298.15 - T) / (273.15 + T))
    DO_sc_manual = DO_c * np.exp(
        S * (B0 + (B1 * T_scaled) + (B2 * T_scaled ** 2) + (B3 * T_scaled ** 3))
        + C0 * S ** 2
    )

    DO_sc = rinko.salinity_correction(DO_c, T, S)
    assert np.allclose(DO_sc, DO_sc_manual)

    #   NaN handling
    S = np.array([35.0, np.nan, 34.5])
    DO_sc = rinko.salinity_correction(DO_c, T, S)
    assert np.isnan(DO_sc[1])

@pytest.mark.parametrize("coefs, inputs", [(coefs, inputs)])
def test_Uchida_DO_eq(coefs, inputs):
    V_r, P, T, S, o2_sol = inputs
    K_sv = coefs[0] + (coefs[1] * T) + (coefs[2] * T ** 2)
    V0 = (1 + coefs[3] * T)
    Vc = (coefs[4] + coefs[5] * V_r)
    o2_sat = ((V0 / Vc) - 1) / K_sv
    DO = o2_sat * o2_sol
    DO_c = DO * (1 + coefs[6] * P / 1000) ** (1 / 3)
    DO_sc_manual = rinko.salinity_correction(DO_c, T, S)

    DO_sc = rinko._Uchida_DO_eq(coefs, inputs)
    assert np.allclose(DO_sc, DO_sc_manual)

@pytest.mark.parametrize("coefs, pressure, inputs, refoxy", [(coefs, pressure, inputs, refoxy)])
def oxy_weighted_residual(coefs, weights, inputs, refoxy, pressure):

    weights = oxy_fitting.calculate_weights(pressure)
    residuals_manual = np.sum(
        (weights * (refoxy - rinko._Uchida_DO_eq(coefs, inputs)) ** 2)
    ) / np.sum(weights ** 2)

    residuals = oxy_weighted_residual(coefs, weights, inputs, refoxy, pressure)

    assert np.allclose(residuals, residuals_manual)
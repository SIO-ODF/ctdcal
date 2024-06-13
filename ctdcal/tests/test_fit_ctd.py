import numpy as np
import pandas as pd
import pytest
import yaml

from ctdcal import fit_ctd


@pytest.mark.parametrize("xN, yN", [(1, 0), (0, 1), (1, 1), (2, 1), (1, 2)])
def test_multivariate_fit(xN, yN):
    data = [0.0, 0.0]
    coef_names = ["x", "y"]
    x_coefs = [f"x{n}" for n in np.arange(1, xN + 1)] if xN > 0 else []
    y_coefs = [f"y{n}" for n in np.arange(1, yN + 1)] if yN > 0 else []

    # ensure number of coefs produced is correct (including constant offset)
    np.testing.assert_allclose(
        fit_ctd.multivariate_fit(data, (data, xN), (data, yN)), np.zeros(xN + yN + 1)
    )

    # check that coefs are named properly
    fit = fit_ctd.multivariate_fit(data, (data, xN), (data, yN), coef_names=coef_names)
    assert all(coef in fit.keys() for coef in x_coefs + y_coefs)
    assert "c0" in fit.keys()

    # check error if len of inputs and coefs differ
    with pytest.raises(ValueError):
        fit_ctd.multivariate_fit(data, (data, xN), (data, yN), coef_names=["x"])

    # check error if input is not tuple
    with pytest.raises(TypeError):
        fit_ctd.multivariate_fit(data, (data, xN), [data, yN])


def test_apply_polyfit():
    y = np.array([1, 2, 3])

    # check correction is applied correctly (no dependent variables)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0,)), y)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (1, 0)), y + 1)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0, 1)), y + y)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0, 0.5)), y + 0.5 * y)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0, 0, 1)), y + y ** 2)

    # check correction is applied correctly (with dependent variables)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0,), (y, (0,))), y)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0,), (y, (1,))), y + y)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0,), (y, (1.0,))), y + y)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0.0,), (y, (1,))), y + y)
    np.testing.assert_array_equal(fit_ctd.apply_polyfit(y, (0.0,), (y, (1.0,))), y + y)
    np.testing.assert_array_equal(
        fit_ctd.apply_polyfit(y, (0,), (y, (0, 1))), y + y ** 2
    )

    # check error if input is not tuple
    with pytest.raises(TypeError):
        fit_ctd.apply_polyfit(y, (0,), [y, (0,)])

def test_generate_yaml(tmp_path):
    fname = str(tmp_path) + "filename.yaml"
    fit_ctd.generate_yaml("filename.yaml", str(tmp_path))
    with open(fname, 'r') as f:
        generated_data = yaml.safe_load(f)

    assert generated_data['t2']['ssscc_t1']['T_order'] == 0
    assert generated_data['c1']['ssscc_c1']['zRange'] == '1000:6000'

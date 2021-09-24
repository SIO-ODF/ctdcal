from ctdcal import fit_ctd
import numpy as np
import pandas as pd
import pytest


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

"""
Classes, definitions and utilities for use across ctdcal fitting modules.
"""
import json

import numpy as np
from munch import Munch


class BottleFlags(Munch):
    """
    A dictionary class with the attribute-style access of Munch, plus methods
    for adding nodes and flag data, and loading or saving to/from a JSON file.
    TODO: Move this to ctdcal.flagging.common after reorg
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_node(self, label, keys):
        """
        Add a new empty node.

        Parameters
        ----------
        label : str
            name of node
        keys : list
            names of node keys
        """
        key_dict = {k: [] for k in keys}
        node = ({label: BottleFlags(key_dict)})
        self.update(node)

    def update_node(self, **kwargs):
        """
        Add a row of values to a node. Requires an existing node with keys
        that are already defined.

        Parameters
        ----------
        kwargs
            Named keyword arguments as required for a node.

        """
        for k, v in kwargs.items():
            self[k].append(v)

    def save(self, fname):
        """
        Export the flags to a JSON file. Existing contents of the file
        will be overwritten. The file must already exist.

        Pretty printing may be removed or made optional in future updates.

        Parameters
        ----------
        fname : str or Path-like
            filename
        """
        with open(fname, 'w') as f:
            f.write(self.toJSON(indent=4))


class NodeNotFoundError(Exception):
    pass


# Function definitions
# --------------------

# BottleFlag wrangling
# TODO: Move to ctdcal.flagging.common after reorg
def df_node_to_BottleFlags(df):
    """
    Convert a flag node from a DataFrame to a formatted BottleFlags object

    Parameters
    ----------
    df : DataFrame
        Node data to convert.

    Returns
    -------
    BottleFlags object
    """
    node_dict = df.to_dict()
    for k, v in node_dict.items():
        node_dict[k] = [vv for kk, vv in v.items()]
    return BottleFlags(node_dict)


def get_node(fname, label):
    """
    Return a node from a BottleFlags file.

    Parameters
    ----------
    fname : str or Path-like
        Name of the BottleFlags file.
    label : str
        Name of the node to return.

    Returns
    -------
    BottleFlags object
    """
    with open(fname, 'r') as f:
        flags = json.load(f)
        if label in flags:
            return BottleFlags(flags[label])
        else:
            raise NodeNotFoundError


def save_node(fname, node, label, create_new=False):
    """
    Save an updated node to a BottleFlags file. The file must already exist.
    Optionally create a new node if the named node does not exist.

    Parameters
    ----------
    fname : str or Path-like
        Name of the BottleFlags file.
    node : BottleFlags object
        The formatted node data.
    label : str
        Name of the node to save.
    create_new : Bool
        If true, create a new node if it does not exist. Default is false.

    """
    with open(fname, 'r') as f:
        buf = f.read()
    if buf == '':
        # File exists but is empty
        flags = BottleFlags()
    else:
        flags = BottleFlags.fromJSON(buf)
    if label in flags or create_new is True:
        flags[label] = node
    else:
        raise NodeNotFoundError("The node '%s' was not found in %s" % (label, fname))
    flags.save(fname)


def multivariate_fit(y, *args, coef_names=None, const_name="c0"):
    """
    Least-squares fit data using multiple dependent variables. Dependent variables must
    be provided in tuple pairs of (data, order) as positional arguments.

    If coef_names are defined, coefficients will be returned as a dict. Otherwise,
    coefficients are return as an array in the order of the dependent variables, sorted
    by decreasing powers.

    Parameters
    ----------
    y : array-like
        Indepedent variable to be fit
    args : tuple
        Pairs of dependent variable data and fit order (i.e., (data, order))
    coef_names : list-like, optional
        Base names for coefficients (i.e., "a" for 2nd order yields ["a2", "a1"])
    const_name : str, optional
        Name for constant offset term

    Returns
    -------
    coefs : array-like
        Least-squares fit coefficients in decreasing powers

    Examples
    --------
    Behavior when coef_names is None:

    >>> z = [1, 4, 9]
    >>> x = [1, 3, 5]
    >>> y = [1, 2, 3]
    >>> multivariate_fit(z, (x, 2), (y, 1))
    array([0.25, 0.375, 0.25, 0.125])  # [c1, c2, c3, c4]

    where z = (c1 * x ** 2) + (c2 * x) + (c3 * y) + c4

    Behavior when coef_names is given:

    >>> z = [1, 4, 9]
    >>> x = [1, 3, 5]
    >>> y = [1, 2, 3]
    >>> multivariate_fit(z, (x, 2), (y, 1), coef_names=["a", "b"])
    {"a2": 0.25, "a1": 0.375, "b1": 0.25, "c0": 0.125}

    where z = (a2 * x ** 2) + (a1 * x) + (b1 * y) + c0
    """
    to_dict = True
    if coef_names is None:
        to_dict = False
        coef_names = [""] * len(args)  # needed to prevent zip() error
    elif len(args) != len(coef_names):
        raise ValueError(
            "length of coef_names must match the number of dependent variables"
        )

    # iteratively build fit matrix
    rows, names = [], []
    for arg, coef_root in zip(args, coef_names):
        if type(arg) is not tuple:
            raise TypeError(f"Positional args must be tuples, not {type(arg)}")

        series, order = arg
        for n in np.arange(1, order + 1)[::-1]:
            rows.append(series ** n)  # n is np.int64 so series will cast to np.ndarray
            names.append(f"{coef_root}{n}")

    # add constant offset term
    rows.append(np.ones(len(y)))
    names.append(const_name)

    fit_matrix = np.vstack(rows)
    coefs = np.linalg.lstsq(fit_matrix.T, y, rcond=None)[0]

    return dict(zip(names, coefs)) if to_dict else coefs


def apply_polyfit(y, y_coefs, *args):
    """
    Apply a polynomial correction to series of data. Coefficients should be provided in
    increasing order (i.e., a0, a1, a2 for y_fit = y + a2 * y ** 2 + a1 * y + a0)

    For the independent variables (y), coefficients start from the zero-th order (i.e.,
    constant offset). For dependent variables (args), coefficients start from the first
    order (i.e., linear term).

    Parameters
    ----------
    y : array-like
        Independent variable data to be corrected
    y_coefs : tuple of float
        Independent variable fit coefficients (i.e., (coef0, ..., coefN))
    args : tuple of (array-like, (float, float, ...))
        Dependent variable data and fit coefficients (i.e., (data, (coef1, ..., coefN)))

    Returns
    -------
    fitted_y : array-like
        Independent variable data with polynomial fit correction applied

    Examples
    --------
    Behavior without additional args:

    >>> y = [2, 4, 6]
    >>> apply_polyfit(y, (1, 2, 3))  # y0 = 1; y1 = 2; y2 = 3
    array([ 19.,  61., 127.])

    where fitted_y = y + y0 + (y1 * y) + (y2 * y ** 2)

    Behavior with additional args:

    >>> y = [2, 4, 6]
    >>> x = [1, 2, 3]
    >>> apply_polyfit(y, (1,), (x, (2, 3)))  # y0 = 1; x1 = 2; x2 = 3
    array([ 8., 21., 40.])

    where fitted_y = y + y0 + (x1 * x) + (x2 * x ** 2)
    """
    fitted_y = np.copy(y).astype(float)
    for n, coef in enumerate(y_coefs):
        fitted_y += coef * np.power(y, n)

    for arg in args:
        if type(arg) is not tuple:
            raise TypeError(f"Positional args must be tuples, not {type(arg)}")

        series, coefs = arg
        for n, coef in enumerate(coefs):
            fitted_y += coef * np.power(series, n + 1)

    return fitted_y

"""Holodeck2 utility functions.

"""

import numpy as np
import scipy as sp
import scipy.stats

from .constants import *    # noqa


def midpoints_lin(vals):
    return 0.5*(vals[1:] + vals[:-1])


def midpoints_log(vals):
    return np.power(10.0, 0.5*(np.log10(vals[1:]) + np.log10(vals[:-1])))


def stats(vals, percs=None, prec=2, weights=None) -> str:
    """Return a string giving quantiles of the given input data.

    Parameters
    ----------
    vals : npt.ArrayLike,
        Input values to get quantiles of.
    percs : npt.ArrayLike, optional
        Quantiles to calculate.
    prec : int, optional
        Precision in scientific notation of output.

    Returns
    -------
    rv : str
        Quantiles of input formatted as a string of scientific notation values.

    Raises
    ------
    TypeError: raised if input data is not iterable.

    """
    try:
        if len(vals) == 0:        #### type: ignore
            return "[]"
    except TypeError:
        raise TypeError(f"`vals` (shape={np.shape(vals)}) is not iterable!")

    if percs is None:
        percs = [sp.stats.norm.cdf(1), 0.95, 1.0]
        percs = np.array(percs)
        percs = np.concatenate([1-percs[::-1], [0.5], percs])

    # stats = np.percentile(vals, percs*100)
    stats = quantiles(vals, percs, weights=weights)
    _rv = ["{val:.{prec}e}".format(prec=prec, val=ss) for ss in stats]
    rv = ", ".join(_rv)
    return rv


def quantiles(values, percs=None, sigmas=None, weights=None, axis=None,
              values_sorted=False, filter=None):
    """Compute weighted percentiles.

    NOTE: if `values` is a masked array, then only unmasked values are used!

    Parameters
    ----------
    values: (N,)
        input data
    percs: (M,) scalar
        Desired quantiles of the data.  Within range of [0.0, 1.0].
    weights: (N,) or `None`
        Weights for each input data point in `values`.
    axis: int or `None`,
        Axis over which to calculate quantiles.
    values_sorted: bool
        If True, then input values are assumed to already be sorted.
        Otherwise they are sorted before calculating quantiles (for efficiency).

    Returns
    -------
    percs : (M,) float
        Array of quantiles of the input data.

    """
    if not isinstance(values, np.ma.MaskedArray):
        values = np.asarray(values)

    if (percs is None) == (sigmas is None):
        err = "either `percs` or `sigmas`, and not both, must be given!"
        log.error(err)
        raise ValueError(err)

    if percs is None:
        percs = sp.stats.norm.cdf(sigmas)

    if np.ndim(values) > 1:
        if axis is None:
            values = values.flatten()
    elif (axis is not None):
        raise ValueError("Cannot act along axis '{}' for 1D data!".format(axis))

    percs = np.array(percs)
    if weights is None:
        ww = np.ones_like(values)
    else:
        ww = np.array(weights)

    try:
        ww = np.ma.masked_array(ww, mask=values.mask)  # type: ignore
    except AttributeError:
        pass

    assert np.all(percs >= 0.0) and np.all(percs <= 1.0), 'percentiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values, axis=axis)
        values = np.take_along_axis(values, sorter, axis=axis)
        ww = np.take_along_axis(ww, sorter, axis=axis)

    if axis is None:
        weighted_quantiles = np.cumsum(ww) - 0.5 * ww
        weighted_quantiles /= np.sum(ww)
        percs = np.interp(percs, weighted_quantiles, values)
        return percs

    ww = np.moveaxis(ww, axis, -1)
    values = np.moveaxis(values, axis, -1)

    weighted_quantiles = np.cumsum(ww, axis=-1) - 0.5 * ww
    weighted_quantiles /= np.sum(ww, axis=-1)[..., np.newaxis]
    percs = [np.interp(percs, weighted_quantiles[idx], values[idx])
             for idx in np.ndindex(values.shape[:-1])]
    percs = np.array(percs)
    return percs
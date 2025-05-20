/**
 *
 */

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "hdf5.h"

using namespace std;

namespace utils {

    bool is_almost_equal(double a, double b, double atol = 1E-8, double rtol = 1E-6) {
        double rval;
        if ((fabs(a) > 0) && (fabs(b) > 0)) {
            rval = fmin(fabs(a), fabs(b));
        } else {
            rval = fmax(fabs(a), fabs(b));
        }
        // #ifdef DEBUG
        //     printf(
        //         "is_almost_equal(): diff=%.8e, atol=%.8e, rtol*rval=%.8e*%.8e=%.8e\n",
        //         fabs(a-b), atol, rtol, rval, rtol * rval
        //     );
        // #endif
        return (fabs(a - b) <= (atol + rtol * rval));
    }

    void index_2d_to_1d(int i, int j, int dim1, int dim2, int* idx) {
        *idx = (i * dim2) + j;
    }

    void index_3d_to_1d(int i, int j, int k, int dim1, int dim2, int dim3, int* idx) {
        *idx = (i * dim2 * dim3) + (j * dim3) + k;
    }

    void index_1d_to_3d(int idx, int dim1, int dim2, int dim3, int* i, int* j, int* k) {
        *i = idx / (dim2 * dim3);
        *j = (idx % (dim2 * dim3)) / dim3;
        *k = idx % dim3;
    }

    void index_1d_to_2d(int idx, int dim1, int dim2, int* i, int* j) {
        *i = idx / dim2;
        *j = idx % dim2;
    }

    /*
    !NOT WORKING!
    double* quantiles(
        double* values, int num_vals, double* percs, int num_percs,
        double* weights = nullptr, bool values_sorted = false
    ) {

        int* indices;
        if (values_sorted) {
            indices = new int[num_vals];
            for (int i = 0; i < num_vals; ++i) {
                indices[i] = i;
            }
        } else {
            indices = argsort(values, num_vals);
        }

        // printf("quantiles(): created sorted indices.\n");

        // Calculate cumulative weights
        double* cum_weights = new double[num_vals];
        cum_weights[0] = (weights != nullptr) ? weights[0] : 1.0;
        for (int i = 1; i < num_vals; i++) {
            // printf("indices[%d]=%d, cum_weights[%d-1]=%.8e\n", i, indices[i], i, cum_weights[i-1]);
            cum_weights[i] = cum_weights[indices[i - 1]] + ((weights != nullptr) ? weights[indices[i]] : 1.0);
            #ifdef DEBUG
            if(values[indices[i]] < values[indices[i - 1]]) {
                std::string err_msg = std::format(
                    "Error `quantiles()` - values are not sorted!: values[{:d}]={:.8e} < values[{:d}]={:.8e}\n",
                    indices[i], values[indices[i]], indices[i-1], values[indices[i-1]]
                );
                throw std::runtime_error(err_msg);
            }
            #endif
        }

        // printf("quantiles(): created cumulative weights.\n");

        double total_weight = cum_weights[num_vals - 1];
        double* quantiles = new double[num_percs];
        double target;
        int idx;
        for (int i = 0; i < num_percs; i++) {
            target = percs[i] * total_weight;
            idx = 0;
            while (idx < num_vals && cum_weights[idx] < target) {
                idx++;
            }
            if (idx == 0) {
                quantiles[i] = values[indices[0]];
            } else if (idx == num_vals) {
                quantiles[i] = values[indices[num_vals - 1]];
            } else {
                quantiles[i] = (
                    values[indices[idx - 1]] + (
                        (target - cum_weights[idx - 1]) *
                        (values[idx] - values[idx - 1]) /
                        (cum_weights[idx] - cum_weights[idx - 1])
                    )
                );
            }
        }

        // printf("quantiles(): calculated quantiles.\n");

        delete[] indices;
        delete[] cum_weights;
        return quantiles;
    };
    */


    // =========================================================================
    // ====     HDF5 functions    ====
    // =========================================================================


    hid_t open_or_create_group(hid_t parent, const char* group_name) {
        // Suppress HDF5 automatic error printing
        H5E_auto2_t old_func;
        void* old_client_data;
        H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
        H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

        // Try to open the group
        hid_t group = H5Gopen(parent, group_name, H5P_DEFAULT);

        // Restore original error printing
        H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

        // If group didn't exist, create it
        if (group < 0) {
            group = H5Gcreate(parent, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }

        return group;  // Caller must close with H5Gclose()
    }

}

/*
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
    percs : (M,) double
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


def trapz_loglog( yy, xx, axis=-1, dlogx=None, lntol=1e-2):
    """Calculate integral, given `y = dA/dx` or `y = dA/dlogx` w/ trapezoid rule in log-log space.

    We are calculating the integral `A` given sets of values for `y` and `x`.
    To associate `yy` with `dA/dx` then `dlogx = None` [default], otherwise,
    to associate `yy` with `dA/dlogx` then `dlogx = True` for natural-logarithm, or `dlogx = b`
    for a logarithm of base `b`.

    For each interval (x[i+1], x[i]), calculate the integral assuming that y is of the form,
        `y = a * x^gamma`

    Parameters
    ----------
    yy : ndarray
    xx : (X,) array_like of scalar,
    bounds : (2,) array_like of scalar,
    axis : int,
    dlogx : scalar or None,
    lntol : scalar,

    Returns
    -------
    integ

    """

    if np.ndim(yy) != np.ndim(xx):
        if (np.ndim(xx) != 1) or (np.size(xx) != np.shape(yy)[axis]):
            raise ValueError(f"Mismatch between shapes of `xx` and `yy` ({np.shape(xx)=}, {np.shape(yy)=})")

        # convert `xx` from shape (N,) to (1, ... N, ..., 1) where all
        # dimensions besides `axis` have length one
        cut = [np.newaxis for ii in range(np.ndim(yy))]
        cut[axis] = slice(None)
        xx = xx[tuple(cut)]

    log_base = np.e
    # If `dlogx` is True, then we're using log-base-e (i.e. natural-log)
    # Otherwise, set the log-base to the given value
    if dlogx not in [None, True]:
        log_base = dlogx
    elif dlogx is None:
        dlogx = False

    # Numerically calculate the local power-law index
    delta_logx = np.diff(np.log(xx), axis=axis)
    gamma = np.diff(np.log(yy), axis=axis) / delta_logx
    xx = np.moveaxis(xx, axis, 0)
    yy = np.moveaxis(yy, axis, 0)
    # aa = np.mean([xx[:-1] * yy[:-1], xx[1:] * yy[1:]], axis=0)
    # aa = np.moveaxis(aa, 0, axis)
    # xx = np.moveaxis(xx, 0, axis)
    # yy = np.moveaxis(yy, 0, axis)

    # Integrate dA/dx   ::   A = (x1*y1 - x0*y0) / (gamma + 1)
    if (dlogx is False):
        dz = np.diff(yy * xx, axis=0)
        integ = dz / (gamma + 1)
        # when the power-law is (near) '-1' then, `A = a * log(x1/x0)`
        idx = np.isclose(gamma, -1.0, atol=lntol, rtol=lntol)

    # Integrate dA/dlogx    ::    A = (y1 - y0) / gamma
    else:
        dy = np.diff(yy, axis=0)
        integ = dy / gamma
        # when the power-law is (near) '-1' then, `A = a * log(x1/x0)`
        idx = np.isclose(gamma, 0.0, atol=lntol, rtol=lntol)

    if np.any(idx):
        aa = np.mean([xx[:-1] * yy[:-1], xx[1:] * yy[1:]], axis=0)
        # if `xx.shape != yy.shape` then `delta_logx` should be shaped (N-1, 1, 1, 1...)
        # broadcast `delta_logx` to the same shape as `idx` in this case
        if np.shape(xx) != np.shape(yy):
            delta_logx = delta_logx * np.ones_like(aa)
        integ[idx] = aa[idx] * delta_logx[idx]

    integ = integ / np.log(log_base)
    xx = np.moveaxis(xx, 0, axis)
    yy = np.moveaxis(yy, 0, axis)
    integ = np.moveaxis(integ, 0, axis)
    return integ


*/


#endif  // UTILS_H

/**
 *
 */

#pragma once

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "hdf5.h"

#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"

#include "quill/sinks/FileSink.h"
#include "quill/sinks/ConsoleSink.h"

using namespace std;



#define PATH_LOG_OUTPUT "logs/main.txt"


inline quill::Logger* get_logger()
{
    static quill::Logger* logger = []() -> quill::Logger* {
        // Start backend thread (only once, thread-safe)
        quill::Backend::start();

        // --- Console Sink ---
        auto console_sink = quill::Frontend::create_or_get_sink<quill::ConsoleSink>("console_sink");

        // Optional: restrict console to warnings+
        console_sink->set_log_level_filter(quill::LogLevel::Info);

        // --- File Sink ---
        auto file_sink = quill::Frontend::create_or_get_sink<quill::FileSink>(
            PATH_LOG_OUTPUT,
            [] {
                quill::FileSinkConfig cfg;
                cfg.set_open_mode('w');  // or 'a' for append
                // cfg.set_filename_append_option(quill::FilenameAppendOption::StartDateTime);
                cfg.set_filename_append_option(quill::FilenameAppendOption::None);
                return cfg;
            }(),
            quill::FileEventNotifier{}
        );

        file_sink->set_log_level_filter(quill::LogLevel::Debug);

        // --- Combine sinks into logger ---
        std::vector<std::shared_ptr<quill::Sink>> sinks = {console_sink, file_sink};

        auto* logger = quill::Frontend::create_or_get_logger(
            "global_logger", std::move(sinks),
            quill::PatternFormatterOptions{
                "%(time) [%(thread_id)] %(short_source_location:<28) "
                "%(log_level:<7) %(message)",
                "%H:%M:%S.%Qus"
            }
        );

        // NOTE: `logger` level needs to be lower than any desired sink levels!
        logger->set_log_level(quill::LogLevel::Debug);

        //! FIX: NOT WORKING
        // LOG_INFO(logger, "Logging to file '%s'.", file_sink->get_filename());

        return logger;
    }();

    return logger;
}






namespace utils {

    template <typename T>
    int* argsort(const T* array, size_t size) {
        int* indices = new int[size];
        for (int i = 0; i < size; ++i)
            indices[i] = i;

        std::sort(
            indices, indices + size,
            [&](size_t i, size_t j) { return array[i] < array[j]; }
        );

        for (int i = 0; i < size; ++i) {
            printf("%d: array[%d]=%.2e\n", i, indices[i], array[indices[i]]);
        }

        return indices;
    }


    void index_2d_to_1d(int i, int j, int dim1, int dim2, int* index);

    void index_3d_to_1d(int i, int j, int k, int dim1, int dim2, int dim3, int* index);

    void index_1d_to_3d(int index, int dim1, int dim2, int dim3, int* i, int* j, int* k);

    void index_1d_to_2d(int index, int dim1, int dim2, int* i, int* j);

    bool is_almost_equal(double a, double b, double atol = 1E-8, double rtol = 1E-6);

    // !NOT WORKING!
    double* quantiles(
        double* values, int num_vals, double* percs, int num_percs,
        double* weights = nullptr, bool values_sorted = false
    );


    // =========================================================================
    // ====     HDF5 functions    ====
    // =========================================================================


    template <typename T> hid_t get_hdf5_type();
    template <> inline hid_t get_hdf5_type<int>() { return H5T_NATIVE_INT; }
    template <> inline hid_t get_hdf5_type<float>() { return H5T_NATIVE_FLOAT; }
    template <> inline hid_t get_hdf5_type<double>() { return H5T_NATIVE_DOUBLE; }
    template <> inline hid_t get_hdf5_type<long>() { return H5T_NATIVE_LONG; }

    hid_t hdf5_open_or_create_group(hid_t parent, const char* group_name);


    template <typename T>
    void hdf5_write_scalar(hid_t h5_file, const char* group_name, const char* name, const T& value) {
        // printf("hdf5_write_scalar(%s/%s)\n", group_name, name);
        static_assert(
            std::is_arithmetic<T>::value,
            "hdf5_write_scalar only supports primitive numeric types"
        );

        hid_t dtype = get_hdf5_type<T>();
        hid_t group = hdf5_open_or_create_group(h5_file, group_name);
        hid_t space = H5Screate(H5S_SCALAR);
        hid_t dset = H5Dcreate(
            group, name, dtype, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );

        // printf("\thdf5_write_scalar(%s/%s): writing...\n", group_name, name);
        H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);

        H5Dclose(dset);
        H5Sclose(space);
        H5Gclose(group);
    }


    template <typename T>
    void hdf5_write_array1d(
        hid_t h5_file, const char* group_name, const char* dataset_name,
        const T* data, const int data_size
    ) {
        // printf("hdf5_write_array1d(%s/%s)\n", group_name, dataset_name);
        static_assert(
            std::is_arithmetic<T>::value,
            "hdf5_write_scalar only supports primitive numeric types"
        );

        hid_t dtype = get_hdf5_type<T>();
        hid_t group = hdf5_open_or_create_group(h5_file, group_name);
        hsize_t dims[1] = {static_cast<hsize_t>(data_size)};
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t dset = H5Dcreate(
            group, dataset_name, dtype, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );

        // printf("\thdf5_write_array1d(%s/%s): writing...\n", group_name, dataset_name);
        H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

        // Close the dataset and group
        H5Dclose(dset);
        H5Sclose(space);
        H5Gclose(group);
    }


    template <typename T>
    void hdf5_write_array2d(
        hid_t h5_file, const char* group_name, const char* dataset_name,
        T** data, const int xsize, const int ysize
    ) {
        // printf("hdf5_write_array2d(%s/%s)\n", group_name, dataset_name);
        static_assert(
            std::is_arithmetic<T>::value,
            "hdf5_write_array2d only supports primitive numeric types"
        );

        hid_t dtype = get_hdf5_type<T>();
        hid_t group;
        if (group_name != nullptr) {
            group = hdf5_open_or_create_group(h5_file, group_name);
        } else {
            group = h5_file;
        }
        hsize_t dims[2] = {static_cast<hsize_t>(xsize), static_cast<hsize_t>(ysize)};
        hid_t space = H5Screate_simple(2, dims, nullptr);
        hid_t dset = H5Dcreate(
            group, dataset_name, dtype, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );

        // Flatten the 2D array into a 1D array
        T* data_flat = new T[xsize * ysize];
        int idx;
        for (int i = 0; i < xsize; ++i) {
            for (int j = 0; j < ysize; ++j) {
                index_2d_to_1d(i, j, xsize, ysize, &idx);
                data_flat[idx] = data[i][j];
            }
        }

        // Write the flattened data to the dataset
        // printf("\thdf5_write_array2d(%s/%s): writing...\n", group_name, dataset_name);
        H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_flat);

        // Close the dataset and group
        H5Dclose(dset);
        H5Sclose(space);
        if (group_name != nullptr) {
            H5Gclose(group);
        }
        delete[] data_flat;

    }


    template <typename T>
    void hdf5_write_array3d(
        hid_t h5_file, const char* group_name, const char* dataset_name,
        T*** data, const int xsize, const int ysize, const int zsize
    ) {
        // printf("hdf5_write_array3d(%s/%s) : shape %d,%d,%d\n", group_name, dataset_name, xsize, ysize, zsize);

        static_assert(
            std::is_arithmetic<T>::value,
            "hdf5_write_array3d only supports primitive numeric types"
        );

        hid_t dtype = get_hdf5_type<T>();
        hid_t group;
        if (group_name != nullptr) {
            group = hdf5_open_or_create_group(h5_file, group_name);
        } else {
            group = h5_file;
        }
        hsize_t dims[3] = {static_cast<hsize_t>(xsize), static_cast<hsize_t>(ysize), static_cast<hsize_t>(zsize)};
        hid_t space = H5Screate_simple(3, dims, nullptr);
        hid_t dset = H5Dcreate(
            group, dataset_name, dtype, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );

        // Flatten the 3D array into a 1D array
        T* data_flat = new T[xsize * ysize * zsize];
        int idx;
        for (int i = 0; i < xsize; ++i) {
            for (int j = 0; j < ysize; ++j) {
                for (int k = 0; k < zsize; ++k) {
                    index_3d_to_1d(i, j, k, xsize, ysize, zsize, &idx);
                    data_flat[idx] = data[i][j][k];
                }
            }
        }

        // Write the flattened data to the dataset
        // printf("\thdf5_write_array3d(%s/%s): writing...\n", group_name, dataset_name);
        H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_flat);

        // Close the dataset and group
        H5Dclose(dset);
        H5Sclose(space);
        if (group_name != nullptr) {
            H5Gclose(group);
        }
        delete[] data_flat;
    }

}


class H5Slice_4D {
    public:
    // hid_t file;
    const char* data_name;
    hsize_t dims[4];
    hid_t dspace;
    hid_t dset;
    int slice_size;

    H5Slice_4D(hid_t h5_file, const char* data_name, int size0, int size1, int size2, int size3) {
        // this->file = h5_file;
        this->data_name = data_name;
        dims[0] = static_cast<hsize_t>(size0);
        dims[1] = static_cast<hsize_t>(size1);
        dims[2] = static_cast<hsize_t>(size2);
        dims[3] = static_cast<hsize_t>(size3);
        slice_size = size1 * size2 * size3;

        dspace = H5Screate_simple(4, dims, nullptr);
        dset = H5Dcreate(
            h5_file, data_name, H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
        );
        // H5Sclose(dspace);  // dataset now owns its space, so the dataspace can be closed
    }
    ~H5Slice_4D() {
        H5Sclose(dspace);
        H5Dclose(dset);
    }

    void write_slice_at(int idx0, double*** data_slice, int size1, int size2, int size3) {

        // printf("H5Slice_4D(%s).write_slice_at(%d)\n", data_name, idx0);

        hsize_t h5_offset[4] = {static_cast<hsize_t>(idx0), 0, 0, 0};

        if (size1*size2*size3 != slice_size)
            throw std::runtime_error("INCONSISTENT SLICE!  total size does not match `slice_size`!");
        if ((size1 != dims[1]) || (size2 != dims[2]) || (size3 != dims[3]))
            throw std::runtime_error("INCONSISTENT SLICE!  `size` does not match `dims`!");
        if ((idx0 >= dims[0]))
            throw std::runtime_error("Cannot write at slice beyond dims[0]!");

        double* data_flat = (double *)malloc(slice_size * sizeof(double));
        int idx;
        for (int i = 0; i < size1; i++) {
            for (int j = 0; j < size2; j++) {
                for (int k = 0; k < size3; k++) {
                    utils::index_3d_to_1d(i, j, k, size1, size2, size3, &idx);
                    data_flat[idx] = data_slice[i][j][k];
                }
            }
        }

        hsize_t h5_count[4] = { 1, dims[1], dims[2], dims[3] };

        hid_t filespace = H5Dget_space(dset);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, h5_offset, nullptr, h5_count, nullptr);

        hid_t memspace = H5Screate_simple(4, h5_count, nullptr);
        // printf("\tH5Slice_4D(%s).write_slice_at(%d): writing...\n", data_name, idx0);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, data_flat);

        H5Sclose(memspace);
        H5Sclose(filespace);
        free(data_flat);
    }

};


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

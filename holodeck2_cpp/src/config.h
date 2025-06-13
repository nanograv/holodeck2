/**
 *
 */

#pragma once


// ---- Required variables

#define FNAME_OUTPUT_HDF5 "output.hdf5"
#define PATH_LOG_OUTPUT "logs/main.txt"

// skip SAM bins with expected-numbers of binaries less than this value divided by the number of bin-realizations
// e.g. if there are 1e3 bins, 1e2 reals, and SAM_NUM_EXPECT_FLOOR=1e-2, then skip values below 1e-7
constexpr double SAM_NUM_EXPECT_FLOOR = 1.0E-2;


// ---- Optional behavior

// Maximum cutoff value for `tauf` (in years) for number-density calculations.
#define MAX_TAUF_YR 1.0E10    // [yr]


// ---- Performance / optimization


// ---- Output

#define HDF5_OUTPUT
#ifdef    HDF5_OUTPUT
#define HDF5_OUTPUT_DETAILS  // save additional arrays to HDF5 output file
#endif // HDF5_OUTPUT


// ---- Debugging / diagnostic

#define DEBUG

#ifdef DEBUG
#define DEBUG_FREQ_STATS   // calculate & print per-frequency-bin statistics
#endif // DEBUG


/**
 *
 */

#pragma once

// #define DEBUG

#define HDF5_OUTPUT
#ifdef HDF5_OUTPUT
// #define HDF5_OUTPUT_DETAILS  // save additional arrays to HDF5 output file
#endif  // HDF5_OUTPUT

#ifdef DEBUG
// #define DEBUG_FREQ_STATS   // calculate & print per-frequency-bin statistics
#endif // DEBUG

#define FNAME_OUTPUT_HDF5 "output.hdf5"
#define PATH_LOG_OUTPUT "logs/main.txt"

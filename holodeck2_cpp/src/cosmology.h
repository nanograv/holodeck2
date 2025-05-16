/**
 *
 */

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <format>
#include <fstream>
#include <filesystem>
#include <sstream>   // for std::istringstream
#include <string>

#include "utils.h"
#include "constants.h"


class Cosmology {
public:
    Cosmology();

    ~Cosmology();

    // Comoving Distance [Mpc]
    void dcom_from_redz(double redz, double* dcom) { value_from_redz(redz, dcom_log10, dcom); }

    // E-function E(z)=H(z)/H0  [unitless]
    void efunc_from_redz(double redz, double* efunc) { value_from_redz(redz, efunc_log10, efunc); }

    // Comoving Volume [Mpc^3]
    void vcom_from_redz(double redz, double* vcom) { value_from_redz(redz, vcom_log10, vcom); }

    // Comoving Volume [Myr]
    void tlook_from_redz(double redz, double* tlook) { value_from_redz(redz, tlook_log10, tlook); }

    // ---- Derived parameters

    void dtdz_from_redz(double redz, double* dtdz) {
        // dt/dz = -1 / (H0 * E(z))
        double efunc;
        efunc_from_redz(redz, &efunc);
        *dtdz = H0_inv / (efunc * (1.0 + redz));
    }


private:
    // Cosmology parameters (loaded from input file at initialization)
    double H0;            // Hubble constant at z=0   [km/s/Mpc]
    double Omega_m0;      // Matter density parameter (z=0)
    double Omega_l0;      // Dark energy density parameter (z=0)
    double Omega_b0;      // Baryon density parameter (z=0)

    // Pre-calculated arrays loaded in from file at initialization
    int grid_size;        // Size of the interpolation grid
    double *redz_log10;   // Redshift (on regular log-spaced grid)
    double *scafa_log10;  // Scale-factor a = 1/(1+z)
    double *dcom_log10;   // Comoving distance [Mpc]
    double *vcom_log10;   // Comoving volume [Mpc^3]
    double *tlook_log10;  // Lookback time [Myr]
    double *efunc_log10;  // E(z) = H(z)/H0

    // Derived parameters
    double H0_cgs;
    double H0_inv;
    double dlog10z;

    int find_lower_index_uniform_redz_log10(double zlog10);

    double lin_interp_at_index(double target, int idx, double* xx, double* yy);

    void load_cosmo_data(const std::string& path);

    // int brute_force_lower_index(double target, double* xx) {
    //     int idx = 0;
    //     while (idx < grid_size && xx[idx] < target) {
    //         idx++;
    //     }
    //     return idx - 1;
    // }

    void value_from_redz(double redz, double* values_log10, double* out_val);

};


#endif  // COSMOLOGY_H

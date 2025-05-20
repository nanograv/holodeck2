/**
 *
 */

#include "cosmology.h"


namespace fs = std::filesystem;

namespace {
    // path to cosmology data file (with pre-calculated grid of cosmographic parameters)
    fs::path full_path = fs::path(PATH_DATA_DIR) / "cosmology.txt";
}


Cosmology::Cosmology() {
    printf("Loading cosmology data from %s\n", full_path.c_str());

    load_cosmo_data(full_path.string());
    H0_cgs = H0 * KM_S_MPC;    // [1/sec]
    H0_inv = 1.0 / H0_cgs;     // [sec]

    printf(
        "Loaded grid with %d points, redshifts between [%.2e, %.2e]\n",
        grid_size, pow(10.0, redz_log10[0]), pow(10.0, redz_log10[grid_size-1])
    );
};


Cosmology::~Cosmology() {
    delete[] redz_log10;
    delete[] scafa_log10;
    delete[] dcom_log10;
    delete[] vcom_log10;
    delete[] tlook_log10;
    delete[] efunc_log10;
};


int Cosmology::find_lower_index_uniform_redz_log10(double zlog10) {
    int idx = (int)floor((zlog10 - redz_log10[0]) / dlog10z);
    #ifdef DEBUG
    if (idx < 0 || idx >= grid_size) {
        std::string msg = std::format(
            "Index out of range for interpolation! idx={:d} vs. [0, {:d}]",
            idx, grid_size-1
        );
        throw std::out_of_range(msg);
    }
    #endif
    return idx;
}

double Cosmology::lin_interp_at_index(double target, int idx, double* xx, double* yy) {
    // Linear interpolation in log10 space
    double x1 = xx[idx];
    double y1 = yy[idx];
    #ifdef DEBUG
    if (idx < 0 || idx >= grid_size - 1) {
        std::string msg = std::format(
            "Index out of range for interpolation! idx={:d} vs. [0, {:d}]",
            idx, grid_size-1
        );
        throw std::out_of_range(msg);
    }
    if (target < x1 || target > xx[idx + 1]) {
        std::string msg = std::format(
            "Target value out of range for interpolation! x=[{:.8e}, {:.8e}] vs. {:.8e}",
            x1, xx[idx+1], target
        );
        throw std::out_of_range(msg);
    }
    #endif
    return y1 + (yy[idx + 1] - y1) * (target - x1) / (xx[idx + 1] - x1);
}

void Cosmology::load_cosmo_data(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Failed to open file: " + path);

    std::string line;

    // --- First line: metadata
    std::getline(file, line);
    std::istringstream meta(line);
    std::string token;

    while (meta >> token) {
        if (token.find("H0=") == 0)       H0  = std::stof(token.substr(3));
        else if (token.find("Om0=") == 0) Omega_m0 = std::stof(token.substr(4));
        else if (token.find("Ob0=") == 0) Omega_b0 = std::stof(token.substr(4));
        else if (token.find("Ol0=") == 0) Omega_l0 = std::stof(token.substr(4));
        else if (token.find("N=") == 0)   grid_size = std::stoi(token.substr(2));
    }

    if (grid_size <= 0) throw std::runtime_error("Failed to parse N from metadata");

    // Allocate arrays
    redz_log10 = new double[grid_size];
    scafa_log10    = new double[grid_size];
    dcom_log10     = new double[grid_size];
    vcom_log10     = new double[grid_size];
    tlook_log10    = new double[grid_size];
    efunc_log10    = new double[grid_size];

    // --- Skip column header line
    std::getline(file, line);

    // --- Read data lines
    double temp;
    for (int i = 0; i < grid_size; i++) {
        if (!std::getline(file, line)) throw std::runtime_error("Unexpected EOF before N lines");

        std::istringstream ss(line);
        double z, a, d, v, t, e;
        if (!(ss >> z >> a >> d >> v >> t >> e)) {
            throw std::runtime_error("Failed to parse line: " + line);
        }

        redz_log10[i]  = log10(z);
        scafa_log10[i] = log10(a);
        dcom_log10[i]  = log10(d);
        vcom_log10[i]  = log10(v);
        tlook_log10[i] = log10(t);
        efunc_log10[i] = log10(e);

        // Check for evenly spaced grid in log10(z)
        if (i > 0) {
            temp = redz_log10[i] - redz_log10[i-1];
            if (i == 1) {
                dlog10z = temp;
            } else if (i > 1){
                if (!utils::is_almost_equal(temp, dlog10z)) {
                    std::string msg = std::format(
                        "Grid is not evenly spaced in log10(z)!  {:.8e} vs. {:.8e} at i={:d}\n",
                        dlog10z, temp, i
                    );
                    throw std::runtime_error(msg);
                }
            }
        }

    } // i

}

void Cosmology::value_from_redz(double redz, double* values_log10, double* out_val) {
    double zlog10 = log10(redz);
    int idx = find_lower_index_uniform_redz_log10(zlog10);
    *out_val = pow(10.0, lin_interp_at_index(zlog10, idx, redz_log10, values_log10));
}


/**
 *
 */

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib>

namespace fs = std::filesystem;

fs::path full_path = fs::path(PATH_DATA_DIR) / "cosmology.txt";

// const std::string data_file_path = std::string(PATH_DATA_DIR) + "/my_data.txt";


class Cosmology {
public:
    Cosmology() {
        printf("Loading cosmology data from %s\n", full_path.c_str());
        load_cosmo_data(full_path.string());
        printf(
            "Loaded grid with %d points, redshifts between [%.2e, %.2e]\n",
            grid_size, redshift[0], redshift[grid_size-1]
        );
    };




    ~Cosmology() {
        delete[] redshift;
        delete[] scafa;
        delete[] dcom;
        delete[] vcom;
        delete[] tlook;
        delete[] efunc;
    }

private:
    // Cosmology parameters
    float H0;       // Hubble constant at z=0   [km/s/Mpc]
    float Omega_m0; // Matter density parameter
    float Omega_l0; // Dark energy density parameter
    float Omega_b0; // Baryon density parameter
    int grid_size;  // Size of the interpolation grid

    // Derived parameters
    // float H0_cgs;
    // float H0_inv;
    float dlog10z;
    float *redshift;
    float *scafa;
    float *dcom;
    float *vcom;
    float *tlook;
    float *efunc;

    void load_cosmo_data(const std::string& path) {
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
        redshift = new float[grid_size];
        scafa    = new float[grid_size];
        dcom     = new float[grid_size];
        vcom     = new float[grid_size];
        tlook    = new float[grid_size];
        efunc    = new float[grid_size];

        // --- Skip column header line
        std::getline(file, line);

        // --- Read data lines
        for (int i = 0; i < grid_size; i++) {
            if (!std::getline(file, line)) throw std::runtime_error("Unexpected EOF before N lines");

            std::istringstream ss(line);
            float z, a, d, v, t, e;
            if (!(ss >> z >> a >> d >> v >> t >> e)) {
                throw std::runtime_error("Failed to parse line: " + line);
            }

            redshift[i] = z;
            scafa[i]    = a;
            dcom[i]     = d;
            vcom[i]     = v;
            tlook[i]    = t;
            efunc[i]    = e;

            if (i == 1) {
                dlog10z = log10(redshift[i]) - log10(redshift[i-1]);
            } else if (i > 1) {

            }

        } // i

    }
};


#endif  // COSMOLOGY_H

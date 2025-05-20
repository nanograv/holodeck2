/**
 *
 */


#ifndef SAM_H
#define SAM_H

#include <cstdlib> // For malloc and free
#include <cmath>
#include <random>

#include "constants.h"
#include "physics.h"
#include "cosmology.h"

constexpr double FOUR_PI_C_MPC = 4 * PI * SPLC / MPC;    // 4*pi*c   [Mpc/s]

inline std::default_random_engine rng(42);          // RNG with fixed seed


class PTA {

public:
    double obs_dur;
    int num_freq_cents;
    double* fobs_cents;
    double* fobs_edges;

    PTA(double dur=20.0, double nfreqs=30) : obs_dur(dur), num_freq_cents(nfreqs) {
        fobs_cents = (double*)malloc(num_freq_cents * sizeof(double));
        fobs_edges = (double*)malloc((num_freq_cents+1) * sizeof(double));
        double df = 1.0 / obs_dur;
        fobs_edges[0] = df * 0.5;
        for (int ii = 0; ii < num_freq_cents; ii++) {
            fobs_cents[ii] = df * (ii + 1);
            fobs_edges[ii+1] = df * (ii + 1.5);
        }
    };

    ~PTA() {
        free(fobs_cents);
        free(fobs_edges);
    };

};


class GravWaves {
public:
    int num_freq_cents;
    int num_reals;
    int num_louds;

    double** gwb;
    double*** cws;

    GravWaves(PTA *pta, int nreals = 100, int nloud = 5) : pta(pta), num_reals(nreals), num_louds(nloud) {
        this->pta = pta;
        num_freq_cents = pta->num_freq_cents;
        gwb = (double** )calloc(num_freq_cents, sizeof(double*));
        cws = (double*** )calloc(num_freq_cents, sizeof(double**));
        for (int ii = 0; ii < num_freq_cents; ii++) {
            gwb[ii] = (double* )calloc(num_reals, sizeof(double));
            cws[ii] = (double** )calloc(num_reals, sizeof(double*));
            for (int jj = 0; jj < num_reals; jj++) {
                cws[ii][jj] = (double* )calloc(num_louds, sizeof(double));
            }
        }
    }

    ~GravWaves() {
        for (int ii = 0; ii < pta->num_freq_cents; ii++) {
            for (int jj = 0; jj < num_reals; jj++) {
                free(cws[ii][jj]);
            }
            free(cws[ii]);
            free(gwb[ii]);
        }
        free(cws);
        free(gwb);
    };

private:
    PTA* pta;

};


class SAM {
public:
    int num_mass_edges;
    int num_redz_edges;
    double *mass_edges;
    double *redz_edges;
    double *mass_cents;
    double *redz_cents;

//! FIX: make these private
// private:
    // Grid
    // double mass_pars[3] = {1E6, 1E12, 101};
    // double redz_pars[3] = {1E-2, 1E1, 91};
    double mass_pars[3] = {1E6, 1E12, 21};
    double redz_pars[3] = {1E-2, 1E1, 11};

    // GSMF - Double Schechter
    double gsmf_log10_phi1[3] = {-2.383, -0.264, -0.107};
    double gsmf_log10_phi2[3] = {-2.818, -0.368, +0.046};
    double gsmf_log10_mchar[3] = {+10.767, +0.124, -0.033};
    double gsmf_alpha1 = -0.28;
    double gsmf_alpha2 = -1.48;

    // GMR - Illustris; Galaxy Merger Rate
    double gmr_norm0_log10 = -2.2287;          // -2.2287 ± 0.0045    A0 [log10(A*Gyr)]
    double gmr_normz       = +2.4644;          // +2.4644 ± 0.0128    eta
    double gmr_malpha0     = +0.2241;          // +0.2241 ± 0.0038    alpha0
    double gmr_malphaz     = -1.1759;          // -1.1759 ± 0.0316    alpha1
    double gmr_mdelta0     = +0.7668;          // +0.7668 ± 0.0202    delta0
    double gmr_mdeltaz     = -0.4695;          // -0.4695 ± 0.0440    delta1
    double gmr_qgamma0     = -1.2595;          // -1.2595 ± 0.0026    beta0
    double gmr_qgammaz     = +0.0611;          // +0.0611 ± 0.0021    beta1
    double gmr_qgammam     = -0.0477;          // -0.0477 ± 0.0013    gamma
    double gmr_mref_delta  = 2.0E11;           // fixed value [Msol]
    double gmr_mref_gamma  = 1.0E10;           // fixed value [Msol]

    // M-MBulge
    double mmbulge_mass_amp_log10 = 8.69;      // 8.69 ± 0.05  [log10(M/Msol)]  approximate uncertainties!
    double mmbulge_mass_ref = 1E11;            // 1E11 Msol
    double mmbulge_plaw = 1.17;                // 1.17 ± 0.08
    double mmbulge_scatter_dex = 0.28;         // scatter stdev in dex
    double mmbulge_bulge_frac = 0.7;

    SAM(Cosmology* cosmo) {
        this->cosmo = cosmo;

        num_mass_edges = mass_pars[2];
        num_redz_edges = redz_pars[2];
        mass_cents = (double*)malloc(num_mass_edges * sizeof(double));
        redz_cents = (double*)malloc(num_redz_edges * sizeof(double));
        mass_edges = (double*)malloc((num_mass_edges+1) * sizeof(double));
        redz_edges = (double*)malloc((num_redz_edges+1) * sizeof(double));

        // ---- Initialize mass edges and centers

        double dlog10m = (log10(mass_pars[1]) - log10(mass_pars[0])) / num_mass_edges;
        double next, prev = log10(mass_pars[0]);
        mass_edges[0] = pow(10.0, prev);
        // printf("mass_edges[0]=%.2e, dlog10m=%.2e\n", mass_edges[0], dlog10m);
        for (int ii = 0; ii < num_mass_edges; ii++) {
            next = prev + dlog10m;
            mass_edges[ii+1] = pow(10.0, next);
            mass_cents[ii] = pow(10.0, 0.5 * (prev + next));
            // printf("mass_edges[%d]=%.2e, mass_cents[%d]=%.2e (next=%.2e)\n", ii, mass_edges[ii], ii, mass_cents[ii], next);
            prev = next;
        }

        // ---- Initialize redshift edges and centers

        double dlog10z = (log10(redz_pars[1]) - log10(redz_pars[0])) / num_redz_edges;
        prev = log10(redz_pars[0]);
        redz_edges[0] = pow(10.0, prev);
        // printf("redz_edges[0]=%.2e, dlog10z=%.2e\n", redz_edges[0], dlog10z);
        for (int ii = 0; ii < num_redz_edges; ii++) {
            next = prev + dlog10z;
            redz_edges[ii+1] = pow(10.0, next);
            redz_cents[ii] = pow(10.0, 0.5 * (prev + next));
            // printf("redz_edges[%d]=%.2e, redz_cents[%d]=%.2e (next=%.2e)\n", ii, redz_edges[ii], ii, redz_cents[ii], next);
            prev = next;
        }

    }

    double dngal_to_dnmbh(double mbh1, double mbh2, double mstar_tot, double mstar_rat, double rz);

    double gmr_illustris(double mstar_tot, double mstar_rat, double rz);

    void grav_waves(PTA &pta, GravWaves &gw);

    double gsmf_double_schechter(double mstar, double redz);

    double gsmf_single_schechter(double mstar, double redz, double phi_terms[3], double mchar_terms[3], double alpha);

    double mmbulge_mstar_from_mbh(double mbh);

    double number_density(double mbh1, double mbh2, double rz, bool return_galgal = false);

    ~SAM() {
        free(mass_cents);
        free(redz_cents);
        free(mass_edges);
        free(redz_edges);
    };

private:
    Cosmology* cosmo;

};


#endif  // SAM_H

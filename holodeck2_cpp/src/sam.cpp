/**
 *
 */

#include "sam.h"
#include "physics.h"
#include "cosmology.h"


double SAM::dngal_to_dnmbh(double mbh1, double mbh2, double ms1, double ms2, double rz) {
    // # we want ``dn_mbhb / [dlog10(M_bh) dq_bh qz]``
    // # so far we have ``dn_gal / [dlog10(M_gal) dq_gal dz]``

    // # dn / [dM dq dz] = (dn_gal / [dM_gal dq_gal dz]) * (dM_gal/dM_bh) * (dq_gal / dq_bh)
    double dqbh_dqgal = mmbulge_plaw * pow(ms2/ms1, mmbulge_plaw - 1.0);
    // # (dMstar-pri / dMbh-pri) * (dMbh-pri/dMbh-tot) = (dMstar-pri / dMstar-tot) * (dMstar-tot/dMbh-tot)
    // # ==> (dMstar-tot/dMbh-tot) = (dMstar-pri / dMbh-pri) * (dMbh-pri/dMbh-tot) / (dMstar-pri / dMstar-tot)
    // #                           = (dMstar-pri / dMbh-pri) * (1 / (1+q_bh)) / (1 / (1+q_star))
    // #                           = (dMstar-pri / dMbh-pri) * ((1+q_star) / (1+q_bh))
    double dmstar_dmbh_pri = (ms1 / (mmbulge_plaw * mbh1));
    // #! ==== FIX: this allows for q>1 ! ====
    double qterm = (1.0 + ms2/ms1) / (1.0 + (mbh2/mbh1));
    return ((mbh1 + mbh2) / (ms1 + ms2)) * (dmstar_dmbh_pri * qterm / dqbh_dqgal);
}


double SAM::gmr_illustris(double mstar_tot, double mstar_rat, double rz) {
    double zplus1 = rz + 1.0;
    double norm = (pow(10.0, gmr_norm0_log10) / GYR) * pow(zplus1, gmr_normz);
    double malpha = gmr_malpha0 * pow(zplus1, gmr_malphaz);
    double mdelta = gmr_mdelta0 * pow(zplus1, gmr_mdeltaz);
    double qgamma = gmr_qgamma0 * pow(zplus1, gmr_qgammaz) + gmr_qgammam * log10(mstar_tot/gmr_mref_gamma);

    return (
        norm *
        pow(mstar_tot/gmr_mref_gamma, malpha) *
        pow(1.0 + mstar_tot/gmr_mref_delta, mdelta) *
        pow(mstar_rat, qgamma)
    );
}


double SAM::gsmf_single_schechter(double mstar, double redz, double phi_terms[3], double mchar_terms[3], double alpha) {
    double rz2 = pow(redz, 2);
    double mchar = MSOL * pow(10.0, mchar_terms[0] + mchar_terms[1] * redz + mchar_terms[2] * rz2);
    double mm = mstar / mchar;
    return (
        LOG10 *
        pow(10.0, phi_terms[0] + phi_terms[1] * redz + phi_terms[2] * rz2) *
        pow(mm, 1.0 + alpha) *
        exp(-mm)
    );
}


double SAM::gsmf_double_schechter(double mstar, double redz) {
    double phi = gsmf_single_schechter(mstar, redz, this->gsmf_log10_phi1, this->gsmf_log10_mchar, this->gsmf_alpha1);
    phi = phi +  gsmf_single_schechter(mstar, redz, this->gsmf_log10_phi2, this->gsmf_log10_mchar, this->gsmf_alpha2);
    return phi;
}


double SAM::mmbulge_mstar_from_mbh(double mbh) {
    // get bulge mass
    double mstar = (log10(mbh) - this->mmbulge_mass_amp_log10) / this->mmbulge_plaw;
    mstar = this->mmbulge_mass_ref * pow(10.0, mstar);
    // convert from bulge-mass to stellar-mass
    mstar = mstar / this->mmbulge_bulge_frac;
    return mstar;
}


double SAM::number_density(double mbh1, double mbh2, double rz, bool return_galgal) {

    // ---- Get galaxy stellar masses

    // NOTE: assume primary MBH ==> primary galaxy mass ==> GSMF
    double ms1 = mmbulge_mstar_from_mbh(mbh1);

    // ---- Get galaxy stellar mass function

    //! UNITS?
    double phi = gsmf_double_schechter(ms1, rz);

    // ---- Get galaxy-galaxy merger rate

    double ms2 = mmbulge_mstar_from_mbh(mbh1);
    double gmr = gmr_illustris(ms1 + ms2, ms2/ms1, rz);

    // ---- Get galaxy-galaxy merger, number density

    double dtdz;
    cosmo->dtdz_from_redz(rz, &dtdz);

    double ndens = phi * gmr * dtdz;

    // Return only galaxy-galaxy merger rate
    if (return_galgal) { return ndens; }

    // ---- Convert from galaxy-galaxy to MBH-MBH merger rate

    // we want ``dn_mbhb / [dlog10(M_bh) dq_bh qz]``
    // so far we have ``dn_gal / [dlog10(M_gal) dq_gal dz]``

    ndens = ndens * dngal_to_dnmbh(mbh1, mbh2, ms1, ms2, rz);

    // ---- Add scatter from the M-Mbulge relation

    //! ========================================================================
    //! FIX: still need to add scatter
    // ndens = add_scatter_to_masses(self.mtot, self.mrat, ndens, scatter, log=log)
    //! ========================================================================

    return ndens;
}


void SAM::grav_waves(PTA &pta, GravWaves &gw) {
    // Grid/domain sizes
    int num_mass = this->num_mass;
    int num_redz = this->num_redz;
    int num_freqs = pta.num_freqs;
    int num_reals = gw.num_reals;
    int num_louds = gw.num_louds;

    // Grid edges & centers
    double* fobs_edges = pta.fobs_edges;
    double* fobs_cents = pta.fobs_cents;
    double* mass_edges = this->mass_edges;
    double* mass_cents = this->mass_cents;
    double* redz_edges = this->redz_edges;
    double* redz_cents = this->redz_cents;

    // Values to be calculated/filled
    double** gwb = gw.gwb;   // [F][R]     : frequencies, realizations
    double*** cws = gw.cws;  // [F][R][L]  : frequencies, realizations, loudest

    // Temporary variables
    double hs2;                // spectral-strain-squared of an individual binary
    double tauf;               // GW Hardening timescale in frequency (f / df/dt)  [s]
    double numb_expect;        // number of binaries (expectation-value) in each bin (bin centers) []
    double mbh_mass_grams[num_mass];  // MBH mass [grams]
    double dlog10_mbh[num_mass-1];  // log10(mbh) bin-width
    double dz[num_redz-1];           // redshift bin-width
    double dlnf[num_freqs-1];
    double frst_orb[num_redz][num_freqs];   // Rest-frame orbital frequency [1/s]
    double dist_com_Mpc[num_redz]; // comoving distance [MPC]
    double dist_com_cm[num_redz];  // comoving distance [cm]
    double cosmo_fact[num_redz];   // comoving distance [cm]
    int numb_in_real;

    int m1, m2, z, f, r;           // loop iteration variables

    double** mchirp;               // [M1][M2] : chirp mass (rest-frame) [grams]
    double*** num_dens;            // [M1][M2][Z] : differential number-density of binaries [Mpc^-3]
    // double*** numb_expect;         // [M1][M2][Z] : number of binaries (expectation-value) in each bin (bin centers) []
    double*** diff_numb;          // [M1][M2][Z] : differential-number of binaries at grid edges [] (unitless, but differential!)

    // ---- Perform initializations & pre-calculations

    mchirp = (double**)malloc(num_mass * sizeof(double*));
    num_dens = (double***)malloc(num_mass * sizeof(double**));
    // numb_expect = (double***)malloc(num_mass * sizeof(double**));
    diff_numb = (double***)malloc(num_mass * sizeof(double**));

    for(m1 = 0; m1 < num_mass; m1++) {
        if (m1 > 0) { dlog10_mbh[m1-1] = log10(mass_edges[m1]) - log10(mass_edges[m1-1]); }

        mchirp[m1] = (double*)malloc(num_mass * sizeof(double));
        num_dens[m1] = (double**)malloc(num_mass * sizeof(double*));
        // numb_expect[m1] = (double**)malloc(num_mass * sizeof(double*));
        diff_numb[m1] = (double**)malloc(num_mass * sizeof(double*));

        for(m2 = 0; m2 < num_mass; m2++) {
            if (m1 == 0) { mbh_mass_grams[m2] = mass_edges[m2] * MSOL; }
            // `mchirp` [gram], Chirp-Mass (rest-frame)
            physics::chirp_mass_from_m1m2(mbh_mass_grams[m1], mbh_mass_grams[m2], &mchirp[m1][m2]);

            num_dens[m1][m2] = (double*)malloc(num_redz * sizeof(double));
            // numb_expect[m1][m2] = (double*)malloc(num_redz * sizeof(double));
            diff_numb[m1][m2] = (double*)malloc(num_redz * sizeof(double));
        }
    }

    for(z = 0; z < num_redz; z++) {
        if (z > 0) { dz[z-1] = redz_edges[z] - redz_edges[z-1]; }

        this->cosmo->dcom_from_redz(redz_edges[z], &dist_com_Mpc[z]);
        // convert from Mpc to cm
        dist_com_cm[z] = dist_com_Mpc[z] * MPC;

        // [Mpc^3/s] : this is (dVc/dz) * (dz/dt)
        cosmo_fact[z] = FOUR_PI_C_MPC * (1.0 + redz_edges[z]) * pow(dist_com_Mpc[z], 2);

        for(f = 0; f < num_freqs; f++) {
            if (z == 0 && f > 1) { dlnf[f] = log(fobs_edges[f]) - log(fobs_edges[f-1]); }
            physics::frst_from_fobs(fobs_cents[f], redz_cents[z], &frst_orb[z][f]);
        }
    }

    // ---- GWB Calculation

    for(f = 0; f < num_freqs; f++) {

        for(m1 = 0; m1 < num_mass; m1++) {
            for(m2 = 0; m2 < num_mass; m2++) {
                for(z = 0; z < num_redz; z++) {

                    // calculation differential number-densitiy of binary merger events
                    if (f == 0) {
                        num_dens[m1][m2][z] = number_density(mbh_mass_grams[m1], mbh_mass_grams[m2], redz_edges[z]);
                    }

                    //! FIX: consider implementing `num_dens` cutoff: if too low, skip
                    //! NOTE: make sure skipped values are zero

                    // ---- Calculate (expectation-value) number of binaries in the universe at this bin

                    // `tauf` [sec], GW Hardening timescale in frequency
                    physics::gw_hardening_time_freq(mbh_mass_grams[m1], mbh_mass_grams[m2], frst_orb[z][f], &tauf);

                    // expectation value for differential number of binaries
                    // this is:  $d^3 n / [dlog10(M1) dlog10(M2) dz dlog_e(f)]$
                    //! FIX: make sure this is d4n/[dlog10(m1) dlog10(m2) dz dlnf]  --- i.e. the log10 of both masses.
                    diff_numb[m1][m2][z] = num_dens[m1][m2][z] * cosmo_fact[z] * tauf;

                    // ---- Integrate over differential number

                    if ( m1 == 0 || m2 == 0 || z == 0 ) continue;

                    numb_expect = 0.0;
                    // stencil over corners of the bin in 3D
                    for (int ii = 0; ii < 2; ii++) {
                        for (int jj = 0; jj < 2; jj++) {
                            for (int kk = 0; kk < 2; kk++) {
                                numb_expect += (
                                    diff_numb[m1-ii][m2-jj][z-kk] *
                                    dlnf[f] *
                                    dlog10_mbh[m1-1] *
                                    dlog10_mbh[m2-1]
                                );
                            }
                        }
                    }
                    // we've added together 2-sides on 3 dimensions, so divide by 2^3 = 8
                    numb_expect /= 8.0;

                    // construct Poisson distribution centered on the expectation value
                    std::poisson_distribution<int> dist(numb_expect);

                    // ---- Calculate GW amplitude for each binary

                    // `hs2` [], Spectral strain-squared of the binary
                    physics::gw_strain_source_sq(mchirp[m1][m2], dist_com_cm[z], frst_orb[z][f], &hs2);

                    for(r = 0; r < num_reals; r++) {
                        numb_in_real = dist(rng);
                        if (numb_in_real == 0) continue;
                        gwb[f][r] += hs2 * numb_in_real / dlnf[f];
                    } // r

                } // z
            } // m2
        } // m1

        for(r = 0; r < num_reals; r++) {
            gwb[f][r] = sqrt(gwb[f][r]);
        }

    } // f

    // ---- Clean up

    for(m1 = 0; m1 < num_mass; m1++) {
        for(m2 = 0; m2 < num_mass; m2++) {
            free(num_dens[m1][m2]);
            free(diff_numb[m1][m2]);
        }
        free(mchirp[m1]);
        free(num_dens[m1]);
        free(diff_numb[m1]);
    }
    free(mchirp);
    free(num_dens);
    free(diff_numb);

};


/**
 *
 */

#include <boost/math/distributions/poisson.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <memory>   // for std::unique_ptr

#include "sam.h"

#include "cosmology.h"
#include "physics.h"
#include "utils.h"

using RNGType   = boost::random::mt19937;
using DistType  = boost::random::poisson_distribution<>;
using VGType    = boost::random::variate_generator<RNGType&, DistType>;

constexpr double FLOOR_NUMB_EXPECT = 1.0E-10;
constexpr double FLOOR_NUM_DENS = 1.0E-20;


RNGType rng(42);


void hdf5_write_meta(hid_t h5_file, SAM &sam, PTA &pta, GravWaves &gw) {
    // Write metadata to the HDF5 file

    LOG_DEBUG(get_logger(), "Writing metadata to HDF5...");

    //! FIX: WRITE SETTINGS/CONFIG TO HDF5 FILE ALSO
    //! ......

    // Write grid information to HDF5 file
    LOG_DEBUG(get_logger(), "Writing grid metadata to HDF5 file...");
    const char* group_name = "grid";
    utils::hdf5_write_array1d(h5_file, group_name, "mass_edges", sam.mass_edges, sam.num_mass_edges);
    utils::hdf5_write_array1d(h5_file, group_name, "mass_cents", sam.mass_cents, sam.num_mass_edges-1);
    utils::hdf5_write_array1d(h5_file, group_name, "redz_edges", sam.redz_edges, sam.num_redz_edges);
    utils::hdf5_write_array1d(h5_file, group_name, "redz_cents", sam.redz_cents, sam.num_redz_edges-1);
    utils::hdf5_write_array1d(h5_file, group_name, "fobs_edges", pta.fobs_edges, pta.num_freq_cents+1);
    utils::hdf5_write_array1d(h5_file, group_name, "fobs_cents", pta.fobs_cents, pta.num_freq_cents);
    utils::hdf5_write_scalar<int>(h5_file, group_name, "num_reals", gw.num_reals);
    utils::hdf5_write_scalar<int>(h5_file, group_name, "num_louds", gw.num_louds);
}


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


/**
 * Evaluate a Schechter Galaxy Stellar-Mass Function (GSMF), to find galaxy number-density.
 *
 * @param mstar        Galaxy stellar-mass [Msol].
 * @param redz         Redshift [].
 * @param phi_terms
 * @param mchar_terms
 * @param alpha
 *
 * @return             Galaxy number density for the given stellar mass [Mpc^-3].
 *
 */
double SAM::gsmf_single_schechter(double mstar, double redz, double phi_terms[3], double mchar_terms[3], double alpha) {
    double rz2 = pow(redz, 2);
    double mchar = pow(10.0, mchar_terms[0] + mchar_terms[1] * redz + mchar_terms[2] * rz2);
    double mm = mstar / mchar;
    return (
        LOG10 *
        pow(10.0, phi_terms[0] + phi_terms[1] * redz + phi_terms[2] * rz2) *
        pow(mm, 1.0 + alpha) *
        exp(-mm)
    );
}


/**
 * Evaluate a Double-Schechter Galaxy Stellar-Mass Function (GSMF) to find galaxy number-density.
 *
 * @param mstar  Galaxy stellar-mass [Msol].
 * @param redz   Redshift [].
 *
 * @return       Galaxy number density [Mpc^-3].
 *
 */
double SAM::gsmf_double_schechter(double mstar, double redz) {
    double phi = gsmf_single_schechter(mstar, redz, this->gsmf_log10_phi1, this->gsmf_log10_mchar, this->gsmf_alpha1);
    phi +=       gsmf_single_schechter(mstar, redz, this->gsmf_log10_phi2, this->gsmf_log10_mchar, this->gsmf_alpha2);
    return phi;
}


/**
 * Convert from BH mass to stellar mass using Mbh–-Mbulge relation.
 *
 * @param mbh   Black-hole (MBH) mass [Msol].
 * @return      Stellar-mass of host galaxy [Msol].
 *
 */
double SAM::mmbulge_mstar_from_mbh(double mbh) {
    // get bulge mass
    double mstar = (log10(mbh) - this->mmbulge_mass_amp_log10) / this->mmbulge_plaw;
    mstar = this->mmbulge_mass_ref * pow(10.0, mstar);
    // convert from bulge-mass to stellar-mass
    mstar = mstar / this->mmbulge_bulge_frac;
    return mstar;
}


/**
 * Calculate the number-density of MBH binary mergers with the given masses and redshift.
 *
 * @param mbh1           Mass of the primary   MBH [Msol].
 * @param mbh2           Mass of the secondary MBH [Msol].
 * @param redz           Redshift of the merger [].
 * @param return_galgal  If `true`, return the number-density of galaxy mergers (not MBH mergers).
 *
 * @return               Number-density of binary mergers [Mpc^-3].
 *
 */
double SAM::number_density(double mbh1, double mbh2, double rz, bool return_galgal) {

    // ---- Get galaxy stellar masses

    // NOTE: assume primary MBH ==> primary galaxy mass ==> GSMF
    // primary-galaxy stellar-mass [Msol]
    double ms1 = mmbulge_mstar_from_mbh(mbh1);

    // ---- Get galaxy stellar mass function

    // number-density of galaxies, [Mpc^-3]
    double phi = gsmf_double_schechter(ms1, rz);

    // ---- Get galaxy-galaxy merger rate

    // secondary-galaxy stellar-mass [Msol]
    double ms2 = mmbulge_mstar_from_mbh(mbh2);
    // galaxy merger rate (GMR)
    double gmr = gmr_illustris(ms1 + ms2, ms2/ms1, rz);

    // ---- Get galaxy-galaxy merger, number density

    double dtdz;
    cosmo->dtdz_from_redz(rz, &dtdz);

    double ndens = phi * gmr * dtdz;

    // printf(
    //     "mbh1=%.2e, ms1=%.2e | mbh2=%.2e, ms2=%.2e | phi=%.2e, gmr=%.2e, dtdz=%.2e\n",
    //     mbh1, ms1, mbh2, ms2, phi, gmr, dtdz
    // );

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
    LOG_DEBUG(get_logger(), "SAM::grav_waves(): allocating...");

    // Grid/domain sizes
    int num_mass_edges = this->num_mass_edges;
    int num_mass_cents = num_mass_edges - 1;
    int num_redz_edges = this->num_redz_edges;
    int num_redz_cents = num_redz_edges - 1;
    int num_freq_cents = pta.num_freq_cents;
    int num_reals = gw.num_reals;
    int num_louds = gw.num_louds;

    // Grid edges & centers
    double* fobs_edges = pta.fobs_edges;
    double* fobs_cents = pta.fobs_cents;
    double* mbh_mass_edges = this->mass_edges;
    double* mbh_mass_cents = this->mass_cents;
    double* redz_edges = this->redz_edges;
    double* redz_cents = this->redz_cents;

    // Values to be calculated/filled
    double** gwb = gw.gwb;   // [F][R]     : frequencies, realizations
    // double*** cws = gw.cws;  // [F][R][L]  : frequencies, realizations, loudest

    // Temporary variables
    double hs2;                // spectral-strain-squared of an individual binary
    // double tauf;               // GW Hardening timescale in frequency (f / df/dt)  [s]
    int numb_in_real;
    int m1, m2, z, f, r;      // loop iteration variables
    int m1c, m2c, zc;         // loop iteration variables for bin-centers
    int ii, jj, kk;           // loop iteration variables for integration stencil
    int idx;
    std::string msg;

    double mbh_mass_edges_grams[num_mass_edges];    // MBH mass, at grid edges   [grams]
    double mbh_mass_cents_grams[num_mass_cents];  // MBH mass, at grid centers [grams]
    double frst_orb_zedges[num_redz_edges][num_freq_cents];   // Rest-frame orbital freq (at redshift-edges)   [1/s]
    double frst_orb_zcents[num_redz_cents][num_freq_cents]; // Rest-frame orbital freq (at redshift-centers) [1/s]
    double dcom_edges_Mpc[num_redz_edges]; // comoving distance, at redshift bin edges [MPC]
    double dcom_edges_cm[num_redz_edges];    // comoving distance, at redshift bin edges [cm]
    double dcom_cents_cm[num_redz_cents];  // comoving distance, at redshift bin *centers* [cm]
    double cosmo_fact[num_redz_edges];     // dVc/dz  [Mpc^3/s]
    // bin-widths
    double dlog10_mbh[num_mass_cents];   // log10(mbh) bin-width
    double dz[num_redz_cents];           // redshift bin-width
    double dlnf[num_freq_cents];

    double** mchirp_cents_grams;    // [M1-1][M2-1] : chirp mass (rest-frame) at grid centers [grams]
    double*** num_dens;       // [M1][M2][Z] : differential number-density of binaries at grid edges [Mpc^-3]
    double*** diff_numb;      // [M1][M2][Z] : differential-number of binaries at grid edges []

    //! FIX: this doesn't need to be stored always, add preprocessor directive to turn-on and turn-off as desired
    double*** numb_expect;    // [M1-1][M2-1][Z-1] : number of binaries (expectation-value) in each bin (bin centers) []
    double*** tauf;           // [M1][M2][Z] :
    // double numb_expect;        // number of binaries (expectation-value) in each bin (bin centers) []

    boost::random::poisson_distribution<> dist_poisson;
    // boost::random::variate_generator<boost::random::mt19937&, boost::random::poisson_distribution<>> draw_poisson;
    std::unique_ptr<VGType> draw_poisson;

    // ---- Perform initializations & pre-calculations

    mchirp_cents_grams = (double**)malloc((num_mass_cents) * sizeof(double*));
    num_dens = (double***)malloc(num_mass_edges * sizeof(double**));
    diff_numb = (double***)malloc(num_mass_edges * sizeof(double**));
    tauf = (double***)malloc((num_mass_edges) * sizeof(double**));
    numb_expect = (double***)malloc((num_mass_cents) * sizeof(double**));     // bin-centers!

    #ifdef DEBUG_FREQ_STATS
    // int size4d_cedges = num_freq_cents * num_mass_edges * num_mass_edges * num_redz_edges;
    // int size4d_cents = num_freq_cents * num_mass_cents * num_mass_cents * num_redz_cents;
    int size3d_edges = num_mass_edges * num_mass_edges * num_redz_edges;
    int size3d_cents = num_mass_cents * num_mass_cents * num_redz_cents;
    double* num_dens_flat = (double*)malloc(size3d_edges * sizeof(double));
    double* tauf_flat = (double*)malloc(size3d_edges * sizeof(double));
    double* diff_numb_flat = (double*)malloc(size3d_edges * sizeof(double));
    double* numb_expect_flat = (double*)malloc(size3d_cents * sizeof(double));
    #endif  // DEBUG_FREQ_STATS

    LOG_DEBUG(get_logger(), "SAM::grav_waves(): allocation completed.  Initializing...\n");

    for(m1 = 0; m1 < num_mass_edges; m1++) {
        mbh_mass_edges_grams[m1] = mbh_mass_edges[m1] * MSOL;
        if (m1 > 0) {
            m1c = m1 - 1;
            dlog10_mbh[m1c] = log10(mbh_mass_edges[m1]) - log10(mbh_mass_edges[m1-1]);
            mbh_mass_cents_grams[m1c] = mbh_mass_cents[m1] * MSOL;
        }
    }

    for(m1 = 0; m1 < num_mass_edges; m1++) {

        // initialize arrays shaped corresponding to bin *edges*
        num_dens[m1] = (double**)malloc(num_mass_edges * sizeof(double*));
        diff_numb[m1] = (double**)malloc(num_mass_edges * sizeof(double*));
        tauf[m1] = (double**)malloc(num_mass_edges * sizeof(double*));

        if (m1 > 0) {
            m1c = m1 - 1;
            mchirp_cents_grams[m1c] = (double*)malloc((num_mass_cents) * sizeof(double));
            numb_expect[m1c] = (double**)malloc((num_mass_cents) * sizeof(double*));
        }

        for(m2 = 0; m2 < num_mass_edges; m2++) {
            num_dens[m1][m2] = (double*)malloc(num_redz_edges * sizeof(double));
            diff_numb[m1][m2] = (double*)malloc(num_redz_edges * sizeof(double));
            tauf[m1][m2] = (double*)malloc(num_redz_edges * sizeof(double));

            // `mchirp_cents_grams` [gram], Chirp-Mass (rest-frame)
            if (m1 > 0 && m2 > 0) {
                m2c = m2 - 1;
                numb_expect[m1c][m2c] = (double*)malloc((num_redz_cents) * sizeof(double));  // bin-centers!
                physics::chirp_mass_from_m1m2(
                    mbh_mass_cents_grams[m1c], mbh_mass_cents_grams[m2c], &mchirp_cents_grams[m1c][m2c]
                );
            }
        } // m2
    } // m1

    for(z = 0; z < num_redz_edges; z++) {

        this->cosmo->dcom_from_redz(redz_edges[z], &dcom_edges_Mpc[z]);

        // convert from Mpc to cm
        dcom_edges_cm[z] = dcom_edges_Mpc[z] * MPC;

        // [Mpc^3/s] : this is (dVc/dz) * (dz/dt)
        cosmo_fact[z] = FOUR_PI_C_MPC * (1.0 + redz_edges[z]) * pow(dcom_edges_Mpc[z], 2);

        if (z > 0) {
            zc = z - 1;
            dz[zc] = redz_edges[z] - redz_edges[z-1];
            // Get comoving distance at bin-centers in Mpc
            this->cosmo->dcom_from_redz(redz_cents[zc], &dcom_cents_cm[zc]);
            // convert from Mpc to cm
            dcom_cents_cm[zc] = dcom_cents_cm[zc] * MPC;
        }

        for(f = 0; f < num_freq_cents; f++) {
            physics::frst_from_fobs(fobs_cents[f], redz_edges[z], &frst_orb_zedges[z][f]);
            if (z == 0) {
                dlnf[f] = log(fobs_edges[f+1]) - log(fobs_edges[f]);
            } else {
                physics::frst_from_fobs(fobs_cents[f], redz_cents[zc], &frst_orb_zcents[zc][f]);
            }
        }
    }

    // ---- Setup HDF5 file and outputs

    #ifdef HDF5_OUTPUT
    LOG_DEBUG(get_logger(), "SAM::grav_waves(): intializing hdf5 file '{}' for output...", FNAME_OUTPUT_HDF5);
    hid_t h5_file = H5Fcreate(FNAME_OUTPUT_HDF5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    #ifdef HDF5_OUTPUT_DETAILS
    H5Slice_4D h5_numb_expect = H5Slice_4D(
        h5_file, "numb_expect",
        num_freq_cents, num_mass_cents, num_mass_cents, num_redz_cents
    );
    H5Slice_4D h5_tauf = H5Slice_4D(
        h5_file, "tauf",
        num_freq_cents, num_mass_edges, num_mass_edges, num_redz_edges
    );
    #endif  // HDF5_OUTPUT_DETAILS
    #endif  // HDF5_OUTPUT

    // ---- GWB Calculation

    LOG_DEBUG(get_logger(), "SAM::grav_waves(): initialization completed.  Calculating GWB...");

    // int num_floor_numb_expect = 0;
    // int num_floor_num_dens = 0;

    for(f = 0; f < num_freq_cents; f++) {

        for(m1 = 0; m1 < num_mass_edges; m1++) {
            for(m2 = 0; m2 < num_mass_edges; m2++) {
                for(z = 0; z < num_redz_edges; z++) {

                    // calculation differential number-densitiy of binary merger events
                    if (f == 0) {
                        num_dens[m1][m2][z] = number_density(
                            mbh_mass_edges[m1], mbh_mass_edges[m2], redz_edges[z]
                        );
                    }

                    //! FIX: consider implementing `num_dens` cutoff: if too low, skip
                    //! NOTE: make sure skipped values are set to zero as needed (i.e. if initialized w/ malloc)

                    // ---- Calculate (expectation-value) number of binaries in the universe at this bin

                    // `tauf` [sec], GW Hardening timescale in frequency
                    physics::gw_hardening_time_freq(
                        mbh_mass_edges_grams[m1], mbh_mass_edges_grams[m2], frst_orb_zedges[z][f], &tauf[m1][m2][z]
                    );

                    //! FIX: ENABLE FLOOR, MAKE SURE TO ZERO UNSET VALUES CORREECTLY
                    // if ((FLOOR_NUM_DENS > 0.0) && (num_dens[m1][m2][z] < FLOOR_NUM_DENS)) {
                    //     num_dens[m1][m2][z] = 0.0;
                    //     diff_numb[m1][m2][z] = 0.0;
                    //     num_floor_num_dens++;
                    //     continue;
                    // }

                    // expectation value for differential number of binaries
                    // this is:  $d^3 n / [dlog10(M1) dlog10(M2) dz dlog_e(f)]$
                    //! FIX: make sure this is d4n/[dlog10(m1) dlog10(m2) dz dlnf]  --- i.e. the log10 of both masses.
                    diff_numb[m1][m2][z] = num_dens[m1][m2][z] * cosmo_fact[z] * tauf[m1][m2][z];

                    // ---- Integrate over differential number

                    // We perform a 3D trapezoid-rule integration over M1, M2, and z
                    // to get the expectation value for the number of binaries in the bin.
                    // The 'centers' have one less element than the 'edges', so we skip the 0th
                    // element in the 3D grid.  Centers are indexed as 1-less than the edges:
                    // i.e. 'm1-1', 'm2-1', 'z-1'
                    // Remember that the frequencies we're interating over are already bin centers!

                    if ( m1 == 0 || m2 == 0 || z == 0 ) continue;

                    m1c = m1 - 1;
                    m2c = m2 - 1;
                    zc = z - 1;

                    numb_expect[m1c][m2c][zc] = 0.0;

                    // stencil over corners of the bin in 3D
                    // We subtract 0,1 to ensure we get the 0th elements that we skipped above
                    for (ii = 0; ii < 2; ii++) {
                        for (jj = 0; jj < 2; jj++) {
                            for (kk = 0; kk < 2; kk++) {
                                numb_expect[m1c][m2c][zc] += (
                                    diff_numb[m1-ii][m2-jj][z-kk] *
                                    dlnf[f] *
                                    dlog10_mbh[m1c] *
                                    dlog10_mbh[m2c]
                                );
                            }
                        }
                    }

                    if (numb_expect[m1c][m2c][zc] == 0.0) continue;

                    // we've added together 2-sides on 3 dimensions, so divide by 2^3 = 8
                    numb_expect[m1c][m2c][zc] /= 8.0;
                    //! FIX: ENABLE FLOOR, MAKE SURE TO ZERO UNSET VALUES CORREECTLY
                    // if ((FLOOR_NUMB_EXPECT > 0.0) && (numb_expect[m1c][m2c][zc] < FLOOR_NUMB_EXPECT)) {
                    //     numb_expect[m1c][m2c][zc] = 0.0;
                    //     num_floor_numb_expect++;
                    //     continue;
                    // }

                    // construct Poisson distribution centered on the expectation value
                    // boost::random::poisson_distribution<> dist_poisson(numb_expect[m1c][m2c][zc]);
                    dist_poisson = boost::random::poisson_distribution(numb_expect[m1c][m2c][zc]);
                    // draw_poisson = boost::random::variate_generator<boost::random::mt19937&, boost::random::poisson_distribution<>>(rng, dist_poisson);
                    draw_poisson = std::make_unique<VGType>(rng, dist_poisson);

                    // ---- Calculate GW amplitude for each binary

                    // `hs2` [], Spectral strain-squared of the binary
                    physics::gw_strain_source_sq(
                        mchirp_cents_grams[m1c][m2c], dcom_cents_cm[zc], frst_orb_zcents[zc][f], &hs2
                    );

                    for(r = 0; r < num_reals; r++) {
                        //! FIX: consider implementing a cutoff to approximate the Poisson distribution
                        //!      with a normal distribution.
                        numb_in_real = (*draw_poisson)();
                        if (numb_in_real == 0) continue;
                        gwb[f][r] += hs2 * numb_in_real / dlnf[f];
                    } // r

                    // Prepare flatted output for printing statistics
                    #ifdef DEBUG_FREQ_STATS
                    utils::index_3d_to_1d(m1c, m2c, zc, num_mass_cents, num_mass_cents, num_redz_cents, &idx);
                    numb_expect_flat[idx] = numb_expect[m1c][m2c][zc];
                    utils::index_3d_to_1d(m1, m2, z, num_mass_edges, num_mass_edges, num_redz_edges, &idx);
                    diff_numb_flat[idx] = diff_numb[m1][m2][z];
                    tauf_flat[idx] = tauf[m1][m2][z];
                    if (f == 0) {
                        num_dens_flat[idx] = num_dens[m1][m2][z];
                    }
                    #endif  // DEBUG_FREQ_STATS

                } // z
            } // m2
        } // m1

        for(r = 0; r < num_reals; r++) {
            gwb[f][r] = sqrt(gwb[f][r]);
        }

        #ifdef DEBUG_FREQ_STATS
        LOG_INFO(get_logger(), "frequency {:03d}\n", f);
        if (f == 0) {
            msg = utils::quantiles_string(num_dens_flat, size3d_edges, {0.10, 0.25, 0.5, 0.75, 0.90});
            LOG_INFO(get_logger(), "`num_dens`       : {}\n", msg);
        }
        msg = utils::quantiles_string(tauf_flat, size3d_edges, {0.10, 0.25, 0.5, 0.75, 0.90});
        LOG_INFO(get_logger(), "`tauf[f]`        : {}\n", msg);
        msg = utils::quantiles_string(diff_numb_flat, size3d_edges, {0.10, 0.25, 0.5, 0.75, 0.90});
        LOG_INFO(get_logger(), "`diff_numb[f]`   : {}\n", msg);
        msg = utils::quantiles_string(numb_expect_flat, size3d_cents, {0.10, 0.25, 0.5, 0.75, 0.90});
        LOG_INFO(get_logger(), "`numb_expect[f]` : {}\n", msg);
        msg = utils::quantiles_string(gwb[f], num_reals, {0.10, 0.25, 0.5, 0.75, 0.90});
        LOG_INFO(get_logger(), "`gwb[f]`         : {}\n", msg);
        #endif  // DEBUG_FREQ_STATS

        // ---- Save output to HDF5 file

        // Write hyperslab at [f, :, :, :]
        #ifdef HDF5_OUTPUT_DETAILS
        h5_numb_expect.write_slice_at(f, numb_expect, num_mass_cents, num_mass_cents, num_redz_cents);
        h5_tauf.write_slice_at(f, tauf, num_mass_edges, num_mass_edges, num_redz_edges);
        #endif  // HDF5_OUTPUT_DETAILS

        LOG_DEBUG(get_logger(), "freq {:03d}/{:03d} done.", f, num_freq_cents);

    } // f

    LOG_DEBUG(get_logger(), "SAM::grav_waves(): Calculating GWB completed.");

    // ---- Check/report diagnostics

    // if (FLOOR_NUM_DENS > 0.0) {
    //     int tot_num = num_freq_cents * num_mass_edges * num_mass_edges * num_redz_edges;
    //     LOG_INFO(
    //         get_logger(), "`num_dens`    below floor ({:.2e}): {}/{}={:.2e}",
    //         FLOOR_NUM_DENS, num_floor_num_dens, tot_num, (double)num_floor_num_dens/tot_num
    //     );
    // }

    // if (FLOOR_NUMB_EXPECT > 0.0) {
    //     int tot_num = num_freq_cents * num_mass_cents * num_mass_cents * num_redz_cents;
    //     LOG_INFO(
    //         get_logger(), "`numb_expect` below floor ({:.2e}): {}/{}={:.2e}",
    //         FLOOR_NUMB_EXPECT, num_floor_numb_expect, tot_num, (double)num_floor_numb_expect/tot_num
    //     );
    // }

    #ifdef HDF5_OUTPUT
    LOG_DEBUG(get_logger(), "Writing data to HDF5...");

    // ---- Write additional data to HDF5 file

    utils::hdf5_write_array2d<double>(
        h5_file, NULL, "gwb", gwb,
        num_freq_cents, num_reals
    );

    #ifdef HDF5_OUTPUT_DETAILS
    utils::hdf5_write_array3d<double>(
        h5_file, NULL, "num_dens",
        num_dens, num_mass_edges, num_mass_edges, num_redz_edges
    );
    #endif  // HDF5_OUTPUT_DETAILS

    hdf5_write_meta(h5_file, *this, pta, gw);

    LOG_DEBUG(get_logger(), "SAM::grav_waves(): Wrote data to HDF5.");
    #endif  // HDF5_OUTPUT

    // ---- Clean up

    LOG_DEBUG(get_logger(), "SAM::grav_waves(): Cleaning up...");

    for(m1 = 0; m1 < num_mass_edges; m1++) {
        for(m2 = 0; m2 < num_mass_edges; m2++) {
            free(num_dens[m1][m2]);
            free(diff_numb[m1][m2]);
            free(tauf[m1][m2]);
            if (m1 > 0 && m2 > 0) {
                m1c = m1 - 1;
                m2c = m2 - 1;
                free(numb_expect[m1c][m2c]);
            }
        }
        free(num_dens[m1]);
        free(diff_numb[m1]);
        free(tauf[m1]);
        if (m1 > 0) {
            m1c = m1 - 1;
            free(mchirp_cents_grams[m1c]);
            free(numb_expect[m1c]);
        }
    }
    free(num_dens);
    free(diff_numb);
    free(mchirp_cents_grams);
    free(tauf);

    #ifdef DEBUG_FREQ_STATS
    free(num_dens_flat);
    free(tauf_flat);
    free(diff_numb_flat);
    free(numb_expect_flat);
    #endif  // DEBUG_FREQ_STATS

    // Close HDF5 dataset and file
    #ifdef HDF5_OUTPUT
    H5Fclose(h5_file);
    #endif  // HDF5_OUTPUT

    LOG_DEBUG(get_logger(), "SAM::grav_waves(): completed.");

};



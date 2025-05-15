/**
 *
 */


#ifndef SAM_H
#define SAM_H

#include <cstdlib> // For malloc and free
#include <cmath>

#include "constants.h"
#include "physics.h"


class PTA {

public:
    double obs_dur;
    int num_freqs;
    double* fobs_cents;
    double* fobs_edges;

    PTA(double dur=20.0, double nfreqs=30) : obs_dur(dur), num_freqs(nfreqs) {
        fobs_cents = (double*)malloc(num_freqs * sizeof(double));
        fobs_edges = (double*)malloc((num_freqs+1) * sizeof(double));
        double df = 1.0 / obs_dur;
        fobs_edges[0] = df * 0.5;
        for (int ii = 0; ii < num_freqs; ii++) {
            fobs_cents[ii] = df * (ii + 1);
            fobs_edges[ii+1] = df * (ii + 1.5);
        }
    };

    // double* fobs_cents() { return fobs_cents; };
    // double* fobs_edges() { return fobs_edges;};

    ~PTA() {
        free(fobs_cents);
        free(fobs_edges);
    };

// private:
//     double* fobs_cents;
//     double* fobs_edges;

};


class GravWaves {
public:
    int num_freqs;
    int num_reals;
    int num_louds;

    double** gwb;
    double*** cws;

    GravWaves(PTA *pta, int nreals = 100, int nloud = 5) : pta(pta), num_reals(nreals), num_louds(nloud) {
        num_freqs = pta->num_freqs;
        gwb = (double** )malloc(num_freqs * sizeof(double*));
        cws = (double*** )malloc(num_freqs * sizeof(double**));
        for (int ii = 0; ii < num_freqs; ii++) {
            gwb[ii] = (double* )malloc(num_reals * sizeof(double));
            cws[ii] = (double** )malloc(num_reals * sizeof(double*));
            for (int jj = 0; jj < num_reals; jj++) {
                cws[ii][jj] = (double* )malloc(num_louds * sizeof(double));
            }
        }
    }

    ~GravWaves() {
        for (int ii = 0; ii < pta->num_freqs; ii++) {
            free(gwb[ii]);
            for (int jj = 0; jj < num_reals; jj++) {
                free(cws[ii][jj]);
            }
            free(cws[ii]);
        }
        free(gwb);
        free(cws);
    };

private:
    PTA* pta;

};


class SAM {
public:
    int num_mass;
    int num_redz;
    double *mass_edges;
    double *redz_edges;
    double *mass_cents;
    double *redz_cents;

//! FIX: make these private
// private:
    // Grid
    double mass_pars[3] = {1E6, 1E12, 101};
    double redz_pars[3] = {1E-2, 1E1, 91};

    // GSMF - Double Schechter
    double gsmf_log10_phi1_z0  = -2.383;
    double gsmf_log10_phi1_z1  = -0.264;
    double gsmf_log10_phi1_z2  = -0.107;
    double gsmf_log10_phi2_z0  = -2.818;
    double gsmf_log10_phi2_z1  = -0.368;
    double gsmf_log10_phi2_z2  = +0.046;
    double gsmf_log10_mchar_z0 = +10.767;
    double gsmf_log10_mchar_z1 = +0.124;
    double gsmf_log10_mchar_z2 = -0.033;
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

    SAM() {
        num_mass = mass_pars[2];
        num_redz = redz_pars[2];
        mass_cents = (double*)malloc(num_mass * sizeof(double));
        redz_cents = (double*)malloc(num_redz * sizeof(double));
        mass_edges = (double*)malloc((num_mass+1) * sizeof(double));
        redz_edges = (double*)malloc((num_redz+1) * sizeof(double));

        // ---- Initialize mass edges and centers

        double dlog10m = (log10(mass_pars[1]) - log10(mass_pars[1])) / num_mass;
        double next, prev = log10(mass_pars[0]);
        mass_edges[0] = pow(10.0, prev);
        for (int ii = 0; ii < num_mass; ii++) {
            next = prev + dlog10m;
            mass_edges[ii+1] = pow(10.0, next);
            mass_cents[ii] = pow(10.0, 0.5 * (prev + next));
        }

        // ---- Initialize redshift edges and centers
        double dlog10z = (log10(redz_pars[1]) - log10(redz_pars[0])) / num_redz;
        prev = log10(redz_pars[0]);
        redz_edges[0] = pow(10.0, prev);
        for (int ii = 0; ii < num_redz; ii++) {
            next = prev + dlog10z;
            redz_edges[ii+1] = pow(10.0, next);
            redz_cents[ii] = pow(10.0, 0.5 * (prev + next));
        }

    }

    ~SAM() {
        free(mass_cents);
        free(redz_cents);
        free(mass_edges);
        free(redz_edges);
    };

    /*

    @classmethod
    def _gsmf_single_schechter(self, mstar, redz, phi_terms, mchar_terms, alpha):
        rz2 = redz**2

        phi = np.power(10.0, phi_terms[0] + phi_terms[1] * redz + phi_terms[2] * rz2)
        mchar = MSOL * np.power(10.0, mchar_terms[0] + mchar_terms[1] * redz + mchar_terms[2] * rz2)
        alpha = alpha

        xx = mstar / mchar
        rv = np.log(10.0) * phi * np.power(xx, 1.0 + alpha) * np.exp(-xx)
        return rv

    def _mstar_from_mbh(self, mbh, mass_amp_log10, mass_ref, plaw, bulge_mfrac):

        # get bulge mass
        mstar = (np.log10(mbh) - mass_amp_log10) / plaw
        mstar = mass_ref * np.power(10.0, mstar)
        # convert from bulge-mass to stellar-mass
        mstar = mstar / bulge_mfrac
        return mstar

    def number_density_3d(self, return_galgal=False):
        edges = self.edges_3d
        mbh1, mbh2, rz = self.grid3d

        # ---- Get galaxy stellar masses

        # NOTE: assume primary MBH ==> primary galaxy mass ==> GSMF
        ms1 = self._mstar_from_mbh(
            mbh1, self.mmbulge_mass_amp_log10, self.mmbulge_mass_ref, self.mmbulge_plaw, self.mmbulge_bulge_frac
        )

        # ---- Get galaxy stellar mass function

        phi = self._gsmf_single_schechter(ms1, rz, self.gsmf_log10_phi1, self.gsmf_log10_mchar, self.gsmf_alpha1)
        phi += self._gsmf_single_schechter(ms1, rz, self.gsmf_log10_phi2, self.gsmf_log10_mchar, self.gsmf_alpha2)

        # ---- Get galaxy-galaxy merger rate

        ms2 = self._mstar_from_mbh(
            mbh1, self.mmbulge_mass_amp_log10, self.mmbulge_mass_ref, self.mmbulge_plaw, self.mmbulge_bulge_frac
        )
        mstar_tot = ms1 + ms2
        mstar_rat = ms2 / ms1
        zplus1 = 1.0 + rz

        norm = ((10.0 ** self.gmr_norm0_log10) / GYR) * np.power(zplus1, self.gmr_normz)
        malpha = self.gmr_malpha0 * np.power(zplus1, self.gmr_malphaz)
        mdelta = self.gmr_mdelta0 * np.power(zplus1, self.gmr_mdeltaz)
        qgamma = self.gmr_qgamma0 * np.power(zplus1, self.gmr_qgammaz)
        qgamma = qgamma + self.gmr_qgammam * np.log10(mstar_tot/self.gmr_mref_gamma)

        gmr = (
            norm *
            np.power(mstar_tot/self.gmr_mref_gamma, malpha) *
            np.power(1.0 + mstar_tot/self.gmr_mref_delta, mdelta) *
            np.power(mstar_rat, qgamma)
        )

        # ---- Get galaxy-galaxy merger, number density

        ndens = phi * gmr * cosmo.dtdz(self.redz)[np.newaxis, np.newaxis, :]

        # Return only galaxy-galaxy merger rate
        if return_galgal:
            edges = [ms1, ms2, rz]
            return edges, ndens

        # ---- Convert from galaxy-galaxy to MBH-MBH merger rate

        # we want ``dn_mbhb / [dlog10(M_bh) dq_bh qz]``
        # so far we have ``dn_gal / [dlog10(M_gal) dq_gal dz]``

        # dn / [dM dq dz] = (dn_gal / [dM_gal dq_gal dz]) * (dM_gal/dM_bh) * (dq_gal / dq_bh)
        dqbh_dqgal = self.mmbulge_plaw * np.power(mstar_rat, self.mmbulge_plaw - 1.0)
        # (dMstar-pri / dMbh-pri) * (dMbh-pri/dMbh-tot) = (dMstar-pri / dMstar-tot) * (dMstar-tot/dMbh-tot)
        # ==> (dMstar-tot/dMbh-tot) = (dMstar-pri / dMbh-pri) * (dMbh-pri/dMbh-tot) / (dMstar-pri / dMstar-tot)
        #                           = (dMstar-pri / dMbh-pri) * (1 / (1+q_bh)) / (1 / (1+q_star))
        #                           = (dMstar-pri / dMbh-pri) * ((1+q_star) / (1+q_bh))
        dmstar_dmbh_pri = (ms1 / (self.mmbulge_plaw * mbh1))
        #! ==== FIX: this allows for q>1 ! ====
        qterm = (1.0 + mstar_rat) / (1.0 + (mbh2/mbh1))
        ndens *= ((mbh1 + mbh2) / mstar_tot) * (dmstar_dmbh_pri * qterm / dqbh_dqgal)

        # ---- Add scatter from the M-Mbulge relation

        #! =====================================================================
        #! FIX: still need to add scatter
        print("\n!!NOT ADDING SCATTER!!\n")
        # ndens = add_scatter_to_masses(self.mtot, self.mrat, ndens, scatter, log=log)
        #! =====================================================================

        return edges, ndens
    */
};


void grav_waves(PTA* pta, GravWaves *gw, SAM* sam);


#endif  // SAM_H

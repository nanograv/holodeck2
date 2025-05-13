/**
 *
 */


#ifndef SAM_H
#define SAM_H

#include <cstdlib> // For malloc and free
#include <cmath>

#include "constants.h"


class PTA {

public:
    float obsDur;
    int numFreqs;
    float* fobsCents;
    float* fobsEdges;

    PTA(float dur=20.0, float nfreqs=30) : obsDur(dur), numFreqs(nfreqs) {
        fobsCents = (float*)malloc(numFreqs * sizeof(float));
        fobsEdges = (float*)malloc((numFreqs+1) * sizeof(float));
        float df = 1.0 / obsDur;
        fobsEdges[0] = df * 0.5;
        for (int ii = 0; ii < numFreqs; ii++) {
            fobsCents[ii] = df * (ii + 1);
            fobsEdges[ii+1] = df * (ii + 1.5);
        }
    };

    // float* fobsCents() { return fobsCents; };
    // float* fobsEdges() { return fobsEdges;};

    ~PTA() {
        free(fobsCents);
        free(fobsEdges);
    };

// private:
//     float* fobsCents;
//     float* fobsEdges;

};


class GravWaves {
public:
    int numFreqs;
    int numReals;
    int numLouds;

    float** gwb;
    float*** cws;

    GravWaves(PTA *pta, int nreals = 100, int nloud = 5) : pta(pta), numReals(nreals), numLouds(nloud) {
        numFreqs = pta->numFreqs;
        gwb = (float** )malloc(numFreqs * sizeof(float*));
        cws = (float*** )malloc(numFreqs * sizeof(float**));
        for (int ii = 0; ii < numFreqs; ii++) {
            gwb[ii] = (float* )malloc(numReals * sizeof(float));
            cws[ii] = (float** )malloc(numReals * sizeof(float*));
            for (int jj = 0; jj < numReals; jj++) {
                cws[ii][jj] = (float* )malloc(numLouds * sizeof(float));
            }
        }
    }

    ~GravWaves() {
        for (int ii = 0; ii < pta->numFreqs; ii++) {
            free(gwb[ii]);
            for (int jj = 0; jj < numReals; jj++) {
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
    int numMass;
    int numRedz;
    float *massEdges;
    float *redzEdges;
    float *massCents;
    float *redzCents;

//! FIX: make these private
// private:
    // Grid
    float mass_pars[3] = {1E6, 1E12, 101};
    float redz_pars[3] = {1E-2, 1E1, 91};

    // GSMF - Double Schechter
    float gsmf_log10_phi1_z0  = -2.383;
    float gsmf_log10_phi1_z1  = -0.264;
    float gsmf_log10_phi1_z2  = -0.107;
    float gsmf_log10_phi2_z0  = -2.818;
    float gsmf_log10_phi2_z1  = -0.368;
    float gsmf_log10_phi2_z2  = +0.046;
    float gsmf_log10_mchar_z0 = +10.767;
    float gsmf_log10_mchar_z1 = +0.124;
    float gsmf_log10_mchar_z2 = -0.033;
    float gsmf_alpha1 = -0.28;
    float gsmf_alpha2 = -1.48;

    // GMR - Illustris; Galaxy Merger Rate
    float gmr_norm0_log10 = -2.2287;          // -2.2287 ± 0.0045    A0 [log10(A*Gyr)]
    float gmr_normz       = +2.4644;          // +2.4644 ± 0.0128    eta
    float gmr_malpha0     = +0.2241;          // +0.2241 ± 0.0038    alpha0
    float gmr_malphaz     = -1.1759;          // -1.1759 ± 0.0316    alpha1
    float gmr_mdelta0     = +0.7668;          // +0.7668 ± 0.0202    delta0
    float gmr_mdeltaz     = -0.4695;          // -0.4695 ± 0.0440    delta1
    float gmr_qgamma0     = -1.2595;          // -1.2595 ± 0.0026    beta0
    float gmr_qgammaz     = +0.0611;          // +0.0611 ± 0.0021    beta1
    float gmr_qgammam     = -0.0477;          // -0.0477 ± 0.0013    gamma
    float gmr_mref_delta  = 2.0E11;           // fixed value [Msol]
    float gmr_mref_gamma  = 1.0E10;           // fixed value [Msol]

    // M-MBulge
    float mmbulge_mass_amp_log10 = 8.69;      // 8.69 ± 0.05  [log10(M/Msol)]  approximate uncertainties!
    float mmbulge_mass_ref = 1E11;            // 1E11 Msol
    float mmbulge_plaw = 1.17;                // 1.17 ± 0.08
    float mmbulge_scatter_dex = 0.28;         // scatter stdev in dex
    float mmbulge_bulge_frac = 0.7;

    SAM() {
        numMass = mass_pars[2];
        numRedz = redz_pars[2];
        massCents = (float*)malloc(numMass * sizeof(float));
        redzCents = (float*)malloc(numRedz * sizeof(float));
        massEdges = (float*)malloc((numMass+1) * sizeof(float));
        redzEdges = (float*)malloc((numRedz+1) * sizeof(float));

        // ---- Initialize mass edges and centers

        float dlog10m = (log10(mass_pars[1]) - log10(mass_pars[1])) / numMass;
        float next, prev = log10(mass_pars[0]);
        massEdges[0] = pow(10.0, prev);
        for (int ii = 0; ii < numMass; ii++) {
            next = prev + dlog10m;
            massEdges[ii+1] = pow(10.0, next);
            massCents[ii] = pow(10.0, 0.5 * (prev + next));
        }

        // ---- Initialize redshift edges and centers
        float dlog10z = (log10(redz_pars[1]) - log10(redz_pars[0])) / numRedz;
        prev = log10(redz_pars[0]);
        redzEdges[0] = pow(10.0, prev);
        for (int ii = 0; ii < numRedz; ii++) {
            next = prev + dlog10z;
            redzEdges[ii+1] = pow(10.0, next);
            redzCents[ii] = pow(10.0, 0.5 * (prev + next));
        }

    }

    ~SAM() {
        free(massCents);
        free(redzCents);
        free(massEdges);
        free(redzEdges);
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



#endif  // SAM_H

/**
 * Holodeck2 physics module.
 *
 */

#pragma once

#include <cmath>

#include "constants.h"


namespace physics {

    void chirp_mass_from_m1m2(double m1, double m2, double *mc);

    void mtmr_from_m1m2(double m1, double m2, double* mtot, double* mrat);

    void m1m2_from_mtmr(double mt, double mr, double *m1, double* m2);

    /**
     * Convert from rest-frame frequency to observer-frame frequency.
     *
     * The output units will match the input units, so any units can be used, but the convention is
     * to use [1/sec] for frequency.
     *
     * NOTE:   If `nharm` is not 1, then `fobs` is interpretted as the *harmonic* frequency while `frst` is interpretted
     *         as the *orbital* frequency, so the output (`fobs`) will be multiplied by `nharm`.
     *
     * @param frst       Rest-frame frequency [1/sec].
     * @param redz       Redshift [].
     * @param fobs       [output] Pointer for the observer-frame frequency [1/sec].
     * @param nharm      Harmonic number *corresponding to observer-frame* frequency (default is 1).
     *
     */
    void fobs_from_frst(double frst, double redz, double* fobs, int nharm=1);

    /**
     * Convert from observer-frame frequency to rest-frame frequency.
     *
     * The output units will match the input units, so any units can be used, but the convention is
     * to use [1/sec] for frequency.
     *
     * NOTE:   If `nharm` is not 1, then `fobs` is interpretted as the *harmonic* frequency while `frst` is interpretted
     *         as the *orbital* frequency, so the output (`frst`) will be divided by `nharm`.
     *
     * @param fobs       Observer-frame frequency [1/sec].
     * @param redz       Redshift [].
     * @param frst       [output] Pointer for the rest-frame frequency [1/sec].
     * @param nharm      Harmonic number *corresponding to observer-frame* frequency (default is 1).
     *
     */
    void frst_from_fobs(double fobs, double redz, double* frst, int nharm=1);

    void gw_hardening_rate_dadt(double m1, double m2, double sepa, double* dadt);

    /**
     * Calculate the hardening time in terms of frequency, $\tau_f = f/(df/dt) = dt/dln(f)$.
     *
     * @param m1          Mass of the primary [gram].
     * @param m2          Mass of the secondary [gram].
     * @param frst_orb    Rest-frame orbital frequency [1/sec].
     * @return tauf       Output pointer for the hardening time [sec].
     *
     */
    void gw_hardening_time_freq(double m1, double m2, double frst_orb, double *tauf);

    void gw_strain_source(double mchirp, double dcom, double freq_rest_orb, double* hs);

    void gw_strain_source_sq(double mchirp, double dcom, double freq_rest_orb, double* hs2);

    void kepler_forb_from_sepa(double mass, double sepa, double* forb);

    void kepler_sepa_from_freq(double mass, double freq, double* sepa);

} // namespace physics


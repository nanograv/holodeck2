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

    void fobs_from_frst(double frst, double redz, double* fobs, int nharm=1);

    void frst_from_fobs(double fobs, double redz, double* frst, int nharm=1);

    void gw_hardening_rate_dadt(double m1, double m2, double sepa, double* dadt);

    void gw_hardening_time_freq(double m1, double m2, double frst_orb, double *tauf);

    void gw_strain_source(double mchirp, double dcom, double freq_rest_orb, double* hs);

    void gw_strain_source_sq(double mchirp, double dcom, double freq_rest_orb, double* hs2);

    void kepler_forb_from_sepa(double mass, double sepa, double* forb);

    void kepler_sepa_from_freq(double mass, double freq, double* sepa);

} // namespace physics


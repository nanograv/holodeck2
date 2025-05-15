/**
 *
 */

//  import numbers

// import numpy as np
// from .constants import MSOL, GYR, SPLC, MPC

// from holodeck2 import cosmo, physics, utils
#include "sam.h"


void grav_waves(PTA* pta, GravWaves *gw, SAM* sam) {
    int num_freqs = pta->num_freqs;
    int num_reals = gw->num_reals;
    int num_louds = gw->num_louds;
    int num_mass = sam->num_mass;
    int num_redz = sam->num_redz;

    double* fobs_edges = pta->fobs_edges;
    double* fobs_cents = pta->fobs_cents;
    double* mass_edges = sam->mass_edges;
    double* mass_cents = sam->mass_cents;
    double* redz_edges = sam->redz_edges;
    double* redz_cents = sam->redz_cents;

    double** gwb = gw->gwb;
    double*** cws = gw->cws;

    double hs2;    // spectral strain of an individual binary
    double mchirp; // chirp mass
    double dist_com; // comoving distance

    for(int m1 = 0; m1 < num_mass; m1++) {
        for(int m2 = 0; m2 < num_mass; m2++) {

            // mchirp = physics.chirp_mass(m1=mass_cents[m1], m2=mass_cents[m2]);
            // dist_com = cosmo.comoving_distance(redz_cents).cgs.value;

            for(int z = 0; z < num_redz; z++) {
                for(int f = 0; f < num_freqs; f++) {

                    // hs2 = pow(physics.gw_strain_source(mchirp*MSOL, dist_com, frst_orb), 2);
                    // gwb_sq = (numb_4d / d_lnf) * hs2


                    for(int r = 0; r < num_reals; r++) {

                    }

                }
            }

        }
    }



}

/*
def number_expect_4d_gwonly_instant(fobs_gw_edges, edges_3d, ndens_3d):
    """Total number of binaries in the universe, assuming GW-Only & instantaneous evolution.

    Assumptions:
    * GW-Only: this means that the residence time at each frequency is determined by the rate of GW emission.
    * Instantaneous: this means that there is no delay between galaxy merger and reaching each GW frequency.

    """

    # NOTE: this is a 'linear'-spaced center-point
    fobs_gw_cents = utils.midpoints_lin(fobs_gw_edges)
    dln_fobs = np.diff(np.log(fobs_gw_edges))

    mbh1, mbh2, redz = edges_3d
    cents_4d = [fobs_gw_cents, utils.midpoints_log(mbh1), utils.midpoints_log(mbh2), utils.midpoints_lin(redz)]

    # ---- Residence Time : time spent in the given frequency-band

    # (F, Z) rest-frame orbital frequencies
    frst_orb = physics.frst_from_fobs(fobs_gw_cents[:, np.newaxis], redz[np.newaxis, :], nharm=2)
    # (F, M, M, Z)
    tauf = physics.gw_hardening_time_freq(
        mbh1[np.newaxis, :, np.newaxis, np.newaxis] * MSOL,
        mbh2[np.newaxis, np.newaxis, :, np.newaxis] * MSOL,
        frst_orb[:, np.newaxis, np.newaxis, :],
    )

    # ---- Cosmological Factor : convert from number-density to number-in-Universe

    # (Z,) this is `(dVc/dz) * (dz/dt)` in units of [Mpc^3/s]
    cosmo_fact = 4 * np.pi * (SPLC/MPC) * (1.0 + redz) * (cosmo.comoving_distance(redz).to('Mpc').value ** 2)

    # (F, M, M, Z) expectation value for differential number of binaries
    # d^3 n / dlog10(M) dq dz dlogf
    #! FIX: make sure this is d4n/[dlog10(m1) dlog10(m2) dz dlnf]  --- i.e. the log10 of both masses.
    diff_numb = ndens_3d[np.newaxis, ...] * cosmo_fact * tauf

    # ==== Integrate over differential number ==============================

    # Frequency ($ln f$), axis 0   :   (F,M,M,Z) ==> (F,M,M,Z)
    numb = diff_numb[:, ...] * dln_fobs[:, np.newaxis, np.newaxis, np.newaxis]

    # Mass 1 ($log_{10} m_1$), axis 1   :   (F,M,M,Z) ==> (F,M-1,M,Z)
    dlog10_mbh = np.diff(np.log10(mbh1))
    numb = (numb[:, 1:, ...] + numb[:, :-1, ...]) * dlog10_mbh[np.newaxis, :, np.newaxis, np.newaxis]

    # Mass 2 ($log_{10} m_1$), axis 2   :   (F,M-1,M,Z) ==> (F,M-1,M-1,Z)
    numb = (numb[:, :, 1:, :] + numb[:, :, :-1, :]) * dlog10_mbh[np.newaxis, np.newaxis, :, np.newaxis]

    # Redshift ($z$), axis 3   :   (F,M-1,M-1,Z) ==> (F,M-1,M-1,Z-1)
    d_redz = np.diff(redz)
    numb = (numb[:, :, :, 1:] + numb[:, :, :, :-1]) * d_redz[np.newaxis, np.newaxis, np.newaxis, :]

    # we've added together 2-sides on 3 dimensions, so 2^3 = 8
    numb /= 8.0

    return cents_4d, numb


def gws_from_number_expect_instant(fobs_gw_edges, cents_4d, numb_4d, realize=None):

    # (F,)
    d_lnf = np.diff(np.log(fobs_gw_edges))
    # (F, Z)
    frst_orb = physics.frst_from_fobs(cents_4d[0][:, np.newaxis], cents_4d[3][np.newaxis, :], nharm=2)
    # (M1, M2)
    mchirp = physics.chirp_mass(m1=cents_4d[1][:, np.newaxis], m2=cents_4d[1][np.newaxis, :])

    # (Z,)
    dist_com = cosmo.comoving_distance(cents_4d[3]).cgs.value

    d_lnf = d_lnf[:, np.newaxis, np.newaxis, np.newaxis]
    mchirp = mchirp[np.newaxis, :, :, np.newaxis]
    dist_com = dist_com[np.newaxis, np.newaxis, np.newaxis, :]
    frst_orb = frst_orb[:, np.newaxis, np.newaxis, :]

    # `True` :: Single realization
    if (realize is True):
        numb_4d = np.random.poisson(numb_4d)

    # `None`|`False` :: No realizations (GWB expectation value)
    elif (realize in [None, False]):
        pass

    elif isinstance(realize, numbers.Integral):
        # New shape will add realizations as a new dimension, (F, M1, M2, Z) ==> (F, M1, M2, Z, R)
        reals_shape = numb_4d.shape + (realize,)

        d_lnf = d_lnf[..., np.newaxis]
        mchirp = mchirp[..., np.newaxis]
        dist_com = dist_com[..., np.newaxis]
        frst_orb = frst_orb[..., np.newaxis]
        numb_4d = np.random.poisson(numb_4d[..., np.newaxis], size=reals_shape)

    else:
        raise ValueError(f"Unrecognized type/value for `realize`!  {realize=}")

    gwb_sq = (numb_4d / d_lnf) * physics.gw_strain_source(mchirp*MSOL, dist_com, frst_orb) ** 2

    return gwb_sq
*/
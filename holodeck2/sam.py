"""
"""

import numpy as np
from .constants import MSOL, GYR, SPLC, MPC

from holodeck2 import cosmo, physics, utils


class SAM:

    # Grid
    mass_pars = (1e6, 1e12, 21)
    # mrat_pars = (1e-4, 1.0, 71)
    redz_pars = (1e-2, 1e1, 11)

    # GSMF - Double Schechter
    gsmf_log10_phi1 = [-2.383, -0.264, -0.107]
    gsmf_log10_phi2 = [-2.818, -0.368, +0.046]
    gsmf_log10_mchar = [+10.767, +0.124, -0.033]
    gsmf_alpha1 = -0.28
    gsmf_alpha2 = -1.48

    # GMR - Illustris Galaxy Merger Rate
    gmr_norm0_log10 = -2.2287      # -2.2287 ± 0.0045    A0 [log10(A*Gyr)]
    gmr_normz = +2.4644            # +2.4644 ± 0.0128    eta
    gmr_malpha0 = +0.2241          # +0.2241 ± 0.0038    alpha0
    gmr_malphaz = -1.1759          # -1.1759 ± 0.0316    alpha1
    gmr_mdelta0 = +0.7668          # +0.7668 ± 0.0202    delta0
    gmr_mdeltaz = -0.4695          # -0.4695 ± 0.0440    delta1
    gmr_qgamma0 = -1.2595          # -1.2595 ± 0.0026    beta0
    gmr_qgammaz = +0.0611          # +0.0611 ± 0.0021    beta1
    gmr_qgammam = -0.0477          # -0.0477 ± 0.0013    gamma
    gmr_mref_delta = 2.0e11 * MSOL   # fixed value
    gmr_mref_gamma = 1.0e10 * MSOL   # fixed value

    # M-MBulge
    mmbulge_mass_amp_log10 = 8.69             # 8.69 ± 0.05  [log10(M/Msol)]  approximate uncertainties!
    mmbulge_mass_ref = MSOL * 1e11            # 1e11 Msol
    mmbulge_plaw = 1.17                       # 1.17 ± 0.08
    mmbulge_scatter_dex = 0.28                # scatter stdev in dex
    mmbulge_bulge_frac = 0.7

    def __init__(self):
        self.mass = np.logspace(*np.log10([self.mass_pars[0], self.mass_pars[1]]), self.mass_pars[2])
        # self.mrat = np.logspace(self.mrat_pars[0], self.mrat_pars[1], self.mrat_pars[2])
        self.redz = np.logspace(*np.log10([self.redz_pars[0], self.redz_pars[1]]), self.redz_pars[2])
        self.grid3d = np.meshgrid(self.mass, self.mass, self.redz, indexing='ij')
        return

    @property
    def edges_3d(self):
        return (self.mass, self.mass, self.redz)

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


    def number_expect_4d_gwonly_instant(self, fobs_gw_edges, edges_3d, ndens_3d):
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

        # print(f"{np.shape(diff_numb)=} : {utils.stats(diff_numb)}")
        # print(f"{np.shape(ndens_3d)=} : {utils.stats(ndens_3d)}")
        # print(f"{np.shape(cosmo_fact)=} : {utils.stats(cosmo_fact)}")
        # print(f"{np.shape(tauf)=} : {utils.stats(tauf/GYR)} [Gyr]")

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


    # def



'''
    # ---- Choose the binary separations over which to integrate the binary evolution.

    # Start at large separations (galaxy merger) and evolve to small separations (coalescense).

    # start from the hardening model's initial separation
    rmax = hard._sepa_init
    # end at the ISCO
    # (M,)
    rmin = utils.rad_isco(self.mtot)
    # Choose steps for each binary, log-spaced between rmin and rmax
    extr = np.log10([rmax * np.ones_like(rmin), rmin])     # (2,M,)
    rads = np.linspace(0.0, 1.0, steps+1)[np.newaxis, :]     # (1,S)
    # (M, S)  <==  (M,1) * (1,S)
    rads = extr[0][:, np.newaxis] + (extr[1] - extr[0])[:, np.newaxis] * rads
    rads = 10.0 ** rads

    # ---- Calculate binary hardening rate (da/dt) at each separation, for each grid point

    # broadcast arrays to a consistent shape
    # (M, Q, S)
    mt, mr, rads, norm = np.broadcast_arrays(
        self.mtot[:, np.newaxis, np.newaxis],
        self.mrat[np.newaxis, :, np.newaxis],
        rads[:, np.newaxis, :],
        hard._norm[:, :, np.newaxis],
    )
    # calculate hardening rate (negative values, in units of [cm/s])
    dadt_evo = hard.dadt(mt, mr, rads, norm=norm)

    # ---- Integrate evolution
    # to find times and redshifts at which binaries reach each separation

    # (M, Q, S-1)
    # Integrate (inverse) hardening rates to calculate total lifetime to each separation
    times_evo = -utils.trapz_loglog(-1.0 / dadt_evo, rads, axis=-1, cumsum=True)
    # ~~~~ RIEMANN integration ~~~~
    # times_evo = 2.0 * np.diff(rads, axis=-1) / (dadt_evo[..., 1:] + dadt_evo[..., :-1])
    # times_evo = np.cumsum(times_evo, axis=-1)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # add array of zero time-delays at starting point (i.e. before the first step)
    # with same shape as a slice at a single step
    zpad = np.zeros_like(times_evo[..., 0])
    times_evo = np.concatenate([zpad[..., np.newaxis], times_evo], axis=-1)

    # ---- Convert from time to redshift

    # initial redshift (of galaxy merger)
    rz = self.redz[np.newaxis, np.newaxis, :, np.newaxis]    # (1, 1, Z, 1)

    tlbk_init = cosmo.z_to_tlbk(rz)
    tlbk = tlbk_init - times_evo[:, :, np.newaxis, :]
    # Combine the binary-evolution time, with the galaxy-merger time (if it is defined)
    if self._gmt_time is not None:
        tlbk -= self._gmt_time[:, :, :, np.newaxis]

    # (M, Q, Z, S)
    redz_evo = cosmo.tlbk_to_z(tlbk)

    #! age of the universe version of calculation is MUCH less accurate !#
    # Use age-of-the-universe
    # times_tot = times_evo[:, :, np.newaxis, :]
    # # Combine the binary-evolution time, with the galaxy-merger time (if it is defined)
    # if self._gmt_time is not None:
    #     times_tot += self._gmt_time[:, :, :, np.newaxis]
    # redz_evo = utils.redz_after(times_tot, redz=rz)
    #! ---------------------------------------------------------------- !#

    # ---- interpolate to target frequencies

    # convert from separations to rest-frame orbital frequencies
    # (M, Q, S)
    frst_orb_evo = utils.kepler_freq_from_sepa(mt, rads)
    # (M, Q, Z, S)
    fobs_orb_evo = frst_orb_evo[:, :, np.newaxis, :] / (1.0 + redz_evo)

    # (M, Q, Z, S)  ==>  (M*Q*Z, S)
    fobs_orb_evo, redz_evo = [tt.reshape(-1, steps+1) for tt in [fobs_orb_evo, redz_evo]]
    # `ndinterp` interpolates over 1th dimension
    # (M*Q*Z, X)
    redz_final = utils.ndinterp(fobs_orb, fobs_orb_evo, redz_evo, xlog=True, ylog=False)

    # (M, Q, Z, X)  <===  (M*Q*Z, X)
    redz_final = redz_final.reshape(self.shape + (fobs_orb.size,))

    coal = (redz_final > 0.0)
    frst_orb = fobs_orb * (1.0 + redz_final)
    frst_orb[frst_orb < 0.0] = 0.0
    redz_final[~coal] = -1.0

    # (M, Q, Z, X) comoving-distance in [Mpc]
    dc = np.zeros_like(redz_final)
    dc[coal] = cosmo.comoving_distance(redz_final[coal]).to('Mpc').value

    # (M, Q, Z, X) this is `(dVc/dz) * (dz/dt)` in units of [Mpc^3/s]
    cosmo_fact = np.zeros_like(redz_final)
    cosmo_fact[coal] = 4 * np.pi * (SPLC/MPC) * np.square(dc[coal]) * (1.0 + redz_final[coal])

    # (M, Q) calculate chirp-mass
    mt = self.mtot[:, np.newaxis, np.newaxis, np.newaxis]
    mr = self.mrat[np.newaxis, :, np.newaxis, np.newaxis]

    # Convert from observer-frame orbital freq, to rest-frame orbital freq
    sa = utils.kepler_sepa_from_freq(mt, frst_orb)
    mt, mr, sa, norm = np.broadcast_arrays(mt, mr, sa, hard._norm[:, :, np.newaxis, np.newaxis])
    # hardening rate, negative values, units of [cm/sec]
    dadt = hard.dadt(mt, mr, sa, norm=norm)
    # Calculate `tau = dt/dlnf_r = f_r / (df_r/dt)`
    # dfdt is positive (increasing frequency)
    dfdt, frst_orb = utils.dfdt_from_dadt(dadt, sa, frst_orb=frst_orb)
    tau = frst_orb / dfdt

    # (M, Q, Z, X) units: [1/s] i.e. number per second
    dnum = dens[..., np.newaxis] * cosmo_fact * tau
    dnum[~coal] = 0.0

    if details:
        tau[~coal] = 0.0
        dadt[~coal] = 0.0
        sa[~coal] = 0.0
        cosmo_fact[~coal] = 0.0
        # (M, Q, X)  ==>  (M, Q, Z, X)
        dets = dict(tau=tau, cosmo_fact=cosmo_fact, dadt=dadt, fobs=fobs_orb, sepa=sa)
        return edges, dnum, redz_final, dets

    return edges, dnum, redz_final
'''
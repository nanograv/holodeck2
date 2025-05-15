/**
 * Holodeck2 physics module.
 *
 */

#pragma once
#include "constants.h"

namespace physics {

    inline double _C5 = pow(SPLC, 5);
    inline double _G3 = pow(NWTG, 3);

    inline double _GW_CONST = (64.0 / 5.0) * _G3 / _C5;
    // [Sesana2004]_ Eq.36
    inline double _GW_SRC_CONST = (8.0 / sqrt(10.0)) * pow(NWTG, 5.0/3.0) * pow(PI, 2.0/3.0) / pow(SPLC, 4);
    inline double _GW_DADT_SEP_CONST = - (64.0 / 5.0) * _G3 / _C5;
    inline double _GW_DEDT_ECC_CONST = - (304.0 / 15.0) * _G3 / _C5;
    // [EN2007]_, Eq.2.2
    inline double _GW_LUM_CONST = (32.0 / 5.0) * pow(NWTG, 7.0/3.0) / _C5;
    // [EN2007]_, Eq.2.9
    inline double _GW_TAU_CONST = (5.0 / 96.0) * pow(NWTG, -5.0/3.0) * _C5 * pow(TWO_PI, -8.0/3.0);

    // double _AGE_UNIVERSE_GYR = cosmo.age(0.0).to('Gyr').value;  // [Gyr]  ~ 13.78



    /*
    # ==============================================================================
    # ====    Basics    ====
    # ==============================================================================


    def dfdt_from_dadt(dadt, sepa, mtot=None, frst_orb=None):
        """Convert from hardening rate in separation to hardening rate in frequency.

        Parameters
        ----------
        dadt : array_like
            Hardening rate in terms of binary separation.
        sepa : array_like
            Binary separations.
        mtot : None or array_like
            Binary total-mass in units of [gram].
            Either `mtot` or `frst_orb` must be provided.
        frst_orb : None or array_like
            Binary rest-frame orbital-frequency in units of [1/sec].
            Either `mtot` or `frst_orb` must be provided.

        Returns
        -------
        dfdt :
            Hardening rate in terms of rest-frame frequency.  [1/sec^2]
            NOTE: Has the opposite sign as `dadt`.
        frst_orb :
            Orbital frequency, in the rest-frame.  [1/sec]

        """
        if (mtot is None) and (frst_orb is None):
            err = "Either `mtot` or `frst_orb` must be provided!"
            raise ValueError(err)
        if frst_orb is None:
            frst_orb = kepler_freq_from_sepa(mtot, sepa)

        dfdt = - 1.5 * (frst_orb / sepa) * dadt
        return dfdt, frst_orb


    def mtmr_from_m1m2(m1, m2=None):
        """Convert from primary and secondary masses into total-mass and mass-ratio.

        NOTE: it doesn't matter if `m1` or `m2` is the primary or secondary.

        Parameters
        ----------
        m1 : array_like,
            Mass.  If this is a single value, or a 1D array, it denotes the mass of one component of
            a binary.  It can also be shaped, (N,2) where the two elements are the two component masses.
        m2 : None or array_like,
            If array_like, it must match the shape of `m1`, and corresponds to the companion mass.

        Returns
        -------
        (2,N) ndarray
            Total mass and mass-ratio.  If the input values are floats, this is just shaped (2,).

        """
        if m2 is not None:
            masses = np.stack([m1, m2], axis=-1)
        else:
            assert np.shape(m1)[-1] == 2, "If only `m1` is given, last dimension must be 2!"
            masses = np.asarray(m1)

        mtot = masses.sum(axis=-1)
        mrat = masses.min(axis=-1) / masses.max(axis=-1)
        return np.array([mtot, mrat])


    def m1m2_from_mtmr(mt, mr):
        """Convert from total-mass and mass-ratio to individual masses.

        Parameters
        ----------
        mt : array_like
            Total mass of the binary.
        mr : array_like
            Mass ratio of the binary.

        Returns
        -------
        (2,N) ndarray
            Primary and secondary masses respectively.
            0-primary (more massive component),
            1-secondary (less massive component)

        """
        mt = np.asarray(mt)
        mr = np.asarray(mr)
        m1 = mt/(1.0 + mr)
        m2 = mt - m1
        return np.array([m1, m2])


    def frst_from_fobs(fobs, redz, nharm=1):
        """Calculate rest-frame frequency from observed frequency and redshift.

        Parameters
        ----------
        fobs : array_like
            Observer-frame frequencies.
        redz : array_like
            Redshifts.
        nharm : array_like
            Harmonic of orbital frequency to calculate.
            If `nharm != 1`, then convert from the given harmonic to the orbital frequency.
            e.g. `nharm = 2` will convert from observer-frame GW-frequency to rest-frame orbital frequency,
            assuming circular orbits.

        Returns
        -------
        fobs : array_like
            Rest-frame frequencies.

        """
        frst = fobs * (1.0 + redz) / nharm
        return frst


    def fobs_from_frst(frst, redz, nharm=1):
        """Calculate observed frequency from rest-frame frequency and redshift.

        Parameters
        ----------
        frst : array_like
            Rest-frame frequencies.
        redz : array_like
            Redshifts.
        nharm : array_like
            Harmonic of orbital frequency to calculate.
            If `nharm != 1`, then convert from orbital frequency to given harmonic.
            e.g. `nharm = 2` will convert from rest-frame orbital frequency to observer-frame GW-frequency,
            assuming circular orbits.

        Returns
        -------
        fobs : array_like
            Observer-frame frequencies.

        """
        fobs = frst * nharm / (1.0 + redz)
        return fobs


    def kepler_freq_from_sepa(mass, sepa):
        """Calculate binary orbital frequency using Kepler's law.

        Parameters
        ----------
        mass : array_like
            Binary total mass [grams].
        sepa : array_like
            Binary semi-major axis or separation [cm].

        Returns
        -------
        freq : array_like
            Binary orbital frequency [1/s].

        """
        freq = np.sqrt(NWTG*mass) / np.power(sepa, 1.5) / (TWO_PI)
        return freq


    def kepler_sepa_from_freq(mass, freq):
        """Calculate binary separation using Kepler's law.

        Parameters
        ----------
        mass : array_like
            Binary total mass [grams]
        freq : array_like
            Binary orbital frequency [1/s].

        Returns
        -------
        sepa : array_like
            Binary semi-major axis (i.e. separation) [cm].

        """
        mass = np.asarray(mass)
        freq = np.asarray(freq)
        sepa = np.power(NWTG*mass/np.square(TWO_PI*freq), 1.0/3.0)
        return sepa


    def frst_isco(m1, m2=0.0, **kwargs):
        """Get rest-frame orbital frequency of ISCO orbit.

        Arguments
        ---------
        m1 : array_like, units of [gram]
            Total mass, or mass of the primary.  Added together with `m2` to get total mass.
        m2 : array_like, units of [gram]  or  None
            Mass of secondary, or None if `m1` is already total mass.

        Returns
        -------
        fisco : array_like, units of [Hz]

        """
        risco = rad_isco(m1, m2, **kwargs)
        fisco = kepler_freq_from_sepa(m1+m2, risco)
        return fisco


    '''
    def redz_after(time, redz=None, age=None):
        """Calculate the redshift after the given amount of time has passed.

        Parameters
        ----------
        time : array_like, [s]
            Amount of time to pass, in units of seconds.
        redz : None  or  array_like, []
            Redshift of starting point after which `time` is added.  Unitless.
        age : None  or  array_like, [s]
            Age of the Universe at the starting point, after which `time` is added.  Units of seconds.

        Returns
        -------
        new_redz : array_like, []
            Redshift of the Universe after the given amount of time.  Unitless

        """
        if (redz is None) == (age is None):
            raise ValueError("One of `redz` and `age` must be provided (and not both)!")

        if redz is not None:
            age = cosmo.age(redz).to('s').value
        new_age = age + time

        if np.isscalar(new_age):
            if new_age < _AGE_UNIVERSE_GYR * GYR:
                new_redz = cosmo.tage_to_z(new_age)
            else:
                new_redz = -1.0

        else:
            new_redz = -1.0 * np.ones_like(new_age)
            idx = (new_age < _AGE_UNIVERSE_GYR * GYR)
            new_redz[idx] = cosmo.tage_to_z(new_age[idx])

        return new_redz
    '''


    def pta_freqs(dur=16.03*YR, num=40, cad=None):
        """Get Fourier frequency bin specifications for the given parameters.

        Arguments
        ---------
        dur : double,
            Total observing duration, which determines the minimum sensitive frequency, ``1/dur``.
            Typically `dur` should be given in units of [sec], such that the returned frequencies are
            in units of [1/sec] = [Hz]
        num : int,
            Number of frequency bins.  If `cad` is not None, then the number of frequency bins is
            determined by `cad` and the `num` value is disregarded.
        cad : double or `None`,
            Cadence of observations, which determines the maximum sensitive frequency (i.e. the Nyquist
            frequency).  If `cad` is not given, then `num` frequency bins are constructed.

        Returns
        -------
        cents : (F,) ndarray
            Bin-center frequencies for `F` bins.  The frequency bin centers are at:
            ``F_i = (i + 1.5) / dur`` for i between 0 and `num-1`.
            The number of frequency bins, `F` is the argument `num`,
            or determined by `cad` if it is given.
        edges : (F+1,) ndarray
            Bin-edge frequencies for `F` bins, i.e. `F+1` bin edges.  The frequency bin edges are at:
            ``F_i = (i + 1) / dur`` for i between 0 and `num`.
            The number of frequency bins, `F` is the argument `num`,
            or determined by `cad` if it is given.

        """
        fmin = 1.0 / dur
        if cad is not None:
            num = dur / (2.0 * cad)
            num = int(np.floor(num))

        cents = np.arange(1, num+2) * fmin

        edges = cents - fmin / 2.0
        cents = cents[:-1]
        return cents, edges


    def rad_isco(m1, m2=0.0, factor=3.0):
        """Inner-most Stable Circular Orbit, radius at which binaries 'merge'.

        ENH: allow single (total) mass argument.
        ENH: add function to calculate factor as a function of BH spin.

        Parameters
        ----------
        m1 : array_like,
            Mass of first (either) component of binary [grams].
        m2 : array_like,
            Mass of second (other) component of binary [grams].
        factor : double,
            Factor by which to multiple the Schwarzschild radius to define the ISCO.
            3.0 for a non-spinning black-hole.

        Returns
        -------
        rs : array_like,
            Radius of the inner-most stable circular orbit [cm].

        """
        return factor * schwarzschild_radius(m1+m2)


    def schwarzschild_radius(mass):
        """Return the Schwarschild radius [cm] for the given mass [grams].

        Parameters
        ----------
        m1 : array_like
            Mass [grams]

        Returns
        -------
        rs : array_like,
            Schwzrschild radius for this mass.

        """
        rs = SCHW * mass
        return rs


    def velocity_orbital(mt, mr, per=None, sepa=None):
        sepa, per = _get_sepa_freq(mt, sepa, per)
        v2 = np.power(NWTG*mt/sepa, 1.0/2.0) / (1 + mr)
        # v2 = np.power(2*np.pi*NWTG*mt/per, 1.0/3.0) / (1 + mr)
        v1 = v2 * mr
        vels = np.moveaxis([v1, v2], 0, -1)
        return vels


    def lambda_factor_dlnf(frst, dfdt, redz, dcom=None):
        """Account for the universe's differential space-time volume for a given hardening rate.

        For each binary, calculate the factor: $$\\Lambda \\equiv (dVc/dz) * (dz/dt) * [dt/dln(f)]$$,
        which has units of [Mpc^3].  When multiplied by a number-density [Mpc^-3], it gives the number
        of binaries in the Universe *per log-frequency interval*.  This value must still be multiplied
        by $\\Delta \\ln(f)$ to get a number of binaries across a frequency in.

        Parameters
        ----------
        frst : ArrayLike
            Binary frequency (typically rest-frame orbital frequency; but it just needs to match what's
            provided in the `dfdt` term.  Units of [1/sec].
        dfdt : ArrayLike
            Binary hardening rate in terms of frequency (typically rest-frame orbital frequency, but it
            just needs to match what's provided in `frst`).  Units of [1/sec^2].
        redz : ArrayLike
            Binary redshift.  Dimensionless.
        dcom : ArrayLike
            Comoving distance to binaries (for the corresponding redshift, `redz`).  Units of [cm].
            If not provided, calculated from given `redz`.

        Returns
        -------
        lambda_fact : ArrayLike
            The differential comoving volume of the universe per log interval of binary frequency.

        """
        zp1 = redz + 1
        if dcom is None:
            dcom = cosmo.z_to_dcom(redz)

        # Volume-factor
        # this is `(dVc/dz) * (dz/dt)`,  units of [Mpc^3/s]
        vfac = 4.0 * np.pi * SPLC * zp1 * (dcom**2)
        # Time-factor
        # this is `f / (df/dt) = dt/d ln(f)`,  units of [sec]
        tfac = frst / dfdt

        # Calculate weighting
        lambda_fact = vfac * tfac
        return lambda_fact


    # =================================================================================================
    # ====    Gravitational Waves    ====
    # =================================================================================================


    def chirp_mass(m1, m2=None):
        """Calculate the chirp-mass of a binary.

        Parameters
        ----------
        m1 : array_like,
            Mass [grams]
            This can either be the mass of the primary component, if scalar or 1D array_like,
            or the mass of both components, if 2D array_like, shaped (N, 2).
        m2 : None or array_like,
            Mass [grams] of the other component of the binary.  If given, the shape must be
            broadcastable against `m1`.

        Returns
        -------
        mc : array_like,
            Chirp mass [grams] of the binary.

        """
        # (N, 2)  ==>  (N,), (N,)
        if m2 is None:
            assert (np.ndim(m1) == 2) and (np.shape(m2)[-1] == 2)
            m1, m2 = np.moveaxis(m1, -1, 0)
        mc = np.power(m1 * m2, 3.0/5.0)/np.power(m1 + m2, 1.0/5.0)
        return mc


    def chirp_mass_mtmr(mt, mr):
        """Calculate the chirp-mass of a binary.

        Parameters
        ----------
        mt : array_like,
            Total mass [grams].  This is ``M = m1+m2``.
        mr : array_like,
            Mass ratio.  ``q = m2/m1 <= 1``.
            This is defined as the secondary (smaller) divided by the primary (larger) mass.

        Returns
        -------
        mc : array_like,
            Chirp mass [grams] of the binary.

        """
        mt, mr = _array_args(mt, mr)
        mc = mt * np.power(mr, 3.0/5.0) / np.power(1 + mr, 6.0/5.0)
        return mc


    def gw_freq_dist_func(nn, ee=0.0, recursive=True):
        """GW frequency distribution function.

        See [EN2007]_ Eq. 2.4; this function gives g(n,e).

        NOTE: recursive relation fails for zero eccentricities!
        TODO: could choose to use non-recursive when zero eccentricities are found?

        TODO: replace `ee` variable with `eccen`

        Parameters
        ----------
        nn : int,
            Number of frequency harmonic to calculate.
        ee : array_like,
            Binary eccentricity.

        Returns
        -------
        gg : array_like,
            GW Frequency distribution function g(n,e).

        """

        # Calculate with non-zero eccentrictiy
        bessel = sp.special.jn
        ne = nn*ee
        n2 = np.square(nn)
        jn_m2 = bessel(nn-2, ne)
        jn_m1 = bessel(nn-1, ne)

        # Use recursion relation:
        if recursive:
            jn = (2*(nn-1) / ne) * jn_m1 - jn_m2
            jn_p1 = (2*nn / ne) * jn - jn_m1
            jn_p2 = (2*(nn+1) / ne) * jn_p1 - jn
        else:
            jn = bessel(nn, ne)
            jn_p1 = bessel(nn+1, ne)
            jn_p2 = bessel(nn+2, ne)

        aa = np.square(jn_m2 - 2.0*ee*jn_m1 + (2/nn)*jn + 2*ee*jn_p1 - jn_p2)
        bb = (1 - ee*ee)*np.square(jn_m2 - 2*ee*jn + jn_p2)
        cc = (4.0/(3.0*n2)) * np.square(jn)
        gg = (n2*n2/32) * (aa + bb + cc)
        return gg


    def gw_hardening_rate_dadt(m1, m2, sepa, eccen=None):
        """GW Hardening rate in separation (da/dt).

        NOTE: returned value is negative.

        See [Peters1964]_, Eq. 5.6

        Parameters
        ----------
        m1 : array_like,
            Mass of one component of the binary [grams].
        m2 : array_like,
            Mass of other component of the binary [grams].
        sepa : array_like,
            Binary semi-major axis (i.e. separation) [cm].
        eccen : None or array_like,
            Binary orbital eccentricity.  Treated as zero if `None`.

        Returns
        -------
        dadt : array_like,
            Binary hardening rate [cm/s] due to GW emission.

        """
        dadt = _GW_DADT_SEP_CONST * m1 * m2 * (m1 + m2) / np.power(sepa, 3)
        if eccen is not None:
            fe = _gw_ecc_func(eccen)
            dadt *= fe
        return dadt


    '''
    def gw_hardening_rate_dfdt(m1, m2, frst_orb, eccen=None):
        """GW Hardening rate in frequency (df/dt).

        Parameters
        ----------
        m1 : array_like
            Mass of one component of each binary [grams].
        m2 : array_like
            Mass of other component of each binary [grams].
        freq_orb : array_like
            Rest frame orbital frequency of each binary [1/s].
        eccen : array_like, optional
            Eccentricity of each binary.

        Returns
        -------
        dfdt : array_like,
            Hardening rate in terms of frequency for each binary [1/s^2].

        """
        sepa = kepler_sepa_from_freq(m1+m2, frst_orb)
        dfdt = gw_hardening_rate_dadt(m1, m2, sepa, eccen=eccen)
        # dfdt, _ = dfdt_from_dadt(dfdt, sepa, mtot=m1+m2)
        dfdt, _ = dfdt_from_dadt(dfdt, sepa, frst_orb=frst_orb)
        return dfdt, frst_orb
    '''


    def gw_hardening_time_freq(m1, m2, frst_orb, eccen=None):
        """GW Hardening timescale in frequency: $f / [df/dt]$.

        [EN2007] Eq. 2.9

        Parameters
        ----------
        m1 : array_like
            Mass of one component of each binary [grams].
        m2 : array_like
            Mass of other component of each binary [grams].
        freq_orb : array_like
            Rest frame orbital frequency of each binary [1/s].
        eccen : array_like, optional
            Eccentricity of each binary.

        Returns
        -------
        tauf : array_like,
            Hardening timescale in terms of frequency for each binary [sec].

        """
        tauf = _GW_TAU_CONST * np.power(chirp_mass(m1, m2), -5.0/3.0) * np.power(frst_orb, -8.0/3.0)
        if eccen is not None:
            tauf /= _gw_ecc_func(eccen)
        return tauf


    def gw_lum_circ(mchirp, freq_orb_rest):
        """Calculate the GW luminosity of a circular binary.

        [EN2007]_ Eq. 2.2

        Parameters
        ----------
        mchirp : array_like,
            Binary chirp mass [grams].
        freq_orb_rest : array_like,
            Rest-frame binary orbital frequency [1/s].

        Returns
        -------
        lgw_circ : array_like,
            GW Luminosity [erg/s].

        """
        lgw_circ = _GW_LUM_CONST * np.power(2.0*np.pi*freq_orb_rest*mchirp, 10.0/3.0)
        return lgw_circ


    def gw_strain_source(mchirp, dcom, freq_rest_orb):
        """GW Strain from a single source in a circular orbit.

        For reference, see:
        *   [Sesana2004]_ Eq.36 : they use `f_r` to denote rest-frame GW-frequency.
        *   [Enoki2004]_ Eq.5.

        Parameters
        ----------
        mchirp : array_like,
            Binary chirp mass [grams].
        dcom : array_like,
            Comoving distance to source [cm].
        freq_orb_rest : array_like,
            Rest-frame binary orbital frequency [1/s].

        Returns
        -------
        hs : array_like,
            GW Strain (*not* characteristic strain).

        """
        # The factor of 2 below is to convert from orbital-frequency to GW-frequency
        hs = _GW_SRC_CONST * mchirp * np.power(2*mchirp*freq_rest_orb, 2/3) / dcom
        return hs


    def sep_to_merge_in_time(m1, m2, time):
        """The initial separation required to merge within the given time.

        See: [Peters1964]_

        Parameters
        ----------
        m1 : array_like,
            Mass of one component of the binary [grams].
        m2 : array_like,
            Mass of other component of the binary [grams].
        time : array_like,
            The duration of time of interest [sec].

        Returns
        -------
        array_like
            Initial binary separation [cm].

        """
        a1 = rad_isco(m1, m2)
        return np.power(_GW_CONST*m1*m2*(m1+m2)*time - np.power(a1, 4.0), 1./4.)


    def time_to_merge_at_sep(m1, m2, sepa):
        """The time required to merge starting from the given initial separation.

        See: [Peters1964]_.

        Parameters
        ----------
        m1 : array_like,
            Mass of one component of the binary [grams].
        m2 : array_like,
            Mass of other component of the binary [grams].
        sepa : array_like,
            Binary semi-major axis (i.e. separation) [cm].

        Returns
        -------
        array_like
            Duration of time for binary to coalesce [sec].

        """
        a1 = rad_isco(m1, m2)
        delta_sep = np.power(sepa, 4.0) - np.power(a1, 4.0)
        return delta_sep/(_GW_CONST*m1*m2*(m1+m2))


    def gamma_psd_to_strain(gamma_psd):
        gamma_strain = (gamma_psd + 3.0) / 2.0
        return gamma_strain


    def gamma_strain_to_psd(gamma_strain):
        gamma_psd = 2*gamma_strain - 3.0
        return gamma_psd


    def gamma_strain_to_omega(gamma_strain):
        gamma_omega = (gamma_strain - 2.0) / 2.0
        return gamma_omega


    def char_strain_to_psd(freqs, hc):
        """

        Arguments
        ---------
        freqs : array_like
            Frequencies of interest in [1/sec].
            Note: these should NOT be in units of reference frequency, but in units of [Hz] = [1/sec].
        hc : array_like
            Characteristic strain.

        Returns
        -------
        psd : array_like
            Power spectral density of gravitational waves.

        """
        psd = hc**2 / (12*np.pi**2)
        psd = psd * np.power(freqs, -3)
        # psd = psd * np.power(freqs/fref, -3) * np.power(fref, -3)
        return psd


    def psd_to_char_strain(freqs, psd):
        hc = np.sqrt(psd * (12*np.pi**2 * freqs**3))
        return hc


    def char_strain_to_rho(freqs, hc, tspan):
        psd = char_strain_to_psd(freqs, hc)
        rho = np.sqrt(psd/tspan)
        return rho


    def rho_to_char_strain(freqs, rho, tspan):
        psd = tspan * rho**2
        hc = psd_to_char_strain(freqs, psd)
        return hc


    def _gw_ecc_func(eccen):
        """GW Hardening rate eccentricitiy dependence F(e).

        See [Peters1964]_ Eq. 5.6, or [EN2007]_ Eq. 2.3

        Parameters
        ----------
        eccen : array_like,
            Binary orbital eccentricity [].

        Returns
        -------
        fe : array_like
            Eccentricity-dependence term of GW emission [].

        """
        e2 = eccen*eccen
        num = 1.0 + (73.0/24.0)*e2 + (37.0/96.0)*e2*e2
        den = np.power(1.0 - e2, 7.0/2.0)
        fe = num / den
        return fe

    */

} // namespace physics
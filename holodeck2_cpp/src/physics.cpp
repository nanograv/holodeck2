/**
 * Holodeck2 physics module.
 *
 */

#include "physics.h"


namespace {
    inline double C5 = pow(SPLC, 5);
    inline double G3 = pow(NWTG, 3);

    inline double GW_CONST = (64.0 / 5.0) * G3 / C5;
    // [Sesana2004] Eq.36
    inline double GW_SRC_CONST = (8.0 / sqrt(10.0)) * pow(NWTG, 5.0/3.0) * pow(PI, 2.0/3.0) / pow(SPLC, 4);
    inline double GW_SRC_CONST_SQ = pow(GW_SRC_CONST, 2);
    inline double GW_DADT_SEP_CONST = - (64.0 / 5.0) * G3 / C5;
    inline double GW_DEDT_ECC_CONST = - (304.0 / 15.0) * G3 / C5;
    // [EN2007], Eq.2.2
    inline double GW_LUM_CONST = (32.0 / 5.0) * pow(NWTG, 7.0/3.0) / C5;
    // [EN2007], Eq.2.9
    inline double GW_TAU_CONST = (5.0 / 96.0) * pow(NWTG, -5.0/3.0) * C5 * pow(TWO_PI, -8.0/3.0);

    // double AGE_UNIVERSE_GYR = cosmo.age(0.0).to('Gyr').value;  // [Gyr]  ~ 13.78
}


namespace physics {

    void chirp_mass_from_m1m2(double m1, double m2, double *mc) {
        *mc = pow(m1 * m2, 3.0/5.0) / pow(m1 + m2, 1.0/5.0);
    }

    void mtmr_from_m1m2(double m1, double m2, double* mtot, double* mrat) {
        // Convert from primary and secondary masses into total-mass and mass-ratio.
        *mtot = m1 + m2;
        *mrat = m2 / m1;
    }

    void m1m2_from_mtmr(double mt, double mr, double *m1, double* m2) {
        *m1 = mt/(1.0 + mr);
        *m2 = mt - (*m1);
    }

    void fobs_from_frst(double frst, double redz, double* fobs, int nharm) {
        *fobs = frst * nharm / (1.0 + redz);
    }

    void frst_from_fobs(double fobs, double redz, double* frst, int nharm) {
        *frst = fobs * (1.0 + redz) / nharm;
    }

    void gw_hardening_rate_dadt(double m1, double m2, double sepa, double* dadt) {
        *dadt = GW_DADT_SEP_CONST * m1 * m2 * (m1 + m2) / pow(sepa, 3);
    }

    /**
     * Calculate the hardening time in terms of frequency, $\tau_f = f/(df/dt) = dt/dln(f)$.
     *
     * @param m1          Mass of the primary [gram].
     * @param m2          Mass of the secondary [gram].
     * @param frst_orb    Rest-frame orbital frequency [1/sec].
     * @return tauf       Output pointer for the hardening time [sec].
     *
     */
    void gw_hardening_time_freq(double m1, double m2, double frst_orb, double *tauf) {
        double mc;
        chirp_mass_from_m1m2(m1, m2, &mc);
        *tauf = GW_TAU_CONST * pow(mc, -5.0/3.0) * pow(frst_orb, -8.0/3.0);
    }

    void gw_strain_source(double mchirp, double dcom, double freq_rest_orb, double* hs) {
        // The factor of 2 below is to convert from orbital-frequency to GW-frequency
        *hs = GW_SRC_CONST * pow(mchirp, 5.0/3.0) * pow(2*freq_rest_orb, 2.0/3.0) / dcom;
    }

    void gw_strain_source_sq(double mchirp, double dcom, double freq_rest_orb, double* hs2) {
        // The factor of 2 below is to convert from orbital-frequency to GW-frequency
        *hs2 = GW_SRC_CONST_SQ * pow(mchirp, 10.0/3.0) * pow(2*freq_rest_orb, 4.0/3.0) / pow(dcom, 2);
    }

    void kepler_forb_from_sepa(double mass, double sepa, double* forb) {
        *forb = sqrt(NWTG*mass) / pow(sepa, 1.5) / TWO_PI;
    }

    void kepler_sepa_from_freq(double mass, double freq, double* sepa) {
        *sepa = pow(NWTG*mass/pow(TWO_PI*freq, 2), 1.0/3.0);
    }

} // namespace physics






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
        if new_age < AGE_UNIVERSE_GYR * GYR:
            new_redz = cosmo.tage_to_z(new_age)
        else:
            new_redz = -1.0

    else:
        new_redz = -1.0 * np.ones_like(new_age)
        idx = (new_age < AGE_UNIVERSE_GYR * GYR)
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
    sepa, per = get_sepa_freq(mt, sepa, per)
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


def gw_freq_dist_func(nn, ee=0.0, recursive=True):
    """GW frequency distribution function.

    See [EN2007] Eq. 2.4; this function gives g(n,e).

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


def gw_lum_circ(mchirp, freq_orb_rest):
    """Calculate the GW luminosity of a circular binary.

    [EN2007] Eq. 2.2

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
    lgw_circ = GW_LUM_CONST * np.power(2.0*np.pi*freq_orb_rest*mchirp, 10.0/3.0)
    return lgw_circ


def sep_to_merge_in_time(m1, m2, time):
    """The initial separation required to merge within the given time.

    See: [Peters1964]

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
    return np.power(GW_CONST*m1*m2*(m1+m2)*time - np.power(a1, 4.0), 1./4.)


def time_to_merge_at_sep(m1, m2, sepa):
    """The time required to merge starting from the given initial separation.

    See: [Peters1964].

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
    return delta_sep/(GW_CONST*m1*m2*(m1+m2))


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


def gw_ecc_func(eccen):
    """GW Hardening rate eccentricitiy dependence F(e).

    See [Peters1964] Eq. 5.6, or [EN2007] Eq. 2.3

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


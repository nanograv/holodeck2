# Holodeck2: C++ - Overview

## Gravitational Wave Background (GWB)

The stochastic **Gravitational Wave Background (GWB)** is the superposition of GWs from individual **(super-)massive black-hole binaries (MBHBs)** throughout the universe.  In the binary population literature, the GWB is typically characterized in terms of the **characteristic strain ($h_c$)**.  The strain or **spectral strain ($h$)** is a dimensionless measure of the fractional perturbations to distances and times produced by GWs, such that $h \approx \Delta l / l \approx \Delta t / t$.  The characteristic strain is defined to be a quantity such the measured **signal-to-noise ratio (SNR)** will be equal to the signal characteristic-strain divided by the noise characteristic-strain, over some appropriate frequency band-pass (i.e. integrated over some range of frequencies).  Unfortunately, this means that characteristic strain must be calculated in different ways depending on the regime of the GW signal, and the type of GW detector.

For PTA sources that are monochromatic---where the GW frequency $f$ is effectively constant---the characteristic strain is equal to the spectral strain multiplied by the square-root of the number of cycles over which the signal is observed.  PTA signals are typically analyzed in the frequency domain, using Fourier transforms of time-domain signals.  Thus the natural frequency bin-width is the Rayleigh frequency---the inverse of the total observing time, $\Delta f = 1/T$.  Thus we can write, $h_c^2 = h^2 \cdot N = h^2 \cdot f \cdot T = h^2 \cdot f / \Delta f$.

The sky- and polarization- averaged spectral strain amplitude from a circular binary can be calculated as,
$$
h(f) = \frac{8}{10^{1/2}} \frac{\left(G \mathcal{M} \right)^{5/3}}{c^4 \, d_c} \left(2 \pi f_\mathrm{r,orb} \right)^{2/3},
$$
where $\mathcal{M}$ is the (rest-frame) chirp-mass, and $f_\mathrm{r,orb}$ is the rest-frame orbital frequency such that $f = 2 f_\mathrm{r,orb} / (1+z)$ for a circular orbit.

In the MBHB population literature, it is typical to calculate the GWB neglecting the effects of wave interference such that the squared-strain is the quadratic sum of strains from $N$ individual binaries.  In actuality, this is the expectation value for the GWB amplitude.  The GWB is then calculated as,
$$
h_c^2(f) = \int \left[ h^2(f_r) \, \frac{d^2 N}{dz \, d \ln\!f_r} \right]_{f_r = (1+z) f} \, dz.
$$
In practice, we are able to directly estimate the *(comoving-)number-density* of MBHBs $n_c$, and then extrapolate to the total number of binaries in the comoving-volume of the observer's past light-cone such that $dN = n \cdot dV_c$.  We can also calculate the frequency distribution of binaries by assuming an rate of binary evolution with respect to frequency, called a binary `**hardening**' rate, typically parameterized with a **hardening time** (with respect to frequency) $\tau_f \equiv dt/d\ln f = f / (df/dt)$.

Putting these together, we can convert between binary number and binary number-density as,
$$
\frac{d^2 N}{dz \, d\ln\!f}
    = \frac{d n}{dz} \frac{dt}{d\ln\!f} \frac{d V_c}{dz} \frac{dz}{dt} 
    = \frac{d n}{dz} \cdot \tau_f \cdot \left[4 \pi \, c \, d_c^2 \, (1+z)\right].
$$
The second relation in the preceding equation comes from standard cosmographic relations for the comoving-volume of the Universe as a function of time ($t$).

To calculate the GWB, we need to calculate the comoving number density of binaries $dn_c/dz$, and also the binary hardening rate $\tau_f$.  We need the binary number-density for every value of GW strain $h(f)$, which is determined not only by redshift (i.e. distance to source), but also the other binary parameters: particularly chirp-mass, frequency, and eccentricity.  In this discussion we will assume circular orbits, with zero eccentricity: $e = 0$.  Because the binary evolution can depend on combinations of the two binary masses differently than the chirp mass, we actually need to track two mass parameters.  Here we choose to track each MBH mass: $m_1, m_2$ (although total mass and mass ratio, for example, are equally valid).  Thus we need to define $\frac{d^4 N}{dm_1 \, dm_2 \, dz \, d \ln\!f}$.

## Semi-Analytic Model (SAM) Populations

**Semi-Analytic Models (SAMs)** use simple analytic prescriptions/relations for galaxies and galaxy mergers to then calculate the number of MBH binaries.  In general, the relationships chosen can be arbitrarily simple or complex---a function of as many or as few parameters as desired.  For numerical reasons it will often be far more accurate to describe quantities/distributions in log-space, i.e. in terms of,
$$
\frac{d^4 N}{d \log_{10}(m_1) \, d\log_{10}(m_2) \, dz \, d \ln\!f}.
$$

We choose to perform the calculation by calculating a galaxy merger-rate density and assume a one-to-one mapping (***??BH OCCUPATION FRACTION??***) from galaxies (and thus galaxy mergers) to MBHs (and thus MBH mergers) where the MBH masses are given by a MMBulge relation s.t. $M = M_\mathrm{BH}(M_\star)$.  The galaxy merger-rate density is calculated as the product of a (primary-)galaxy stellar-mass function (GSMF; $\Psi$) multiplied by a galaxy merger-rate ($\mathcal{R}$).  We describe the different components and their parametric dependencies in more detail below.

$$
\frac{d^3 n}{d\log_{10}(m_1) \, d\log_{10}(m_2) \, dz}
    = \frac{d^3 n_{\star\star}}{d\log_{10}(M_{\star 1}) \, d\log_{10}(M_{\star 2}) \, d t}
    \cdot \frac{dn}{dn_{\star \star}}
    \cdot \frac{d\log_{10}(M_{\star 1})}{d\log_{10}(M_1)}
    \cdot \frac{d\log_{10}(M_{\star 2})}{d\log_{10}(M_2)}
    \cdot \frac{dt}{dz} \\
    = \Psi \cdot \mathcal{R} \cdot \frac{d\log_{10}(M_{\star 1})}{d\log_{10}(M_1)}
    \cdot \frac{d\log_{10}(M_{\star 2})}{d\log_{10}(M_2)}
    \cdot \frac{dt}{dz}.
$$

### Galaxy Stellar-Mass Function (GSMF; $\Psi$)

The GSMF is implemented as a double-Schechter function, the sum of two individual Schechter functions, $\Phi(M_\star, z) = \Phi_1(M_\star, z) + \Phi_2(M_\star, z)$.  Each Schechter function is parameterized as,
$$
\Psi_i(M_{\star}, z) \equiv \frac{d n_c}{d \log_{10}(M_\star)} = \ln(10) \cdot \Phi_i(M_\star) \cdot \left(\frac{M_{\star}}{M_{\star,c}}\right)^{1 + \alpha_i}, \\
\log_{10}(\Phi_i[M_\star]) = \phi_{i,0} + \phi_{i,1} \cdot z + \phi_{i,2} \cdot z^2, \\
\log_{10}(M_{\star,c}) = \mu_0 + \mu_1 \cdot z + \mu_2 \cdot z^2.
$$
Here the subscript index $i \in [0, 1]$ denotes the two Schechter functions.  Note that the $\Phi_i(M_\star)$ and $\alpha_i$ are different for each, while the characteristics stellar-mass $M_{\star,c}$ is the same for both functions.  Thus there are 11 parameters in total to describe this GSMF.

### Galaxy Merger Rate (GMR; $\mathcal{R}$)

We define the galaxy merger rate as, for a given primary galaxy (with stellar mass $M_{\star,1}$, at a redshift $z_\star$), the number of mergers per unit time with secondary galaxies of mass $M_{\star 2}$, i.e.
$$
\mathcal{R} \equiv \frac{d^3 n_{\star \star}}{dn_{\star 1}(M_{\star 1}) \, d\log_{10}(M_{\star 2}) \, dt}.
$$
We use a parameterized fit for the galaxy merger rate from the Illustris simulations [Rodriguez-Gomez+2015], which are parameterized in terms of a differential with respect to the mass-ratio, i.e. $q_\star = M_{\star 2} / M_{\star 1}$, so we require an additional conversion as,
$$\mathcal{R} = \mathcal{R}_\mathrm{ill} \cdot \frac{dq_\star}{d\log_{10}(M_{\star 2})} = \mathcal{R}_\mathrm{ill} \cdot q_\star.$$



### MBH-Mass – Galaxy-Mass ($M_\mathrm{BH}$–$M_\mathrm{bulge}$) Relation

The MMBulge (pronounced, 'em-em-bulge') relation maps from a galaxy stellar-bulge mass to an MBH mass.  The stellar-bulge mass is typically calculated as a fraction of the total stellar-mass, and thus we often use 'MMBulge' even when referring to a conversion to/from the total stellar-mass.  We adopt the standard power-law form for the MMBulge relation,
$$
\log_{10}\left( \frac{M_\mathrm{BH}}{M_\odot}\right)
    = \mu + \alpha_\mu \cdot \log_{10}\left( \frac{M_\mathrm{bulge}}{10^{11} \, M_\odot} \right).
$$
Thus the differential relation can be written as,
$$
\frac{d\log_{10}(M_\mathrm{BH})}{d\log_{10}(M_\mathrm{bulge})} = \alpha_\mu,
$$
or equivalently,
$$
\frac{dM_\mathrm{BH}}{ dM_\mathrm{bulge}} = \alpha_\mu \cdot \frac{M_\mathrm{BH}}{M_\mathrm{bulge}},
$$



### Binary Evolution Rate ('Hardening Rate'; $\tau_f$)


## SAM GWB Calculation

$$
h_c^2(f) = \int \left[ h^2(f_r) \, \frac{d^4 N}{d\log_{10}(M_1) \; d\log_{10}(M_2) \; dz \; d \ln\!f_r} \right]_{f_r = (1+z) f} \, d\log_{10}(M_1) \;\, d\log_{10}(M_2) \;\, dz .
$$

$$
\frac{d^4 N}{d\log_{10}(M_1) \; d\log_{10}(M_2) \; dz \; d \ln\!f}
    = \Psi(M_{\star 1}, z_\star) \cdot \mathcal{R}_\mathrm{ill}(M_{\star 1}, q_\star, z_\star)
    \cdot q_\star \cdot \alpha_\mu^2 
    \cdot \tau_f \cdot \left[4 \pi \, c \, d_c^2 / H(z) \right].
$$

***!!Difference in times ($z$ vs $z_\star$) has not been addressed and is not written consistently!!***

## Reference Equations & Relations

Chirp mass, $\mathcal{M} \equiv \frac{\left(m_1 m_2\right)^{3/5}}{M^{1/5}} = M \frac{q^{3/5}}{\left(1 + q\right)^{6/5}}$.

Mass ratio, $q \equiv m_2/m_1 \leq 1$.  And thus: $dq/d\log_{10}(m_2) = m_2 / m_1 = q$.

Total mass $M \equiv m_1 + m_2$.

Component masses, primary mass: $m_1 = M / (1 + q)$, secondary mass: $m_2 = M \cdot q / (1 + q)$.

Kepler's Law: $(2\pi/T)^2 = (G \cdot M / a^3)$.

Cosmography: $dz/dt = (1+z) \cdot H(z)$, and $dV_c/dz = 4 \pi c \, d_c^2 / H(z)$.


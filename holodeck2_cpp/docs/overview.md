# Holodeck2: C++ - Overview

## Gravitational Wave Background (GWB)

The stochastic **Gravitational Wave Background (GWB)** is the superposition of GWs from individual **(super-)massive black-hole binaries (MBHBs)** throughout the universe.  In the binary population literature, the GWB is typically characterized in terms of the **characteristic strain ($h_c$)**.  The strain or **spectral strain ($h$)** is a dimensionless measure of the fractional perturbations to distances and times produced by GWs, such that $h \approx \Delta l / l \approx \Delta t / t$.  The characteristic strain is defined to be a quantity such the measured **signal-to-noise ratio (SNR)** will be equal to the signal characteristic-strain divided by the noise characteristic-strain, over some appropriate frequency band-pass (i.e. integrated over some range of frequencies).  Unfortunately, this means that characteristic strain must be calculated in different ways depending on the regime of the GW signal, and the type of GW detector.

For PTA sources that are monochromatic---where the GW frequency $f$ is effectively constant---the characteristic strain is equal to the spectral strain multiplied by the square-root of the number of cycles over which the signal is observed.  PTA signals are typically analyzed in the frequency domain, using Fourier transforms of time-domain signals.  Thus the natural frequency bin-width is the Rayleigh frequency---the inverse of the total observing time, $\Delta f = 1/T$.  Thus we can write, $h_c^2 = h^2 \cdot N = h^2 \cdot f \cdot T = h^2 \cdot f / \Delta f$.

In the MBHB population literature, it is typical to calculate the GWB neglecting the effects of wave interference such that the squared-strain is the quadratic sum of strains from $N$ individual binaries.  In actuality, this is the expectation value for the GWB amplitude.  The GWB is then calculated as,
$$
h_c^2(f) = \int \left[ h^2(f_r) \, \frac{d^2 N}{dz \, d \ln f_r} \right]_{f_r = (1+z) f} \, dz.
$$
In practice, we are able to directly estimate the *(comoving-)number-density* of MBHBs $n_c$, and then extrapolate to the total number of binaries in the comoving-volume of the observer's past light-cone such that $dN = n_c \cdot dV_c$.  

## Semi-Analytic Model (SAM) Populations

### Galaxy Stellar-Mass Function (GSMF)

The GSMF is implemented as a double-Schechter function, the sum of two individual Schechter functions, $Phi_\star(M_\star, z) = \Phi_1(M_\star, z) + \Phi_2(M_\star, z)$.  Each Schechter function is parameterized as,
$$
\Psi_i(M_{\star}, z) \equiv \frac{d n_c}{d \log_{10}(M_\star)} = \ln(10) \cdot \Phi_i(M_\star) \cdot \left(\frac{M_{\star}}{M_{\star,c}}\right)^{1 + \alpha_i}, \\
\log_{10}(\Phi_i[M_\star]) = \phi_{i,0} + \phi_{i,1} \cdot z + \phi_{i,2} \cdot z^2, \\
\log_{10}(M_{\star,c}) = \mu_0 + \mu_1 \cdot z + \mu_2 \cdot z^2.
$$
Here the subscript index $i \in [0, 1]$ denotes the two Schechter functions.  Note that the $\Phi_i(M_\star)$ and $\alpha_i$ are different for each, while the characteristics stellar-mass $M_{\star,c}$ is the same for both functions.  Thus there are 11 parameters in total to describe this GSMF.

sd
# Design of `QEDData`

## LHAPDF grid for LDFs

Lepton distribution functions (LDFs, and fragmentation functions LFFs) peak at both small-$x$ and large-$x$. In DIS experiments like COMPASS, we usually don't have access to the small-$x$ region, so it is the large-$x$ region that gives dominant contributions.

LHAPDF interpolates in the $\log Q^2\text{-}\log x$ space by default, which is a suitable choice for PDFs that peak at small-$x$. To store an LDF $f(x)$ in an LHAPDF grid, we then choose to store:
$$
g(x) \equiv (1-x) f(1-x).
$$
The inverse relation is
$$
f(x) = \frac{1}{x} g(1-x).
$$

## LHAPDF grid for integrals of LDFs

LDFs actually have "distributional values" at $x=1$, and to calculate their integrals with other functions, we use the subtraction trick:
$$
\int_{x_\mathrm{min}}^1 dx~ f(x) H(x) = \int_{x_\mathrm{min}}^1 dx~ f(x) \big(H(x)-H(1)\big) + H(1) \int_{x_\mathrm{min}}^1 dx~ f(x).
$$
The integral in first term of this formula in suppressed at $x=1$, and we can treat $f(x)$ as an ordinary function, that is, the "distributional value" at $x=1$ is thrown away. To calculate the second term, we treat $f(x)$ as a function on $\mathbb{R}$ with vanishing values outside the physical region, and use Mellin transformation:
$$\begin{aligned}
\int_{x_\mathrm{min}}^1 dx~ f(x)
&= \int_{x_\mathrm{min}}^\infty dx~ \mathcal{M}^{-1}[\mathcal{M}[f]](x) \\
&= \frac{1}{2\pi i} \int ds~ \mathcal{M}[f](s) \int_{x_\mathrm{min}}^\infty dx~ x^{-s} \\
&= \frac{x_\mathrm{min}}{2\pi i} \int ds~ x_\mathrm{min}^{-s} \frac{\mathcal{M}[f](s)}{s-1}.
\end{aligned}$$

It will be convenient to also store these integrals of LDFs in an LHAPDF grid. Since we also want higher precision at large-$x$, we choose to store
$$
g(x) \equiv 1 - \int_{1-x}^1 dx'f(x').
$$
The inverse relation is
$$
\int_{x_\mathrm{min}}^1 dx~ f(x) = 1 - g(1-x_\mathrm{min}).
$$

# Formula summary for `QEDEvolve.jl`

This file summarizes the formulas implemented in `QEDEvolve.jl`.

## Running QED coupling

The code uses the one-loop running coupling

```math
\frac{d\alpha}{d\log \mu^2} = -\beta_0 \alpha^2,
\qquad
\beta_0(n_f) = -\frac{n_f}{3\pi}.
```

With no threshold matching, the solution is

```math
\alpha(\mu^2)
= \left[
    \frac{1}{\alpha(\mu_0^2)}
    + \beta_0(n_f)\log\frac{\mu^2}{\mu_0^2}
  \right]^{-1}.
```

## Mellin-space splitting kernels

Hats denote Mellin moments:

```math
\hat f(N) = \mathcal M[f(x)](N) = \int_0^1 dx\,x^{N-1}f(x).
```

$\gamma_E = -\psi(1)$, $\psi$ is the digamma function, and $\psi_1$ is the trigamma function.

### Unpolarized kernels

```math
\hat P_{ll}(N)
= \frac{3}{2}
  - \frac{1}{N}
  - \frac{1}{N+1}
  - 2\left[\psi(N)-\psi(1)\right],
```

```math
\hat P_{\gamma l}(N)
= \frac{N^2+N+2}{(N-1)N(N+1)},
\qquad
\hat P_{l\gamma}(N)
= \frac{N^2+N+2}{N(N+1)(N+2)},
```

```math
\hat P_{\gamma\gamma}(N) = -\frac{2}{3}.
```

The space-like singlet matrix is

```math
\hat P_s(N)
=
\begin{pmatrix}
\hat P_{ll}(N) & 2\hat P_{l\gamma}(N) \\
\hat P_{\gamma l}(N) & \hat P_{\gamma\gamma}(N)
\end{pmatrix},
```

and the time-like singlet matrix is

```math
\hat{\mathbb P}_s(N)
=
\begin{pmatrix}
\hat P_{ll}(N) & 2\hat P_{\gamma l}(N) \\
\hat P_{l\gamma}(N) & \hat P_{\gamma\gamma}(N)
\end{pmatrix}.
```

### Helicity kernels

```math
\Delta\hat P_{ll}(N)=\hat P_{ll}(N),
\qquad
\Delta\hat P_{\gamma l}(N)=\frac{N+2}{N(N+1)},
```

```math
\Delta\hat P_{l\gamma}(N)=\frac{N-1}{N(N+1)},
\qquad
\Delta\hat P_{\gamma\gamma}(N)=\hat P_{\gamma\gamma}(N).
```

The helicity singlet matrices are

```math
\Delta\hat P_s(N)
=
\begin{pmatrix}
\Delta\hat P_{ll}(N) & 2\Delta\hat P_{l\gamma}(N) \\
\Delta\hat P_{\gamma l}(N) & \Delta\hat P_{\gamma\gamma}(N)
\end{pmatrix},
```

```math
\Delta\hat{\mathbb P}_s(N)
=
\begin{pmatrix}
\Delta\hat P_{ll}(N) & 2\Delta\hat P_{\gamma l}(N) \\
\Delta\hat P_{l\gamma}(N) & \Delta\hat P_{\gamma\gamma}(N)
\end{pmatrix}.
```

## Fixed-order distributions

### Unpolarized distributions

The NLO lepton-to-lepton distribution is

```math
f_{ll}(\alpha,x,m^2,\mu^2)
= \delta(x-1) + \frac{\alpha}{2\pi}
  \left[\frac{1+x^2}{1-x}
  \left(
    \log\frac{\mu^2}{m^2}
    -2\log(1-x)-1
  \right)\right]_+.
```

Its Mellin moment at $\mu = m$ is

```math
\hat f_{ll,0}(\alpha,N)
= 1 + \frac{\alpha}{2\pi}\left[-\hat P_{ll}(N)+R_{ll}(N)\right],
```

with

```math
\begin{align}
R_{ll}(N)
&= -\frac{1}{6N^2(N+1)^2}
\Big\{
12
+ N(N+1)
  \left[
    36 + 12\gamma_E^2N(N+1)
    +(-21+2\pi^2)N(N+1)
    +12\gamma_E(1+2N)
  \right]
-12N(N+1) \\
  &\qquad \times \left[
    -\psi(N)
      \left(
        1+2N(1+\gamma_E+\gamma_E N)
        +N(N+1)\psi(N)
      \right)
    +N(N+1)\psi_1(N)
  \right]
\Big\}.
\end{align}
```

The NLO lepton-to-photon distribution is

```math
f_{\gamma l}(\alpha,x,m^2,\mu^2)
= \frac{\alpha}{2\pi}
  \frac{1+(1-x)^2}{x}
  \left[
    \log\frac{\mu^2}{m^2}
    -2\log x -1
  \right].
```

Its Mellin moment at $\mu = m$ is

```math
\hat f_{\gamma l,0}(\alpha,N)
= \frac{\alpha}{2\pi}
  \left[
    -\hat P_{\gamma l}(N)
    + \frac{4}{(N-1)^2}
    - \frac{4}{N^2}
    + \frac{2}{(N+1)^2}
  \right].
```

### Helicity distributions

The NLO lepton-to-lepton helicity distribution is the same as the unpolarized one:

```math
g_{ll}=f_{ll},
\qquad
\hat g_{ll,0}=\hat f_{ll,0}.
```

The NLO lepton-to-photon helicity distribution is

```math
g_{\gamma l}(\alpha,x,m^2,\mu^2)
= \frac{\alpha}{2\pi}(2-x)
  \left[
    \log\frac{\mu^2}{m^2}
    -2\log x
    -\frac{1-x}{2-x}
  \right].
```

Its Mellin moment at $\mu = m$ is

```math
\hat g_{\gamma l,0}(\alpha,N)
= \frac{\alpha}{2\pi}
  \frac{N^2+7N+4}{N^2(N+1)^2}.
```

### Fragmentation functions

The NLO lepton-to-lepton fragmentation function is the same unpolarized distribution:

```math
D_{ll}=f_{ll},
\qquad
\hat D_{ll,0}=\hat f_{ll,0}.
```

The NLO photon-to-lepton fragmentation function is

```math
D_{l\gamma}(\alpha,x,m^2,\mu^2)
= \frac{\alpha}{2\pi}
  \left[x^2+(1-x)^2\right]
  \log\frac{\mu^2}{m^2}.
```

Its initial Mellin moment at $\mu = m$ is

```math
\hat D_{l\gamma,0}(\alpha,N)=0.
```

## Mellin-space DGLAP evolution

The starting scale is chosen to be $\mu_0 = m$. Define

```math
\rho = \frac{\alpha(\mu_0^2)}{\alpha(\mu^2)}.
```

If $n_f = 0$, the code leaves the Mellin moments unchanged.

### Singlet and nonsinglet basis

The evolution is carried out in the singlet/nonsinglet basis. For distribution functions,

```math
\hat f_v
= \hat f_{ll}-\hat f_{\bar l l},
\qquad
\hat{\boldsymbol f}_s
=
\begin{pmatrix}
\hat f_{ll}+\hat f_{\bar l l} \\
\hat f_{\gamma l}
\end{pmatrix}.
```

The nonsinglet component $\hat f_v$ evolves with $\hat P_{ll}$, while the singlet vector $\hat{\boldsymbol f}_s$ evolves with $\hat P_s$. The physical components are recovered by

```math
\hat f_{ll}
= \frac{1}{2}\left(\hat f_{s,1}+\hat f_v\right),
\qquad
\hat f_{\bar l l}
= \frac{1}{2}\left(\hat f_{s,1}-\hat f_v\right),
\qquad
\hat f_{\gamma l}
= \hat f_{s,2}.
```

For fragmentation functions,

```math
\hat D_v
= \hat D_{ll}-\hat D_{l\bar l},
\qquad
\hat{\boldsymbol D}_s
=
\begin{pmatrix}
\hat D_{ll}+\hat D_{l\bar l} \\
\hat D_{l\gamma}
\end{pmatrix}.
```

It evolves with $\hat P_{ll}$ and the time-like singlet kernel $\hat{\mathbb P}_s$, with

```math
\hat D_{ll}
= \frac{1}{2}\left(\hat D_{s,1}+\hat D_v\right),
\qquad
\hat D_{l\bar l}
= \frac{1}{2}\left(\hat D_{s,1}-\hat D_v\right),
\qquad
\hat D_{l\gamma}
= \hat D_{s,2}.
```

### Nonsinglet solution

The evolved nonsinglet function is

```math
\hat f_v(N,\mu^2)
=
\rho^{\hat P_{ll}(N)/(2\pi\beta_0)}
\hat f_v(N,\mu_0^2).
```

### Singlet matrix solution

For the `2 x 2` Mellin-space singlet kernel matrix $\hat{P}_s(N)$, the eigenvalues are

```math
\lambda_\pm
= \frac{1}{2}
  \left[
    \hat P_{s,11}+\hat P_{s,22}
    \pm
    \sqrt{
      (\hat P_{s,11}-\hat P_{s,22})^2
      +4\hat P_{s,12}\hat P_{s,21}
    }
  \right].
```

The evolution matrix is ($I$ is the `2 x 2` identity matrix)

```math
U_s(N)
=
\frac{\hat P_s-\lambda_- I}{\lambda_+-\lambda_-}
\rho^{\lambda_+/(2\pi\beta_0)}
+
\frac{\hat P_s-\lambda_+ I}{\lambda_- - \lambda_+}
\rho^{\lambda_-/(2\pi\beta_0)}.
```

The evolved singlet vector is

```math
\hat{\boldsymbol f}_s(N,\mu^2)
= U_s(N)\hat{\boldsymbol f}_s(N,\mu_0^2).
```

## Mellin inversion

The Mellin inversion contour is parameterized by

```math
N(r)=c+r e^{i\phi},
\qquad c=1.9,\qquad \phi=\frac{3\pi}{4}.
```

The real one-dimensional integrand is

```math
I(r;x)
= \frac{1}{\pi}
  x^{-(c+r\cos\phi)}
  \operatorname{Im}
  \left[
    e^{i(\phi-r\sin\phi\log x)}
    \hat f(N(r))
  \right],
```

and the inverse Mellin transform is evaluated as

```math
f(x) = \mathcal M^{-1} [\hat f(N)](x) = \int_0^\infty I(r;x)\,dr.
```

The helper `integrate_distribution` computes

```math
\int_{x_{\min}}^1 f(x)\,dx
=
x_{\min}\,
\mathcal M^{-1}
\left[
  \frac{\hat f(N)}{N-1}
\right](x_{\min}).
```

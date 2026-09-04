# Physics model

Status: current
Date: 2026-09-03
Code version: see manifest `commit`

`stiffgwpy_fast` models a flat LCDM cosmology with:

* **radiation + massive neutrinos** (Fermi-Dirac distribution),
* **relativistic particles** including extra radiation `Delta N_eff`,
* **stiff matter** parameterised by `kappa10 = rho_stiff / rho_photon` at 10 MeV,
* a **primordial tensor background** from inflation with amplitude `A_t = A_s * r`,
  tilt `n_t`, and (optionally, `cr > 0`) the single-field consistency relation.

The time variable is `N = ln a`.  The key background quantity is
`sigma(N) = d ln H / d N` (the equation-of-state-weighted factor that sources the
tensor-mode equation) and the horizon position `f_hor(N) = log10(aH/(2 pi))`.

## Tensor-mode equation

Each frequency channel is evolved in the original variables
(`z = ln(k/aH)`, plus the two tensor polarisation combinations `x, y`):

```
z'  = 1.5 sigma - 1
x'  = -3 x + 1.5 sigma x - e^z y
y'  = -y + 1.5 sigma y + e^z x
```

The source is the primordial spectrum
`P_t(k) = A_t (k / k_piv)^{n_t}`.  Today's `Omega_GW(f)` is assembled from the
`Ogw`, `Oj`, `Opgw` combinations of `(x, y, z)` at `N = N_inf` (today), and the
integrated `Delta N_eff` is the bolometric frequency integral of `Ogw - Oj`.

## `Delta N_eff` closure

The SGWB contributes extra radiation, which changes the background, which
changes the SGWB.  This is solved by an outer bisection on `Delta N_eff` until
the successive relative change is below the outer tolerance (`1e-7` for
`production`, `1e-6` for plain-grid).  A physical guard rejects points where
the total `N_eff` exceeds `5` (too much radiation); that is a **physical**
rejection, not a numerical failure.

## Reheating

The instantaneous-reheating kink in `sigma(N)` (matter-like
`sigma = 1` for `N < N_re`, then radiation/neutrino/stiff evolution) is the
hardest feature.  `production` treats it as an exact ODE breakpoint;
the plain-grid profile smears it across the fixed grid, which is the dominant
source of its bias.

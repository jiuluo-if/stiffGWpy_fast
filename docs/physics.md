# Physics model

Status: current model description; parameter-space validation remains scoped

Date: 2026-09-24

Code version: see the run commit in each validation artifact. The validation
manifest is dated 2026-09-17 and is not a current-HEAD certification.

`stiffgwpy_fast` models a flat LCDM cosmology with radiation, massive neutrinos
(Fermi-Dirac distribution), extra relativistic species, stiff matter, and a
primordial tensor background. The stiff component is parameterized by
`kappa10 = rho_stiff / rho_photon` at 10 MeV. Tensor power is
`P_t(k) = A_t (k / k_piv)**n_t`, where `A_t = A_s * r`; for `cr > 0` the
single-field consistency relation determines the tensor tilt and reheating
duration.

The expansion variable is `N = ln(a)`. The background quantity
`sigma(N) = d ln(H) / dN` sources the tensor-mode equation, and
`f_hor(N) = log10(aH/(2 pi))` tracks horizon crossing.

## Tensor-mode equation

For each frequency channel, the solver evolves `z = ln(k/aH)` and the two
tensor combinations `x` and `y`:

```text
z' = 1.5 sigma - 1
x' = -3 x + 1.5 sigma x - exp(z) y
y' = -y + 1.5 sigma y + exp(z) x
```

The present-day spectrum is assembled from `Ogw`, `Oj`, and `Opgw` at
`N = N_inf`. The bolometric `Delta N_eff` contribution is integrated from
`Ogw - Oj` over frequency.

## Self-consistency and physical guards

The SGWB contributes extra radiation, which changes the background and hence
the SGWB. The single user-facing `fast` mode iterates this closure with the
`1e-6` preset tolerance. The original LSODA path has its own default
convergence settings; `engine='reference'` selects the independent
continuous-sigma DOP853 pipeline. A shared guard rejects configurations whose
total `N_eff` exceeds 5. This is an explicit physical rejection, not a
numerical failure.

## Reheating transition

Instantaneous reheating creates a kink in `sigma(N)`: the matter-like
reheating segment ends at `N_re`, after which the radiation/neutrino/stiff
background evolves. The formal `fast` preset splits tensor transfer exactly
at `N_re`, so a transfer step does not cross the discontinuity. Historical
plain-grid experiments smeared this feature and remain useful only as dated
validation evidence; they are not a current user-facing solver mode.

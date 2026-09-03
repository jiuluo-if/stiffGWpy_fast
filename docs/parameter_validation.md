# Full parameter validation

Status: current
Date: 2026-09-03
Code version: see manifest `commit`

## Parameter schema

Every parameter that enters the physics (from the code, not just the README):

| name | meaning | unit | prior range |
|---|---|---|---|
| `Omega_bh2` | baryon density | — | [0.01, 0.04] |
| `Omega_ch2` | CDM density | — | [0.05, 0.25] |
| `H0` | Hubble constant | km/s/Mpc | [50, 90] |
| `DN_eff` | constant extra radiation | — | [0, 3] |
| `A_s` | scalar amplitude | — | log10 in [-9.2, -8.2] |
| `r` | tensor-to-scalar ratio | — | log10 in [-4, -1] |
| `n_t` | tensor tilt | — | [-0.5, 0.5] |
| `cr` | consistency-relation switch | — | [0, 1] |
| `T_re` | reheating temperature | GeV | log10 in [-1, 7] |
| `DN_re` | matter-like reheating e-folds | — | [0, 40] |
| `kappa10` | stiff/photon ratio at 10 MeV | — | log10 in [-3, 3] |

Physical constraints: `r > 0` (no tensor source otherwise), `N_inf` finite
(cut-off set), and total `N_eff <= 5` (the shared guard).

## Coverage

* **Axis edges:** `docs/paramsweep_z8b/` — 16 points on the
  `r`/`n_t`/`cr`/`T_re`/`DN_re`/`kappa10` axis edges + transition interiors.
* **Sobol (production):** `docs/paramsweep_ref/fast_sweep.jsonl` — 240 points:
  **212 ok / 28 explicit shared-`Delta_Neff` guard**.
* **LHS (plain-grid screen):** `docs/validation/param_sweep_plain.json` —
  400 points: **255 success / 145 guard / 0 numerical failure**.

## Classification

Rejections are categorised (never blanket-counted as solver failures):

* `PHYSICAL_INVALID` — `r <= 0`, cut-off not set.
* `PHYSICAL_GUARD` — shared `Delta_Neff > 5` (deterministic, not retried through LSODA).
* `NUMERICAL_FAILURE` — exception / non-finite / iteration failure.
* `FAST_ERROR` / `ORACLE_ERROR` — per-engine solve failures.

The 36% guard fraction in the LHS wide box is a physical rejection of too-stiff /
too-loud corners, reported explicitly, never hidden.

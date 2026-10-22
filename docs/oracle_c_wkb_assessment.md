# Oracle C: analytic higher-order WKB tail correction

## Motivation

The independent reference continues the tensor mode with DOP853 up to a
handoff depth `z_tail` and then uses a frozen-amplitude analytic tail

    Th(N) = coeff * exp(-z_handoff + event_N - N).

The `z_tail = 5/6/7/8/10` sweep (`docs/oracle_tail_convergence_default_full.json`)
reported an observed systematic of `3.63e-3` with a non-monotone sequence.
That sweep alone cannot tell whether the non-monotonicity is physical or an
artefact of the frozen-amplitude approximation.  Oracle C supplies an
independent, fully analytic estimate of the tail term.

## Derivation

Writing the Prüfer amplitude/phase as `x = A sin(theta)`, `y = A cos(theta)`,
the amplitude equation is

    d ln A / dN = 1.5*sigma - 2 + cos(2*theta),
    d theta / dN = -exp(z) - sin(2*theta).

The frozen tail already carries the leading-order WKB amplitude evolution
`d ln A / dN = 1.5*sigma - 2` (equivalently, the physical amplitude
`h = A exp(-z + N)` is frozen at leading order).  The residual is the fast
oscillating `cos(2*theta)` term.  With `d theta / dN ~ -omega` and
`omega = exp(z)`, the stationary-phase boundary term integrates to

    integral cos(2*theta) dN  ~  sin(2*theta_f) / omega_f,

so the transfer-squared (which carries `h^2`) is corrected by

    transfer^2_wkb = transfer^2_frozen * (1 + sin(2*theta_f) / omega_f).

This correction needs only the handoff phase `theta_f` and `omega_f`; it is
analytic and free of any additional ODE solve.

## Protocol

`scripts/benchmark_oracle_c_wkb.py` compares three today estimates on the same
full native grid (76 frequencies) for `default`, `lowT`, `highT`, and `stiff`:

1. `frozen(z5)`: current pipeline handoff at `z=5`;
2. `wkb(z5)`: same handoff with the analytic correction above;
3. `deep(z10)`: handoff at `z=10` (`eps ~ 1e-4`), frozen thereafter.

The deep run is used as a numerical reference.  Resources: workers=1,
Numba=2, BLAS=1.

## Results

| point | DN frozen(z5) | DN wkb(z5) | DN deep(z10) | frozen-vs-deep | wkb-vs-deep | improvement |
|---|---:|---:|---:|---:|---:|---:|
| default | 2.2718752671e-03 | 2.2636716959e-03 | 2.2636593444e-03 | 3.63e-3 | 5.46e-6 | 665x |
| lowT | 5.3034329776e-08 | 5.3105808957e-08 | 5.3105479625e-08 | 1.34e-3 | 6.20e-6 | 216x |
| highT | 5.6637033806e-02 | 5.6433641889e-02 | 5.6434345815e-02 | 3.59e-3 | 1.25e-5 | 288x |
| stiff | 2.3651952547e-02 | 2.3610412290e-02 | 2.3609922928e-02 | 3.65e-3 | 8.09e-6 | 451x |

Per-mode `Ogw` residual against the deep reference drops from a median of
`2.7e-3..4.7e-3` (frozen) to `2.1e-5..2.9e-5` (WKB), with a maximum of
`8.8e-5`.  The handoff adiabaticity is `eps_median ~ 6.7e-3..1.1e-2`, and the
corrected residual is consistent with the expected second-order `eps^2` size.

The frozen `z_tail=5` and deep `z_tail=10` values reproduce the Stage B sweep
exactly (`0.0022718753` and `0.0022636593` for default), which validates the
new script against the existing oracle artifacts.

## Decision

**PASS / VERIFIED.** The frozen-amplitude tail defect is a first-order
adiabatic (`O(eps)`) artefact, not a physical non-monotonicity, and it is
removed by a closed-form correction to `~1e-5` relative on `DN_gw`.  This
narrows the reference tail systematic by more than two orders of magnitude
without any additional ODE solve.  Oracle C is accepted as an analytic tail
oracle and as the explanation for the non-monotone Stage B tail sequence.

## Residual uncertainty

The deep `z=10` reference still carries its own `O(eps(10)) ~ 9e-5` frozen
tail, so the reported `~1e-5` residual is an upper bound on the corrected
tail error.  Promoting the correction into the reference tail itself is a
separate, higher-risk change and is not proposed here.

## 2026-10-22 edge and Sobol follow-up

The full Oracle C comparison covers two low/high parameter edges, one fixed
Sobol point, and two `DN_re` edges. Each point uses its fast self-consistent
`DN_eff` and native frequency grid; 62–66 modes can be compared across the
`z=5` and `z=10` handoffs.

| Point | Frequencies / compared tail modes | WKB-corrected `DN_gw` vs deep `z=10` | Maximum per-mode WKB/deep difference | Fast vs WKB-corrected `DN_gw` | Fast vs deep `DN_gw` |
|---|---:|---:|---:|---:|---:|
| `edge_r_lo` | `76 / 65` | `5.20e-6` | `8.77e-5` | `1.18e-5` | `6.62e-6` |
| `edge_kap_hi` | `77 / 66` | `1.37e-5` | `8.56e-5` | `1.17e-5` | `2.54e-5` |
| `sobol_000` | `76 / 65` | `3.10e-5` | `8.57e-5` | `1.24e-5` | `1.87e-5` |
| `edge_dnre_lo` | `77 / 62` | `2.35e-6` | `8.07e-5` | `2.74e-3` | `2.74e-3` |
| `edge_dnre_hi` | `76 / 65` | `5.66e-6` | `8.23e-5` | `2.49e-3` | `2.50e-3` |

The WKB correction remains independently accurate against the deep handoff on
both `DN_re` edges, but the formal fast result misses the `<2e-4` gate there by
about 12–14x. The discrepancy is not explained by the WKB reference: its
integrated residual is below `6e-6`, with maximum per-mode residual below
`8.3e-5`. Changing the fast frequency quadrature does not bring either point
inside the gate.

A scoped configuration sweep associates the remaining discrepancy with fast
background resolution. At the standard `h=0.005`, the two fast residuals are
`2.74e-3` and `2.50e-3`; deeper tail handoff or smaller `phase_max` barely
changes them. A matched re-run at `h=0.0003125` uses Oracle anchors generated
from the same configuration and `DN_eff`; both the `DN_eff` and native grid
are checked exactly by the comparison script. The fast results are within the
`2e-4` gate, at `1.65e-4` and `1.63e-4` versus WKB and `1.68e-4` versus deep.
The Oracle WKB/deep residuals are `2.89e-6` and `4.56e-6`, with per-mode
maxima below `8.3e-5`.

The accuracy improvement has a substantial latency cost: the `h=0.0003125`
medians are `62.7` and `71.7 ms`, compared with `5.69` and `4.53 ms` for the
standard `h=0.005` fast path (about `11.0x` and `15.8x` slower). Both
measurements used two warmup calls and ten timed fresh-model solves with two
Numba `workqueue` threads and one BLAS thread. This candidate has not been
promoted; optimizing its speed remains open.

A second matched candidate uses the validation `transition-refine` path at
`h=0.0025`, with the vectorized `H2_vec` implementation. It also passes the
gate: fast-vs-deep errors are `3.97e-6` and `1.64e-4`; warm medians are `46.5`
and `73.7 ms` (about `8.2x` and `16.3x` the standard fast timings). It is
faster than the uniform fine-step candidate at `edge_dnre_lo`, but slightly
slower at `edge_dnre_hi`; neither candidate meets the normal fast latency
target, so the formal preset remains unchanged.

A smaller local zone was then exposed as immutable per-call configuration:
`transition_refine_radius=0.25` e-folds and `transition_refine_factor=8` at
the same `h=0.0025`. Regenerated Oracle anchors and exact `DN_eff`/native-grid
matching give fast-vs-deep errors of `5.17e-6` and `1.62e-4`, both within the
scoped gate. The Oracle WKB/deep residuals are `2.90e-6` and `5.72e-6`. With
two warmups and ten timed fresh-model calls, medians are `35.6` and `51.7 ms`;
this is `1.29x`/`1.42x` faster than the radius-2.0 transition-refine candidate,
but still about `6.3x`/`11.4x` slower than standard fast. This remains a
two-edge validation candidate and has not been promoted.

The matched Oracle residuals for this candidate are `1.08e-6` and `1.59e-4`
versus WKB-corrected `DN_gw`; the WKB/deep reference residuals are
`2.89e-6` and `5.65e-6`. Both `DN_eff` and frequency grids are checked exactly
before reporting these values.

Raw Oracle C runs and five-point comparisons are in
[`oracle_c_wkb_corner_20261022/`](oracle_c_wkb_corner_20261022/). The matched
fine-step anchors and comparison are in
[`oracle_c_wkb_corner_20261022/step_h0003125/`](oracle_c_wkb_corner_20261022/step_h0003125/).
The matched transition-refine anchors and comparison are in
[`oracle_c_wkb_corner_20261022/transition_h00025/`](oracle_c_wkb_corner_20261022/transition_h00025/).
The artifact metadata records source commit and resources for each run; the
standard-step `DN_re` Oracle anchors were produced from commit
`5e6e78875cf9998fb64b72c6b216c32c7848abe8`, and the fine-step anchors from
`0729306f4d4d828f7c487a1d613b7588c4ee50f2`.
The matched transition-refine anchors were regenerated together from commit
`e12fb4aae5e574ecb5bd3f61b808a547cbc20700`.
The radius-0.25 anchors and timing record are in
[`oracle_c_wkb_corner_20261022/transition_h00025_r025/`](oracle_c_wkb_corner_20261022/transition_h00025_r025/).

Reproduce the standard, fine-step, and transition-refined anchors and comparisons:

```powershell
python scripts/benchmark_oracle_c_wkb.py --points edge_r_lo edge_kap_hi sobol_000 edge_dnre_lo edge_dnre_hi --out-dir docs/oracle_c_wkb_corner_20261022
python scripts/benchmark_fast_true_error.py --points edge_r_lo edge_kap_hi sobol_000 edge_dnre_lo edge_dnre_hi --anchor-dir docs/oracle_c_wkb_corner_20261022 --warmup 2 --repeats 10 --out docs/oracle_c_wkb_corner_20261022/fast_true_error.json
python scripts/benchmark_oracle_c_wkb.py --points edge_dnre_lo edge_dnre_hi --h 0.0003125 --out-dir docs/oracle_c_wkb_corner_20261022/step_h0003125
python scripts/benchmark_fast_true_error.py --points edge_dnre_lo edge_dnre_hi --anchor-dir docs/oracle_c_wkb_corner_20261022/step_h0003125 --h 0.0003125 --warmup 2 --repeats 10 --out docs/oracle_c_wkb_corner_20261022/step_h0003125/fast_true_error.json
python scripts/benchmark_oracle_c_wkb.py --points edge_dnre_lo edge_dnre_hi --h 0.0025 --transition-refine --out-dir docs/oracle_c_wkb_corner_20261022/transition_h00025
python scripts/benchmark_fast_true_error.py --points edge_dnre_lo edge_dnre_hi --anchor-dir docs/oracle_c_wkb_corner_20261022/transition_h00025 --h 0.0025 --transition-refine --warmup 2 --repeats 10 --out docs/oracle_c_wkb_corner_20261022/transition_h00025/fast_true_error.json
python scripts/benchmark_oracle_c_wkb.py --points edge_dnre_lo edge_dnre_hi --h 0.0025 --transition-refine --transition-refine-radius 0.25 --transition-refine-factor 8 --out-dir docs/oracle_c_wkb_corner_20261022/transition_h00025_r025
python scripts/benchmark_fast_true_error.py --points edge_dnre_lo edge_dnre_hi --anchor-dir docs/oracle_c_wkb_corner_20261022/transition_h00025_r025 --h 0.0025 --transition-refine --transition-refine-radius 0.25 --transition-refine-factor 8 --warmup 2 --repeats 10 --out docs/oracle_c_wkb_corner_20261022/transition_h00025_r025/fast_true_error.json
```

## Artifacts

- `docs/oracle_c_wkb_default.json`
- `docs/oracle_c_wkb_lowT.json`
- `docs/oracle_c_wkb_highT.json`
- `docs/oracle_c_wkb_stiff.json`
- `scripts/benchmark_oracle_c_wkb.py`

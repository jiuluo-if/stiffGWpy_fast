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

## Artifacts

- `docs/oracle_c_wkb_default.json`
- `docs/oracle_c_wkb_lowT.json`
- `docs/oracle_c_wkb_highT.json`
- `docs/oracle_c_wkb_stiff.json`
- `scripts/benchmark_oracle_c_wkb.py`

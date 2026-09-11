# True fast DN_gw error against the Oracle C WKB anchor

## Motivation

`fast` and the frozen-amplitude reference share the same `z_tail` frozen-tail
convention, so a fast-vs-reference comparison cancels the shared tail defect.
After Oracle C established an analytic WKB correction of that defect
(`docs/oracle_c_wkb_assessment.md`), the *true* fast error can be measured as
the difference between the self-consistent fast `DN_gw` and the
WKB-corrected reference `DN_gw` at the same `DN_eff`.

## Protocol

`scripts/benchmark_fast_true_error.py` re-solves fast at the current HEAD and
reads the committed reference anchors from `docs/oracle_c_wkb_*.json` (same
native grid, same `DN_eff`).  Resources: workers=1, Numba=2, BLAS=1.

## Results

| point | fast DN_gw | vs frozen(z5) | vs WKB(z5) | vs deep(z10) |
|---|---:|---:|---:|---:|
| default | 2.2626969518e-03 | 4.04e-3 | **4.31e-4** | 4.25e-4 |
| lowT | 5.2446047710e-08 | 1.11e-2 | **1.24e-2** | 1.24e-2 |
| highT | 5.6391339588e-02 | 4.34e-3 | **7.50e-4** | 7.62e-4 |
| stiff | 1.4965149678e-02 | 4.91e-3 | **1.29e-3** | 1.28e-3 |

## Findings

- The reported fast-vs-reference `DN_gw` (PCHIP `2.94e-4` on the fixed
  `z_tail=8` same-grid check) is a *shared-tail-convention* number.  Against
  the WKB-corrected anchor the true fast error is `4.3e-4..1.3e-3` for
  `default`/`highT`/`stiff`, i.e. 1.5x-4x larger, and `1.24e-2` for `lowT`.
- For `default`/`highT`/`stiff`, fast sits close to the WKB value (and to the
  deep `z=10` value), while the frozen `z=5` reference sits `3.6e-3` away.
  This is consistent with `_tail_match_gamma` in `fast_sgwb.py`, whose
  `T ~ a^-1` amplitude match is already closer to the true tail than the
  plain frozen continuation.
- `lowT` is the outlier: fast differs from both the frozen and the WKB
  references by `~1.1e-2..1.2e-2`, far larger than the `1.3e-3` frozen-to-WKB
  tail correction.  Its `DN_gw` is `5.2e-8`, so this is a genuinely
  different, non-tail error source (quadrature / grid / low-amplitude
  cancellation) rather than a tail artefact.

## Decision

**ACCEPTED as a calibration finding.** The true fast `DN_gw` error is larger
than the previously reported same-convention number and is not yet inside the
`2e-4` release gate.  The dominant remaining term for `default`/`highT`/
`stiff` is the residual fast-vs-WKB difference (`4e-4..1.3e-3`); for `lowT`
the error is dominated by a separate low-amplitude effect.

## Next experiment

Decompose the fast-vs-WKB residual per frequency (tail vs deep-subhorizon
stepping vs frequency quadrature) and, if the residual is tail-dominated,
apply the Oracle C `1 + sin(2*theta_f)/omega_f` correction directly inside the
fast tail assembly under acceptance criteria (`spectrum max < 1e-3`,
`DN rel < 5e-4`, no new failure, determinism pass).

## Artifacts

- `docs/fast_true_error.json`
- `scripts/benchmark_fast_true_error.py`

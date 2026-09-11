# Fast DN_gw error decomposition against the Oracle C anchor

## Motivation

`fast` and the frozen-amplitude reference share the same `z_tail` frozen-tail
convention, so a fast-vs-reference number cancels the shared tail defect.  The
previous calibration (`docs/fast_true_error_assessment.md`) measured the true
fast error against the Oracle C WKB anchor, but that number still mixed three
theoretical sources: frequency quadrature, tail amplitude, and the
deep-subhorizon propagation.  This script separates them.

## Protocol

`scripts/benchmark_fast_residual_decomposition.py` re-solves fast at the
current HEAD and reuses the committed `docs/oracle_c_wkb_*.json` anchors on the
same native grid.  For every point it reports

1. `fast_vs_wkb_*_rel`: the DN difference at the *same* quadrature rule;
2. `fast_simpson_vs_pchip_rel`: the quadrature-only difference of the fast
   per-node values;
3. the Simpson-weighted per-node contribution split into analytic-tail and
   non-tail nodes, plus the node that dominates it.

The script also verifies its own quadrature: `dn.fast.simpson` reproduces the
fast-reported `DN_gw` bitwise (e.g. default `2.2626969518e-03`).

## Results

| point | fast DN (Simpson) | fast DN (PCHIP) | WKB anchor (PCHIP) | fast-vs-WKB (Simpson) | fast-vs-WKB (PCHIP) | Simpson-vs-PCHIP | tail \|abs\| share | dominant node |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| default | 2.2626969518e-03 | 2.2636474151e-03 | 2.2636716959e-03 | 9.49e-06 | **1.07e-05** | 4.20e-04 | 1.000 | f=-0.2232 (17.0%) |
| lowT | 5.2446047710e-08 | 5.3114602834e-08 | 5.3105808957e-08 | 1.41e-04 | **1.66e-04** | 1.26e-02 | 0.110 | f=-18.4594 (52.7%) |
| highT | 5.6391339588e-02 | 5.6433084271e-02 | 5.6433641889e-02 | 9.32e-06 | **9.88e-06** | 7.40e-04 | 1.000 | f=1.1751 (18.3%) |
| stiff | 1.4965149678e-02 | 1.4984335009e-02 | 1.4984507082e-02 | 1.10e-05 | **1.15e-05** | 1.28e-03 | 1.000 | f=0.0268 (20.5%) |

## Findings

- **The fast propagation kernel is ~1e-5 accurate.** At a fixed quadrature
  rule the fast-vs-WKB per-node DN difference is `9.3e-06..1.4e-04`, i.e. one
  to two orders of magnitude smaller than the `4.2e-04..1.3e-02` that the
  Simpson-vs-PCHIP rule alone produces.  The true fast error is therefore
  quadrature-dominated, not tail- or propagation-dominated.
- For `default`/`highT`/`stiff` every bit of the residual lives on analytic
  tail nodes (`tail |abs| share = 1.000`), but only at the `1e-5` level; this
  is the residual `O(eps^2)` difference between fast's `_tail_match_gamma`
  tail and the WKB anchor, not the `O(eps)` tail defect Oracle C removed.
- `lowT` is different again: only `11.0%` of the residual magnitude is on tail
  nodes, and a single node, the lowest frequency `f=-18.4594`, carries `52.7%`
  of it.  That node sits at the low-amplitude end of the spectrum
  (`DN_gw = 5.2e-08`), so it is a cancellation/quadrature artefact rather than
  a physical tail effect.
- Consequence (EMPIRICALLY VALIDATED for these four named points): switching
  the default frequency quadrature from Simpson to PCHIP would leave a true
  DN error of `1.07e-05 / 9.88e-06 / 1.15e-05 / 1.66e-04`
  (default/highT/stiff/lowT), all inside the `2e-04` release gate and well
  inside the `5e-04` target.
- This also explains the long-standing "non-monotone node-count convergence"
  and "PCHIP is opt-in" observations: node refinement and PCHIP both reduce
  the same quadrature term, and the Simpson default masks the underlying
  `1e-5` kernel accuracy.

## Decision

**ACCEPTED as a diagnosis finding.**  No formal kernel is changed.  The fast
propagation/tail implementation is confirmed accurate to `~1e-5`; the
dominant remaining accuracy term is the default Simpson frequency quadrature
(`4.2e-04..1.3e-02`), with a secondary cancellation node at the low-frequency
end of `lowT`.

## Next experiment

Hypothesis: making PCHIP the default frequency quadrature (Numba-backed if
the scipy call is too slow) cuts the true DN error by one to two orders of
magnitude at a runtime cost below the accepted budget.

Acceptance criteria (fixed before running):

- true DN rel vs the Oracle C WKB anchor `< 2e-04` on default/highT/stiff;
- `spectrum max` rel unchanged or better than the current default;
- warm runtime median increase `< 10%`;
- no new failure and an unchanged guard/fallback set;
- determinism replay identical.

## Artifacts

- `docs/fast_residual_decomposition.json`
- `scripts/benchmark_fast_residual_decomposition.py`

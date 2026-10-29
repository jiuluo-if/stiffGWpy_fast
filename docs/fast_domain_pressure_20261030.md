# Fast parameter-domain pressure screen (2026-10-30)

This is a dated screening record for `fast_v0.2`, not a certification of the
full parameter space or a claim that the largest reliable domain is known.
The reproducible scanner is [`scripts/pressure_fast_domain.py`](../scripts/pressure_fast_domain.py).

## Scope and method

- Boundary suite: 514 unique points covering zero/near-zero and signed/log
  sweeps over the cosmological inputs, both `cr` branches, transition-adjacent
  and coupled cases. The scanner classifies derivation, preflight, background,
  propagation, and output outcomes separately. `N_inf > 120` is recorded as a
  resource preflight skip, not as an input clamp or a solver failure.
- Stiff corpus replay: 20 retained cases from the 2026-10-29 random solver
  study, replayed once at `c2b1b54`.
- Random pilots: 16-point Sobol and 16-point Latin-hypercube pilots. Their
  signed/log mapping and resource preflight produced no successful points;
  these are harness smoke pilots, not estimates of domain-wide success rate.
- Failure minimization: coordinatewise probes on selected nonfinite cases.
  This does not establish a global nearest boundary.

The relevant source snapshots are boundary suite `11b2eaf` and stiff replay
`c2b1b54`. Timing records use fixed two-core affinity on Windows and one
process; values are local measurements, not portable performance guarantees.
The full cold-start, thread-scaling, and batch-throughput matrix remains
incomplete; the pressure scanner records those fields when each suite is run.

## Results

| Screen | Outcome |
|---|---|
| 514 boundary points | 328 finite positive spectra; 136 classified solver failures (94 nonfinite numerical state, 39 nonpositive background H², 2 invalid outer-iterate cutoffs, 1 explicit guard); 24 derivation/solver exceptions; 14 invalid-derived-cutoff preflight rejects; 10 resource preflight skips; 2 resource timeouts |
| 20 retained stiff cases | 16 finite positive spectra and 4 nonfinite failures (IDs 2, 5, 11, 18); no uncaught exceptions in this replay |
| Named-point paired Fast profile | 6 cases × 25 alternating pairs against `9f1a70a`; outputs matched exactly. Worst warm P50 increase was 3.66%, worst P95 increase 2.28%; all stayed within the 5% budget |

The boundary “finite positive” count is only a numerical output classification;
it is not an accuracy pass. The two timed-out cases were boundary IDs 129
(`A_s=100`, `cr=1`) and 316 (coupled extreme `n_t`, `kappa10`, and `T_re`).
They remain unresolved resource/time cases. A one-step fixed-point map for ID
129 had positive residual at sampled `DN_eff` values through `1e50`; this is
not a proof that no solution exists. Several other cases grow to nonfinite
state in fewer than ten outer iterations, so increasing the 60-iteration cap
alone would not repair those failures.

## Accuracy and guarded recovery

The retained difficult sample 9 has `DN_gw=2.926979295585161e-4` on the raw
Fast result. Against a fixed-background, 353-node converged Reference check,
raw Fast differs by `2.4645e-3` in `DN_gw` and about `2.925e-3 dex` at the worst
native spectrum node. It therefore fails the project accuracy limits
(`2e-4` and `1e-3`). Tightening tested Fast step/phase settings and refining
only the Fast frequency grid did not close that gap.

With the opt-in `auto_escalate=True, escalate_to_reference=True` path, the
same point returns `FAST_ESCALATED`. A separate fixed-background pass through
the Reference solver at the returned `DN_eff` reproduces `DN_gw` and the
spectrum exactly. The raw Fast-to-Reference accuracy comparison remains the
353-node result above. The current tighter Reference escalation took about
26.3 seconds; this is a guarded recovery for this point, not normal Fast-path
performance.

A second minimized counterexample (`r=0.19649`, `n_t=1.06225`,
`kappa10=0.01160`, `T_re=1761.31 GeV`, `cr=0`) exposed a more severe
frequency-grid failure. At `DN_eff=0`, the formal `goal` grid produced a Fast
one-step map of about `2.91e13`; changing `freq_res` from 1 to 8 did not change
it. The same Fast map with `construct`, `grid_independent`, or `adaptive`
sampling was about `2.24`, consistent with the fixed-background Reference map
`2.23857`. The high-level Reference solve found a self-consistent result near
`DN_gw=2.2966`, and a fixed-background Reference evaluation at that value
reproduced the map to `5e-11` absolute residual.

The alternate Fast grids avoided the catastrophic map but did not meet the
accuracy gates: `construct` at `freq_res=1` differed from its fixed-background
Reference by `3.47e-4` in `DN_gw` and `3.06e-3 dex` in the spectrum. Increasing
`freq_res` to 4 or reducing propagation step `h` from `0.005` to `0.000625`
did not improve these errors. Therefore the alternate grid is not promoted as
a Fast accuracy fix.

The explicitly requested `auto_escalate=True, escalate_to_reference=True`
path now also retries a returned Fast numerical failure in Reference, using
`rtol <= 1e-9` while preserving caller tail and grid density. It
continues to exclude input/guard rejections and raised Fast
exceptions. A regression test confirms this point returns a finite Reference
result with `FAST_ESCALATED`; the run took about 61 seconds. This expands the
opt-in reliable evaluation path for this point, not the Fast-only domain.

## Accepted fixes in this work

- Grid changes now invalidate outer propagation reuse, preventing a false
  successful all-zero spectrum. A/B results match the full-recompute path and
  remain within the warm-path timing budget.
- Nonpositive background H² and invalid outer-iterate cutoffs now fail with
  explicit reasons instead of leaking arithmetic/type exceptions.
- Repeated requested evaluation frequencies are deduplicated before matrix
  construction, avoiding a zero frequency interval.
- An explicit Reference escalation request now retries recoverable Fast
  numerical failures, with a strict reference tolerance and without changing
  ordinary Fast or physical-rejection behavior.
- The pressure scanner records staged classifications, resource telemetry,
  resumable chunks, and optional coordinatewise minimization.

These fixes improve correctness and observability. They do not yet certify a
larger Fast-only accuracy domain. Remaining failures, exception cases, and
timeouts need separate root-cause work; no timeout or nonfinite result is
classified as a mathematical no-solution case.

## Reproduction

```powershell
python scripts/pressure_fast_domain.py --suite boundaries --threads 2 --max-n-inf 120 --output .planning/pressure_boundaries.jsonl
python scripts/pressure_fast_domain.py --suite sobol --power 4 --seed 20261030 --threads 2 --max-n-inf 120 --output .planning/pressure_sobol.jsonl
python scripts/pressure_fast_domain.py --suite latin --power 4 --seed 20261030 --threads 2 --max-n-inf 120 --output .planning/pressure_latin.jsonl
pytest -q tests/test_physics_limits.py::test_fast_auto_escalates_uncertified_hard_point_to_reference -o addopts=''
pytest -q tests/test_physics_limits.py::test_fast_auto_escalates_nonfinite_solver_failure_to_reference -o addopts=''
```

The date-stamped per-point records and paired timing JSON used for this
summary were local ignored research artifacts. Re-running the commands creates
fresh records; the counts above must be interpreted with the stated source
commits and screen definitions.

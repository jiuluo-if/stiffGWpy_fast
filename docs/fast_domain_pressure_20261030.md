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
`kappa10=0.01160`, `T_re=1761.31 GeV`, `cr=0`) exposed a frequency-grid
boundary error. Before the fix, the blue-tilt amplitude cutoff was
`fmax=-7.277`, while the reheating frequency was `f_re=0.0735`. The goal-grid
builder still added an `f_re±1` reserve, extending to `-0.926`, outside the
valid band. At `DN_eff=0`, Fast and Reference on those same 80 nodes both
integrated to about `2.9e13`; `freq_res=1..8` could not change the extra
feature nodes. This was grid construction, not a tensor-propagation mismatch.

The builder now adds the transition reserve only where it overlaps the valid
frequency interval. On this case its 64 seed nodes stay within
`[-18.459,-7.277]`. Fast returns a finite raw result, but it is not accurate
enough: against fixed-background Reference on the same nodes, `DN_gw` differs
by `3.92e-3` and the spectrum by `4.09e-3 dex`. Alternate `construct` and
`grid_independent` grids also miss the gates (`3.47e-4` DN and `3.06e-3 dex`);
increasing `freq_res` or reducing `h` to `0.000625` did not fix that residual.
No Fast-only accuracy candidate is promoted.

With `auto_escalate=True, escalate_to_reference=True` and `freq_res=1`, the
local budget marks this point `uncertified-path-configuration` with estimated
`DN_gw_error=0.0539` (the default escalation threshold is `5e-3`). It returns
`FAST_ESCALATED`, but its Reference root near `2.296589` is not frequency
converged: a fixed-background map at `freq_res=2` differs by `1.92e-3` from
`freq_res=1`, and at `freq_res=4` differs by `2.18e-3`. This default Reference
fallback is a recovery path, not an accuracy pass.

At `freq_res=8`, the full self-consistent Reference solve returns
`DN_gw=2.2914118670` on 800 frequencies in 481.6 seconds. At that root, a
fixed-background `freq_res=16` map gives `2.2914094306`, a relative change of
`1.06e-6`; the full `freq_res=8` solve has a self-consistency residual of
`5.0e-10`. This supports convergence of the integrated `DN_gw`. The spectrum
is still not certified to `1e-3 dex`:
8/16 uniform grids are nonnested, so their interpolated max difference
(`0.013 dex`) is not an exact-node comparison. Separate curvature-adaptive
Reference grids at `1e-3` and `2.5e-4 dex` changed the integrated result by
`1.04e-3` and the common-node spectrum by `1.24e-3 dex`. The curvature-only
estimator was false-safe on this case: at 6000 points its local estimate was
`1.84e-3 dex`, while a dense probe still found `2.74e-3 dex` change versus the
3213-point grid.

The adaptive helper now checks actual midpoint residuals as well as curvature
and enforces `max_points` by prioritizing intervals with the largest measured
or estimated error. A Gaussian peak centered between seed nodes reproduced the
old false-safe result (dense interpolation error `1.0`) and now passes the
direct-midpoint regression. On sample 2, this method remains expensive and
uncertified: the `1e-3 dex` target used 1822 retained nodes but still had
`9.09e-3 dex` sampled error after 14 rounds; the `2.5e-4 dex` target used 3467
nodes and 655 seconds, while its DN proxy changed by `3.03e-4` from the coarse
result and its spectrum changed by `9.43e-3 dex` on a dense probe.

The Reference `freq_res=16` attempt initially hit an artificial
`quad(limit=800)` ceiling with 1600 breakpoints. The limit now scales with the
breakpoint count and retains 50 adaptive subdivisions; an exact 1600-node
regression covers the failure. The high-resolution rerun then completed.

## Accepted fixes in this work

- Grid changes now invalidate outer propagation reuse, preventing a false
  successful all-zero spectrum. A/B results match the full-recompute path and
  remain within the warm-path timing budget.
- Nonpositive background H² and invalid outer-iterate cutoffs now fail with
  explicit reasons instead of leaking arithmetic/type exceptions.
- Repeated requested evaluation frequencies are deduplicated before matrix
  construction, avoiding a zero frequency interval.
- The goal-grid transition reserve stays inside the valid frequency band when
  blue-tilt amplitude cutoff places `fmax` below the reheating feature.
- An explicit Reference escalation request now retries recoverable Fast
  numerical failures, with a strict reference tolerance and without changing
  ordinary Fast or physical-rejection behavior.
- The pressure scanner records staged classifications, resource telemetry,
  resumable chunks, and optional coordinatewise minimization.
- Reference frequency quadrature now sizes QUADPACK's interval limit to dense
  breakpoint grids, avoiding the old 800-subinterval failure.
- Adaptive frequency refinement now honors `max_points` as a hard cap and
  prioritizes intervals with the largest sampled or estimated interpolation
  error; this detects missed peaks but does not certify sample 2's spectrum.

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
pytest -q tests/test_physics_limits.py::test_fast_auto_escalates_uncertified_blue_tilt_cutoff_case_to_reference -o addopts=''
pytest -q tests/test_freq_adaptive.py::test_goal_oriented_freqs_never_crosses_blue_tilt_amplitude_cutoff -o addopts=''
pytest -q tests/test_freq_adaptive.py::test_adaptive_detects_narrow_feature_between_seed_nodes -o addopts=''
pytest -q tests/test_reference.py::test_integrate_spectrum_scales_quad_budget_with_breakpoint_count -o addopts=''
```

For the high-resolution DN check, run `LCDM_SG(...).SGWB_iter(engine="reference", z_tail=5.0, rtol=1e-9, freq_res=8.0)` with the sample-2 parameters above. It takes about eight minutes on the recorded Windows machine.

The date-stamped per-point records and paired timing JSON used for this
summary were local ignored research artifacts. Re-running the commands creates
fresh records; the counts above must be interpreted with the stated source
commits and screen definitions.

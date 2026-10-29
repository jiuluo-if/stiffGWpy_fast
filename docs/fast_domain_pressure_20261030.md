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

With `auto_escalate=True, escalate_to_reference=True`, the local budget marks
this point uncertified and the wrapper runs Reference at `rtol<=1e-9`. It
returns `FAST_ESCALATED` near `DN_gw=2.296589`; a separate fixed-background
Reference pass at that value has `5e-11` absolute map residual. This is a
same-solver consistency check, not a separate oracle. The run took about 60
seconds. This expands the opt-in reliable evaluation path for the point, not
the Fast-only accuracy domain.

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
```

The date-stamped per-point records and paired timing JSON used for this
summary were local ignored research artifacts. Re-running the commands creates
fresh records; the counts above must be interpreted with the stated source
commits and screen definitions.

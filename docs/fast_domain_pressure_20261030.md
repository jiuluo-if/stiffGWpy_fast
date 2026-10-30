# Fast parameter-domain pressure screen (2026-10-30)

This is a dated screening record for `fast_v0.2`, not a certification of the
full parameter space or a claim that the largest reliable domain is known.
The reproducible scanner is [`scripts/pressure_fast_domain.py`](../scripts/pressure_fast_domain.py).

## Scope and method

- Historical boundary-suite record: 514 points covering zero/near-zero and signed/log
  sweeps over the cosmological inputs, both `cr` branches, transition-adjacent
  and coupled cases. The scanner classifies derivation, preflight, background,
  propagation, and output outcomes separately. `N_inf > 120` is recorded as a
  resource preflight skip, not as an input clamp or a solver failure.
- The current scanner removes `n_t` as an effective dimension wherever
  `cr>0`, since `input_nt()` derives it from `r`. The reproducible boundary
  suite is now 402 cases and the coupled sub-suite is 140 cases; the 514/240
  counts and results below remain historical measurements and are not
  silently recomputed.
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

## Current-head coupled recheck and accuracy triage

The deterministic `coupled` suite is now directly runnable and emits each
Fast result's local accuracy budget without conflating finite output with an
independent accuracy pass:

```powershell
python scripts/pressure_fast_domain.py --suite coupled --threads 2 --max-n-inf 120 --output .planning/pressure_coupled.jsonl
```

At source commit `b607750d43df67cc6b132c8ab3806d5026eee4f2`, all 240
`cr × n_t × kappa10 × T_re` points return finite positive spectra. Runtime P50
is `6.33 ms`, P95 `11.65 ms`, batch throughput `115.8 cases/s`, and peak
observed RSS `180.3 MiB` for this single-process, two-thread screen. The CLI
summary also aggregates local accuracy statuses and `DN_gw_error` estimates.
This revives the old nonfinite cluster as
computable outputs, but every point's Fast budget is
`uncertified-path-configuration`: estimated `DN_gw_error` ranges from
`4.05e-3` (median `2.56e-2`) to `1.02`, and none is below `2e-4`. These budget
values are not measured errors; this screen does not establish a larger
reliable Fast region.

For `cr=0, n_t=2, T_re=1e4 GeV, kappa10=0.01, DN_re=10`, the formal Fast result
at `freq_res=1` is `DN_gw=0.133901825` on 64 nodes. Independent Reference
self-consistency gives `0.133262472` at `freq_res=1` (65 nodes),
`0.133042970` at `freq_res=4` (260 nodes), and `0.133041909` at `freq_res=8`
(521 nodes). The res4/res8 roots agree to `8.0e-6` relative, while res1 is
`1.65e-3` away; the low-density Reference answer is not converged.

Frequency densification alone does not repair Fast: `grid_independent` at
res4/res8 gives roots `0.133707459`/`0.133706435`, which stabilize against each
other but remain about `5.0e-3` from the converged Reference root. At Fast's
res8 root and the identical 521 frequencies, a fixed-background Reference map
is `0.133055770` (`-4.87e-3` relative); the spectrum differs by max
`3.20e-3 dex`, P95 `3.18e-3`, median `1.90e-3`. Varying Fast `z_tail` 5/8/10,
`h` 0.005/0.0025, phase cap 0.25/0.1, `sigma_exact`, and transition
refinement did not meet either accuracy gate. No such Fast candidate is
promoted.

LSODA also needs a careful comparison here. At `freq_res=1`, tightening
`rtol` from `1e-9` to `1e-11` changes its root only from `0.133247367` to
`0.133247361`; it agrees with Reference res1 within `1.14e-4`, but both are
under-resolved relative to Reference res4/res8. At `freq_res=4`, LSODA gives
`0.132952454`, `6.80e-4` below Reference res4. A fixed-background Reference
map on the exact same 467 LSODA frequencies has a `1.43e-3` relative root
residual and spectrum differences max `5.95e-3 dex` / P95 `4.64e-3`. Refining
the LSODA background around reheating by 16× (local spacing about
`6.25e-4`, 11,340 background nodes) leaves its root unchanged to `1.3e-9`
absolute, so this discrepancy is not explained by the reheating grid alone.
At the stated `z_tail=5`, the Reference frequency sequence is converged, but
this point is not certified across independent propagation algorithms; both
the Fast and LSODA paths need further root-cause work before claiming the point
reliable. The follow-up below shows why the tail setting itself also needs
convergence checking for this sample.

### Propagation-tail recheck and rejected Fast candidate

A follow-up checked whether the earlier `z_tail=5` Reference comparisons were
deep enough for this blue-tilt point. On a shared mode at the Fast
`grid_independent/res8` root, extending the Reference propagation from
`z_tail=5` to `z_tail=10` moves the Reference log spectrum by about
`2.88e-3 dex`; changing only the Reference start from `z_start=-12` to Fast's
actual `z_start=-6.911807` changes it by only `2.71e-4 dex`. At the low
frequency end, extending the tail from 5 to 8 changes the Fast mode by about
`2.60e-3 dex`; halving its step from `h=0.005` to `0.0025` then changes the
mode by `1.60e-4 dex`, and further halving to `0.000625` changes it by less
than `1.4e-6 dex`. This isolates the earlier mismatch mainly to tail handoff
and low-frequency step error, not the initial-time cutoff. It also means the
previous `z_tail=5` root comparison is not a deep-tail oracle for this sample.

The best screened Fast candidate was
`grid_independent`, `freq_res=4`, `h=0.0025`, `z_tail=8`. Its root is
`DN_gw=0.133708325663`. On the same 260 frequencies and fixed background, a
Reference propagation to `z_tail=9` gives a map of `0.133734847235`, a
`1.984e-4` relative residual, and spectrum max/P95/median differences of
`4.554e-4 / 4.320e-4 / 1.946e-4 dex`. The same-background spectrum gate
passes, and the map residual is just inside the DN gate. That Reference map
took `2207.5 s`; a second map was not run, so a fully converged Reference
fixed point at res4 was not independently established. At res1, repeated
tail-9 Reference maps converge near `0.133952642`; the candidate's
`0.133925317` root is about `2.04e-4` relative away, narrowly outside the
DN gate. The evidence therefore does not certify the candidate's
self-consistent root at either density.

A 25-pair alternating warm timing on the same res4 case measured the
`h=0.005,z_tail=5` baseline at P50/P95 `41.0/43.4 ms` and the candidate at
`255.0/263.0 ms` (about `6.2×` slower); output digests differ as expected
from the changed numerical settings. No ordinary Fast defaults or solver code
were changed. Follow-up candidate tests did not close the accuracy gap; the
candidate is rejected for promotion.

### Adaptive-handoff and selective-step experiments (not retained)

The first-order `handoff_eps=|1.5 sigma-1| exp(-z)` is not a sufficient
standalone handoff rule for this sample. A first calculation reported
`0.999998` correlation between the log proxy and tail-depth error, but it
incorrectly included the `-1` sentinel for modes that never use the analytic
tail. Recomputing only the 47 actual handoff modes gives a log-log correlation
of `0.434`. In particular, at `f=-16.292921`, the mode's maximum `z` is only
`4.98977`; as `sigma` approaches `2/3`, the first-order proxy falls below
`6.76e-6` even though the mode never reaches `z=8`. An epsilon threshold alone
can therefore hand off too early. The prototype was removed.

The `freq_res=4`, `z_tail=8`, `h=0.005` spectrum comparison against the saved
`z_tail=9` Reference spectrum (computed at the `h=0.0025` candidate's
background; the `h=0.005` Fast root differs by `6.5e-6` relative) partitions
by maximum mode `z`: all 53 bins over `1e-3 dex` lie among 118 modes with
`z_max<8`; the 142 modes that reach `z=8` have maximum error `4.51e-4 dex`.
This identifies the no-WKB-handoff modes as the affected subset on this
sample, not a global parameter-space rule.

Two selective substep prototypes did not improve those modes. First, halving
only negative-`z` segments left the res4 max spectrum difference at about
`1.98e-3 dex`. Then full propagation modes were selected by `z_max<z_tail` and
substepped on the existing coarse background grid using its endpoint and
midpoint values. The warm 15-pair profile measured P50/P95 `230.5/234.1 ms`,
versus `227.0/231.3 ms` for `h=0.005,z_tail=8`; spectrum max/P95 remained
`1.982e-3/1.699e-3 dex` with 53 bins over threshold. This is evidence that
substepping the tensor state over the coarse `Phi` background does not replace
refining that background grid. The `sigma_exact=True` and
`transition_refine=True` controls at `h=0.005,z_tail=8` also left the compared
spectrum unchanged to the displayed precision.

Uniform `h=0.0025,z_tail=8` does pass the same-background spectrum gate:
max/P95/median `4.554e-4/4.320e-4/1.946e-4 dex`, with zero bins over
`1e-3 dex`. At res4 its warm P50/P95 was `256.8/269.9 ms` versus
`227.0/231.3 ms` for the `h=0.005,z_tail=8` baseline. Its deep-tail Reference
map residual is only `1.984e-4`, and the 2207.5-second Reference pass did not
complete a second map. The res1 self-consistent Reference comparison misses
the DN limit by a narrow margin. Thus even the uniform refinement remains a
research candidate, not a certified Fast optimization.

The per-mode epsilon prototype, both selective substep prototypes, and the
uniform-refinement candidate were not retained in source. The remaining
promising path is a mode-aware refinement of the *background* integration
grid (or an independent continuous-background propagation), followed by the
same deep-tail fixed-point and performance checks; no accuracy gain is
claimed until those checks pass.

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

With `auto_escalate=True, escalate_to_reference=True` and requested
`freq_res=1`, the local budget marks this point `uncertified-path-configuration`
with estimated `DN_gw_error=0.0539` (the default escalation threshold is
`5e-3`). Automatic Reference escalation now raises `freq_res` to at least 8
for uncertified paths (and returned Fast numerical failures), while preserving
caller values above 8. The old direct res1 Reference root near `2.296589` is
not frequency converged: fixed-background maps at res2/res4 differed from
res1 by `1.92e-3`/`2.18e-3`.

At `freq_res=8`, the full self-consistent Reference solve returns
`DN_gw=2.2914118670` on 800 frequencies in 481.6 seconds. At that root, a
fixed-background `freq_res=16` map gives `2.2914094306`, a relative change of
`1.06e-6`; the full `freq_res=8` solve has a self-consistency residual of
`5.0e-10`. This supports convergence of the integrated `DN_gw`. A separate
`z_tail=7` cross-check on all 800 native frequencies had maximum spectrum
difference `9.90096e-4 dex` (even/odd halves), with P95 `6.48e-4` and
`7.10e-4 dex`; at nine nearby frequencies, halving the LSODA step from
`0.00125` to `0.000625` reduced the maximum to `5.40e-4 dex`. That result does
not certify the automatic output, which preserves the requested `z_tail=5`.
For the matching `z_tail=5` convention, the full-grid LSODA comparison at
`h=0.00125` initially exceeded the `1e-3 dex` gate (`1.0522e-3 dex`). Around
its worst point, step refinement gave `1.0522e-3`, `5.3351e-4`, and
`2.6545e-4 dex` at `h=0.00125`, `0.000625`, and `0.0003125`. The complete
800-frequency `h=0.0003125` comparison then had maximum difference
`3.57898e-4 dex`, P95 `2.31351e-4 dex`, and median absolute difference
`7.00581e-5 dex` (Reference 67.6 s, LSODA 287.4 s). This passes the sampled
native-bin spectral gate for the requested tail convention; it does not
certify interpolation between the 800 output bins.

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
result and its spectrum changed by `9.43e-3 dex` on a dense probe. These are
continuous-interpolant convergence diagnostics, distinct from the native-bin
cross-solver gate above. A further
`6.25e-5 dex` run hit the 6000-node hard cap after 7982 Reference evaluations
and 681 seconds. Its DN proxy was within `1.20e-4` of uniform `freq_res=8`,
but the remaining local estimate was `4.35e-3 dex` and its spectrum differed
by `4.38e-3 dex` from the 3467-node grid. This reaches the practical per-point
Reference budget without certifying between-node continuous interpolation; the
sampled 800-bin output passes the refined independent LSODA comparison above.

The real automatic-escalation regressions verify this resolution floor: sample
2 returns the 800-bin `freq_res=8` Reference root (`DN_gw=2.2914118670`) in
483.7 seconds; sample 9 still passes its prior DN anchor in 207.2 seconds.
These costs occur only after explicit
`auto_escalate=True, escalate_to_reference=True` use.

The Reference `freq_res=16` attempt initially hit an artificial
`quad(limit=800)` ceiling with 1600 breakpoints. The limit now scales with the
breakpoint count and retains 50 adaptive subdivisions; an exact 1600-node
regression covers the failure. The high-resolution rerun then completed.

## Reclassified steep-blue cutoff case

The existing `test_fast_guard_off_reports_nonfinite_background_instead_of_index_error`
was introduced at `2cb7a4d`, before the goal-grid transition reserve was clipped
to the valid blue-tilt band. At current HEAD `0814be1`, that sample now returns
a finite 64-bin raw Fast result (`DN_gw=0.108301603`); the old expected
`nonfinite_background` failure is stale. The result is not accurate: its local
budget is `uncertified-path-configuration` with `DN_gw_error=0.01896`; a dense
same-background Reference map differs by `5.18e-4` relative and the native
spectrum differs by `3.05e-3 dex` max. The regression now asserts the finite
result remains uncertified and exercises explicit Reference escalation.

Automatic escalation returns `FAST_ESCALATED` with 554 bins in 254.0 s at
`DN_eff=0.1082445735`, `DN_gw=0.1082445636`. Full-grid LSODA at
`h=0.000625` initially put the DN map `4.21e-4` below the Reference root, just
outside the DN gate; halving to `h=0.0003125` moved it to `-1.91e-4`, narrowly
inside. The final `h=0.00015625` map is `-8.36e-5` relative, with max spectrum
difference `5.38e-4 dex` (P95 `4.09e-4`, median `5.47e-5 dex`). The full LSODA
pass took 455.7 s. This reclassifies the old background failure as a recoverable
but uncertified raw-Fast point; only the explicit Reference path passes both
independent gates.

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
  numerical failures, with `rtol<=1e-9` and a minimum `freq_res=8` on
  uncertified paths, without changing ordinary Fast or physical-rejection behavior.
- The pressure scanner records staged classifications, resource telemetry,
  resumable chunks, and optional coordinatewise minimization.
- Reference frequency quadrature now sizes QUADPACK's interval limit to dense
  breakpoint grids, avoiding the old 800-subinterval failure.
- Adaptive frequency refinement now honors `max_points` as a hard cap and
  prioritizes intervals with the largest sampled or estimated interpolation
  error; this detects missed peaks but does not certify sample 2's spectrum.
- The independent Reference engine now honors the existing `physical_guard`
  option instead of unconditionally rejecting total `Delta N_eff > 5`;
  explicit Reference escalation forwards the same opt-in guard.
- The exact zero-tensor limit (`r*A_s=0`, `cr<=0`) now returns zero energy and
  `DN_gw=0`, with `log10OmegaGW=-inf`; the single-field `cr>0` branch continues
  to require positive `r`. The pressure scanner records this as
  `exact_zero_tensor_spectrum`, separate from invalid/nonfinite output. A
  nonzero `r*A_s` product that underflows is reported as
  `tensor_amplitude_underflow`, never reclassified as the exact-zero limit.
- The Fast entry path adds one exact-zero eligibility check. A two-core,
  25-pair alternating profile compared the real helper with a stub that still
  makes the same function call but returns `False`; all five output digests
  matched for each of six ordinary cases. Across cases, helper-vs-stub warm
  P50 changes ranged from `-1.33%` to `+3.57%`, and P95 changes from `-0.42%`
  to `+4.65%`. This isolates the helper-body overhead conservatively; it is
  not a whole-solver before/after comparison. Candidate absolute warm P50 was
  `6.65–11.30 ms`, P95 `7.24–12.17 ms` in this run.

## High-`DN_gw` Reference recovery

The 514-case screen at `11b2eaf` contains 125 finite Fast outputs with
`DN_gw > 10` and no `max_iter` failures. This does not establish accuracy. A
fresh replay of its `baseline_cr_0` point at `9418e3a`
(`r=0.01`, `n_t=0.5`, `cr=0`, `T_re=2000 GeV`, `kappa10=0.01`) returns
`DN_gw=252043.75` in 0.41 s, but its raw Fast result
misses the dense same-background Reference map by `2.49e-3` relative and the
spectrum by `2.95e-3 dex`. A high-resolution Fast solve (`grid_independent`,
`freq_res=8`) remains inaccurate (`2.28e-3` DN-map difference and `2.95e-3 dex`
max spectrum difference); halving `h` to `0.00125` does not materially change
the discrepancy. Thus `DN_gw_max=10` is not the cause of this point's accuracy
failure. The scan also provides no real `max_iter` counterexample proving the
fixed bisection ceiling is the active failure cause, so the ceiling was not
changed.

The real blocker for robust recovery was an unconditional total-`Delta N_eff`
guard inside `run_reference`, even though the public `physical_guard` option
defaults off. The Reference path now makes that check conditional: with the
guard off, values above 5 are evaluated; with it on, the solver raises an
explicit physical-guard error. Fast→Reference escalation forwards the
caller's setting. A synthetic `DN_gw=12` fixed-point regression covers both
branches.

On this high-blue-tilt sample, explicit automatic escalation now returns a
finite 1568-bin Reference result in 633.3 s at
`DN_eff=251417.467588`, `DN_gw=251417.467439`; relative self-consistency
residual is `5.93e-10`. At the same background, a `freq_res=16` Reference map
differs from the `freq_res=8` map by `1.70e-6` relative. Full-grid LSODA at
`h=0.000625` independently gives a fixed-background `DN_gw` map `1.144e-4`
relative to the Reference root and a sampled spectrum max difference of
`7.760e-4 dex` (P95 `4.167e-4`, median `4.153e-5 dex`), passing both acceptance
gates. LSODA used 242.6 s. The coarser `h=0.00125` run failed both gates,
which motivated the step-refinement check. This validates one difficult point
through the explicit Reference recovery path; it does not promote raw Fast at
this point or establish a general high-`DN_gw` boundary.

A second screen point just beyond the scanner's `N_inf=120` resource preflight
(`DN_re=100`, `N_inf=137.7723`, with the other `baseline_cr_0` parameters)
also returns a finite raw Fast spectrum in 0.31 s, but its `DN_gw=251883.50`
is `2.635e-3` away from the independent Reference root. Explicit escalation
returns a finite 1802-bin Reference result in 752.8 s at
`DN_eff=251221.539082`, `DN_gw=251221.539201` (self-consistency residual
`4.75e-10` relative). Full-grid LSODA at `h=0.000625` takes 884.8 s and gives
`DN_gw=251255.308391`, `1.344e-4` relative to the Reference root; maximum
spectrum difference is `7.623e-4 dex` (P95 `4.502e-4`, median `5.025e-5 dex`).
At its worst frequency, halving the LSODA step reduces the local difference
from `7.623e-4` to `2.980e-4 dex`. This certifies this point beyond the scanner
preflight boundary; it does not certify all `N_inf>120` points or the raw Fast
result.

## Extreme-amplitude runaway classification

The `A_s=1e100`, `cr=1` boundary point (`N_inf=152.227` at input
`DN_eff=0`) is not merely hidden by the `N_inf=120` scanner preflight. A fresh
Fast replay exits with `invalid_iterate_cutoff`: tracing the first outer maps
gives `F(0)=1.07845e106`, then `F(1.07845e106)=1.50199e204`, after which the
derived interval is invalid. With a one-map diagnostic, 20 inputs from
`DN_eff=1` through `1e158` all produce finite maps. On the logarithmic probes
`1e10` through `1e158`, `F(x)/x` ranges from `1.28096e98` to `1.58571e98`,
with no observed fixed-point sign change. The derived `Delta_N>=0` boundary is near
`log10(DN_eff)=158.44162414575`; `N_inf` is `37.77234` just below it and
becomes `None` just above it because `Delta_N<0`.

This is classified as runaway with no fixed point observed inside the sampled
derived-valid interval, not as a `DN_gw_max=10` or `MAX_ITER=60` failure. The
finite log sweep is evidence, not a proof over the continuum, so this case is
not labeled mathematically unsolvable and no approximate root is returned.

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
pytest -q tests/test_physics_limits.py::test_blue_tilt_cutoff_point_is_finite_but_uncertified -o addopts=''
pytest -q tests/test_physics_limits.py::test_blue_tilt_cutoff_point_escalates_to_reference -o addopts=''
pytest -q tests/test_freq_adaptive.py::test_goal_oriented_freqs_never_crosses_blue_tilt_amplitude_cutoff -o addopts=''
pytest -q tests/test_freq_adaptive.py::test_adaptive_detects_narrow_feature_between_seed_nodes -o addopts=''
pytest -q tests/test_reference.py::test_integrate_spectrum_scales_quad_budget_with_breakpoint_count -o addopts=''
pytest -q tests/test_reference.py::test_reference_physical_guard_is_optional_for_large_fixed_point -o addopts=''
pytest -q tests/test_physics_limits.py::test_fast_zero_tensor_limit_returns_exact_zero_spectrum -o addopts=''
pytest -q tests/test_physics_limits.py::test_lsoda_zero_tensor_limit_returns_exact_zero_spectrum -o addopts=''
pytest -q tests/test_reference.py::test_reference_zero_tensor_limit_is_exact_and_grid_safe -o addopts=''
pytest -q tests/test_pressure_fast_domain.py::test_zero_tensor_limit_is_success_with_separate_classification -o addopts=''
pytest -q tests/test_pressure_fast_domain.py::test_underflowed_nonzero_tensor_amplitude_is_not_zero_success -o addopts=''
```

For the high-resolution DN check, run `LCDM_SG(...).SGWB_iter(engine="reference", z_tail=5.0, rtol=1e-9, freq_res=8.0)` with the sample-2 parameters above. It takes about eight minutes on the recorded Windows machine.

The date-stamped per-point records and paired timing JSON used for this
summary were local ignored research artifacts. Re-running the commands creates
fresh records; the counts above must be interpreted with the stated source
commits and screen definitions.

## Fresh branch-priority audit (2026-10-30)

At the subsequent numeric baseline `12c920e`, code inspection and red tests
found that `transition_refine=True,kink_split=True` built a nonuniform grid,
then selected the `kink_split` primitive branch and passed `h_arr=None` to the
tensor kernel. The kernel therefore stepped by the uniform configured `h` on
refined intervals. The same branch ordering made `sigma_exact=True` inactive
whenever `kink_split=True`, while telemetry still recorded it as requested.
The implementation now routes transition-refined grids through their local
interval widths. It uses the node-based fast primitive by default and the
continuous primitive when `sigma_exact=True`; `kink_split` remains covered by
the exact reheating node. The historical A/B sentence above saying those
options did not meet the accuracy gates is not valid evidence for the masked
combination; fresh reference comparisons are required.

A preliminary res4 A/B compared against a fixed-background `z_tail=9`
spectrum evaluated at another Fast background, so it could not certify the
candidate. The later same-background and tail-convergence checks are recorded
below; transition refinement remains an opt-in validation path and the global
Fast domain remains open.
Separately, a fresh `A_s=1e20` trace localizes the first nonfinite value to the
last `S2` background-primitive node on the second outer iterate: `sigma` and
`Phi` are huge but finite before `exp(Psi)` loses representability. A local
H2-log panel correction plus an unbounded root bracket produced a finite Fast
fixed point, but the root changed by `4.99e-3` relative when `h` halved and by
much more as frequency density increased; sampled spectra also shifted by
`2.17e-3 dex`. The independent Reference overflowed in `math.exp(z)`. That
candidate was rolled back and does not establish an accurate root or a larger
reliable region.

### Shared-grid res4 result

The mode decomposition was rerun at the Fast node-primitive candidate
background `DN_eff=0.13374433666555344` on 260 native frequency nodes. The
Reference spectrum/map uses exactly this background and these frequencies.
The uniform
`h=0.005` Fast path differs from the `z_tail=9` Reference by median/P95/max
`2.897e-4 / 2.231e-3 / 2.366e-3 dex`; 66 modes exceed `1e-3 dex`. The
uniform-grid `sigma_exact+kink_split` path differs by
`2.897e-4 / 2.232e-3 / 2.367e-3 dex` and still has 66 modes over `1e-3`, so
continuous-sigma primitive interpolation is not the leading source for this
case. The transition-refined path with `radius=0.005` differs by
`3.969e-5 / 3.581e-4 / 4.198e-4 dex`, with no modes above `5e-4 dex`.
All 66 uniform-path outliers have `zmax<8` and no WKB handoff, isolating the
tested error from the analytic tail. At this exact background the Reference
`z_tail=9` map is `0.133735606198`, `-6.528e-5` relative to the Fast
node-primitive root. A full `z_tail=10` map at the nearby `h=0.0025` Fast root
`0.133745206036` gives `0.133736995454` (`-6.139e-5` relative). On the same
260 frequencies, the Reference spectrum changes from `z_tail=9` to10 by
median/P95/max `5.69e-8 / 7.28e-5 / 2.06e-4 dex`; the radius-`0.005` h=.005
candidate differs from that tail-10 output by max `3.888e-4 dex`, with no
frequency above `5e-4`. The two nearby `z_tail=9` maps give local slope
`dF/dDN_eff≈0.0214`; extrapolating that slope with the tail-10 map estimates
a Reference root near `0.133736816`, `-5.62e-5` relative to the node-primitive
Fast root. This local secant estimate is within the `2e-4` root gate; it is not
a full outer Reference iteration. The direct tail9 map residual and tail-depth
check both stay inside the requested map/root tolerance. The global Fast
domain remains open.

Using the node-based fast primitive instead of the full-grid exact primitive
on the transition-refined grid changes the blue-tilt spectrum by at most
`3.26e-6 dex` and the root by `1.40e-6` relative. On sample11 it changes the
spectrum by at most `3.13e-6 dex` and the root by `2.08e-7` relative. In
25-pair warm profiles it is 1.3% faster at blue res4 P50 and 2.7% at P95; at
sample11 it is 6.1% faster at P50 and 6.4% at P95. Both remain inside their
same-frequency Reference gates; the Reference background differs by only
`1.4e-6` relative.

The second retained stiff failure, random-corpus sample 11 (`r=6.41847e-4`,
`n_t=1.44008`, `kappa10=1.1694`, `T_re=3.10117 GeV`, `cr=0`), was also
replayed with the radius-`0.005` transition-refined path on its 94 native
grid-independent frequencies. The node-fast Fast root is `DN_eff=4.857618635`;
the Reference map `4.857475102` was computed `2.08e-7` relative away at the
full-exact-primitive root `4.857617627`, giving an approximately `-2.955e-5`
map residual at the node-fast value. The
sampled spectrum error median/P95/max is
`1.288e-5 / 4.513e-5 / 7.917e-5 dex`; Reference quadrature and interpolation
estimates are `6.68e-19` and `4.42e-7` in `DN_gw`, and no frequency exceeds
`5e-4 dex`. At that nearly identical background and frequency set, uniform
Fast has max error `7.134e-4 dex`, P95 `6.249e-4 dex`, and 13 frequencies above
`5e-4 dex`. This
recovers a second previously nonfinite/out-of-tolerance counterexample. A
second Reference map at the first map value `4.857475102` returns
`4.857472632` (`-5.084e-7` relative), confirming the self-consistent Reference
root on these 94 native nodes. Other parameter directions remain open.

### Fast-path performance check

Paired warm runs against numeric baseline `12c920e` cover six ordinary Fast
cases, two warmups and 25 alternating pairs per case. Outputs are bitwise
identical. The worst P50 ratio is `1.013` and the worst P95 ratio is `1.017`
(the high-T P95 was separately repeated for 101 pairs and measured `0.999`);
the normal-path change stays within the 5% budget. A same-point res4 profile
with 25 pairs measures uniform Fast P50/P95 `239.1/243.2 ms` and the
transition-refined node-fast candidate `265.7/276.4 ms` (about 11.1%/13.7% slower).
This misses
the aspirational 100 ms difficult-point target, so the candidate is not made
the default. Profiles exclude cold JIT startup; thread scaling was not
measured in this focused comparison.

## Current 402-point boundary rescreen

The deduplicated `boundaries` suite was rerun at source commit `2b7316a`; its
402 inputs contain 401 completed rows and one manually stopped resource case.
Counts are: 297 finite successes (295 positive spectra plus two exact-zero
limits), 57 solver failures (39 nonpositive H2, 15 nonfinite numerical states,
3 invalid outer cutoffs), 22 exceptions (15 solver exceptions, 7 parameter
derivation exceptions), 14 invalid-cutoff preflight rejects, 10 resource
preflight skips, one returned nonfinite spectrum, and one resource timeout.
The timeout is boundary case 129 (`A_s=100,cr=1`), stopped after more than 220
CPU seconds; it is not counted as a mathematical no-solution case. Completed
records have runtime P50/P95 `5.24/14.60 ms` and peak observed RSS `183.9 MiB`
for a two-thread, single-process scan.

Of the 297 finite successes, 93 have `DN_gw>10`; the largest is
`7.78e96`. All 295 nonzero spectra remain
`uncertified-path-configuration`; only the two exact-zero spectra are
analytically certified. The scan maps computability and failure classes but
does not certify those extreme values' precision.

The fresh 140-point coupled suite at `8e26916` covers the deduplicated
`cr × effective n_t × kappa10 × T_re` combinations: 140/140 outputs are finite,
with zero solver failures. All 140 remain `uncertified-path-configuration`;
local `DN_gw_error` estimates range from `4.05e-3` (median `3.04e-2`, P95
`1.003`) to `1.023`. On an uncontended single-process/two-thread run, per-case
P50/P95 is `5.86/11.51 ms`, throughput `112.4 cases/s`, and observed peak RSS
`179.5 MiB`. This is computability/sampler evidence only, not an accuracy pass.

The current Sobol and Latin-hypercube pilots (`power=4`, seed `20261030`) each
sample 16 points after canonicalizing the unused `n_t` coordinate for `cr>0`.
Sobol outcomes are 6 resource preflight skips, 6 invalid-derived rejects, 2
derivation/solver exceptions, and 2 nonpositive-H2 failures; Latin-hypercube
outcomes are 5 resource skips, 6 invalid-derived rejects, 4 exceptions, and 1
nonpositive-H2 failure. Neither pilot contains a finite output, so these broad
signed/log pilots do not add a precision-certified region.

For the branch-priority correction, 25-pair warm Fast timing across six named
normal cases gives identical spectrum arrays. P50 ratios to numeric baseline
are `0.996/1.013/0.982/0.999/1.001/1.007`; P95 ratios are
`0.940/0.955/1.005/0.974/0.948/1.017`, with the high-T P95 from a separate
101-pair replay. These stay within the 5% normal-path target. At the difficult
blue-tilt res4 point the transition-refined candidate P50/P95 is
`265.7/276.4 ms` versus `239.1/243.2 ms` uniform (11.1%/13.7% slower), so it remains an
explicit validation path and misses the aspirational 100 ms target. Cold JIT
startup and thread scaling are not measured in this focused profile; paired
profile RSS is about `234 MiB` with two-thread affinity.

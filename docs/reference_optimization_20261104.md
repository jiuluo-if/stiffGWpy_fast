# Reference acceleration and algorithm screens — 2026-11-04

## Scope and baseline

The Reference path remains an independent continuous-background DOP853
calculation, separate from the fixed-step Fast kernel. Its exact-input res4
tail10 map at `DN_eff=0.133732083004284` is archived with 260 modes and took
`5691.7 s` with four workers at `rtol=1e-11`.

An instrumented one-worker sample at the same `DN_eff`, `z_tail=10`, and
`rtol=1e-11` attributed the sampled cost to oscillatory propagation:

| `log10(f/Hz)` | Runtime | DOP853 `nfev` | Background evaluations | Tail used |
|---:|---:|---:|---:|:---:|
| -18.4594 | 0.0255 s | 2,093 | 2,701 | no |
| -16.5303 | 0.2075 s | 18,452 | 19,060 | no |
| -14.6013 | 17.1879 s | 1,398,719 | 1,399,327 | no |
| -12.6723 | 19.3333 s | 1,708,643 | 1,710,368 | yes |
| -10.7729 | 11.4255 s | 1,350,578 | 1,353,041 | yes |

The fixed 200-step start search took only about 3–4 ms per mode. It is not a
material acceleration target. The heavy modes use over one million DOP853 RHS
evaluations each.

## Prüfer coordinate backend

From the existing tensor system, set `x=A sin(theta)` and `y=A cos(theta)`.
The transformed equations are

```text
z'          = 1.5 sigma - 1
(log A)'    = 1.5 sigma - 2 + cos(2 theta)
theta'      = -exp(z) - sin(2 theta)
```

This is an exact change of state variables. The tail convention remains the
existing phase-averaged frozen-amplitude approximation. The implementation is
available as an explicit opt-in: `backend="prufer"`. The public default stays
Cartesian until broader parameter and root coverage is collected.

At the independently closed res4 Reference root
`DN_eff=0.13373681817755848`, the production backend ran the same 260 frequencies
as the archived Cartesian tail10 spectrum. With candidate `rtol=1e-9`,
`atol=[1e-12,1e-12,1e-12]`, one worker, it returned
`DN_gw=0.1337368198597089`; the Cartesian root map at that input is
`0.13373681807168084`. The candidate-to-Cartesian map difference is
`1.337e-8` relative, and the candidate fixed-point residual is `1.258e-8`
relative. Same-input spectrum median/P95/maximum differences are
`1.465e-7 / 2.602e-6 / 4.279e-5 dex`. Candidate quadrature/interpolation
estimates are `1.84e-20 / 1.83e-10`. Candidate propagation time was `657.4 s`;
this full-grid run was not paired with a fresh Cartesian run on the same
resource assignment.

For a repeated paired sample, four representative modes were run twice each,
alternating Cartesian and Prüfer order. Candidate `rtol=1e-9` was compared to
Cartesian `rtol=1e-11` at the same `DN_eff`, frequencies, and tail. Aggregate
speedup was `4.269x`; the eight per-mode ratios had P50/P95 `4.144x / 4.704x`.
The maximum spectrum error was `1.401e-5 dex`. The lowest-frequency endpoint
mode was `2.03x` faster, so the measured ratio is not a guarantee for every
mode or regime.

Run records:

- Full production backend, root spectrum and map:
  [`reference_prufer_production_res4_refroot_20261104.json`](reference_prufer_production_res4_refroot_20261104.json)
- Repeated paired subset:
  [`reference_prufer_production_ab_smoke_20261104.json`](reference_prufer_production_ab_smoke_20261104.json)
- Driver: [`benchmark_res4_prufer_candidate.py`](../scripts/benchmark_res4_prufer_candidate.py)
- Original Cartesian root/spectrum anchor:
  [`res4_reference_secant_root_tail10_20261104.json`](res4_reference_secant_root_tail10_20261104.json)

## ARDC screen

The official `riccati==2.0.0` package by Agocs and Barnett was installed into a
temporary target only; the repository gained no package dependency. For the
tensor equation `T'' + (3 - 1.5 sigma) T' + exp(2z) T = 0`, its documented form
uses `gamma=1.5-0.75 sigma` and `omega=exp(z)`. A res4 handoff mode took about
`0.145 s` at `eps=epsh=1e-9`, with `3.585e-4 dex` difference from Cartesian.
Tightening the package tolerances to `1e-10` and `1e-11` did not materially
reduce that observed difference.

The 260-mode ARDC run did not reach completion within the 45-minute CPU budget.
It was stopped while still active; no partial spectrum or precision result was
produced. The timeout record is
[`reference_ardc_res4_timeout_20261104.json`](reference_ardc_res4_timeout_20261104.json).
It is an incomplete resource screen, not a full-grid numerical rejection.

## DOP853 dense-output overhead

The Cartesian Reference solver requested `dense_output=True` and used the
interpolant only once, to recover the tail-event state. SciPy 1.17.1's local
DOP853 implementation builds a dense interpolant at every accepted step; that
construction evaluates three additional RHS stages. With a terminal event,
`solve_ivp` still computes the event interpolant and exposes the event state
through `y_events` when global dense output is disabled. The Cartesian path now
uses that state and requests `dense_output=False`.

The alternating same-process A/B used the same res4 background, frequencies,
`z_tail=10`, `rtol=1e-11`, and Cartesian equations. A no-tail mode and a
tail-handoff mode were run twice per variant. All output triplets matched
bitwise. The no-tail mode improved from `18.54/18.55 s` P50/P95 to
`14.55/14.56 s` (`1.275x`); its `nfev` fell from `1,398,878` to `1,119,218`.
The tail mode improved from `21.15/21.40 s` to `16.37/16.41 s` (`1.292x`);
`nfev` fell from `1,708,436` to `1,366,937`. Sampled peak RSS delta fell from
`81.5` to `29.4 MB` and `94.0` to `35.2 MB`, respectively. These are two
difficult modes, not a full-grid speed guarantee. The paired record and
reproduction driver are
[`reference_dense_output_benchmark_20261104.json`](reference_dense_output_benchmark_20261104.json)
and [`benchmark_reference_dense_output.py`](../scripts/benchmark_reference_dense_output.py).

Because this change also accelerates the Cartesian baseline, the earlier
Prüfer A/B speedup was remeasured against the new default. On the same four
res4 frequencies, two alternating pairs per mode, candidate `rtol=1e-9` versus
Cartesian `rtol=1e-11`, the Prüfer aggregate speedup is now `3.349x`; the
per-mode ratio P50/P95 is `3.085x/3.990x`. Maximum spectrum error against the
archived same-input Cartesian spectrum is `1.401e-5 dex`. The low-frequency
endpoint is only `1.736x` faster. These four modes do not provide a full-grid
speed distribution. Updated samples are in
[`reference_prufer_production_ab_postdense_20261104.json`](reference_prufer_production_ab_postdense_20261104.json).

## Compiled scalar background-spline evaluation

Profiling the Cartesian DOP853 RHS found the Fermi–Dirac density and pressure
`CubicSpline` scalar calls consumed `2.59 s` of a `7.35 s` no-tail solve and
`2.43 s` of a `7.58 s` tail solve at `rtol=1e-9` (about 35% and 32%). A
Numba-compiled evaluator now applies binary interval search and Horner
evaluation to the exact four coefficients of the same SciPy cubic segment.
No spline table, coefficients, or background equations are changed. Across
250,000 random points spanning both spline tables, the maximum absolute
difference from SciPy was `2.22e-16`, including tested extrapolation edges.

Alternating A/B on the same res4 inputs, `DN_eff=0.133732083004284`,
`z_tail=10`, `rtol=1e-9`, three pairs per mode:

| Mode | SciPy P50/P95 | Compiled P50/P95 | Median speedup | Max paired output difference |
|---|---:|---:|---:|---:|
| `log10(f/Hz)=-14.6013`, no tail | `10.22/10.43 s` | `6.19/6.44 s` | `1.65x` | `1.12e-8` relative |
| `log10(f/Hz)=-12.6723`, tail | `11.21/11.29 s` | `6.98/7.10 s` | `1.61x` | `2.29e-7` relative |

The Reference `nfev` was unchanged for the no-tail mode; the tail mode changed
by 84 evaluations out of 767,000 as roundoff changed a few adaptive step
decisions. At `rtol=1e-11`, comparison with the archived same-input Cartesian
res4 spectrum at four frequencies (two tail and two no-tail) had maximum error
`1.93e-8 dex`. First-call JIT compilation took `0.31 s`. Cold dispatch uses the
non-JIT scalar evaluator below the measured frequency crossover and prewarms
once before a spectrum batch containing expensive modes; a cheap low-frequency
single-mode solve took `17.6 ms` and did not compile. A cold 4-worker two-mode
smoke matched serial outputs bitwise; it is not a scaling claim. Exact-zero
tensor spectra bypass JIT. The Fast path does not use this evaluator.
Fresh-process measurements are retained in
[`reference_spline_cold_low_20261104.json`](reference_spline_cold_low_20261104.json)
and [`reference_spline_cold_high_20261104.json`](reference_spline_cold_high_20261104.json),
with reproduction driver
[`benchmark_reference_spline_cold.py`](../scripts/benchmark_reference_spline_cold.py).

Raw paired timing, archived-spectrum checks, source hashes, and resource
metadata are in
[`reference_spline_eval_benchmark_20261104.json`](reference_spline_eval_benchmark_20261104.json);
reproduce with [`benchmark_reference_spline_eval.py`](../scripts/benchmark_reference_spline_eval.py).
After this interpolation change, the production Prüfer backend was run over
the full 260-mode tail10 grid at the independently closed Cartesian Reference
root `DN_eff=0.13373681817755848`, with candidate `rtol=1e-9` and one worker.
Its map differs from the archived Reference map by `4.386e-7` relative, the
fixed-point residual at that input is `-4.394e-7`, and same-input spectrum
P95/max error is `6.651e-6 / 5.246e-5 dex`. Candidate propagation took `569.1 s`
including a `0.325 s` spline JIT warmup. This is a candidate-only run, not a
paired full-grid Cartesian speed ratio; the paired run and its secant-based
root-shift estimate are reported in the following section. The local root
derivative remains unmeasured. Raw per-mode data and provenance are in
[`reference_prufer_spline_res4_refroot_20261104.json`](reference_prufer_spline_res4_refroot_20261104.json).

## Archived Reference map secant sensitivity

The two archived 260-frequency Cartesian tail10 maps at
`DN_eff=0.133732083004284` and `0.13373681817755848` use byte-identical
`reference.py` sources and bitwise-identical frequency arrays. Their secant
slope is `d(DN_gw)/d(DN_eff) ~= 0.0211261`; applying the fixed-point correction
`residual/(1-slope)` to the second map estimates a root shift of
`-1.082e-10`. This helps quantify residual amplification for this res4 pair.
It is an interval-average slope across two stored maps, not a local derivative,
uncertainty bound, or rigorous root bracket. The second map's protocol metadata
is less complete than the first map's, so this result remains a derived
diagnostic rather than a new root certification. Reproduce with
[`analyze_reference_map_sensitivity.py`](../scripts/analyze_reference_map_sensitivity.py);
inputs, source hashes, and derived values are in
[`res4_reference_map_sensitivity_20261104.json`](res4_reference_map_sensitivity_20261104.json).

The crossover is a latency dispatch only. A fresh low-frequency single-mode
call at `log10(f/Hz)=-18.4594` used the Python scalar evaluator, completed in
`17.6 ms`, and left the Numba dispatcher uncompiled. A fresh expensive mode at
`-14.6013` compiled on demand and completed in `5.89 s`. A four-frequency
parallel smoke included the one-time `0.296 s` warmup and returned bitwise the
same spectra as serial, but took `34.79 s` versus `17.50 s`; this confirms
correctness and shows no speedup for that small batch. The parallel sample is
not a thread-scaling estimate. Its record and driver are
[`reference_spline_workers_20261104.json`](reference_spline_workers_20261104.json)
and [`benchmark_reference_spline_workers.py`](../scripts/benchmark_reference_spline_workers.py).

Finally, Prüfer was remeasured against the Numba-accelerated Cartesian path on
the same four res4 frequencies, two alternating pairs each, at candidate
`rtol=1e-9` and Cartesian `rtol=1e-11`. Aggregate speedup is `3.137x`, per-mode
ratio P50/P95 `3.030x/3.539x`, and maximum spectrum error `5.22e-7 dex`; the
low-frequency endpoint is `1.69x` faster. The four-mode result does not
replace the full-grid map result above. Raw samples are in
[`reference_prufer_production_ab_post_spline_20261104.json`](reference_prufer_production_ab_post_spline_20261104.json).

An expanded frequency-spread A/B used indices `[0,65,100,130,160,195,220,259]`
from the same 260-mode archive and two alternating pairs at each index. The
aggregate speedup is `3.626x` excluding the separately logged `0.334 s` JIT
warmup (`3.59x` including it); per-mode speedup P50/P95 is `2.778x/4.174x`.
The maximum same-input spectrum error is `3.94e-6 dex`. This supports a roughly
3.6x gain for this sampled eight-mode workload; the lower P50 and endpoint
speedups show it is not a per-frequency guarantee. Raw samples are in
[`reference_prufer_production_ab_8mode_20261104.json`](reference_prufer_production_ab_8mode_20261104.json).

## Full-grid paired Cartesian A/B

The production Prüfer backend was paired against the spline-accelerated
Cartesian implementation once on each of the same 260 archived frequencies,
alternating execution order by mode. Both used the independent root input
`DN_eff=0.13373681817755848` and `z_tail=10`; Prüfer used `rtol=1e-9`,
Cartesian `rtol=1e-11`. All 260 pairs completed. The paired spectrum error
median/P95/maximum is `1.575e-7 / 6.654e-6 / 5.246e-5 dex`. Total serial time
is `574.13 s` versus `1992.45 s`, an aggregate `3.470x` speedup (`3.468x`
including the separately measured `0.313 s` JIT warmup). Per-mode speedup
P50/P95 is `2.632x/4.332x`; these data cover one parameter point and do not
establish general performance scaling.

Reintegrating the stored per-mode logarithmic spectra reproduces the
benchmark's Prüfer map to `4.2e-16` relative and gives the current Cartesian
map `0.1337368185045801`, within `3.24e-9` relative of the archived Cartesian
map. Prüfer and current Cartesian map outputs differ by `4.42e-7` relative;
the corresponding first-order root shift using the archived interval secant
slope `0.0211261` is `-6.04e-8` in `DN_eff`. This secant is not a local
derivative or rigorous root enclosure. The paired raw data and reconstruction
are
[`reference_prufer_cartesian_paired_fullgrid_20261104.json`](reference_prufer_cartesian_paired_fullgrid_20261104.json)
and
[`reference_prufer_cartesian_paired_fullgrid_analysis_20261104.json`](reference_prufer_cartesian_paired_fullgrid_analysis_20261104.json).

## Newton–Kantorovich phase-transfer single-mode prototype

Read Stojimirovic and Bremer's [original paper](https://arxiv.org/html/2409.18487v2),
sections 2.4, 4, 5.1–5.4 and 6.3. It solves a discretized Riccati equation on
Chebyshev nodes using Newton–Kantorovich updates, with a Liouville–Green
initial guess; this is distinct from the project's earlier Kummer phase
post-fit. The authors require smooth positive coefficients on high-frequency
intervals and warn that phase derivatives become steep near turning points.

The prototype uses the damped tensor equation and removes its first-derivative
term to form `Q=exp(2z)-p'/2-p^2/4`, `p=3-1.5*sigma`. It first integrates the
production `z'=1.5*sigma-1` equation and the original `(T,T')` state with DOP853
to `z=5`; six phase-transfer segments then reach the tail event at `z=10`.
Every segment uses an analytic branch derivative, and the `nu=0.1` background
switch is split explicitly; the right side uses the neutrino-spline branch.
Two tail-using frequencies were screened at this same res4 point:
`log10(f/Hz)=-12.6723` crosses the `nu=0.1` boundary, while `-10.772906`
reaches `z=10` before it. With orders 16/24/32, errors versus current
Cartesian (`rtol=1e-11`) were `1.7461e-8` and `4.1274e-9 dex`; the first mode's
order-to-order spread is `<=5.4e-15 dex`. Dense Riccati residuals were
`<=2.72e-14` and `<=2.25e-14`. Latest total candidate times including the low
DOP853 prefix were `0.138/0.139 s`, versus `12.66/7.57 s` Cartesian. Spline JIT
warmup (`~0.31 s`) was measured separately. This remains two-mode, one-point
evidence; the method is not integrated into the production backend. Reproduce
with [`prototype_nonosc_phase_res4_mode.py`](../scripts/prototype_nonosc_phase_res4_mode.py)
and `--frequency -10.772906114757099`; the records are
[`nonosc_phase_res4_mode_20261104.json`](nonosc_phase_res4_mode_20261104.json)
and [`nonosc_phase_res4_highfreq_mode_20261104.json`](nonosc_phase_res4_highfreq_mode_20261104.json).

## Original papers read so far

The task attachment says a 37-paper appendix exists, but the attached file
contains no bibliography entries. The following named originals were read at
the relevant equations, algorithm descriptions, and assumptions; this list is
not a claim that the missing 37-paper set or a 50-paper survey is complete.

1. Stojimirovic and Bremer, “An accelerated frequency-independent solver for
   oscillatory differential equations,” arXiv:2409.18487v2, published in IMA
   Journal of Numerical Analysis (2026), DOI `10.1093/imanum/draf128`.
   Its positive-coefficient smooth high-frequency theorem and Newton–
   Kantorovich Chebyshev solution of the Riccati phase derivative differ from
   the project’s old retrospective Kummer/Chebyshev fits. It does not directly
   cover a turning point or a reheating derivative kink.
   [Original PDF](https://arxiv.org/pdf/2409.18487)
2. Chow and Bremer, “Airy Phase Functions,” arXiv:2503.02306v2. The method
   assumes a smooth coefficient with a simple turning point and constructs an
   Airy phase via its Airy–Kummer equation; it is not the Round60 Airy+WKB
   splice. Coefficient smoothness at reheating remains an adaptation issue.
   [Original PDF](https://arxiv.org/pdf/2503.02306)
3. Piron and Tacu, “An Explicit Numerical Scheme for Milne’s Phase-Amplitude
   Equations,” arXiv:2411.02621v1, published in *Atoms* 13(6):57 (2025), DOI
   `10.3390/atoms13060057`. It replaces the nonlinear Milne amplitude equation
   by Kiyokawa’s linear third-order equation for squared amplitude; this is
   distinct from the project’s Round62 Ermakov IVP.
   [Original PDF](https://arxiv.org/pdf/2411.02621)
4. Agocs and Barnett, “An adaptive spectral method for oscillatory second-order
   linear ODEs with frequency-independent cost,” arXiv:2212.06924v1, with the
   authors’ `riccati` implementation described in JOSS 8(86):5430 (2023), DOI
   `10.21105/joss.05430`. ARDC’s local defect-correction phase equation accepts
   nonzero damping, but its smooth-coefficient conditions and this screen’s
   full-grid cost limit its immediate use here.
   [Original PDF](https://arxiv.org/pdf/2212.06924),
   [JOSS software paper](https://doi.org/10.21105/joss.05430)
5. Bremer, “Phase function methods for second order linear ordinary
   differential equations with turning points,” arXiv:2209.14561 (2022).
   The construction assumes a smooth coefficient with one turning point,
   combines adaptive Chebyshev ODE solves with coefficient windowing (§4), and
   proves in §3 that an even-order turning point cannot have a single
   nonoscillatory phase across the whole line. The Appell solution can grow rapidly on the
   exponential side. Multiple turns and the reheating kink therefore require
   piecewise phases and matching; this is not a direct production algorithm.
   [Original full text](https://arxiv.org/abs/2209.14561)
6. Bremer, “On the numerical solution of second order differential equations
   in the high-frequency regime,” arXiv:1409.6049 (2014). Theorem 8 assumes a
   smooth, strictly positive coefficient and a rapidly decaying
   Schwarzian-derived correction after a Liouville coordinate transform; the
   numerical method (§4) windows the coefficient, solves logarithmic Kummer
   equations, then spectrally integrates the phase; experiments are in §5. It is foundational for
   the later phase-function work but does not cover a turning point or
   coefficient jump without a separate extension.
   [Original full text](https://arxiv.org/abs/1409.6049)
7. Bremer, “On the numerical calculation of the roots of special functions
   satisfying second order ordinary differential equations,” arXiv:1512.08357
   (2015). Once a nonoscillatory phase is available, roots are obtained by
   inverting the monotone phase (§4), avoiding large-argument trigonometric loss;
   the per-root work is frequency independent under the same positive-
   coefficient assumptions. This targets root locations, not propagation of
   the tensor state, and is not a new propagation backend.
   [Original full text](https://arxiv.org/abs/1512.08357)
8. Chen, Serkh and Bremer, “On the adaptive Levin method,” arXiv:2211.13400
   (2022). Its §3 analysis proves that Chebyshev collocation with truncated SVD
   avoids low-frequency breakdown in the stated Levin setting, including
   stationary points; §7.3 tests stationary-phase integrals. This is an
   oscillatory-integral method rather than a propagation
   method; it may inform future integral work, but its assumptions do not
   directly replace the current spectrum integral.
   [Original full text](https://arxiv.org/abs/2211.13400)
9. Serkh and Bremer, “Phase function methods for second order inhomogeneous
   linear ordinary differential equations,” arXiv:2211.13744 (2022). It
   combines a phase-function basis with adaptive Levin integrals (formula (4),
   §3) for a slowly varying forcing and positive, slowly varying coefficient;
   numerical tests are in §4. The authors note
   that large negative coefficients make one phase-basis function grow rapidly
   and render this representation unstable; the current tensor equation is
   homogeneous and crosses such regimes, so the construction is not adopted.
   [Original full text](https://arxiv.org/abs/2211.13744)
10. Agocs, Handley, Lasenby and Hobson, “An efficient method for solving
    highly oscillatory ordinary differential equations with applications to
    physical systems,” arXiv:1906.01421 (2019, `oscode`). The solver targets
    `x''+2 gamma x'+omega^2 x=0`, can consume numerically tabulated background
    coefficients, and alternates adaptive fifth-order RK with WKB steps using
    an order-difference error estimate. It was applied to the Mukhanov–Sasaki
    equation. The authors caution that tolerances below `1e-6` need a higher
    order RK pair, accumulated WKB phase integrals above about `1e12` exceed
    double-precision phase resolution, and dense arbitrary-point output was
    not available in that version. It is a strong independent candidate for
    single-mode comparison, but reheating kinks require explicit partitioning
    and its documented error ceiling must be tested against the present gates.
    [Original full text](https://arxiv.org/abs/1906.01421)

The exact 37-item appendix remains absent from the supplied attachment and
repository history. The ten entries above are distinct primary works with
targeted full-text reading; they do not reconstruct that list and do not meet
the requested 50-source audit.

## Current limits

- The new backend is opt-in. The Cartesian default and public output meanings
  are unchanged.
- The full 260-mode Prüfer/Cartesian comparison is complete only for the
  res4 Reference root. Broader parameter, turning-point, zero-amplitude,
  strong-coupling, tail-depth, tolerance, batch and thread coverage remains.
- The global Fast reliable parameter domain remains uncertified.
- The missing 37-paper bibliography must be recovered before that requirement
  can be audited honestly; reference performance comparison and full citation
  tracking continue in the active research plan.

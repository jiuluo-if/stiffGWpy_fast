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

## Current limits

- The new backend is opt-in. The Cartesian default and public output meanings
  are unchanged.
- Accuracy is fully checked only at the res4 Reference root and on a small
  paired timing subset; broader parameter, turning-point, zero-amplitude,
  strong-coupling, tail-depth, tolerance, batch and thread coverage remains.
- The global Fast reliable parameter domain remains uncertified.
- The missing 37-paper bibliography must be recovered before that requirement
  can be audited honestly; reference performance comparison and full citation
  tracking continue in the active research plan.

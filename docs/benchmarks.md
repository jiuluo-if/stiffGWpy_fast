# Benchmarks

Status: current summary of the audited scope; results remain commit- and protocol-bound

Date: 2026-09-24

Code version: see `docs/validation/validation_manifest.json` → `commit` and each run artifact's `commit`

## Current formal fast audit

The current formal fast configuration uses `frequency_quadrature=pchip` with
`h=0.005`, `col_step=8`, `z_tail=5`, `phase_max=0.25`, exact kink split, and
the goal frequency grid. The latest scoped measurements are summarized in
[`fast_v02_audit_report.md`](fast_v02_audit_report.md):

| Measure | Observed | Scope / decision |
|---|---:|---|
| Oracle C WKB `DN_gw` relative error, four named points | median `1.1e-5`, max `1.66e-4` | Passes the `<2e-4` gate on these four points |
| Independent reference, same grid, `z_tail=8`, six named points | median `2.93e-4`, max `2.96e-4` | Above the gate; tail-convention limited in this comparison |
| Default point vs independent reference at `z_tail=10` | `4.96e-6` | Passes the `<2e-4` gate for this point |
| Nested native-frequency refinement, six named points | relative changes `7.4e-12`–`2.8e-9` | Native-grid refinement converged in this test |
| Formal default warm runtime | median `4.77 ms`, p95 `5.58 ms` | The stable `<4 ms` target is not met |
| Stability screen | 0 numerical failures; 3 explicit physical guards / 24 points | No silent fallback in the audited sample |
| Production full parameter-space accuracy | — | `NOT VERIFIED` |

The scoped accuracy comparison supports PCHIP for the tested cases. The
independent reference at `z_tail=8` has a frozen-handoff tail systematic of the
same order as its difference from fast; the deeper `z_tail=10` default-point
comparison and Oracle C anchor agree closely. See the audit report for the
full attribution and remaining uncertainty.

## PCHIP and Simpson comparison

The paired 50-repeat comparison in
[`fast_quadrature_default_switch_assessment.md`](fast_quadrature_default_switch_assessment.md)
reports PCHIP / Simpson warm-median ratios of `1.019–1.035` at 2 threads and
`0.987–1.069` at 16 threads for default, low-T, high-T, and stiff cases. PCHIP
is the current default; Simpson remains explicitly selectable. The Oracle C
relative errors for PCHIP are `1.07e-5` (default), `1.66e-4` (low-T),
`9.88e-6` (high-T), and `1.15e-5` (stiff).

## Dated measurements and comparability

The 2026-09-10 fixed-resource benchmark recorded `4.26 ms/point` warm median
at 20 threads. It is a dated snapshot and is not directly comparable with the
latest formal audit, which used a different audited setup. The 2026-09-03
before/after report is also a historical comparison; its protocol and
limitations are documented in
[`performance_comparison_20260903.md`](performance_comparison_20260903.md).

Later `profile_fast_breakdown_*` artifacts include diagnostics through
2026-09-23. They use different thread/resource settings and profiling
protocols, so they do not replace the formal benchmark or establish a new
cross-version speed claim. Browse these by family in
[`experiment_catalog.md`](experiment_catalog.md).

# Engineering audit checklist

Status: baseline review for `main` at commit `8499ed27d2b7a76002bf95072a7e8947ed1df775`.

The working tree already contained the untracked oracle-240 validation artifact
and driver before this audit; they are intentionally preserved.  Baseline
evidence: `python -m pytest -q` = `99 passed, 6 deselected`; sdist and wheel
builds succeed; the repository-wide Ruff check is not yet clean.  “未验证” means
that the repository contains a claim or mechanism, but this audit has not yet
produced sufficient executable evidence for it.

## Current staged evidence

The first implementation stage has now landed locally: `FastSolverConfig` is a
validated immutable per-call record; named fast modes resolve without mutating
legacy globals; the high-level fast default is explicitly `production`; Numba
thread counts are restored on exceptional exits and configured calls are
serialized while that process-wide setting is active; packaged resources use
`importlib.resources`; and an installed-wheel smoke script checks data outside
the checkout. The maintained Ruff surface, the narrow mypy gate, the full
default tests, the Cobaya marker, and the wheel smoke have passed locally.
The Python-version matrix, full legacy Ruff cleanup, concurrent-call stress,
and universal oracle claims remain unverified until their dedicated checks run.

The next CI stage is now wired locally: the main matrix emits a Python 3.11
coverage artifact, and `.github/workflows/slow.yml` runs the slow numerical
marker on demand or weekly. A green local run is recorded below only after the
workflow-equivalent commands have been executed here.

The workflow-equivalent local evidence is now: `python -m pytest -q -m slow` =
6 passed; the coverage-enabled default suite = 103 passed and 64% total line
coverage. Coverage is reported for visibility only; no unsupported percentage
threshold is presented as a scientific correctness gate.

The committed manifest now has an executable structural validator covering its
schema/version, provenance header, oracle semantics, profile configuration,
status counts, and accuracy sections. It intentionally validates structure and
does not promote any existing `NOT VERIFIED` or `FAIL` scientific result.

Repository hygiene and release notes are now documented in `CONTRIBUTING.md`
and `CHANGELOG.md`; generated coverage/build files are ignored, while the two
pre-existing oracle artifacts remain untouched and unstaged.

## P0 — correctness and trust

| ID | Files | Problem / risk | Recommendation | Acceptance |
|---|---|---|---|---|
| P0.1 | `fast_sgwb.py`, `reference.py`, `stiff_SGWB.py`, `LCDM_stiff_Neff.py`, `functions.py`, `global_param.py`, `freq_adaptive.py`, `cobaya/stiffGW.py` | The numerical chain spans legacy globals, mutable model attributes, fast kernels, reference/original engines and Cobaya; the full invariant audit is incomplete. | Map inputs, units, grids, ODE state, outer closure and output contract before each refactor. | Public engine contract and invariant tests cover all three engines; remaining gaps are recorded as 未验证. |
| P0.2 | `fast_sgwb.py`, `stiff_SGWB.py` | Module-level `_THREADS`, `_COL_STEP`, `_FAST_H`, `_Z_TAIL`, `_PHASE_MAX`, `_FREQ_GRID` are mutable process state; calls can affect later calls. | Introduce immutable per-call configuration and restore compatibility globals on exceptional paths. | Two independent solver configurations can be constructed and repeated without cross-call drift; thread setting restoration is tested. |
| P0.3 | `fast_sgwb.py`, `reference.py`, `functions.py` | Failure paths mix `None`, flags and partially written model attributes; finite-value and array-shape checks are not a single contract. | Add explicit validation/result classification at engine boundaries. | NaN/Inf, invalid `r`, invalid cutoff, empty/unsorted grid and non-convergence have deterministic classifications and no stale outputs. |
| P0.4 | `fast_sgwb.py`, `reference.py`, `docs/accuracy.md`, `README.md` | Local telemetry and fiducial-calibrated terms are separated in prose, but the runtime data structure still mixes `kind=local/calibrated` with certification labels. | Make the distinction machine-readable and preserve “not a universal bound” in API/docs. | Tests assert local vs calibrated labels; no doc claims a universal per-point bound. |
| P0.5 | `stiff_SGWB.py`, `cobaya/stiffGW.py`, `stiffGW.yaml` | Default/sentinel resolution and fast/production mapping rely on truthiness and process-global setters. | Centralize resolution: mode → preset → explicit overrides; keep production defaults explicit. | YAML default `engine: fast` with sentinel knobs resolves to production settings; explicit zero/invalid values are rejected or documented. |
| P0.6 | `reference.py`, `exact_background.py`, `docs/accuracy.md` | Oracle frozen-tail and WKB sensitivity is documented, but a repeatable executable gate for oracle self-sensitivity is not yet part of CI/nightly. | Add a slow validation command/schema for `z_tail` and tolerance variants. | Oracle A/B/C output records the measured deltas and status; no oracle result is called exact without the caveat. |

## P1 — tests and CI

| ID | Files | Problem / risk | Recommendation | Acceptance |
|---|---|---|---|---|
| P1.1 | `.github/workflows/ci.yml`, `pyproject.toml` | CI runs only Python 3.11 on Linux. | Add a supported-version matrix for 3.9–3.13 after dependency compatibility is measured; keep an explicit slow/nightly job. | Matrix is green for every supported interpreter, or an incompatibility is pinned and documented. |
| P1.2 | `.github/workflows/ci.yml`, legacy core modules | Five core modules are permanently excluded from Ruff. | Add staged per-file Ruff rules and burn down the legacy debt in small patches. | Exclusion list shrinks each stage and no new violations enter maintained modules. |
| P1.3 | `pyproject.toml`, public API, Cobaya adapter | No type-checking gate exists. | Add a lightweight type checker and annotate public API/config/adapter first. | Type gate runs on the maintained surface with an explicit legacy boundary. |
| P1.4 | `.github/workflows/ci.yml`, `scripts`, `tests` | Wheel is built but no isolated installed-wheel smoke test runs. | Install the wheel into a clean environment and import/load every packaged data family. | Smoke test succeeds outside the checkout and fails if package data is missing. |
| P1.5 | `tests`, `scripts/validate_*` | Unit tests exist, but production/reference regression, edge, invariants, determinism and thread-variation coverage is incomplete. | Add focused tests before each numerical change and retain golden points. | Production/reference golden diffs, edge guards, finite arrays and repeated/thread checks are gated. |
| P1.6 | `scripts/validate_*`, `.github/workflows` | Slow oracle, Sobol, axis-edge and posterior checks are mostly manual. | Define nightly/dispatch validation tiers with fixed seeds and machine-readable results. | Each slow tier can be invoked reproducibly and records rejected-point classes. |
| P1.7 | `scripts/bench_*`, `.github/workflows` | Benchmark claims are artifacts, not a noise-tolerant regression gate. | Add a separate benchmark workflow with broad thresholds and cold/warm separation. | Small timing jitter does not fail CI; material regression does. |
| P1.8 | `.github/workflows`, `pyproject.toml`, `tests` | Coverage is not a required report and fallback/error branches need attention. | Publish coverage and target physics/error/fallback branches. | Coverage artifact is generated; critical branches have named tests. |

## P2 — architecture

| ID | Files | Problem / risk | Recommendation | Acceptance |
|---|---|---|---|---|
| P2.1 | `fast_sgwb.py` | 1187-line module combines presets, kernels, grids, integration, errors and orchestration. | Extract configuration, grids, integration, error model and orchestration behind stable imports; move kernels last. | Public imports remain compatible; each extracted module has focused tests. |
| P2.2 | `fast_sgwb.py`, `stiff_SGWB.py` | Process-global solver configuration is the main control surface. | Use a frozen dataclass; environment variables only seed defaults at import. | Concurrent model objects do not overwrite configuration metadata; compatibility setters are deprecated wrappers. |
| P2.3 | `fast_sgwb.py`, `cobaya/stiffGW.py` | Telemetry/error budgets are untyped dictionaries and ad-hoc attributes. | Introduce typed records with a compatibility mapping at the boundary. | Stable serialized schema and backward-compatible keys are tested. |
| P2.4 | `stiff_SGWB.py`, `cobaya/stiffGW.py` | Cobaya is coupled to internal model attributes and solver setters. | Define a stable result/API boundary and keep adapter translation there. | Adapter tests use only public contract fields. |
| P2.5 | `scripts`, `tests`, `docs/validation` | Validation drivers duplicate parameters, metadata and classification logic. | Extract shared validation schema/utilities. | Two drivers share one classifier/metadata builder; output remains reproducible. |
| P2.6 | all engines | Exceptions and `None` do not consistently distinguish physical invalidity, guard, fast failure and oracle failure. | Establish a small exception/result taxonomy while preserving legacy return behavior during deprecation. | Cobaya maps only recoverable failures to fallback; physical guards never retry. |
| P2.7 | legacy modules and aliases | Legacy compatibility/style debt is broad and undocumented. | Inventory public aliases and add deprecation paths before removing anything. | No unexplained alias remains; deprecations have tests and removal targets. |

## P3 — reproducibility

| ID | Files | Problem / risk | Recommendation | Acceptance |
|---|---|---|---|---|
| P3.1 | `docs/validation`, `scripts` | Manifest is useful but schema/version validation is implicit. | Define and validate a machine-readable schema containing commit, environment, CPU, threads, config, seed, sample and gates. | Invalid manifests fail the build; every artifact identifies its schema version. |
| P3.2 | `README.md`, docs, manifest builder | Human-readable numbers can drift from artifacts. | Generate README/doc snippets or add a consistency checker. | CI detects changed artifact numbers not reflected in docs. |
| P3.3 | `docs/mcmc_posterior`, `.gitignore` | Large chain/data artifacts and canonical evidence are not formally separated. | Classify canonical artifacts vs reproducible temporary outputs and document retention. | Git contains only justified canonical files with provenance. |
| P3.4 | `stiffgwpy/cobaya/likelihoods`, docs | Bundled observational data provenance/license/checksum is incomplete. | Add a data manifest with source, version, license, checksum and regeneration/acquisition method. | Every bundled observational file has a traceable manifest entry. |
| P3.5 | `scripts/validate_*`, posterior scripts | Seeds exist in many scripts but sampling/rejection metadata is not unified. | Record seed, design, domain and rejected categories in the common schema. | Re-running with the same seed reproduces the design and classifications. |
| P3.6 | fast/reference validation artifacts | Golden coverage is concentrated around a small number of points and some stated gates are NOT YET VERIFIED. | Keep a small stable golden set plus independent axis/interior samples; label unverified claims. | Changes require golden diffs and do not silently promote historical or partial evidence. |

## P4 — packaging, API and docs

| ID | Files | Problem / risk | Recommendation | Acceptance |
|---|---|---|---|---|
| P4.1 | `pyproject.toml` | Dependency bounds, metadata and dev tooling are minimal; setuptools emits license deprecation warnings. | Measure supported dependency ranges, add metadata/tool configuration and migrate to SPDX metadata. | Clean build metadata for supported toolchain; compatibility is tested rather than guessed. |
| P4.2 | `pyproject.toml`, Cobaya docs | `mpi4py` is forced by the Cobaya extra although serial operation may not need it. | Split serial Cobaya and MPI extras after import/runtime testing. | Serial Cobaya smoke works without MPI; MPI path is an explicit extra. |
| P4.3 | package data and resource access | Data is listed in package-data and currently enters the wheel, but code uses path joins and no sdist/wheel contract test exists. | Use `importlib.resources` for required resources and test both artifacts. | Clean installed wheel loads all `.dat/.npz/.npy/.yaml/.txt` resources. |
| P4.4 | README, docs | Quick start, production/reference/Cobaya/error/performance docs exist but API surface and release workflow are incomplete. | Add API reference, recommended production config and error-handling examples. | A new user can install, run fast/production/reference and interpret limits. |
| P4.5 | `__init__.py`, docs | Public API is not formally declared; internals such as `_MAX_THREADS` are used by adapter code. | Define `__all__`/public config API and hide internal implementation details. | Tests and adapters do not require private variables. |
| P4.6 | repository root, `.github/workflows` | No explicit changelog/version/release workflow gate. | Add changelog and release workflow including lint, tests, build, installed smoke and numerical gates. | Release workflow is documented and reproducible. |

## P5 — performance

| ID | Files | Problem / risk | Recommendation | Acceptance |
|---|---|---|---|---|
| P5.1 | `scripts/profile_fast_breakdown.py`, fast solver | Optimization evidence exists, but future changes can still be intuition-driven. | Require profile plus accuracy diff before changing kernels. | Every optimization records numerical and timing deltas. |
| P5.2 | benchmark scripts/docs | Cold JIT, warm single point, batch, scaling and end-to-end Cobaya are not one formal schema. | Standardize benchmark dimensions and metadata. | Reports separate cold/warm/thread/batch/Cobaya measurements. |
| P5.3 | fast kernels/grid/integration | Potential allocation/grid/JIT/PCHIP hotspots need current profiling. | Profile before refactoring; inspect cache and repeated background work. | Material optimization is tied to a measured hotspot. |
| P5.4 | production config and benchmarks | Speedup must not be obtained by silently lowering production precision. | Gate accuracy first, then compare broad benchmark thresholds. | Production numerical gate remains unchanged or the change is rejected. |

## P6 — hygiene and maintenance

| ID | Files | Problem / risk | Recommendation | Acceptance |
|---|---|---|---|---|
| P6.1 | legacy core/docs/scripts | Dead code, stale comments, duplicate constants and historical references remain possible. | Remove only after usage search and regression coverage. | No removed symbol is referenced; history/provenance is retained where needed. |
| P6.2 | repository root, docs | No CONTRIBUTING or validation-update workflow is present. | Add contributor, environment, test-tier and validation-refresh instructions. | A maintainer can reproduce gates and update artifacts safely. |
| P6.3 | `.github/workflows`, dependency metadata | Action/dependency supply-chain maintenance is not systematically documented. | Pin stable action major/minor policy and audit dependencies. | CI uses supported pinned actions; audit result is recorded. |
| P6.4 | `LICENSE.md`, bundled data | Code license and observational-data redistribution terms need explicit reconciliation. | Review each data license and document restrictions. | Release gate fails on missing/ incompatible data licensing evidence. |

## Staged implementation order

1. P0.2/P0.5: isolate configuration resolution and add regression tests without
   changing numerical defaults.
2. P0.3/P1.4/P4.3: make failure/resource contracts executable and add installed
   wheel smoke coverage.
3. P1.2/P1.3/P1.5/P1.8: stage lint, type and numerical/invariant coverage.
4. P2/P3: extract typed boundaries and shared validation metadata in small
   compatibility-preserving patches.
5. P4/P5/P6: complete release/docs/benchmark/maintenance gates.

Claims not backed by a current artifact remain explicitly marked `未验证` until
the corresponding command produces evidence.

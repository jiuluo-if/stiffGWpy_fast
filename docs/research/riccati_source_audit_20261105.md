# `riccati` source and reproduction audit — 2026-11-05

## Snapshot and execution

- Official repository: [fruzsinaagocs/riccati](https://github.com/fruzsinaagocs/riccati), source commit `dc9d370a0a51e4c30e8eadd941214e543eb0545e` (2025-05-21). The latest visible release tag is `v2.0.0`.
- License: Apache-2.0. Runtime dependency is NumPy; the official tests use SciPy, mpmath, and pytest. No C/C++ extension is required.
- Source inspected: `riccati/evolve.py`, `step.py`, `stepsize.py`, `solversetup.py`, `chebutils.py`, `tests/test_riccati.py`, and the Airy `docs/quickstart.ipynb` / `docs/examples.ipynb`.
- Ran the official `test_solve_airy` and `test_solve_burst`: **2 passed**, 8 Python invalid-escape deprecation warnings, 0.94 s. Both exercise the package’s mixed Riccati/Chebyshev solver against known Airy/burst solutions.
- Also ran a standalone exact harmonic-oscillator IVP, `u''+40²u=0`, `u(0)=1`, `u'(0)=40i`, on `[0,10]`; the solver took one Riccati step and returned zero reported state and derivative error against `exp(40it)` at its internal nodes. This is a deterministic math smoke test, not a cosmology benchmark.

## Implementation trace

`solversetup` precomputes Chebyshev differentiation, interpolation, and integration matrices. `evolve.solve` switches between `osc_step` (Riccati defect correction) and `nonosc_step` (Chebyshev collocation), and reports accepted step types and work counters. The oscillatory step estimates the Riccati residual and decreases the step if the estimate fails; the step-size chooser also checks interpolation error for `w` and `g` at interstitial nodes. Chebyshev steps double the node count until endpoint solutions from successive grids agree or the configured `nmax` is exhausted.

The implementation has a residual-based local screen, but not a certified global state or final-integral error bound. It uses Python/NumPy orchestration and dense spectral linear algebra, so its single-mode cost cannot be compared directly with the batched Numba kernel without an end-to-end measurement.

## Project fit and decision

The Fast tensor equation can be mapped to the paper’s scalar damped form on smooth coefficient segments. The package’s oscillatory/nonoscillatory switching is relevant to the horizon transition, and exact reheating still has to remain a segment boundary. The algorithm, variables, and residual mechanism are exactly the ARDC/Riccati-Chebyshev family already being assessed in Phase 6c, so this run verifies that upstream method’s implementation and analytic smoke behavior rather than proposing another implementation. It does not test the exact cosmological coefficient, full-frequency integral, self-consistent root, or ≤5% speed gate.

Decision: keep `riccati` as a source/reproduction reference; do not add it as a production dependency. Any continuation should reuse the current Phase 6c candidate and complete its full-spectrum, root, and runtime gates rather than duplicate the package in the Fast path.

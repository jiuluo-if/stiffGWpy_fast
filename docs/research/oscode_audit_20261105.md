# `oscode` source and applicability audit — 2026-11-05

## Source snapshot

- Official repository: [fruzsinaagocs/oscode](https://github.com/fruzsinaagocs/oscode), cloned with its Eigen submodule into a temporary directory outside the project.
- Audited source commit: `4dd643b48596138861aca7fa8ad4b75e07f35c38` (`2025-05-21`, “Add sphinx rtd theme to reqs”). The latest visible release tag is `v1.1.2`; `setup.py` on the audited branch declares package version `1.3.0`, so the source and release metadata are not aligned.
- License: repository `LICENSE` is BSD 3-Clause. The Python extension includes Eigen and compiles as C++17 (`setup.py`).
- Source inspected: `include/oscode/solver.hpp`, `wkbsolver.hpp`, `rksolver.hpp`, `system.hpp`, `pyoscode/__init__.py`, `tests/test_airy.py`, `examples/cosmology.ipynb`, `README.rst`, `setup.py`, and `LICENSE`.

## Algorithm trace

`oscode` solves `y'' + 2 g(t)y' + w(t)^2 y = 0`. On every attempted step the C++ driver computes both an embedded Runge–Kutta step and a WKB step, scales their local error estimates by `rtol`, predicts each method’s next step size, and accepts the method with the larger predicted step. WKB uses an order-1/2/3 asymptotic series and Gauss–Lobatto coefficient samples; truncation and quadrature errors are tracked separately. RK handles intervals where the WKB estimate is poor. Array inputs are interpolated by the `system.hpp` helpers; callable inputs use a separate Python wrapper path.

The package is not an error-controlled black box: the Python API documents `atol` as not implemented, and the WKB validity depends on smooth, slowly varying frequency. The supplied Airy test uses `w(t)=sqrt(t)` on `[1,100]`, so its tested interval starts after the simple turn at zero; it does not test a full turn crossing. The Cosmology notebook is an official example, not an independent guarantee for this project’s exact piecewise background.

## Environment and execution result

- Environment: Python `3.11.9`; `pyoscode` is not installed. No `cl`, `g++`, `clang++`, or `gcc` compiler was discoverable on `PATH`.
- PyPI index reports `1.1.2` as its latest published package, while the repository setup file says `1.3.0`. An attempt to download the source package invoked pip’s isolated build-dependency setup rather than returning promptly; that exact pip process tree was stopped. No package was installed, no repository source was modified, and no official example or upstream test is claimed as run.
- The Eigen submodule was fetched at `e7248b26a1ed53fa030c5c459f7ea095dfd276ac`.

## Project fit and decision

The damped scalar mode equation has a formal mapping with `g=γ` and `w²=q` on a smooth segment. The adaptive RK/WKB comparison is mathematically relevant as a local propagation alternative. Near the horizon turn, however, the oscillatory WKB expansion loses its slow-variation basis and the method must fall back to RK; the reheating discontinuity must remain an exact segment boundary. Array-mode use would also linearly interpolate the supplied coefficients, so feeding a sparse Fast background table would need a separate interpolation-error analysis. Calling this C++ extension once per frequency would not reuse the project’s Numba batched kernel and adds an unmeasured cross-language orchestration cost.

Decision: do not add `oscode` as a production dependency or treat its upstream performance claims as Fast evidence. Its local error split and automatic RK/WKB choice are useful design references, but the project already has substantial WKB/phase experiments. Any new project implementation must first state what changes in the transfer representation and pass exact-turn/kink, full-spectrum, fixed-point-root, resource, and ≤5% trusted-region runtime gates.

## Remaining project audit

This is one source-level project audit. The task’s other named repositories (`riccati`, Chebfun, OrdinaryDiffEq.jl, DifferentialEquations.jl, Tasmanian, python-flint/Arb, AUTO-07p) have not yet been audited or run. The full literature list and domain-level acceptance work remain open.

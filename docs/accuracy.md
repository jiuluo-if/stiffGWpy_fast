# Accuracy boundaries

The solver does not claim uniform accuracy over the full parameter space. Accuracy statements must be tied to the engine, configuration, reference convention, parameter point, and frequency grid that were actually checked.

The default `fast` profile is checked by the core regression and physics-limit tests. In particular, `tests/test_fast_residual_estimator.py` and `tests/test_wkb_singular_mode_completion.py` compare selected outputs against fixed independent reference samples in `tests/data/`. These checks protect known numerical behavior; they do not certify an untested continuous domain.

Use the independent `reference` engine for scoped precision comparisons. Treat `lsoda` as a compatibility and regression path, not as independent truth. For new accuracy claims, retain the input parameters, exact configuration, reference convention, source commit, and machine-readable comparison data in the private research archive before proposing a focused regression fixture.

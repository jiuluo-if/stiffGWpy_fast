# Pre-registered experiment: linear-z Magnus transfer map

## Hypothesis

On each native fast step, the transformed tensor system has a traceless
matrix `A(z)` with `z(N)` varying smoothly.  Replacing the constant-midpoint
transfer by the second-order Magnus map for linearly interpolated `z` may
reduce propagation/model-alignment error without adding real tensor solves.

## Acceptance criteria

The standalone candidate is accepted for further work only if all of the
following hold on default, low-T, high-T, stiff, and high-kappa:

1. handoff amplitude and phase agree with Cartesian DOP853 at least 10x better
   than the current constant-midpoint map on representative native modes;
2. full-grid DN and spectrum improve, with DN relative error below `2e-4`;
3. no physical-guard or numerical-failure classification changes;
4. a 25-repeat runtime A/B shows stable improvement or no more than 5% cost;
5. repeated runs are deterministic under the fixed resource budget.

Failure of any criterion rejects production integration.  The prototype must
remain reference-only until all criteria and the full local gate pass.

## Follow-up: closed Numba implementation

The Python `expm` implementation was replaced by the analytically equivalent
closed exponential of a traceless 2x2 Magnus matrix.  The current-SHA result
is recorded in `docs/transfer_map_phase_numba_round.json`.

This removes the avoidable matrix-library overhead, but it does not satisfy the
experiment gate: the maximum amplitude change versus the midpoint map is only
`1.49e-5`, the maximum phase change is `7.11e-4 rad`, and the 25-repeat
propagation medians are about `4.72-5.01 ms` versus `2.56-2.65 ms` for the
baseline across the five regimes.  The candidate is therefore rejected as a
production algorithm and no full-grid promotion is authorized.

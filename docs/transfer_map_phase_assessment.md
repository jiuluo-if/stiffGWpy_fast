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

# Fast adaptive frequency support with exact mode-row reuse

## Decision

Keep the public `goal` preset on its existing direct propagation path. Add a more accurate, explicitly selected `freq_grid="adaptive"` path for ordinary tensor amplitudes. It starts from the sparse goal support, refines using the existing quadrature-weighted curvature estimator, and reuses already propagated frequency rows when the cosmological background and solver state are identical.

This is an adaptive calculation capability and a correction to prior frequency-sensitive results. It is **not** counted as an expansion of the current HEAD parameter-evaluation domain: current `goal` already returns a result on these samples. No continuous parameter domain is certified by this finite screen.

## Mathematical equivalence of row reuse

For fixed background samples and fixed solver controls, the tensor propagation kernel advances each frequency mode independently. Its row at frequency `f_i` depends on that frequency's initial state and the shared arrays (`Nv`, `Phi_grid`, `Phi_mid`, `S2`, `S2inv`, optional `h_arr`/`Phi_quarters`) but not on any other frequency row.

Therefore, if `DN_eff`, the frequency coordinate, all background arrays, and all propagation controls match exactly, the existing `Ogw`, `Oj`, `Opgw`, handoff and tail-state rows are valid without approximation. The implementation checks those equalities before copying rows into the expanded support. It propagates only the inserted frequencies. If any compatibility condition fails, it recomputes the full support. Log-amplitude and phase-weight paths do not use this cache.

The support error estimate remains the existing PCHIP curvature proxy. In ascending log-frequency coordinates, it estimates local interpolation error as `|y''| h²/8`, weights by each interval's share of the integral, and refines intervals above `freq_grid_target`. This estimate is not a theorem bounding all unresolved features; validation below is finite evidence.

## Paired parameter-slice evidence

The archived CR0 density–tilt rectangle has 9 rectangular samples and 2 independent holdouts. The earlier pre-stabilization baseline recorded `nonfinite_background_primitive` at all 11 samples. The later goal-grid implementation succeeded, but the adaptive results shift the current-goal roots by `1.02%–1.19%`, so same-goal-grid checks were not frequency-convergence evidence.

With `freq_grid_target=1e-6`, adaptive Fast succeeded at all 11 inputs. Current `goal` also succeeded at all 11; adaptive roots differ from those goal roots by `1.02%–1.19%`. Against same-source grid-independent Fast at `freq_res=4.6`, the maximum root difference was `7.22e-5`; at the 265–271 shared native frequency nodes, the maximum spectrum difference was `3.58e-6 dex`. Median adaptive runtime was `0.111 s`; the dense Fast comparison median was `0.257 s`. These dense-grid comparisons are not independent Reference results. The adaptive mode is opt-in and about three times the current goal median on this sample; the default path remains unchanged.

Two exact native-support Reference checks used the unchanged Reference production code (`z_tail=8`, `rtol=1e-9`, two workers):

| Sample | Adaptive nodes | Reference map residual | Max spectrum error | Reference time |
|---|---:|---:|---:|---:|
| Rectangle corner `(Omega_bh2=13.99609375M, n_t=-0.0302734375)` | 265 | `8.45e-5` | `6.11e-4 dex` | `66.3 s` |
| Independent holdout `(Omega_bh2=13.89013671875M, n_t=-0.022705078125)` | 269 | `8.38e-5` after second-map closure | `6.33e-4 dex` | `66.7 s` per map |

The second Reference map on the holdout changed the first mapped root by `3.86e-8` relative. Both Reference checks meet the requested root and spectrum limits. They validate these representatives, not every point in the rectangle.

## Performance A/B

Seven alternating runs at the holdout input, two threads, `freq_grid_target=1e-6`:

| Variant | Warm P50 | Warm P95 sample max | Nodes | Propagated rows | Root |
|---|---:|---:|---:|---:|---:|
| Goal seed + row reuse | `0.107 s` | `0.110 s` | 269 | 613 | `1.3331521683` |
| Goal seed + full repropagation control | `0.193 s` | `0.203 s` | 269 | 1840 | `1.3331521683` |
| Grid-independent seed + row reuse | `0.147 s` | `0.153 s` | 365 | 926 | `1.3330693196` |

The default `goal` path does not enter the subset helper and calls `solve_kernel` directly. A regression test enforces this.

The final alternating HEAD/candidate comparison covered three profiles with 24 warm samples each, two threads, and identical roots. Warm P50/P95 sample-max ratios were `0.994/0.862` (trusted default), `1.013/0.928` (trusted high-T), and `0.984/0.958` (CR0 density edge). All P50 regressions are below 5%, and all sampled P95 values improved. Each process discarded its first solve; the extracted HEAD cache was retained across process rounds.

## Files and verification

- Paired 11-point screen and source hashes: [screen JSON](screen_cr0_density_slice_adaptive_frequency.json), [runner](screen_cr0_density_slice_adaptive_frequency.py).
- Exact mode-row reuse A/B and raw timings: [benchmark JSON](benchmark_adaptive_support_reuse.json), [runner](benchmark_adaptive_support_reuse.py).
- Default-goal HEAD/candidate performance check: [benchmark JSON](benchmark_goal_hot_path_head_ab.json), [runner](benchmark_goal_hot_path_head_ab.py).
- Independent Reference maps and hashes: [rectangle corner](validate_adaptive_frequency_cr0_density_reference.json), [holdout](validate_adaptive_frequency_holdout_reference.json), [holdout closure](close_adaptive_frequency_holdout_reference_root.json).
- After tightening the cache key, the rectangle and holdout's current frequency supports, roots and total `DN_eff` were checked bitwise against the saved Reference inputs; this avoided repeating their expensive Reference runs. See [reuse verification](verify_adaptive_reference_reuse.json).
- Coverage and error figure: [PNG](adaptive_frequency_cr0_density_slice.png).
- Fast source SHA256: `013214fc7063ef73d8e5a91b4274d8b82987aa6d0b66a40faf25f49b182a40fe`.
- Reference source SHA256: `9066adcf69f95398405e77e5a5d47942573041f13d0a0ee1554c4bd119eb17f6` (unchanged).
- Verification: focused adaptive and default-path tests pass; default suite `446 passed, 10 deselected`; Ruff, Mypy, sdist/wheel build and isolated wheel import pass.

## Limits and next work

The change does not alter default-goal outputs and does not yet expand the current HEAD domain. The sampled CR0 rectangle is a re-evaluation of an earlier background-stabilization tranche. A next tranche must start from failure records for current HEAD, establish mathematical validity and cutoff constraints first, then use this adaptive path to check accuracy. The `cr>0` admissible region still needs a separate sampled slice and independent holdout evidence.

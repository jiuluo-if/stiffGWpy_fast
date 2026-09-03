# Plain-grid tier (fast) vs continuous-sigma reference — matched z8

Plain-grid engine (`accuracy_mode="fast"`: h=0.02, col_step=8, no transition_refine, phase_max=0, construct grid) validated against the independent continuous-sigma DOP853 reference (rtol=1e-9) on the SAME `z_tail=8` at the plain-grid's OWN frequency nodes — the residual isolates plain-grid engine error (no frequency-grid interpolation).

| label | n_freq | signal rel max | transition rel max | DN_gw rel | status |
| --- | ---: | ---: | ---: | ---: | --- |
| stiff | 241 | 1.768e-02 | 1.768e-02 | 2.061e-03 | ok |
| default | 246 | 1.867e-02 | 1.620e-02 | -7.345e-03 | ok |
| lowT | 241 | 4.786e-02 | 3.069e-02 | 9.265e-04 | ok |
| cr0_blue | 222 | 1.736e-02 | 1.736e-02 | -9.142e-03 | ok |
| extreme | 236 | 1.718e-02 | 1.718e-02 | -1.017e-02 | ok |
| highT | 242 | 1.660e-02 | 1.612e-02 | -7.681e-03 | ok |
| tiny_r | 243 | 7.019e-02 | 2.609e-02 | -2.725e-02 | ok |
| transition | 245 | 2.440e-02 | 2.440e-02 | -1.705e-02 | ok |
| rad_dominant | 240 | 6.751e-02 | 6.751e-02 | -2.664e-02 | ok |

Aggregates over 9 points: signal-band rel max **7.019e-02** (FAIL <1e-3), transition-band rel max **6.751e-02** (FAIL <1e-3), integrated Delta_Neff rel abs median **9.142e-03** / p95 **2.701e-02** / max **2.725e-02** (<1e-4: FAIL).

Reference runtime per point (workers=3) median 383 s；fast plain-grid 的 `0.76 s`
是本 correctness artifact 生成时的 JIT 前历史测量，不代表当前 runtime。
当前 cold/warm/batch 性能见 `docs/performance_comparison_20260903.md`；本文件的
误差 envelope 仍是有效的 plain-grid validation 结果。

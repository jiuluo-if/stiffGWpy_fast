# 按今天的 `z` 延迟 WKB 交接筛查

日期：2026-11-15

## 候选准则

全局 `z_tail=8` 对本轮最差下溢频点很准，但慢约 1.57 倍。单模态状态显示，在 `z_today` 低于 5 的模式中，Fast z5 可能先越过 `z=5` 并做 WKB 终端外推，但该模态到今天又回到 `z_today≈4.99`；Fast z8 不再触发终端尾，继续数值传播。于是筛查按模式组合两个固定根结果：若 `z_today<5` 用 z8 的行，否则用 z5。这个判据针对已观察到的“交接后没有保持深亚视界至今天”情形；它还不是一条普遍 WKB 余项界。

## 离散结果

阈值 5 在物理振幅 `cr>0` 样本仍有 `0.001429 dex` 最大误差，未过门槛；将延迟范围扩到 `z_today<6` 后该点通过，因此以下统一列出阈值 6：

| 样本 | 网格 | z8 行比例 | 混合谱最大 / P95 差 | 混合映射相对 Reference 差 | 固定根残差 |
|---|---:|---:|---:|---:|---:|
| `n_t=17.38, r=5e-324, T_re=2000, kappa10=.01, cr=0` | 937 | `5.55%` | `0.000517 / 0.000148 dex` | `+4.27e-7` vs z8 | `−2.303e-6` |
| `n_t=17.38, r=1e-316, T_re=2000, kappa10=.01, cr=0` | 937 | `5.55%` | `0.000483 / 0.0000304 dex` | `−3.25e-6` vs z10 | `−3.25e-6` |
| `n_t=17.30, r=5e-324, T_re=2000, kappa10=.01, cr=0` | 235 | `5.53%` | `0.000459 / 0.000148 dex` | `+8.18e-7` vs z8 | `0` |
| `n_t=17.45, r=5e-324, T_re=10, kappa10=.001, cr=0` | 221 | `5.88%` | `0.000578 / 0.000150 dex` | `+1.02e-5` vs z8 | `−2.22e-16` |
| `n_t≈−0.00125, r=.01, T_re=2000, kappa10=.01, cr=1` | 76 | `9.21%` | `0.000374 / 0.000283 dex` | `−4.96e-6` vs z10 | `+4.55e-11` |

五个离散样本的混合最大谱差都低于 `1e-3 dex`，映射/根误差低于 `2e-4`。最小次正规 937 点与 `cr>0` 常规振幅点有直接 Reference z10；其余点用 Reference z8。Reference z8→z10 收敛没有在三个后续点重新测。覆盖了一个 `cr>0` 点，但这不是一致性分支的下溢振幅，因为该极小 `A_t` 会被 `invalid_cutoff` 拒绝。

另尝试 `n_t=17.45,T_re=1e4,kappa10=1,r=5e-324`，关闭物理 guard 的探索性运行给出 `DN_eff≈1.383e9`，超过物理 guard 上限 5。该组合被记为物理 preflight 拒绝，所有谱差结果作废，不计作候选通过。未通过绕过背景 guard 扩展适用域。

## 每模内核原型与 A/B

已增加私有 opt-in 路径：传播核可接收 `z_tail_by_mode`，Fast 默认仍传原标量 `z_tail`；仅设置模型私有属性 `_fast_adaptive_handoff_below=6`、`_fast_adaptive_handoff_to=8` 时启用。具有相位分辨权重的模式不会延迟，避免 phase-tail state 使用不一致的 handoff 定义。TDD 先在旧核上以“不接受 `z_tail_by_mode`”失败，再加入 per-mode 参数；两模测试确认不同交接阈值分别生效，另有测试确认 phase-resolved 活跃行仍使用原阈值。

937 点最小次正规样本做了同根、同网格、单线程、每变体预热后交替 12 组 A/B：

| 路径 | P50 / P95 | RSS P50 / P95 | 最大谱差对 Reference z10 | 映射根残差 |
|---|---:|---:|---:|---:|
| 全局 z5 | `0.2270 / 0.2471 s` | `180.41 / 198.20 MiB` | `0.002441 dex` | `−2.303e-6` |
| 全局 z8 | `0.3079 / 0.3563 s` | `180.92 / 197.36 MiB` | `0.0003722 dex` | `−1.964e-6` |
| per-mode z5 / 延迟交接 | `0.2283 / 0.2341 s` | `178.63 / 196.78 MiB` | `0.0004855 dex` | `−2.303e-6` |

per-mode 相比 z5 的 P50 运行比为 `1.006`（基本持平），P95 约低 `5.2%`；相比全局 z8，P50 快 `1.35×`、P95 快 `1.52×`。三种变体 12 次重复谱的最大离散度为 0。第一次调用另列：z5 `0.483 s`、z8 `0.308 s`、per-mode `0.238 s`，包含各路径首次调用/JIT 成本，不计入预热后统计。该节是一次参数点、单线程的实测，不外推为普遍提速。

同一冻结根、937 点网格、每个线程数预热并交替 12 组三变体，扩到 1/2/8 线程后结果如下：

| 线程 | z5 P50 / P95 (s) | per-mode P50 / P95 (s) | z8 P50 / P95 (s) | per-mode/z5 P50 | per-mode/z8 P50 |
|---:|---:|---:|---:|---:|---:|
| 1 | `0.1400 / 0.1708` | `0.1389 / 0.1940` | `0.2154 / 0.2694` | `1.007×` | `1.551×` |
| 2 | `0.1224 / 0.1380` | `0.1226 / 0.1239` | `0.1745 / 0.1811` | `0.998×` | `1.426×` |
| 8 | `0.0994 / 0.1131` | `0.0986 / 0.1100` | `0.1147 / 0.1223` | `1.008×` | `1.164×` |

三种线程数下 per-mode 的 `max_abs_spectrum_error=0.0004855 dex`、固定根残差 `−2.303e-6`、重复谱离散度 0；z5 超谱差门槛，z8 与 per-mode 通过。RSS P50（z5/per-mode/z8）分别为 1 线程 `179.5/178.2/181.2 MiB`、2 线程 `188.9/189.6/190.4 MiB`、8 线程 `189.7/188.3/191.1 MiB`。因此这一个样本上 per-mode 保持 z5 的稳定成本而接近 z8 的精度；它不是普遍加速结论，需在不同参数、频率密度和 phase-tail 活跃状态复测。

实际 per-mode 内核的直接 Reference z10 留出如下：

| 样本 | 网格 | 最大谱差 / P95 | 映射差 vs Ref z10 | 固定根残差 |
|---|---:|---:|---:|---:|
| `n_t=17.38,r=5e-324,cr=0` | 937 | `0.0004855 / 0.0000317 dex` | `+4.27e-7` | `−2.303e-6` |
| `n_t=17.30,r=5e-324,cr=0` | 235 | `0.0004746 / 0.0000312 dex` | `−4.48e-6` | `0` |
| `n_t=17.45,T_re=10,kappa10=.001,r=5e-324,cr=0` | 221 | `0.0005692 / 0.0000279 dex` | `−2.25e-6` | `−2.22e-16` |
| `n_t≈−.00125,r=.01,cr=1` | 76 | `0.0003741 / 0.0002831 dex` | `−4.96e-6` | `+4.55e-11` |

四个 holdout 的 Reference z8→z10 已直接检查。`n_t=17.30` 点的 z8→z10 最大 / P95 变化为 `0.0001550/0.0001418 dex`，映射变化 `−5.30e-6`；低温 `n_t=17.45` 点为 `0.0001648/0.0001598 dex`、映射 `−1.245e-5`。最小次正规与 `cr=1` 点的收敛见各自研究记录。

per-mode 路径仍是**私有 opt-in 实验路径**，没有改变正式默认或公开接口。接受为生产方案前还需：不同频率网格 Reference z10 留出、更多 `cr>0` 参数点、phase-resolved 活跃/非活跃边界测试及多资源档位 A/B。此之前不得声明连续适用域或修改默认 `z_tail`。

复跑文件：[mask 屏查脚本](../../scripts/underflow_adaptive_handoff_mask_screen_20261115.py)、[倾斜 holdout 脚本](../../scripts/screen_adaptive_handoff_tilt_holdout_20261115.py)、[背景参数 holdout 包装脚本](../../scripts/screen_adaptive_handoff_background_holdout_20261115.py)、[低温 holdout Reference z10](../../scripts/validate_adaptive_handoff_background_z10_20261115.py)、[`n_t=17.30` holdout Reference z10](../../scripts/validate_adaptive_handoff_tilt_z10_20261115.py)、[`cr>0` 分支筛查](../../scripts/screen_adaptive_handoff_consistency_branch_20261115.py)、[单线程 A/B](../../scripts/benchmark_adaptive_handoff_per_mode_20261115.py)、[1/2/8 线程 A/B](../../scripts/benchmark_adaptive_handoff_resources_20261115.py)。逐点 JSON 见[mask](results/underflow_adaptive_handoff_mask_screen_20261115.json)、[单线程 A/B](results/adaptive_handoff_per_mode_benchmark_20261115.json)、[多线程 A/B](results/adaptive_handoff_resource_benchmark_20261115.json)、[倾斜](results/adaptive_handoff_tilt_holdout_20261115.json)、[低温 holdout](results/adaptive_handoff_background_z10_20261115.json)、[`n_t=17.30` holdout](results/adaptive_handoff_tilt_z10_20261115.json)、[`cr>0`](results/adaptive_handoff_consistency_branch_20261115.json)和[物理 preflight 拒绝](results/adaptive_handoff_physical_preflight_rejection_20261115.json)。

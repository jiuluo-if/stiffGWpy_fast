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

实际 per-mode 内核在其他四个离散点也已运行；最大谱差分别为 `0.000517 dex`（最小下溢 937 点）、`0.000459 dex`（`n_t=17.30`）、`0.000578 dex`（`n_t=17.45,T_re=10,kappa10=.001`）、`0.000374 dex`（正常物理振幅 `cr=1`）。这些点仅最后一项有 Reference z10；其余与 Reference z8 比较。它仍是**私有 opt-in 实验路径**，没有改变正式默认或公开接口。接受为生产方案前还需：不同频率网格 Reference z10 留出、三个较新背景点的 z8→z10 收敛、更多 `cr>0` 参数点、phase-resolved 活跃/非活跃边界测试及多资源档位 A/B。此之前不得声明连续适用域或修改默认 `z_tail`。

复跑文件：[mask 屏查脚本](../../scripts/underflow_adaptive_handoff_mask_screen_20261115.py)、[倾斜 holdout 脚本](../../scripts/screen_adaptive_handoff_tilt_holdout_20261115.py)、[背景参数 holdout 包装脚本](../../scripts/screen_adaptive_handoff_background_holdout_20261115.py)、[`cr>0` 分支筛查](../../scripts/screen_adaptive_handoff_consistency_branch_20261115.py)、[per-mode 12 组 A/B](../../scripts/benchmark_adaptive_handoff_per_mode_20261115.py)。逐点结果和源码 SHA-256 见[mask JSON](results/underflow_adaptive_handoff_mask_screen_20261115.json)、[A/B JSON](results/adaptive_handoff_per_mode_benchmark_20261115.json)、[倾斜 JSON](results/adaptive_handoff_tilt_holdout_20261115.json)、[背景 holdout JSON](results/adaptive_handoff_background_holdout_20261115.json)、[`cr>0` JSON](results/adaptive_handoff_consistency_branch_20261115.json)、[preflight 拒绝 JSON](results/adaptive_handoff_physical_preflight_rejection_20261115.json)。

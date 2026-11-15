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

## 取舍

该结果支持继续做**每模 WKB 交接准则**实验，但不支持直接改 Fast：本轮的混合谱是把两个全量 z5/z8 运行输出拼接而成，**没有实现每模内核，也没有测出真实加速或内存成本**。不能用“只选中 5–9% 模式”直接推断相同比例的 runtime 节省，因为模式成本和尾部数值传播长度不相等。尚缺独立频率网格留出、三个后续点的 Reference z10 尾收敛，以及相位分辨尾部状态兼容性。生产 z_tail 未改变。

复跑文件：[mask 屏查脚本](../../scripts/underflow_adaptive_handoff_mask_screen_20261115.py)、[倾斜 holdout 脚本](../../scripts/screen_adaptive_handoff_tilt_holdout_20261115.py)、[背景参数 holdout 包装脚本](../../scripts/screen_adaptive_handoff_background_holdout_20261115.py)、[`cr>0` 分支筛查](../../scripts/screen_adaptive_handoff_consistency_branch_20261115.py)。逐点结果、源码 SHA-256 和物理拒绝记录见[mask JSON](results/underflow_adaptive_handoff_mask_screen_20261115.json)、[倾斜 JSON](results/adaptive_handoff_tilt_holdout_20261115.json)、[背景 holdout JSON](results/adaptive_handoff_background_holdout_20261115.json)、[`cr>0` JSON](results/adaptive_handoff_consistency_branch_20261115.json)、[preflight 拒绝 JSON](results/adaptive_handoff_physical_preflight_rejection_20261115.json)。

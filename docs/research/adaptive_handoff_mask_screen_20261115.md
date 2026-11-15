# 按今天的 `z` 延迟 WKB 交接筛查

日期：2026-11-15

## 候选准则

全局 `z_tail=8` 对本轮最差下溢频点很准，但慢约 1.57 倍。单模态状态显示，在 `z_today` 低于 5 的模式中，Fast z5 可能先越过 `z=5` 并做 WKB 终端外推，但该模态到今天又回到 `z_today≈4.99`；Fast z8 不再触发终端尾，继续数值传播。于是筛查按模式组合两个固定根结果：若 `z_today<5` 用 z8 的行，否则用 z5。这个判据针对已观察到的“交接后没有保持深亚视界至今天”情形；它还不是一条普遍 WKB 余项界。

## 离散结果

| 样本 | 网格 | z8 行比例 | 混合谱最大 / P95 差 | 混合映射相对 Reference 差 | 固定根残差 |
|---|---:|---:|---:|---:|---:|
| `n_t=17.38, r=5e-324, T_re=2000, kappa10=.01` | 937 | `4.70%` | `0.000517 / 0.000148 dex` | `+4.27e-7` vs z8 | `−2.303e-6` |
| `n_t=17.38, r=1e-316, T_re=2000, kappa10=.01` | 937 | `4.70%` | `0.000483 / 0.0000311 dex` | `−3.25e-6` vs z10 | `−3.25e-6` |
| `n_t=17.45, r=5e-324, T_re=10, kappa10=.001` | 221 | `4.98%` | `0.000578 / 0.000151 dex` | `+8.18e-7` vs z8 | `−2.22e-16` |

三点的混合最大谱差都低于 `1e-3 dex`，映射/根误差低于 `2e-4`。最小次正规 937 点的 Reference z10 对照已经直接跑完；另外两个候选点直接比较 Reference z8。Reference z8→z10 收敛没有在后两个点重新测，且三点都在 `cr=0` 自由倾斜支路。

另尝试 `n_t=17.45,T_re=1e4,kappa10=1,r=5e-324`，关闭物理 guard 的探索性运行给出 `DN_eff≈1.383e9`，超过物理 guard 上限 5。该组合被记为物理 preflight 拒绝，所有谱差结果作废，不计作候选通过。未通过绕过背景 guard 扩展适用域。

## 取舍

该结果支持继续做**每模 WKB 交接准则**实验，但不支持直接改 Fast：本轮的混合谱是把两个全量 z5/z8 运行输出拼接而成，**没有实现每模内核，也没有测出真实加速或内存成本**。不能用“只选中 5% 模式”直接推断 95% 的 runtime 节省，因为模式成本和尾部数值传播长度不相等。尚缺 `cr>0` 一致性分支、独立频率网格留出、不同背景的 Reference z10 尾收敛，以及相位分辨尾部状态兼容性。生产 z_tail 未改变。

复跑文件：[mask 屏查脚本](../../scripts/underflow_adaptive_handoff_mask_screen_20261115.py)、[倾斜 holdout 脚本](../../scripts/screen_adaptive_handoff_tilt_holdout_20261115.py)、[背景参数 holdout 包装脚本](../../scripts/screen_adaptive_handoff_background_holdout_20261115.py)。逐点结果、源码 SHA-256 和物理拒绝记录见[mask JSON](results/underflow_adaptive_handoff_mask_screen_20261115.json)、[倾斜 JSON](results/adaptive_handoff_tilt_holdout_20261115.json)、[背景 holdout JSON](results/adaptive_handoff_background_holdout_20261115.json)、[preflight 拒绝 JSON](results/adaptive_handoff_physical_preflight_rejection_20261115.json)。

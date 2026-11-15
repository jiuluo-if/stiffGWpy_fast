# 最小次正规振幅的 937 点尾深收敛

日期：2026-11-15

## 设计

对 `r=5e-324, A_s=2.100549e-9, n_t=17.38, cr=0`，使用其 235 点 Fast 根 `DN_eff=1.3543440396198841e-9` 固定背景，在相同 469/937 直接节点上比较 Fast z5、Fast z8、Reference z8 和 Reference z10。Reference z10 为同一 937 个频点重新直接求解，不由 937 点 z8 插值或尾部缩放近似得到。

## 结果

| 网格 | Fast z5 对 Ref z8 最大谱差 | Fast z8 对 Ref z8 最大 / P95 | Fast z8 对 Ref z10 最大 / P95 | Ref z8→z10 最大 / P95 |
|---:|---:|---:|---:|---:|
| 469 | `0.00055744 dex` | `0.00037591 / 0.00014796 dex` | 未运行 | 未运行 |
| 937 | `0.00244086 dex` | `0.00051668 / 0.00014801 dex` | `0.00037216 / 0.00002325 dex` | `0.00016489 / 0.00014205 dex` |

937 点 z5 超差位于新增中点 `log10(f/Hz)=−16.291944809`，Fast z8 对 Reference z10 的最大差转到 `−12.891944809`。Fast z8 937 映射相对冻结根残差 `−1.964e-6`；Reference z8、z10 映射残差分别为 `−2.731e-6`、`+9.503e-7`。根映射差仍远小于 `2e-4` 筛查门槛。

937 点 Fast 的本轮进程内预热后时间为 z5 `0.1256 s`、z8 `0.1974 s`，比值 `1.57`；这是单组对照，不是完整 P50/P95 性能认证。Reference z10 本次单次耗时 `1673.1 s`，z8 对照为 `243.6 s`，z10 约 6.9 倍；高精度尾部收敛成本显著。Reference z10 仍有 4 条 PCHIP divide warning，49 个归一化积分节点为零。quadrature 相对误差估计 `1.11e-14`，half-density interpolation 估计 `3.22e-7`。

## 判定和限制

此离散点确认 z5 在 937 网格不达 `1e-3 dex`，而 z8 对更深 Reference z10 达 `3.73e-4 dex`，通过本轮门槛。这里的误差变化与尾部渐近深度一致；网格、背景与振幅完全相同，没有以不同根混淆结论。它不证明 z8 在其它倾斜、再热温度、`kappa10` 或 `cr` 参数上通过，也不证明所有相位分辨分支的交接误差；不得据此改生产默认尾深。

复跑文件：[直接 469/937 Fast/Reference z8](../../scripts/underflow_min_subnormal_reference_grid_20261115.py)、[z8 Fast 对照](../../scripts/underflow_min_subnormal_tail_depth_20261115.py)、[937 Reference z10](../../scripts/underflow_min_subnormal_reference_z10_20261115.py)。完整逐点频谱、误差向量、警告、计时、参数和源码哈希分别见[Reference 网格 JSON](results/underflow_min_subnormal_reference_grid_20261115.json)、[Fast z8 JSON](results/underflow_min_subnormal_tail_depth_20261115.json)、[z10 JSON](results/underflow_min_subnormal_reference_z10_20261115.json)。

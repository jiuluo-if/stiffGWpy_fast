# CR0 2M–3M Fast 区域：WKB 基矩阵浮点拒绝门

## 数学修复

Fast 局部 WKB3 将区间首端基本矩阵 `B0` 反演后与末端矩阵相乘，得到原二维状态 `(x,y)` 的传递。旧核无条件用
`det(B0)=a11*a22-a12*a21` 作除数。在实际 3M 密度输入的一个高斜率面板上，WKB3 基底乘积都下溢为有符号零，浮点行列式为零；`solve_kernel` 的并行循环因此发生除零，末端有一批频点没有传播输出。

局部重参数化不变，仍是原来的 `z(N)` 线性段 WKB2/WKB3 映射。现在先算两个行列式乘积，再用二进制 64 位单位舍入 `u=2.220446049250313e-16` 得到保守分辨界

```text
E_det = 4u (|a11*a22| + |a12*a21|).
```

当行列式非有限、落入 subnormal 范围，或 `|det(B0)| <= E_det` 时，返回非有限的候选传递；现有 `wkb_transfer_embedded` 有限性门会拒绝候选，`_phase_segment_wkb3` 随即走原相位子步传播。此条件只依赖浮点可分辨性，与具体参数点无关，也不改物理方程、初态、匹配或再加热转折。

实际触发面板为 `z_start=8.326819262766648`、`z_end=13.006402132438978`、`h=5.853210263495612e-4`，`q≈7994.90`。二阶 WKB 基矩阵行列式仍可分辨（约 `-2.032`），三阶基矩阵的两个乘积都下溢成零。该面板因此回退，而不会终止整个并行频率循环。

方法依据：Higham, *Accuracy and Stability of Numerical Algorithms*, 2e，实际查阅 SIAM 版第 2 章 §2.1–2.2（浮点数与舍入模型）及第 3 章 §3.1（乘积求和误差模型）。本实现把两个乘积和一次减法的项级舍入界保守地放大到 `4u`；遇到 subnormal 时直接拒绝该逆矩阵，避免把正规数模型误用于下溢区。[第 2 章](https://epubs.siam.org/doi/10.1137/1.9780898718027.ch2) · [第 3 章](https://epubs.siam.org/doi/10.1137/1.9780898718027.ch3)。

## 区域证据

固定 `cr=0`、`Omega_ch2=0.1201075`、`H0=67.32117`、`DN_eff=0`、`A_s=2.100549e-9`、`r=.01`、`T_re=2000`、`DN_re=10`、`kappa10=.01`，独立扫描 `n_t`。在 `[2e6,3e6] × [-.5,0]` 的 3×3 矩形格点加两个内部点，旧 Fast `0/11`、新 Fast `11/11`；另外 8 点 scrambled Sobol 留出为 `0/8→8/8`。19 个旧失败全为 `nonfinite_background_primitive`。新解均通过物理 guard，H² 有限正值（样本最小 `0.4532139925`），`N_inf` 有限，75 个输出频点完整，无警告。

Reference 保持原生产源码和阈值。3M、`n_t=-.5` 与独立 Sobol 高密度点均在**相同 Fast 根**上比较频谱，并另做 Reference 根闭合：最大根相对误差 `1.7505e-5`，最大同根谱差 `7.7100e-4 dex`，最大闭合残差 `9.1044e-7`。复用的 3M、`n_t=0` 同根 Reference 为根误差 `3.4198e-6`、最大谱差 `7.5822e-4 dex`。全部通过 `2e-4` 与 `1e-3 dex` 门槛。参考频谱必须在 Fast 同一根比较；将 Fast 根处的 Fast 谱与闭合到稍有差异的 Reference 根谱比较会在窄频段造成约 `5.4e-3 dex` 的假差异，该混根比较不作为精度判据。

原始数据：[3×3 + 内部点](screen_cr0_density_2m_to3m_wkb_guard.json)、[8 点 Sobol](screen_cr0_density_2m_to3m_wkb_guard_sobol.json)、[同根 z8 Reference](validate_cr0_3m_wkb_guard_z8_multi.json)、[3M 零倾斜 Reference](validate_cr0_3m_wkb_guard_z8.json)。

## 性能与图

基线与候选除 `_wkb_transfer.py` 外完全相同，32 组交替完整 Fast 求解，affinity `[0,1]`、Numba 2/workqueue、BLAS 1、goal/PCHIP。可信 750k 与 1M 控制及 2.5M 新域点的 P50/P95 比分别为 `0.999/1.037`、`1.010/1.031`、`0.995/1.015`；根、频率数组和谱逐位相同。首个冷 JIT 的 750k 基线/候选求解为 `18.97/18.54 s`；预热后的 P50 最大回退 `1.01%`、P95 最大回退 `3.67%`。新增区域 18 个非冷启动 screen 样本的完整求解时间 P50/P95 为 `90.9/94.6 ms`。

覆盖、失败类别、同根精度和性能图：[cr0_density_2m_to3m_wkb_guard_domain.png](cr0_density_2m_to3m_wkb_guard_domain.png)。配对 A/B 原始逐次数据：[benchmark_wkb_basis_guard_ab.json](benchmark_wkb_basis_guard_ab.json)；完整来源与 artifact hash：[manifest](cr0_2m_to3m_wkb_basis_guard_manifest.json)。

## 证据边界

这组有限 lattice、内部点、Sobol 留出与三个 Reference 点支持连续矩形内的区域改善，但不构成对所有连续参数点的数学认证。`cr>0` 空间尚未覆盖。WKB 候选条件仍是局部舍入误差拒绝门，不是全局误差证明。

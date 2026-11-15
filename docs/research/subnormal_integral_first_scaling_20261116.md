# 下溢频谱的积分前尺度归一化

日期：2026-11-16

## 数学问题与修正

`cr=0` 时张量方程对正的原初功率尺度是线性的。设单位尺度传播给出的今日正信号为 `G_i`，频点功率的对数为 `p_i=log(P_t(f_i))`，频率求积使用固定节点和现有方法 `Q`。直接构造 `G_i exp(p_i)` 会在每个频点先舍入；多个小于最小次正规数的正节点即使各自变成零，其加权总积分仍可能可表示。

先取 `s=max_i(log(G_i)+p_i)`，构造 `q_i=exp(log(G_i)+p_i-s)`，按原来的节点和所选求积计算 `Q(q)`，最后一次恢复尺度：

`DN_gw = exp(s + log(ln(10) Q(q)) + log(Neff0 / Omega_nu))`。

对固定节点，Simpson、PCHIP、Gauss/Chebyshev 等求积对共同正比例缩放具有齐次性；PCHIP 斜率与区间积分也随正比例因子线性变化。因此该重排不改变实数算术下的节点插值或求积，只避免在积分前分别物化过小的物理频谱节点。输出频谱仍按原来的 float64 线性数组表示；对数频谱保留有限物理值。相位分辨尾部的 `blend_g` 使用它已有的 `log10_base_signal_today`；未被相位修正的行保留原始浮点信号，不做多余的 log/exp 往返。

## 最小正振幅边界分类

独立 Cartesian Reference 在同一 235 点频率网格、`cr=0,n_t=0` 下以单位振幅传播，并在对数域恢复物理尺度：

- `r=5e-324` 的连续映射根为 `0.237` 个最小正 float64 次正规单位，低于 `0.5` 的最近舍入边界，映射为零是正确结果，分类为 `dn_gw_integral_underflow`。
- `r=2*minsub` 仍低于积分可表示边界，也分类为 `dn_gw_integral_underflow`。
- `r=3*minsub` 的连续参考根约为 `0.711` 个单位，float64 映射虽舍入成一个单位，但最近舍入的根相对误差约 41%，不满足 `tol=1e-7`。Fast 拒绝并标记 `dn_gw_root_quantization`，不把量化结果当成精确闭合。
- `r=1e-315` 的根有约 4,800 万个次正规单位，舍入相对上界约 `1.04e-8`，可用 float64 表示。积分前归一化 Fast 根为 `2.3692005e-316`；原按节点恢复尺度的 Fast 根为 `2.3955548e-316`，相对偏高约 `1.08e-2`。参考缩放根为 `2.36921957e-316`（235 点 z10），差 `−8.06e-6`。
- 精确 `r=0` 仍走解析零振幅分支，成功输出根零；其 `log_A_t=-inf` 与非零振幅路径分开。

对外失败码 `tensor_amplitude_underflow` 保持兼容；新 `fast_failure_detail` 区分积分先下溢与根的 float64 量化精度不足。`S2` 原函数溢出仍是背景原函数阶段的 `nonfinite_background_primitive`/`S2`，不与上述两种情况混为一谈。

## 网格、尾部与固定根证据

`r=1e-315` 的直接 Fast 自洽根与单位振幅直接 Reference，在同一背景、各自网格上比较：

| 频率节点 | Fast z5 映射 | Reference z8 映射 | Reference z10 映射 | Fast z5 / Ref z10 映射误差 | Fast z5 vs Ref z10 最大谱差 | Ref 插值误差估计 | Ref z8→z10 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 235 | `2.36920050e-316` | `2.36990227e-316` | `2.36921957e-316` | `−8.06e-6` | 未测 | `9.344e-4` | `−2.881e-4` |
| 457 (`freq_res=2`) | `2.36928970e-316` | `2.36999020e-316` | `2.36930697e-316` | `−7.30e-6` | `4.672e-4 dex` | `1.935e-4` | `−2.883e-4` |

235→457 点 Fast 映射变化约 `3.76e-5`，Reference z10 约 `3.69e-5`。457 点 z12 direct Reference 映射为 `2.36929245e-316`；z10→z12 变化 `−6.14e-6`，z12 quadrature 误差估计 `1.11e-14`，插值误差估计 `1.928e-4`。457 点 Fast 与 z12 Reference 的固定点根相对差为 `−1.17e-6`，低于 `2e-4` 门槛。z8→z10 仍变化 `2.88e-4`，说明 z8 不足以作为尾部收敛结果；z10→z12 已闭合。

457 点逐频谱直接对照的最大误差为 `4.672e-4 dex`，P95 为 `3.073e-5 dex`，最坏节点 `log10(f/Hz)=−15.84708`。因此该点同时通过谱差、根差、网格差和 z10→z12 尾部差门槛。

另外，分别将 235/457 点 Fast 返回根设为新的初始 `DN_eff`，在原频率网格做独立 Fast 固定根复映射。两网格映射残差均为 0，频率哈希与原解相同；`Nv`、`sigma`、`f_hor` 在该根背景与 `DN_eff=0` 背景下逐位相同。这确认了这些次正规 `DN_eff` 在 float64 背景数组中完全被基准 `Neff` 吸收，故 z12 Reference 的 `DN_eff=0` 背景适用于该 fixed-root 对照。

逐点参数、完整计算环境、频率节点哈希和源码哈希保存在[振幅表示/积分下溢 JSON](results/dn_gw_integral_underflow_20261115.json)、[235/457 点 z8/z10 JSON](results/subnormal_map_grid_convergence_20261115.json)、[z12 continuation JSON](results/subnormal_map_tail_z12_20261115.json)、[固定根复映射 JSON](results/subnormal_map_fixed_root_20261116.json)及[457 点逐频谱对照 JSON](results/subnormal_spectrum_reference_z10_20261116.json)。可复跑脚本：[边界分析](../../scripts/analyze_dn_gw_integral_underflow_20261115.py)、[网格/尾部对照](../../scripts/screen_subnormal_map_grid_convergence_20261115.py)、[z12 continuation](../../scripts/screen_subnormal_map_tail_z12_20261115.py)、[固定根复映射](../../scripts/verify_subnormal_map_fixed_root_20261116.py)、[逐频谱 Reference 对照](../../scripts/screen_subnormal_spectrum_reference_z10_20261116.py)。

## 状态与限制

- 生产修改限于 `cr=0` 对数振幅路径的终端积分映射和失败详情；`cr>0` 背景依赖、频率生成、Reference、默认尾部深度和普通振幅路径没有更改。
- 目前的误差对照只有 `cr=0,n_t=0,r=1e-315,T_re=2000,kappa10=.01` 一个参数切片；不表示跨振幅、倾斜、再热温度或 `cr` 连续区域认证。
- 线性数组中的逐点物理谱可能仍舍入为零，只有有限对数谱和终端积分缩放得到保护。低于 float64 根相对精度能力的根会拒绝，而不会恢复成 log-root 公共结果。
- 本轮没有新的性能 A/B；不能据此宣称相对现有 Fast 有提速。

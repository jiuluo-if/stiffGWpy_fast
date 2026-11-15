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

## 独立蓝倾斜留出（2026-11-16）

对独立于 `n_t=0,r=1e-315` 样本的 `r=5e-324,n_t=17.38,cr=0,T_re=2000,kappa10=.01` 参数点，使用仓库原始 JSON 中直接嵌套的 235/469/937 节点和 937 点 Reference z10 完整逐频谱。Reference 文件的 `reference.py` SHA 与当前实现一致，频率哈希也与原 937 点网格完全相同；Reference 以原 235 根 `1.3543440396198841e-9` 固定背景直接求解，未用插值生成新频点。

- 当前 Fast 自洽根为 235/469/937 点 `1.3543440396198835e-9`、`1.3543410383880955e-9`、`1.3543409204417139e-9`；在每个返回根做直接固定根复映射时，相对残差均为 0。
- 固定 Reference 原根时，Fast 235→469、235→937 积分映射分别变化 `−2.2160e-6` 和 `−2.3031e-6`；937 点 Fast z5 映射相对直接 Reference z10 为 `−3.2534e-6`。Reference z10 映射对原 235 根的独立残差为 `+9.5033e-7`。
- 同根、同 937 点网格下，当前 Fast z8 对 Reference z10 的最大谱差为 `3.7216e-4 dex`，P95 为 `2.3253e-5 dex`，通过 `1e-3 dex` 门槛。Reference z8→z10 最大/P95 差为 `1.6489e-4/1.4205e-4 dex`，其 quadrature/interpolation 估计分别为 `1.11e-14/3.22e-7`。
- Reference 原始计算仍记录 4 条 `overflow encountered in divide` 警告和 49 个归一化零节点；本留出读取原始 Reference 数据，不屏蔽警告、不改实现。其直接 `Reference` 源码 SHA 已核验。

逐频谱、三网格根、固定根映射、Reference 原始数据哈希和当前 Fast 源码 SHA 见[原始留出 JSON](results/integral_first_blue_tilt_holdout_20261116.json)，可复跑[脚本](../../scripts/screen_integral_first_blue_tilt_holdout_20261116.py)。这将离散证据增加到两个不同的 `r,n_t` 点，但温度、`kappa10` 和 `cr` 尚未改变，不能认证任何连续参数域。

## 独立低温/高倾斜留出（2026-11-16）

继续使用原始直接 Reference z8/z10 结果，在 `r=5e-324,n_t=17.45,cr=0,T_re=10,kappa10=.001` 的 221 点网格上运行当前 Fast。Reference 源码 SHA 与当前文件一致，网格哈希与原 Reference 完全匹配；未重写 Reference，也没有复用旧 Fast 谱。

- 当前 Fast 根为 `2.832295918247082e-104`，相对 Reference 固定根 `2.832295918247089e-104` 差 `−2.44e-15`；在当前根的 Fast 固定背景复映射残差为 0。
- 固定 Reference 根下 Fast z5 映射与 Reference z10 相差 `−2.251e-6`；Fast z8 映射相差 `−1.338e-6`。
- 当前 Fast z8 对直接 Reference z10 的逐频谱最大差为 `5.6918e-4 dex`，P95 `2.0551e-5 dex`，满足 `1e-3 dex` 门槛。
- Reference z8→z10 最大/P95 谱变化为 `1.6483e-4/1.5985e-4 dex`，映射变化 `−1.2447e-5`；z10 插值/求积估计为 `2.544e-5/1.11e-14`，无警告、无归一化零节点。

逐频谱、根/映射、背景参数和源码哈希见[低温留出 JSON](results/integral_first_low_temperature_holdout_20261116.json)，可复跑[脚本](../../scripts/screen_integral_first_low_temperature_holdout_20261116.py)。现在有三个通过的离散 `cr=0` 点，覆盖 `n_t={0,17.38,17.45}`、`T_re={10,2000}` 与 `kappa10={.001,.01}`；它们共享最小正振幅或同一 `cr`，仍不能构成连续区域认证。

## `kappa10` 十倍变化留出（2026-11-16）

固定低温高倾斜背景 `r=5e-324,n_t=17.45,cr=0,T_re=10`，将 `kappa10` 从已通过的 `.001` 改为 `.01`。当前 Fast 预筛通过物理 guard，存储 `A_t=0` 但 `log_A_t` 有限；在当前 Fast 原生构造的 227 点网格上，直接 Reference z8/z10 均重新求解。没有用 221 点旧样本的谱插值代替新网格。

- Fast 自洽根为 `1.4942521645e-94`；同根 z5 映射残差为 0。Fast z8 对直接 Reference z10 的映射差为 `−3.294e-6`。
- 逐频谱最大差为 `3.1394e-4 dex`，P95 `2.2217e-5 dex`，最坏 `log10(f/Hz)=−12.51653`，通过 `1e-3 dex` 门槛。
- Reference z8→z10 最大/P95 谱变化为 `1.6498e-4/1.5849e-4 dex`，映射变化 `+1.2049e-5`；z10 插值/求积估计 `2.6163e-5/1.11e-14`，无警告且没有积分节点归一化为零。

逐频谱、根/映射和源码哈希见[`kappa10=.01` 留出 JSON](results/integral_first_kappa_decade_holdout_20261116.json)，可复跑[分阶段脚本](../../scripts/screen_integral_first_kappa_decade_holdout_20261116.py)。这一点将同一低温高倾斜切片上的 `kappa10` 扩展到十倍跨度，但没有测试 `cr>0`，不能认证连续的 kappa 区间。

## `cr>0` 一致性分支的当前边界（2026-11-16）

使用已有物理振幅一致性样本 `r=.01,cr=1,T_re=2000,DN_re=20,kappa10=.01` 作当前源码控制。Reference z10 源码哈希和 76 点频率网格与保存结果一致；当前 Fast 根仍为 `0.0022636474124195006`。固定根 z8 谱对 Reference z10 最大差 `3.5795e-4 dex`，根映射相对差 `−4.045e-6`。`A_t=2.100549e-11` 为正规数，直接监视 `_integrate_log_scaled_dn_map` 的调用次数为 0；该候选没有进入物理一致性分支。

另用当前源码重跑 `r=1e-316,cr=1` 的九点 `T_re={10,2000,10000} × kappa10={1e-3,.01,1}` 矩阵。所有点的 `Delta_N<0`、`N_inf=None`，Fast 均在构造传播原语之前按 `invalid_cutoff` 拒绝。独立正常振幅压力控制仍到达 `primitive:S2`，并报 `nonfinite_background_primitive`；没有把 cutoff 拒绝和 `S2` 溢出合并。矩阵当前源码结果见[JSON](results/consistency_underflow_boundary_20261115.json)，物理振幅控制见[JSON](results/integral_first_cr_positive_control_20261116.json)和[复跑脚本](../../scripts/verify_integral_first_cr_positive_control_20261116.py)。

这仍未找到 `cr>0`、cutoff 有效且原初振幅落在次正规边界的可传播点；一致性分支适用域没有扩展。

- 生产修改限于 `cr=0` 对数振幅路径的终端积分映射和失败详情；`cr>0` 背景依赖、频率生成、Reference、默认尾部深度和普通振幅路径没有更改。
- 目前的误差对照只有 `cr=0` 下四个离散点，不表示跨振幅、倾斜、再热温度、`kappa10` 或 `cr` 连续区域认证。
- 线性数组中的逐点物理谱可能仍舍入为零，只有有限对数谱和终端积分缩放得到保护。低于 float64 根相对精度能力的根会拒绝，而不会恢复成 log-root 公共结果。
- 本轮没有新的性能 A/B；不能据此宣称相对现有 Fast 有提速。

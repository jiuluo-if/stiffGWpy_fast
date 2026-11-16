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

## `cr>0` cutoff 与 `A_t` 可表示边界的解析关系

在固定 `T_re,kappa10,DN_eff=0` 的 `cr>0` 支路，`N_re`、`rho_rad,re` 与 `rho_stiff,re` 不依赖 `r`。一致性方程可写为

`Delta_N = 4/3*(N_10-N_re) + 4/3*ln(M_Pl) + 1/3*(ln(45/2)+ln(r)+ln(A_s)-ln(rho_rad,re+rho_stiff,re))`。

因此 `Delta_N` 对 `ln(r)` 的斜率是 `1/3`，cutoff 零点有闭式表达：

`ln(r_cut) = -4*(N_10-N_re) - 4*ln(M_Pl) - ln(45/2) - ln(A_s) + ln(rho_rad,re+rho_stiff,re)`。

我用该式给出候选边界，再在 `r_cut*(1±1e-6)` 直接重建模型验证 `Delta_N` 符号。对 `T_re={10,2000,10000}`、`kappa10={1e-3,.01,1}` 九个组合，边界 `log10(r_cut)` 分别落在 `−57.47..−36.25`；float64 `A_t` 正规/次正规边界是 `log10(r)=-298.975`。所以 cutoff 有效下界比 `A_t` 开始次正规高至少 `241.5` 个 `r` 数量级；在 `A_t` 刚进入次正规的样本中，`Delta_N` 仍为 `−185..−202`，`N_inf=None`。在这些离散背景上，不存在同时“cutoff 有效”和“张量振幅次正规/舍入为零”的交集，不能通过单位振幅传播绕开。

完整九点解析根、两侧符号、正规边界和舍入零控制见[`cr>0` cutoff/振幅边界 JSON](results/cr_positive_cutoff_amplitude_boundary_20261116.json)，复跑[脚本](../../scripts/map_cr_positive_cutoff_amplitude_boundary_20261116.py)。这只证明指定九个 `T_re,kappa10` 点，不外推到所有宇宙学参数。

### cutoff 有效侧的物理振幅验证

选择中位背景 `T_re=2000,kappa10=.01`，取 `r=1.1*r_cut=3.962216254776857e-43`。在 `DN_eff=0` 起始背景上 `N_inf-N_re=0.03177006>0`；当前 Fast 的自洽根是 `1.067074488428421e-43`，固定该根时 `A_t=8.322829391755272e-52`，仍为正规数，`N_inf=37.80411134404018`、`N_re=37.7723412841054`。传播没有绕过 `cr>0` 背景依赖。

当前 Fast 在新生成的 77 点频率网格上做直接自洽；之后在同一根和同一网格做固定根 Fast 与 Reference 对照。Reference `reference.py` 源哈希写入 JSON；Cartesian z8/z10/z12 均为直接求解，未插值生成频点、无警告、无非正/非有限原始谱节点。Fast 固定根积分为 `1.0670751194058553e-43`，相对 Reference z12 映射 `−6.086e-6`；Fast 自洽根相对 Reference z12 映射约 `−6.68e-6`。Reference z10→z12 积分变化 `−1.953e-5`，谱最大变化 `3.647e-4 dex`；Fast z8 对 Reference z12 原始正值谱的最大/P95 差为 `3.776e-4/2.567e-4 dex`，满足本轮 `2e-4` 根参考和 `1e-3 dex` 谱门槛。z8→z10 映射变化为 `−2.733e-4`，但更深的 z10→z12 已闭合到上述值。

需要特别说明 Reference 对外 `log10OmegaGW` 在本点的 77/77 个节点都触及公开下限 `−40`，不能拿它与 Fast 直接比。本对照使用直接 Reference 结果中的 `Ogw-Oj` 原始正值取 `log10`；所有原始值经脚本检查为正且有限，并将被公开下限截断的节点数记录在 JSON。逐点频谱、原始对数谱、网格与所有源码哈希见[cutoff 有效侧 JSON](results/cr_positive_cutoff_valid_sample_20261116.json)，可复跑[直接对照脚本](../../scripts/verify_cr_positive_cutoff_valid_sample_20261116.py)。

这是一个物理振幅、cutoff 有效的 `cr>0` 离散样本，验证的是该点 Fast/Reference 一致性；它不扩展 `cr>0` 次正规振幅的传播域，也不认证 `T_re`、`kappa10` 或其它参数的连续范围。Reference z12 单点运行约 1291 秒，暂无相同样本的性能 A/B 结论。

### 低温背景独立留出

再取边界矩阵低温端 `T_re=10,kappa10=.001`，同样选 `r=1.1*r_cut=3.7536771286729356e-58`，起始 `Delta_N=.03177006>0`。当前 Fast 自洽根为 `3.287879454649206e-64`；固定根 `A_t` 为正规数且 `N_inf=32.42228249191523>N_re`。新生成的 76 点 Fast 网格上直接计算 Reference z8/z10/z12：Fast 固定根积分相对 z12 Reference 为 `+2.571e-6`，Fast 自洽根相对 Reference z12 映射约 `−2.39e-7`；Fast z8 对 Reference z12 原始正值谱最大/P95 差为 `5.693e-4/4.881e-4 dex`。Reference z10→z12 映射变化 `+3.389e-6`、逐点谱最大变化 `5.410e-4 dex`，无警告。Reference z12 耗时约 1494 秒。

与中温样本相同，Reference 对外 `log10OmegaGW` 的 76 个节点全部被 `−40` 下限截断，因而使用直接求解的 `Ogw-Oj` 正有限值取对数；计数、原始谱和源码哈希见[低温 cutoff 有效侧 JSON](results/cr_positive_cutoff_valid_lowT_20261116.json)。同一[复跑脚本](../../scripts/verify_cr_positive_cutoff_valid_sample_20261116.py)通过 `--t-re`、`--kappa10` 和 `--output` 选择这组参数与产物路径。

现有两个 cutoff 有效侧物理振幅点覆盖了低温 (`10,.001`) 和中温 (`2000,.01`) 两组离散背景，均通过根/谱/尾部门槛；这仍不是连续域认证，也不提供 `cr>0` 次正规振幅下的传播证明。z12 参考成本分别约 1494 秒和 1291 秒；没有进行 Fast 性能 A/B。

### 高温/高 `kappa10` 角点留出

第三个对角角点取 `T_re=10000,kappa10=1`，`r=1.1*r_cut=6.209353644613012e-37`。初始 `Delta_N=.03177006>0`，自洽 Fast 根为 `4.186526420733873e-34`。在直接新生成的 77 点 Fast 网格上，Fast 根对 Reference z12 映射差 `−6.074e-6`，固定根 Fast 积分对 Reference z12 差 `−5.402e-6`；Fast z8 对 Reference z12 原始正值谱最大/P95 差 `8.357e-5/2.859e-5 dex`。Reference z10→z12 映射变化 `−3.942e-5`、最大逐点谱变化 `8.712e-5 dex`，无警告；z12 耗时约 1198 秒。

z8→z10 映射变化为 `−2.528e-4`，略高于 `2e-4` 根筛查参考线；故保留这项不收敛记录，并以更深的 z10→z12 变化判断尾部闭合，未用 z8 结果替代 z12。Reference 对外谱的 59/77 节点触及 `−40` 下限，其余节点也统一以原始 `Ogw-Oj` 对数谱做比较。完整数据见[高温角点 JSON](results/cr_positive_cutoff_valid_highT_20261116.json)。

当前三个有效侧验证点是 `(T_re,kappa10)=(10,.001),(2000,.01),(10000,1)`，即低/中/高的对角角点，而不是正交参数扫描；温度与 `kappa10` 同时改变，不能据此区分各自影响或认证连续域。三点均是物理正常振幅，仍没有找到 cutoff 有效的 `cr>0` 次正规振幅传播样本。没有生产代码改变或性能 A/B 结果。

### 固定 `kappa10` 的温度对照

为拆分温度和 `kappa10` 的共同变化，再固定 `kappa10=.01`，取低温 `T_re=10`，参数为 `r=1.1*r_cut=3.7512651663963e-57`。初始 cutoff 裕量仍为 `Delta_N=.03177006`。Fast 自洽根为 `2.2509204142520693e-62`，对 Reference z12 映射差 `+6.326e-6`；固定根 Fast 映射差 `−1.235e-6`。Fast z8 对 Reference z12 原始正值谱最大/P95 差 `3.752e-4/2.980e-4 dex`。Reference z10→z12 映射变化 `+4.972e-6`、最大逐点谱变化 `3.741e-4 dex`，无警告；z12 耗时约 1467 秒。

该 76 点网格的 Reference 对外谱仍全部落在 `−40` 下限，使用直接求得的 `Ogw-Oj` 比较；原始信号为正且有限。产物及哈希见[固定 `kappa10` 温度对照 JSON](results/cr_positive_cutoff_valid_lowT_fixed_kappa_20261116.json)，运行命令可由参数化脚本 `--t-re 10 --kappa10 0.01` 重现。

目前四个通过的离散 cutoff 有效侧样本是三个低/中/高对角角点及一个固定 `kappa10=.01` 的低温对照。它们支持这些点上的 Fast/Reference 一致性，不代表 `T_re` 或 `kappa10` 连续范围认证，也没有改变 `cr>0` 次正规振幅的无效 cutoff 结论。

固定 `kappa10=.01` 的高温端取 `T_re=10000,r=1.1*r_cut=6.209353644953826e-39`，初始 `Delta_N=.03177006>0`。Fast 自洽根对 Reference z12 映射差 `−9.856e-6`，固定根 Fast map 差 `−9.103e-6`；Fast z8 对 Reference z12 原始正值谱最大/P95 差为 `3.845e-4/2.394e-4 dex`。Reference z10→z12 映射变化 `−3.432e-5`、谱最大变化 `3.549e-4 dex`，无警告，z12 耗时约 1243 秒。z8→z10 映射变化 `−2.542e-4` 超过筛查参考线，已保留；更深 z10→z12 对照满足收敛判断。公开 Reference 谱的 77 个节点全部触及 `−40` 下限，比较使用经脚本验证为正有限的原始 `Ogw-Oj`。完整数据见[高温固定 kappa JSON](results/cr_positive_cutoff_valid_highT_fixed_kappa_20261116.json)。

固定 `kappa10=.01` 的 `T_re={10,2000,10000}` 三个离散温度现在都有直接 Reference z12 和固定根谱对照；中温控制与此前中温 cutoff 样本相同。三点根和频谱比较通过，但不认证中间温度的连续区间。总计五个不同的正常振幅 `cr>0` cutoff 有效侧点通过，覆盖三个低/中/高对角角点和两个固定 `kappa10` 温度端点；没有生产代码改变、性能 A/B 或 `cr>0` 次正规振幅传播证据。

- 生产修改包括 `cr=0` 对数振幅路径的终端积分映射/失败详情，以及按模态 WKB 延迟交接；Reference、频率生成规则和频率积分没有更改。
- `cr>0` 一致性支路会用张量一致性关系派生实际 `nt`；传入的 `cosmo_param['n_t']` 不控制该支路的物理谱。`n_t=17.3` 与 `17.45` 输入对照都得到实际 `nt=-4.952770318471071e-44`、相同频谱和网格，因此不能算高倾斜参数留出。
- 离散样本覆盖若干振幅、倾斜、再热温度、`kappa10` 与 `cr` 组合，不表示这些变量的连续区域认证。
- 线性数组中的逐点物理谱可能仍舍入为零，只有有限对数谱和终端积分缩放得到保护。低于 float64 根相对精度能力的根会拒绝，而不会恢复成 log-root 公共结果。
- 按模态交接的 A/B 只对已记录的同根、同网格、单线程参数点成立；不能宣称广域 Fast 提速。

## 解析纹波阈值的逐模 WKB 交接（2026-11-16）

### 根因与阈值推导

在有效 cutoff `cr>0` 物理振幅点，当前默认 z5 Fast 与固定根直接 Reference z12 同网格比较时，最大误差 `1.886e-3 dex` 出现在 `log10(f/Hz)=-16.0917663`。该模式 z5 `handoff_eps=3.3378e-3`，说明在 `z=5` 触发了 WKB 终端交接；将同一模式改为 Fast z8/z10 后 `handoff_eps=-1`，即不触发解析尾部、继续数值传播到今天，该节点误差降至 `3.35e-6 dex`。误差主要来自浅层 WKB 交接的逐频谱偏差，根映射变化仅约 `6.7e-6`。

用当前代码已存在的领先端点相位纹波界

`b(z) = log10(1 + exp(-z) + exp(-2z))`

将允许误差设为 `1e-3 dex`。令 `A=10^(1e-3)-1`、`x=exp(-z)`，满足 `x^2+x=A` 的正根为

`x = 2A / (1 + sqrt(1 + 4A))`，`z_cut=-ln(x)=6.074868627758201`。

因此仅当模态在当前 `z_tail` 会触发交接且 `z_today<z_cut` 时，才把该模态的交接推迟到 `z_tail+3`。低于当前 `z_tail` 的模态原本就会传播到今天，不增加无用逐模数组；超过 `z_cut` 的模态继续走原 z5 交接。相位分辨权重大于零的模态不延迟，继续用原相位尾部修正。该纹波界约束领先相位纹波，不是一般 WKB 高阶余项的严格上界。

### 当前 Fast 的离散验证

固定根、网格、单线程 A/B：

| 物理点 | z5 无延迟最大谱差 | 自适应最大谱差 | z5 P50/P95 | 自适应 P50/P95 | 自适应选中频点 |
|---|---:|---:|---:|---:|---:|
| `cr>0`, `T_re=2000`, `kappa10=.01`, `r=1.1*r_cut`，实际 `nt≈-5e-44` | `1.886e-3 dex` | `5.639e-4 dex` | `4.00/5.51 ms` | `4.45/6.41 ms` | 1/77 |
| `cr>0`, `r=.01`, `T_re=2000`, `kappa10=.01`，实际 `nt=-.00125` | `1.429e-3 dex` | `3.741e-4 dex` | `6.42/6.79 ms` | `6.87/9.17 ms` | 1/76 |
| `cr=0`, `r=5e-324`, `nt=17.45`, `T_re=10`, `kappa10=.001` | `5.695e-4 dex` | `5.695e-4 dex` | `29.42/32.46 ms` | `29.24/34.65 ms` | 0/221 |

三个点自适应 Fast 最大谱误差均低于 `1e-3 dex`，映射/根误差低于 `2e-4`。参考侧分别使用直接 Cartesian Reference z12（前两个点）和源哈希/网格哈希匹配的直接 Reference z10（次正规 `cr=0` 点）。其 z8→z10 变化为 `1.648e-4 dex`，映射变化 `-1.245e-5`。次正规点相位活跃模式没有被延迟；自适应路径和 z5 数值结果逐点相同。

同条件下自适应相对全局 z8 的 P50 加速分别为 `2.93×`、`2.35×`、`1.76×`；相对 z5 的 P50 时间比分别为 `1.11×`、`1.07×`、`0.994×`。自适应 P95 分别为 `6.41 ms`、`9.17 ms`、`34.65 ms`，RSS P50/P95 与对照相近。高精度和省下的全局 z8 成本来自逐模推迟，但该结果限于三个离散点、单线程，尚非跨参数性能保证。

按模态开关之前是私有原型；本轮把由 `1e-3 dex` 纹波界计算的默认阈值接入 Fast 正式路径。原型复用按相位权重保护活动相位行的逻辑。全量 `pytest` 为 `411 passed, 9 deselected`；Ruff、Mypy 和 sdist/wheel 构建通过。产物与复跑脚本：[`cr>0` cutoff 点 Reference](results/cr_positive_consistency_cutoff_20261116.json)、[当前默认 Fast cutoff 复核](results/cr_positive_consistency_adaptive_cutoff_20261116.json)、[`cr>0` 自适应 A/B](results/cr_positive_consistency_tail_handoff_20261116.json)、[常规振幅控制 A/B](results/cr_positive_normal_control_handoff_20261116.json)、[`cr=0` 次正规控制 A/B](results/cr0_underflow_tail_handoff_20261116.json)、[cutoff 预筛脚本](../../scripts/screen_cr_positive_consistency_cutoff_20261116.py)、[当前 Fast cutoff 复核脚本](../../scripts/validate_cr_positive_consistency_adaptive_cutoff_20261116.py)、[逐模 A/B 脚本](../../scripts/diagnose_cr_positive_consistency_tail_handoff_20261116.py)。

当前候选通过这三个离散点，但仍须在不同 `T_re`、`kappa10`、普通 `cr=0` 振幅与更高频率网格上继续独立验证；它没有形成一般 WKB 余项证明，也不能认证连续参数域。

### 937 点 `cr=0` 最小次正规留出与非单调背景修正

在 canonical nested 937 点网格上，以 `r=5e-324, n_t=17.38, T_re=2000, kappa10=.01, cr=0` 的直接 Reference z10 谱和冻结根作同根、同网格、单线程 A/B。该网格逐点直接求解，没有把 235/469 点谱插值到 937 点。z5 最大谱差为 `2.4408647e-3 dex`，z8 为 `3.7215932e-4 dex`；默认自适应方案延迟 44/937 个模式，最大差为 `4.8550909e-4 dex`，满足本轮 `1e-3 dex` 频谱门槛。Reference z8→z10 映射变化 `3.6809e-6`、最大谱变化 `1.6489e-4 dex`；Fast 自洽根相对冻结 z10 Reference 根差 `−2.3031e-6`，满足本轮 `2e-4` 根筛查线。

该留出暴露了首版模式筛选条件的问题：最差 z5 模式今天的 `z_today=4.99072<z_tail=5`，但在非单调背景中曾短暂越过 `z_tail`，所以当前确实会触发 z5 交接；仅筛选 `z_today>z_tail` 会漏掉它。最终规则对所有 `z_today<z_cut` 且相位权重为零的模式延迟到 `z_tail+3`，相位分辨活动模式仍走原路径。逐模阈值变体当前只在 `z_today` 纹波超过解析门槛时应用；这个修正由 937 点直接谱对照支持，不代表 WKB 高阶余项存在统一严格上界。

同一固定根/网格/线程，预热后 12 次交错运行：z5、自适应、z8 的 P50 分别为 `125.39、126.22、199.15 ms`，P95 为 `127.89、132.32、219.92 ms`。自适应较全局 z8 快 `1.578×`，但相对 z5 P50 慢约 `0.67%`；因此这里显示的是精度修复以极小运行成本维持 z5 级速度，不是相对原 z5 的加速。RSS P50 分别 `227.96、229.46、229.48 MiB`。冷启动首调含 JIT，单列于 JSON，不与预热数据混合。

937 点产物[`A/B JSON`](results/cr0_underflow_937_handoff_20261116.json)记录完整逐节点谱、根残差、重复计时、RSS、环境、输入和源码 SHA；[复跑脚本](../../scripts/diagnose_cr0_underflow_937_handoff_20261116.py)会核验 nested grid 与 Reference 源哈希。Reference 的 PCHIP 仍发出 4 次 `overflow encountered in divide`，49 个缩放节点舍入为零；本次直接 Reference 谱完整且有限、网格/尾部收敛筛查满足门槛，但这些警告和 Reference 插值稳定性问题尚未修复，Reference 仍需专门审计。该结论只覆盖一个高倾斜、最小次正规振幅的 `cr=0` 离散点，不认证 469/937 的一般网格收敛或连续参数范围。

### 937 点 Reference PCHIP 警告诊断

复用同一直接 Reference z10 的完整有限 `log10(OmegaGW)` 节点和 canonical 937 点网格，对全局归一化后的线性数组重跑原始 `integrate_spectrum`，捕获到 4 条 `overflow encountered in divide`；937 点中 49 个节点在全局归一化后成为精确零，另有 10 个为正次正规数。节点动态范围为 `362.22 dex`。这与 PCHIP 内部加权调和割线斜率的倒数运算相符；warning 不是 QUADPACK 求积产生的。

以现有 `LogPchip` 从同一 log 节点计算**同一个正值线性 PCHIP Hermite 多项式**的割线、导数和积分，未对 `log(y)` 构造 PCHIP。该计算无 warning，总积分与原始 SciPy PCHIP 相差 `1.42e-14` 相对。49 个触及归一化零节点的插值区间，稳定候选累计积分仅占总量 `10^-322.716`，故这些零值在此离散点对 `DN_gw` 积分可忽略。对 937 个新增区间中点，原 SciPy PCHIP 有 48 个也舍入为零；在仍可表示为正的中点，候选与 SciPy 输出 P95 差约 `2.84e-14 dex`，最大局部差 `0.3009 dex` 出现在极弱尾部 `log10(f/Hz)=-15.51694`。因此总积分稳定不等于每个极小尾点插值稳定；若该频带之后被加权、扩展或用于其它可观测量，仍需独立评估。

数据、源码哈希和可复跑脚本见[PCHIP 937 点诊断 JSON](results/reference_pchip_937_diagnostic_20261116.json)与[脚本](../../scripts/diagnose_reference_pchip_937_20261116.py)。这是对警告及其积分影响的诊断，不改 Reference 正式实现，也不把普通 log-PCHIP 当替代方案。该稳定算法先前只在 235 点振幅阶梯上评估；此次 937 点对照进一步支持总积分数学等价性，但尚未在不同背景、独立频率网格与更广动态范围重复，因此不能据此替换独立 Reference 基准。

### 跨网格和异背景 PCHIP 复核

继续用同一套诊断代码复算直接 Reference 数据：相同 `cr=0,r=minsub,T_re=2000,kappa10=.01,n_t=17.38` 背景上的 469 点 z8、937 点 z8/z10，以及异背景 `T_re=10,kappa10=.001,n_t=17.45` 的 221 点 z10。所有网格哈希及 `reference.py` SHA 都与归档来源匹配；不同频率网格的谱均是各自直接求解结果，不做插值拼接。

| Reference 点 | 对数节点动态范围 | 全局缩放零/正次正规节点 | 原 PCHIP warning | 同多项式稳定积分相对差 | 原中点为零 | 可表示中点差 P95 / 最大值 |
|---|---:|---:|---:|---:|---:|---:|
| 469 点，z8 | `362.216 dex` | `25 / 5` | 4 | 0 | 24 | `2.84e-14 / 0.2713 dex` |
| 937 点，z8 | `362.216 dex` | `49 / 10` | 4 | `1.42e-14` | 48 | `4.69e-14 / 0.3009 dex` |
| 937 点，z10 | `362.216 dex` | `49 / 10` | 4 | `1.42e-14` | 48 | `2.84e-14 / 0.3009 dex` |
| 221 点，低温 z10 | `267.338 dex` | `0 / 0` | 0 | 0 | 0 | `5.68e-14 / 1.14e-13 dex` |

在发生下溢的三个 2000 K 数据中，稳定同多项式法测得零节点相邻区间的总积分分别只占 `10^-321.72`、`10^-322.72`、`10^-322.72`；这解释了本样本总积分对下溢零节点不敏感。PCHIP warning 与尾节点规模同步出现，在没有缩放零/次正规点的 221 点低温背景消失。局部中点对照则显示另外一面：原始 SciPy 在两个下溢网格都产生零值，有限但极小的尾部节点处仍有最高 `0.27–0.30 dex` 的偏差，尽管这些区间对当前积分可忽略。不能把积分一致当成点值插值一致。

跨网格原始结果、参数、来源哈希和运行环境见[多网格 PCHIP JSON](results/reference_pchip_multigrid_diagnostic_20261116.json)，可由[多网格诊断脚本](../../scripts/diagnose_reference_pchip_multigrid_20261116.py)复跑。当前证据支持把原 Reference warning 分类为极弱正值割线上的浮点条件问题，且对这四个积分没有可见影响；它不证明其它背景、加权量或局部弱尾观测量安全。Reference 保持不变；下一步需检查更广正值动态范围、独立新节点及其权重，方可判断是否将稳定同多项式积分作为可选诊断分支。

## 当前自适应 WKB 交接在 469 点独立网格的验证

在已有直接 Cartesian Reference z8 的精确 469 点高倾斜、最小次正规 `cr=0` 网格上，以 Reference 冻结根、相同频率节点、单线程，对照当前 Fast z5、z8 与默认自适应；另外直接运行当前 Fast 的 469 点自洽根。直接 Reference 的 `reference.py` SHA 与当前源码一致，四条 PCHIP warning 和 25 个归一化零节点继续原样记录。

| Fast 分支 | 相对直接 Ref z8 最大 / P95 谱差 | P50 / P95 | RSS P50 / P95 | 延迟模式 |
|---|---:|---:|---:|---:|
| z5 | `5.5744e-4 / 1.4844e-4 dex` | `61.60 / 63.35 ms` | `178.59 / 189.52 MiB` | 0 |
| 自适应 | `4.7462e-4 / 1.4838e-4 dex` | `62.42 / 63.68 ms` | `179.28 / 188.08 MiB` | 22/469 |
| z8 | `3.7591e-4 / 1.4796e-4 dex` | `98.89 / 100.58 ms` | `178.57 / 187.19 MiB` | 0 |

三种分支重复运行逐点确定；三者最大谱差均低于 `1e-3 dex`。本点 z5 已满足门槛，自适应只把误差再降低约 15%，且 P50 比 z5 慢 `1.33%`，所以没有相对 z5 的性能收益；相对全局 z8 的 P50 加速为 `1.58×`。当前 Fast 自洽根为 `1.3543410383880955e-9`，相对 Reference 冻结根差 `−2.2160e-6`，在同一根下 Fast 映射残差为 0；自洽根相对 Reference z8 映射差 `+2.3209e-7`。这些结果支持该离散点上精度与根闭合，不证明更广参数域，也不意味着每个已过 z5 门槛的点都需要自适应延迟。

完整 A/B 节点谱、计时样本、RSS、交接数量、环境和源码/输入哈希见[469 点 A/B JSON](results/cr0_underflow_469_handoff_20261116.json)，可用[复跑脚本](../../scripts/diagnose_cr0_underflow_469_handoff_20261116.py)重现。

## `cr>0` 一致性分支高温交接留出

选用已通过 cutoff 检查的 `T_re=10000,kappa10=.01,r=1.1*r_cut=6.209353644953826e-39,cr=1` 点，沿用同一 77 点直接 Reference z12 网格和 Reference 冻结根，对当前 Fast z5、自适应、z8 作同根/网格/单线程 A/B；同时在该网格运行当前 Fast 自洽根。初始 `Delta_N=0.03177006>0`，Fast 根上 `A_t=1.3043e-47` 为正规数，未绕过一致性背景。输入 `n_t=0` 不控制该分支的物理倾斜；当前模型派生 `n_t=−7.7617e-40`，所以此点不是高倾斜留出。

Reference z8→z10 积分变化为 `−2.5422e-4`，超过筛查参考线；z10→z12 收敛到 `−3.4323e-5`，最大逐点谱变化 `3.5486e-4 dex`。以直接 z12 为谱和根基准，不用 z8/z10 代替更深参考。Reference 对外谱 77 个节点全部触及 `−40` 下限，因此谱比较使用直接 Reference 求解保存的原始正值对数谱；Reference 无 warning。

| Fast 分支 | 相对 Ref z12 最大 / P95 谱差 | P50 / P95 | RSS P50 / P95 | 延迟模式 |
|---|---:|---:|---:|---:|
| z5 | `1.3817e-3 / 3.7713e-4 dex` | `4.15 / 5.64 ms` | `175.34 / 188.06 MiB` | 0 |
| 自适应 | `5.4243e-4 / 2.7457e-4 dex` | `4.50 / 6.15 ms` | `174.10 / 186.84 MiB` | 9/77 |
| z8 | `3.8447e-4 / 2.3939e-4 dex` | `12.69 / 20.15 ms` | `174.80 / 187.24 MiB` | 0 |

当前 Fast 自洽根 `4.184938229560105e-38` 相对 Reference z12 映射差 `−9.8558e-6`，当前根映射自身残差为 0。自适应把 z5 超门槛误差降至 `5.42e-4 dex`；P50 比 z5 慢 `8.5%`，比全局 z8 快 `2.82×`。这是一个 cutoff 有效高温 `cr>0` 离散点，支持该交接策略跨过 `cr=0/cr>0` 分支边界的可用性，但不认证连续温度、`kappa10` 或振幅范围。

完整频谱、逐次计时、RSS、派生倾斜、根、Reference 尾部收敛、环境及源码哈希见[`cr>0` 高温交接 JSON](results/cr_positive_highT_handoff_20261116.json)，可由[复跑脚本](../../scripts/diagnose_cr_positive_highT_handoff_20261116.py)重现。

随后使用同一脚本在固定 `kappa10=.01` 的低温 cutoff 点 `T_re=10,r=1.1*r_cut=3.7512651663963e-57` 复跑，构成低/高温两个离散 `cr>0` 对照。此点初始 `Delta_N=.03177006>0`，输入 `n_t=0` 派生为 `−4.6891e-58`，正规 `A_t=7.8797e-66`；Fast 自洽根对 Reference z12 差 `+5.8515e-6`。Reference z8→z10 根差 `−2.9090e-5`，z10→z12 根/谱变化 `+4.9719e-6 / 3.7411e-4 dex`，无 warning；76 个公开 Reference 对数谱节点都在 `−40` floor，比较用直接求得的原始正值 log 谱。

| Fast 分支 | 相对 Ref z12 最大 / P95 谱差 | P50 / P95 | RSS P50 / P95 | 延迟模式 |
|---|---:|---:|---:|---:|
| z5 | `2.5217e-3 / 3.7191e-4 dex` | `4.14 / 5.93 ms` | `172.36 / 182.56 MiB` | 0 |
| 自适应 | `4.3519e-4 / 3.3615e-4 dex` | `4.69 / 6.11 ms` | `171.46 / 181.60 MiB` | 12/76 |
| z8 | `3.7522e-4 / 2.9797e-4 dex` | `13.99 / 22.31 ms` | `171.94 / 182.04 MiB` | 0 |

因此两个固定 `kappa10=.01` 温度端点的 z5 都超 `1e-3 dex`，默认自适应均通过；低温点 P50 比 z5 慢 `13.2%`、较 z8 快 `2.98×`。结合上一高温点，证据支持这两个离散 `cr>0` cutoff 有效样本，不认证中间温度连续段。低温逐点数据见[低温交接 A/B JSON](results/cr_positive_lowT_handoff_20261116.json)；高温和低温均由[同一参数化复跑脚本](../../scripts/diagnose_cr_positive_highT_handoff_20261116.py)直接从各自的 cutoff/Reference 原始 JSON 读取。

### 固定低温的 `kappa10` 十倍边界对照

在相同 `T_re=10,cr=1` 下将 `kappa10` 从 `.01` 降至 `.001`，两点分别取自身解析 cutoff 的 `1.1` 倍振幅，保持 cutoff 有效。注意 `r` 同时从 `3.7513e-57` 改到 `3.7537e-58`；因此这是沿 cutoff 边界跨一个 `kappa10` decade 的联合样本，不是固定振幅的 kappa 单变量隔离。两点输入 `n_t=0` 都由一致性支路派生出接近零的负倾斜，张量 `A_t` 为正规数。

`kappa10=.001` 点使用直接 76 节点 Reference z12；z8→z10 与 z10→z12 映射变化分别 `+2.5594e-6`、`+3.3891e-6`，z10→z12 最大谱变化 `5.4096e-4 dex`，无警告。Reference 对外 76 个谱节点都达到 `−40` 下限，使用原始正值 log 谱比较。Fast 自洽根对 Reference z12 映射差 `+7.4970e-7`；z5、自适应、z8 最大谱差 `1.9847e-3`、`5.6956e-4`、`5.6931e-4 dex`。自适应延迟 13/76 模式，P50/P95 为 `4.14/4.25 ms`，较 z5 慢 `10.3%`、较 z8 快 `3.41×`。和 `.01` 点相比，两点的 adaptive 误差均通过，且都保留 PCHIP/尾部来源证据。

这扩大到 `cr>0` 低温有效 cutoff 的两个离散 kappa decade 端点，不认证 kappa 连续区间或固定振幅横切。新点逐节点谱、根、A/B 计时、RSS、环境和来源哈希见[`kappa10=.001` 交接 JSON](results/cr_positive_lowT_kappa001_handoff_20261116.json)，由同一[参数化诊断脚本](../../scripts/diagnose_cr_positive_highT_handoff_20261116.py)复跑。

### 高温 cutoff 边界上的高 `kappa10` 留出

再在 `T_re=10000` 固定温度，将 `kappa10` 从 `.01` 扩到 `1`，选高 kappa 点自己的 cutoff 边界上方 `r=1.1*r_cut=6.209353644613012e-37`。由于 r 也随解析 cutoff 调整，这与低温 decade 对照一样是沿 cutoff 边界的组合样本。该点初始 `Delta_N=.03177006>0`；一致性支路派生 `n_t=−7.7617e-38`（传入值仍为 0），`A_t=1.3043e-45` 为正规数。

Reference 的 z8→z10 映射变化 `−2.5281e-4` 高于根筛查参考线；保留该观察，用更深 z10→z12 的 `−3.9420e-5` 映射变化和 `8.7119e-5 dex` 最大逐点谱变化判断收敛。z12 求积/插值估计分别为 `5.7599e-53`、`1.7732e-41`，无 warning；公开 Reference 对数谱 59/77 节点触及 `−40`，比较使用直接求得的原始正值 log 谱。

同一 frozen root/grid/thread 的 z5、自适应、z8 最大/P95 谱差分别为 `1.7319e-3/1.0291e-4`、`7.0358e-4/7.6946e-5`、`8.3570e-5/2.8592e-5 dex`。z5 超 `1e-3 dex`，adaptive 延迟 8/77 后通过；当前 Fast 自洽根相对 Reference z12 映射为 `−6.0744e-6`，自映射残差为零。A/B P50/P95：z5 `3.75/4.30 ms`、adaptive `3.99/4.60 ms`、z8 `11.72/14.54 ms`；adaptive 比 z5 慢约 `6.4%`，比 z8 快 `2.94×`。这把高温 cutoff 边界的离散 `kappa10` 对照扩展到 `{.01,1}`，不认证固定 r 或连续范围。

完整节点谱、根、交接数量、逐次计时、RSS、源码/输入哈希见[高温 `kappa10=1` 交接 JSON](results/cr_positive_highT_kappa1_handoff_20261116.json)，由[同一参数化脚本](../../scripts/diagnose_cr_positive_highT_handoff_20261116.py)复跑。

## `cr>0` cutoff 有效域的 Sobol 独立预筛

为从 3×3 背景矩阵走向独立采样，对 `log10(T_re)∈[1,4]`、`log10(kappa10)∈[-3,0]` 生成固定 seed `20261116` 的 8 点 scrambled Sobol 设计。每点在一致性方程的精确 cutoff 根上取 `r=1.1*r_cut`。固定其余背景时

`Delta_N(r) = (1/3) ln(r/r_cut)`，

故预计 cutoff 裕量 `(1/3)ln(1.1)=0.0317700599`；代码仍对每个真实参数点重新计算 `Delta_N`、`A_t`，并只在物理 guard 下调用 Fast。8 点全在 `Delta_N≈0.03177006` cutoff 有效侧、`A_t` 为正规数、Fast 成功，映射固定点残差均为 0，节点数为 76 或 77，Fast 无 warning。这是小规模设计/脚本预筛，不能认证连续域。

Sobol index 0 `(T_re,kappa10)=(1186.4313,.0743135)` 和随后 index 1 `(186.6991,.0257680)` 用各自 Fast 原生网格和完全相同 Fast 根直接求 Reference z8/z10/z12。index 0 为预筛中预先指定的 Reference 点；index 1 是同一固定 8 点设计中的追加 Reference 留出，没有按它的结果选参数。两个点的 Reference 都是当前源码 SHA 一致的直接 Cartesian 求解，77/77 个信号节点正有限且无 Reference warning。

| Sobol 点 | T_re | kappa10 | Fast 根对 Reference z12 | Fast 默认自适应谱最大 / P95 差 | Reference z8→z10 映射 | z10→z12 映射 / 最大谱变化 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | `1186.43` | `.07431` | `−5.246e-6` | `5.600e-4 / 1.885e-4 dex` | `−2.628e-4` | `−3.041e-5 / 2.022e-4 dex` |
| 1 | `186.70` | `.02577` | `−7.590e-6` | `8.133e-4 / 2.538e-4 dex` | `−2.883e-4` | `+3.847e-6 / 2.417e-4 dex` |
| 2 | `16.73` | `.75557` | `−6.323e-6` | `4.509e-4 / 6.403e-5 dex` | `−9.944e-5` | `−4.311e-6 / 1.209e-4 dex` |

三点 Fast 根误差低于 `2e-4`，全部 76/77 直接比较节点的最大谱差低于 `1e-3 dex`。index 0/1 的 z8→z10 Reference 映射变化超过 `2e-4`，已保留；更深 z10→z12 根变化和逐点谱变化满足本轮收敛判断。index 2 位于低温、高 `kappa10` 端，z8→z10、z10→z12 变化均小于门槛，z12 Reference 用时约 `1333 s`，无 Reference warning。Sobol preflight 原脚本用 `r=minsub` 的仅背景 cutoff 探针读取 `N_re,rho_rad,rho_stiff`，它本身产生 “Delta_N < 0” 提示。单独审计表明：该探针没有送入 Fast 或 Reference；同一背景改用有效安全探针 `r=.01` 后 `N_re` 和两种能量密度逐位相同、无提示，解析 cutoff 与 8 点原结果逐点相同。真实 Sobol 样本全部用 `1.1*r_cut`，故提示归类为无效探针构造噪声，不是有效点的失败或传播警告；完整捕获和对照见[cutoff 探针审计 JSON](results/cr_positive_sobol_cutoff_probe_audit_20261116.json)。

8 点真实参数、Sobol 单位坐标、Fast 根/网格/错误状态、源码哈希和 index 0 原始直接 Reference 见[Sobol 预筛 JSON](results/cr_positive_sobol_cutoff_holdout_20261116.json)；index 1/2 分别见[Reference JSON](results/cr_positive_sobol_index1_reference_20261116.json)和[index 2 Reference JSON](results/cr_positive_sobol_index2_reference_20261116.json)。index 2 的初次运行使用了 index1 的复制脚本，计算输入、root 和 grid 已经交叉核验；[元数据校正脚本](../../scripts/correct_sobol_index2_reference_metadata_20261116.py)只修正 status/scope，并保留原计算脚本快照 SHA 与校正脚本 SHA，科学数据未改。因原 preflight 字段名把默认 adaptive 输出误标成 z5，对当前源码 8 点重新运行默认 adaptive 并明确对三个 Reference 留出重比较，结果见[纠正标签后的 adaptive Sobol JSON](results/cr_positive_sobol_default_adaptive_20261116.json)，可复跑[纠正验证脚本](../../scripts/verify_cr_positive_sobol_default_adaptive_20261116.py)。边界探针告警可用[独立审计脚本](../../scripts/audit_cr_positive_sobol_cutoff_probe_20261116.py)重现。该设计只有 3/8 点完成直接 Reference 收敛/谱检查，另外 5 点是 Fast 预筛；不构成连续参数域认证，后续需继续分层追加 Reference 留出。

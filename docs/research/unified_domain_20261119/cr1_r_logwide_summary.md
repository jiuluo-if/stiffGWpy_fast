# CR1 对数宽域 `r×T_re` 扩展核验

## 验收结论

当前 Fast 在固定 `Omega_bh2=Omega_ch2=7e5`、`kappa10=.01`、`H0=67.32117`、`DN_eff=0`、`A_s=2.100549e-9`、`cr=1`、`DN_re=10` 下，完成了候选区域

`T_re∈[1750,2000] GeV`，`r∈[4e-43,0.005]`。

`r` 方向跨约 40 个数量级，因此格点中点与 scrambled Sobol 均在 `log10(r)` 上均匀抽样。19 个共同输入中旧 Fast `1/19` 成功、当前 Fast `19/19` 成功：18 个旧失败均为正有限 `H²` 下的 `nonfinite_background_primitive`；旧版唯一成功点仍成功。所有当前结果提供完整 75–77 点频谱、有限根、无 warning。

这是有限格点、内部点和独立留出支持的区域证据，不是整个连续矩形的数学认证。生产 Fast 与 Reference 方程均未改；本轮验证当前已发布 Fast 在低张量幅度至通常 `r` 区域内的新求解能力。

## 数学有效域与截止

`cr>0` 实际派生 `n_t=-r/8`，在候选区间内始终为负（约 `-5e-44` 至 `-6.25e-4`），因此正倾斜功率截止不触发。模型派生的再热截止仍由 `N_inf=N_re+Delta_N` 决定，且必须 `Delta_N≥0`：

`Delta_N = 4/3 (N_10-N_re) + 4/3 ln(M_Pl) + 1/3 ln(45 r A_s / (2 (rho_rad,re+rho_stiff,re)))`。

固定其他参数后，真实下界为 `r_min(T)=2(rho_rad,re+rho_stiff,re)/(45 A_s) × exp[-4(N_10-N_re)-4 ln(M_Pl)]`。由等价关系 `r_min=r0 exp(-3 Delta_N(r0))` 求值，并在 `T_re∈[1750,2000]` 上 101 点检查：端点下界约 `1.6161e-43`、`3.6020e-43`，最大值出现在高温端。候选下界 `4e-43` 在所有检查温度上满足 `N_inf>N_re`；最小间隔约 `0.03493` e-fold。所有 19 个配对样本 `H²` 正有限，最小样本值约 `0.453214`，`N_inf` 约在 `[37.8073,68.7164]`；未把非正 `H²`、无效截止或不可表示输出算成新成功。

更高 `r` 的宽盒 `[.001,.04]` 虽然 Fast screen 为 `1/19→19/19`，但不能整盒验收：`(T_re=1750,r=.04)` 的 z8 Reference 根差约 `9.47e-4`、谱差 `0.001270 dex`；`r=.0205` 探点也超过根误差门槛。因此本候选只报告通过精度验证的上界 `.005`，不把高 `r` 探点失败分类隐藏为物理拒绝。

## 精度与 Reference 下限处理

七个区域代表样本（四角、对数中段旧成功点、`r=1e-40` 低幅度点、独立 log-Sobol 留出）均以当前 Fast 同一频率数组做 z8 Reference。最大根误差为 `1.2886e-4`；对最大根误差角继续做两次 Reference map，第二至第三次闭合相对差 `1.2115e-7`，Fast 根相对第三次 map 为 `1.2886e-4`。对 Reference 的原始正线性量 `Ogw-Oj` 直接取 `log10` 后，最大谱差为 `8.7991e-4 dex`。均满足 `2e-4` 与 `1e-3 dex` 门槛。

Reference 返回的 `log10OmegaGW` 对小信号固定钳到 `-40`；在 `r≤1e-40` 的样本中，该字段不适合作谱误差比较。本轮未改 Reference，而是对其保留的正 `Ogw-Oj` 原始线性值取对数，与 Fast 的 `log10OmegaGW` 比较。比如 `r=4e-43,T_re=2000` 的根误差 `2.894e-5`、原始线性谱误差 `8.131e-4 dex`；`r=1e-40,T_re=1875` 分别为 `1.022e-4`、`7.982e-4 dex`。信号本身约 `1e-42`，仍远高于 binary64 下溢，根和完整频谱均非零且有限。

唯一旧 Fast 成功点为 `(T_re=2000,r=4.4721e-23)`。旧根 `6.91584e-22` 在旧根处 Reference map 的相对误差为 `5.304%`；新根 `6.56745e-22` 对同一频率网格的当前源码 Reference map 误差 `1.172e-5`，旧/新频率网格逐位相同。旧 Fast 谱对原始 Reference 线性输出最大差 `4.63 dex`，当前谱差 `8.125e-4 dex`。因此该点从“旧版返回成功”转为更准确的成功解，没有丢失成功输入。

## 500 对完整 Fast 性能 A/B

使用同机 CPU affinity `[0,1]`、2 Numba 线程、TBB、BLAS 1 线程、独立常驻进程和 Numba cache，每侧预热 3 次后交替 500 对完整 Fast 求解。对照源 `d4ee1b9` 与当前生产源码相同；所有样本根、频率数组、频谱逐位相同。

| 样本 | P50 比 | P95 比 |
|---|---:|---:|
| 低 `r` cutoff 角 | `0.9829` | `0.9893` |
| 高 `r` 角 | `1.0231` | `1.0250` |
| 独立 log-Sobol | `0.9906` | `1.0149` |
| 旧 Fast 成功中段点 | `0.9853` | `0.9650` |
| 可信 CR0 控制点 | `0.9993` | `0.9810` |

最大新增区 P50/P95 比 `1.0231/1.0250`，可信区均未出现 5% 回退。首次进程调用（含 JIT）基线/当前为 `20.675/20.404 s`；不将它与预热后的稳态计时混合，也不把微小计时差宣称为算法加速。

## 复现产物

- 对数宽域 paired screen、格点与独立 Sobol 原始数据：`screen_cr1_r_logwide.py`、`screen_cr1_r_logwide_lattice.json`、`screen_cr1_r_logwide_sobol.json`。
- 当前源码 Reference、原始线性谱、cutoff 邻点、低振幅端和旧成功点核验：`validate_cr1_r_logwide_z8.py/json`、`validate_cr1_r_cutoff_probe_z8.py/json`、`validate_cr1_r_log_tail_probe_z8.py/json`、`audit_cr1_r_old_success_reference.py/json`。
- 根闭合、500 对完整 A/B、覆盖与误差图：`close_reference_cr1_r_lowstrip_z8.py/json`、`benchmark_cr1_r_logwide_ab.py/json`、`plot_cr1_r_logwide.py`、`cr1_r_logwide_evidence.png`。
- 更高 `r` 被拒绝的 Fast/Reference 边界记录：`screen_cr1_r_reheat_*`、`validate_cr1_r_reheat_z8.json`、`validate_cr1_r_reheat_upper_probe_z8.json`。
- 源码和产物 SHA-256 汇总见 `cr1_r_logwide_manifest.json`。

# CR1 温度轴相邻区域扩展

## 本次新增的有限采样区域

固定 `cr=1`、`kappa10=.01`、`Omega_bh2=Omega_ch2=7e5`、`H0=67.32117`、`DN_eff=0`、`A_s=2.100549e-9`、`DN_re=10`，在已有高倾斜区域 `r∈[.005,.01]` 的温度上边界继续筛查。CR1 的实际倾斜为 `n_t=-r/8`；本区域 `n_t∈[-.00125,-.000625]`，均为红倾斜。所有屏幕样本 `N_inf` 有限（`68.465–68.814`）、`H²` 正且有限（样本最小值 `0.4532139932`），Fast 物理 guard 均接受，输出完整。

配对格点加独立 Sobol 共 19 点：`T_re∈[2000,2250] GeV`、`r∈[.005,.01]`。旧 Fast 成功 `1/19`，其余 `18/19` 均明确为 `nonfinite_background_primitive`；旧成功 `(2250,.0075)` 单独保留。当前 Fast 为 `19/19` 成功。该片与已验证的 `[1875,2000]×[.005,.01]` 在温度边界相邻，因此有限样本支持把已发布 Fast 的该 CR1 参数片从 `T_re≤2000` 扩到 `T_re≤2250`。

这是有限采样覆盖证据，不是整个连续矩形的数学认证。

## 精度与自洽闭合

Reference 使用未修改的生产代码，在当前 Fast 自洽根和相同频率节点计算固定背景 map；每个新增 Reference 样本保留原始正线性 `Ogw-Oj` 后再取对数。代表点包括新低 `r` 温度边角、旧 Fast 成功点、两个独立 Sobol 留出点。新增样本最大根误差 `1.6741e-4`，最大有效谱误差 `8.6204e-4 dex`，均低于 `2e-4` / `1e-3 dex` 门槛。

最紧 Sobol 留出完成第三次 map：Fast 根与第三次 map 相差 `1.6753e-4`，第二至第三次 map 相对变化 `6.73e-9`。复用的边界核验中，`(T_re=2000,r=.01)` 根误差为 `1.0933e-4`；CR1 新温度高 `r` 边角 `(2250,.01)` 为 `1.5426e-4`。谱误差汇总仅覆盖 Reference 原始线性差分为正、可直接比较的频点；低于线性表示范围的点不伪造对数误差。

较高 `r` 的两个候选矩形 `[1875,2250]×[.01,.015]` 和 `[1875,2250]×[.01,.01125]` 均拒绝：前者的 Reference 高 `r` 边角及一 Sobol 点根误差超过门槛；后者的低温高 `r` 角、旧 Fast 成功内点及一 Sobol 点超门槛。筛查记录、失败分类和 Reference 点保留在相邻 `nextband` / `narrowband` 文件中，未将这些点并入接受区域。

## 性能 A/B

同机、固定两个 Numba 线程、TBB、BLAS 单线程、独立常驻进程、每边预热三次，对七组完整 Fast 输入交替运行 500 对。基准代码为 `d4ee1b9c8f9fad401c5e22c2316094a608074727`，当前为 `c0bde6469293905413a9309a7c3cc04a4d1531cd`。最差稳态分位数比为 P50 `1.0052`、P95 `1.0199`；可信 CR0 控制分别为 `0.9954` / `0.9918`。所有七组 root、频率摘要和频谱完全相同（root delta `0`、频率摘要不匹配 `0`、最大谱差 `0 dex`）。这次只扩展已发布数值机制的证据范围，没有改生产 Fast 或 Reference 源码，因此无新增热路径回退。

单次首调时长与预热后 500 对样本都在原始 JSON 中独立保存。`paired_ratio_p95` 对少数配对测量的瞬时抖动较敏感；正式门槛按两侧各自的稳定运行 P50/P95 分位数比较，同时保留每对原始时间，未用冷启动数据替代稳态数据。

## 可复核文件

- 配对原始 screen：`screen_cr1_r_temperature_extension.py`、`screen_cr1_r_temperature_extension_lattice.json`、`screen_cr1_r_temperature_extension_sobol.json`。
- Reference 与 map 闭合：`validate_cr1_r_temperature_extension_z8.py/.json`、`close_reference_cr1_r_temperature_extension_sobol_z8.py/.json`。
- 性能 A/B：`benchmark_cr1_r_temperature_extension_ab.py`、`benchmark_cr1_r_temperature_extension_ab_500pairs.json`。
- 被拒高 `r` 候选：对应 `screen_cr1_r_nextband_*`、`screen_cr1_r_narrowband_*`、`validate_cr1_r_nextband_z8.*`、`validate_cr1_r_narrowband_z8.*`。
- 覆盖、失败类别与精度图：`cr1_r_temperature_extension_evidence.png`。
- 源码和全部原始数据 SHA-256：`cr1_r_temperature_extension_manifest.json`。

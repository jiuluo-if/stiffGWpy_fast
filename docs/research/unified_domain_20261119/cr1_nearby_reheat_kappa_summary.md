# CR1 低温相邻矩形扩域记录

## 本轮结论

在 `cr=1` 下，验证了与既有 `T_re∈[1750,2000]` 区域相邻的有限采样矩形 `T_re∈[1687.5,1750]`、`kappa10∈[0.0125,0.015]`。同一组 19 个输入上，旧 Fast `0/19`、当前 Fast `19/19`；旧版失败全部归为 `nonfinite_background_primitive`，新 Fast 均返回完整 76 点频谱且无 warning。结论只覆盖所采样格点、内部点和独立 Sobol 留出，不构成整个矩形的连续认证。

这次没有改生产求解器或 Reference。可求值能力来自当前已发布的统一背景表示；本轮新增的是低温相邻参数块的旧/新配对、独立精度、性能和拒绝边界证据。

## 数学有效性与 CR1 约束

固定 `Omega_bh2=Omega_ch2=7e5`、`H0=67.32117`、`DN_eff=0`、`A_s=2.100549e-9`、`r=0.01`、`DN_re=10`。项目派生关系为 `n_t=-r/8`，所以本块实际 `n_t=-0.00125`，不使用输入 `n_t=0` 代替。全体样本 `N_inf` 有限且 `H²` 正有限，均通过 `physical_guard=True`。由于实际倾斜为负，正倾斜功率截止条件不触发；有限 `N_inf` 截止仍有效。没有把数学无效点或非正 `H²` 点计为 Fast 失败扩域。

旧 Fast 失败记录来自 `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45`；性能对照为发布源基线 `d4ee1b9c8f9fad401c5e22c2316094a608074727`。本轮 `stiffgwpy_fast` 相对性能基线没有源码差异。

## 精度核验

当前源码在低温边两角和一个独立 Sobol 留出上，以相同 Fast 频率数组调用未修改的 z8 Reference：

| 样本 | `T_re` | `kappa10` | Fast 根相对 Reference map | 76 点最大谱差 |
|---|---:|---:|---:|---:|
| 低温、低 κ 角 | 1687.5 | 0.0125 | `1.517930e-4` | `9.223931e-4 dex` |
| 低温、高 κ 角 | 1687.5 | 0.015 | `1.347906e-4` | `8.277459e-4 dex` |
| 独立 Sobol 留出 | 1708.1084 | 0.0128840 | `1.941642e-4` | `8.466765e-4 dex` |

Sobol 留出继续迭代 Reference map 至第三次，第二到第三次闭合相对差 `8.504492e-8`，Fast 根相对第三次 map 为 `1.942825e-4`。三项代表样本均满足根误差 `≤2e-4` 和频谱误差 `≤1e-3 dex`；Sobol 根余量较小，因此不向更低温或更宽 κ 范围外推。

更低温试探矩形 `[1500,1750]×[0.0125,0.015]` 未获接受：`(1531.25,0.013125)` 与 `(1500,0.0125)` 的 Fast 根相对第三次 Reference map 分别约 `2.533240e-4`、`2.693125e-4`，超出根误差门槛；其谱误差虽通过，也不足以挽救根误差失败。对应首轮 z8 谱差为 `8.538028e-4` 和 `8.649403e-4 dex`。这是本轮明确保留的失败边界证据。

## 性能 A/B

Windows 同机、CPU affinity `[0,1]`、2 个 Numba 线程、TBB、BLAS 1 线程、独立常驻进程和 Numba cache。每个样本先做 3 次预热，再交替运行 500 对完整 Fast；包含频谱装配和自洽根。对照源 `d4ee1b9` 与当前 `9681c0c` 的生产目录相同。

| 样本 | P50 比 | P95 比 | 根/频率/谱差 |
|---|---:|---:|---|
| 新块低 κ 角 | `0.9656` | `0.9771` | 根 0，频率逐位相同，谱 0 |
| 新块高 κ 角 | `0.9668` | `0.9744` | 根 0，频率逐位相同，谱 0 |
| 新块独立 Sobol | `0.9676` | `0.9610` | 根 0，频率逐位相同，谱 0 |
| 可信 CR0 控制点 | `0.9963` | `0.9901` | 根 0，频率逐位相同，谱 0 |

首次进程调用（含 JIT）基线/当前为 `20.595/20.417 s`；该次编译计时不与后续热调用混合。由于生产源码完全相同，这些比例用于核实当前发布 Fast 在新区域没有可测性能回退，不将计时差值宣称为本轮算法加速。

## 产物与复现

- 区域配对与留出：`screen_cr1_nearby_reheat_kappa.py` 及两个 `screen_...json`。
- Reference 与 map closure：`validate_cr1_nearby_reheat_kappa_z8.py/json`、`close_reference_cr1_nearby_reheat_kappa_z8.py/json`。
- 完整 Fast A/B：`benchmark_cr1_nearby_reheat_kappa_ab.py/json`。
- 覆盖、失败分类与精度图：`plot_cr1_nearby_reheat_kappa.py`、`cr1_nearby_reheat_kappa_evidence.png`。
- 更低温拒绝样本及闭合数据分别保存在 `screen_cr1_low_reheat_kappa_*`、`validate_cr1_low_reheat_kappa_z8.json`、`close_reference_cr1_low_reheat_kappa_z8.json`。
- 源码哈希、数据哈希和有限样本声明见 `cr1_nearby_reheat_kappa_manifest.json`。

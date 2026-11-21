# CR0 极端密度闭合相消扩域核验

## 结论范围

本轮候选只修复一个已定位的数学合法失败：当重子密度远大于今天的有限闭合值时，旧式 Friedmann 闭合展开把多个约 `10^16` 的项直接相加，舍入后得到 `H²=0`。候选改用严格等价的闭合恒等式，保留 H² 的正负号和非有限状态检测。它没有把数学无效 H² 改成有限值，也没有更改再加热转折、张量方程、WKB 交接或物理参数含义。

有限样本显示，CR0 的两个二维区域均从 HEAD `160527d` 的 `0/19` 成功提高到候选的 `19/19`：

- `Omega_bh2 ∈ [1e16, 2e16]`、`H0 ∈ [60, 75]`，其余参数固定。
- `Omega_bh2 ∈ [1e16, 2e16]`、`T_re ∈ [1500, 2500] GeV`，其余参数固定。

每个区域均有 3×3 矩形格点、2 个预先指定内部留点和 seed 固定的 8 点 scrambled Sobol 样本。HEAD 的 38/38 失败均为 `nonpositive_h2`；候选的 38/38 均返回完整有限频谱、正有限 H² 且没有警告。候选区域内频谱根求解保持 Fast 计算量级。以上是有限覆盖证据，不构成连续域认证。

## 数学改写与选择条件

记 `x=N_inf−N`、`E=exp(x)`、`q=exp(−3x)`。原闭合项包含

```text
Omega_m + Omega_sh E^3 + Omega_L q
Omega_L = h^2 - Omega_m - Omega_mnu - (2/3) Omega_nu - Omega_ph - Omega_sh - Omega_er
```

将相消的今天闭合值保留为独立锚点

```text
A = h^2 - Omega_mnu - (2/3) Omega_nu - Omega_ph
B(x) = Omega_m (1-q) + Omega_sh (E^3-q) + A q
      = -Omega_m expm1(-3x) + 2 Omega_sh sinh(3x) + A exp(-3x)
```

不同中微子分支继续添加其原方程中的辐射、热中微子和额外辐射项；`Omega_er` 项写成 `Omega_er (E-q)`，所以 `x=0` 时额外辐射对闭合的净贡献严格为零。该恒等式来自对原 `Omega_L` 代入并逐项合并，没有改变背景模型。

旧式和式的舍入绝对误差尺度约为 `u (|Omega_m|+|Omega_sh|)`，其中 `u` 是 float64 unit roundoff。只有当该误差尺度超过 `sqrt(u)|A|`（等价于闭合条件数大于约 `1/sqrt(u)`）时才启用稳定式；这是统一的浮点条件数门限，不依赖参数名或采样点。正常可信区保持旧热路径。稳定式用 `expm1` 保留小 `x` 下的 `1−exp(−3x)` 差值，用 `sinh` 直接形成 `exp(3x)−exp(−3x)`。`H²≤0` 仍由原传播失败路径报告；中间量不可表示和最终输出不可表示也未被钳位。

背景 `sigma` 的分子及所有分支条件不变，只将分母 `H²` 换为等价稳定值。`fast_sgwb.gen_kernel` 与 `exact_background.H2_vec/sigma_vec` 使用同一恒等式，防止准备阶段和传播阶段采用不一致的闭合值。

## 方法来源与独立核验

- Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed.：依照项目已有论文审查记录，查阅第 2 章 §2.1–2.2 的浮点舍入模型和第 3 章 §3.1 的误差传播分析。适用假设是 IEEE 双精度的正常舍入模型；它支持条件数/舍入误差门控，不提供宇宙学模型结论。项目已有对应记录见 `cr0_2m_to3m_wkb_basis_guard_summary.md`。
- NumPy `expm1` 官方 API 说明其在小输入时比直接 `exp(x)-1` 精度更高；本候选将其用于同一解析差值。来源：[NumPy `expm1`](https://numpy.org/doc/stable/reference/generated/numpy.expm1.html)。
- 极端密度处未调用生产 Reference，因为其绝对 `N` 坐标的背景闭合会经历同一种浮点相消。本轮只在研究验证器中使用 `x=N_inf−N` 作为独立积分坐标；生产 `reference.py` 未修改，哈希与 HEAD 完全相同：`9066adcf69f95398405e77e5a5d47942573041f13d0a0ee1554c4bd119eb17f6`。
- 验证器重写现有 Prüfer 方程，状态为相位 `theta` 与归一化对数振幅 `b=log(A)+0.5 log(H²)+0.5N`。利用背景恒等式 `d log(H²)/dN=3(1−sigma)` 可得 `db/dx=−cos(2theta)`；另有 `dtheta/dx=exp(z)+sin(2theta)`，`z=(log10 f−f_hor(x))ln(10)` 直接由同一背景算出。再加热前继续使用 Reference 的解析物质主导视界延拓，并在真实 `x=DN_re` 边界处分段。初态、状态方程、尾部交接阈值 `z_tail=8` 及输出换算与 Prüfer 定义一致。
- 状态变换可逆：Reference 初态 `(z0, log(A0), theta0)=(z0,z0,0)` 对应 `b0=z0+0.5log(H²(x_start))+0.5(N_inf−x_start)`；任意点恢复 `log(A)=b−0.5log(H²)−0.5(N_inf−x)`。因此相位、振幅、WKB 交接系数及今天的 `Ogw/Oj/Opgw` 仍由原方程输出。再加热前视界延拓使用恒定 `H²_re` 锚点，且在 `x=DN_re` 两侧分段，不平滑真实转折。
- 验证器用正常密度控制点与未修改 Reference Prüfer 单模对照：低频差 `5.02e-5 dex`，再加热前高频差 `7.44e-9 dex`。极端密度处对 Fast 返回的完整自适应频率节点逐点传播；根映射使用项目现有 PCHIP + Gauss–Kronrod 积分，并在独立根估计处重算一次映射。该结果核验 Fast 张量输出和自洽根，不声称验证未采样频率之间的连续谱。

## 精度与成本

| 区域代表点 | Fast 自洽根 | 独立闭合根 | 根相对误差 | 最大频谱误差 | 频率节点 | Fast 暖运行 | 独立映射闭合 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 密度×H0：`(1e16, 67.32117)` | `1.0295676908793673e-5` | `1.0296897130666834e-5` | `1.1852e-4` | `3.6131e-4 dex` | 265 | `1.306 s` | 相对变化 `4.46e-14` |
| 密度×T_re：`(1.5e16, 2500 GeV)` | `1.544515250888175e-5` | `1.5445682151870935e-5` | `3.4292e-5` | `1.4572e-4 dex` | 277 | `1.167 s` | 相对变化 `2.00e-14` |

两个代表点均满足根误差 `≤2e-4`、最大有效频谱误差 `≤1e-3 dex`。独立映射的自洽复算差远小于根门槛。每个新区域 19 点的正式 goal 配置候选暖运行 P50/P95 分别为 `0.383/0.617 s` 和 `0.390/0.615 s`；自适应精度核验的完整 Fast 调用约 `1.17–1.31 s`。

旧可信配置 A/B 在 3 个参数组各 48 次交替暖样本，候选与 HEAD 根和谱逐项相同；P50/P95 候选/HEAD 比率：

| 参数组 | P50 比率 | P95 比率 |
|---|---:|---:|
| 默认 CR1 | `0.985` | `0.974` |
| 高温 CR1 | `1.016` | `1.041` |
| 已接受 CR0 高密度边缘 | `1.017` | `1.011` |

三组可信样本的暖运行 P50/P95 均未超过 5% 回退上限。三个独立空 Numba 缓存冷启动样本：HEAD 首次默认调用为 `[11.732, 11.102, 11.810] s`，候选为 `[11.259, 13.553, 11.303] s`；两侧中位数分别为 `11.732 s` 和 `11.303 s`。候选的一次冷样本较慢，保留在原始数据中。新区域 goal 配置 Fast 暖 P50/P95 为 `0.383/0.617 s` 与 `0.390/0.615 s`；完整自适应 Fast 调用在两个代表点分别用时 `1.306 s` 和 `1.167 s`。

## 代码、数据与重现

- 生产差异：`stiffgwpy_fast/exact_background.py`、`stiffgwpy_fast/fast_sgwb.py`；红绿回归测试：`tests/test_physics_limits.py`。
- 旧/新 38 个同输入结果：`screen_cr0_extreme_density_h0.json`、`screen_cr0_extreme_density_reheat.json`；屏幕运行器：`screen_cr0_extreme_density_slices.py`。
- 两个代表点的逐频原始数组、每模 ODE 调用数、源码哈希和根闭合数据：`cr0_extreme_density_x_oracle.json`、`cr0_extreme_density_reheat_x_oracle.json`；验证器：`verify_cr0_extreme_density_x_oracle.py`。
- 可信区性能 A/B：`benchmark_extreme_density_closure_trusted_ab.json`；图：`cr0_extreme_density_expansion.png`；文件哈希与汇总：`cr0_extreme_density_expansion_manifest.json`；绘图脚本：`plot_cr0_extreme_density_expansion.py`。
- HEAD 基线为 `160527dc9ac7ed5a315f7d6a2087883a8ae71b1a`。候选源码 SHA256：`exact_background.py=d07ead49e420823882abc9c471739aa2e3b0b9d74023afb249213ebbe6ed53b0`，`fast_sgwb.py=011582d4ffd4b6efbf2b3ed97cecd2965ea33b2cc79183586dc88ddb0eb24553`。
- 本地交付门禁：全量默认 pytest `447 passed, 10 deselected`；相关回归 `37 passed, 5 deselected`；Ruff、Mypy、compileall、manifest、兼容 smoke、sdist/wheel、发布归档边界及已安装 wheel smoke 全部通过。GitHub CI 需在推送后复核。

## 限制

1. 两个切片的 Sobol/格点是有限采样；不可外推为完整连续域数学认证。
2. 代表点的独立 Prüfer 验证器是为避免绝对坐标闭合失败所写的研究核验工具，不是生产 Reference 的替代品，也没有覆盖所有 CR0/CR1 参数族。
3. 频谱验收覆盖 Fast 自适应支撑上的全部原生节点；没有把节点之间未采样频率称为严格连续误差界。
4. 本候选只扩展高密度闭合相消区域，不解决更大 look-back 下 `exp(N_inf−N)`、`S2`/`S2inv` 及最终不可表示输出的一般动态范围问题。此类异常仍须按真实可表示性分类。

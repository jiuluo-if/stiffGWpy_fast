# CR0 高密度负倾斜的温度轴扩域

## 区域与约束

固定 `cr=0`、`Omega_bh2=1e6`、`Omega_ch2=0.1201075`、`H0=67.32117`、`DN_eff=0`、`A_s=2.100549e-9`、`r=.01`、`DN_re=10`、`kappa10=.01`。筛查区域为
`T_re∈[4e3,1e7] GeV × n_t∈[-.5,-.25]`。

`cr=0` 时 `n_t` 是独立输入；该区域全为负倾斜，不触发正倾斜 UV cutoff。温度、倾斜、`r`、`A_s`、`kappa10` 均处于 Cobaya 对应范围，`Omega_bh2=1e6` 超出其采样先验，但代码方程可求值；因此这里报告 Fast 的数学/数值有效域，不声称它属于 Cobaya posterior 范围。所有接受样本的 `H²` 正且有限、`N_inf` 有限，并通过 `physical_guard=True`（共享总 `N_eff≤5`）。

## 配对 Fast 区域证据

按 `log10(T_re)` 与 `n_t` 构造 3×3 格点、2 个内部点及独立 8 点 Sobol 留出。共 19 个同输入配对样本：旧 Fast `2/19` 成功、`17/19` 因 `nonfinite_background_primitive` 失败；当前 Fast `19/19` 成功。两个旧版原成功点都保留，并分别通过当前 Reference 核验。所有样本 `N_inf` 有限（`48.466–56.291`），样本最小 H² 为 `0.4532139930`，输出频率为 75–76 点。

一个更宽的 `n_t∈[-.5,0]` 探索框有 3 点触发 `shared_Neff_guard`，所以没有把整框纳入接受区；正倾斜高温候选 `T_re∈[1e6,1.5e6]`、`n_t∈[.25,.5]` 的 19 点也全部被 guard 拒绝（旧/新 Fast 相同）。接受区域因此限定为 `n_t≤-.25`。这些 guard 拒绝仍保留在原始筛查数据中，没有算作数值失败或求值成功。另一个 `log10(kappa10)×n_t` 候选的 19 点旧/新 Fast 都成功，没有产生旧失败到新成功的扩展，不计入成果。

## Reference 精度与自洽闭合

当前源码 z8 Reference 核验覆盖四个新矩形边角、两个旧版原成功点和一个独立 Sobol 留出。Reference 使用相同 Fast 频率网格、`rtol=1e-9`、未修改生产源码，并对其原始正线性 `Ogw-Oj` 取对数。最大首次 map 根差为 `1.88163e-4`，最大有效频谱差为 `8.81429e-4 dex`。

最紧的右上角 `(T_re=1e7 GeV,n_t=-.25)` 经第三次 map 核验：二至三次 map 相对变化 `3.27e-8`，Fast 根相对第三次 map 差 `1.88135e-4`。所有代表样本都通过 `2e-4` 根与 `1e-3 dex` 谱门槛；最紧角的余量较小，不能据有限样本向矩形外推。

## 完整求解性能 A/B

对照 Fast 源码提交 `d4ee1b9c8f9fad401c5e22c2316094a608074727` 与当前 `67a02040f85180341fbada64d6bba6672fb2b0c3`。`exact_background.py`、`fast_sgwb.py`、`freq_adaptive.py`、`_wkb_transfer.py` 的 Git blob 相同。同机 TBB、2 Numba 线程、BLAS 单线程、独立持久进程与缓存、每边预热三次，每个区域/留出/控制配置交替进行 500 次完整求解。

接受区最大稳态 P50/P95 比为 `1.0153/1.0211`；高密度 CR0 控制为 `0.9590/0.9513`。普通 CR0 控制最初 500 对的 P95 比接近门槛，因此另做 1000 对复测；合并 1500 对后的 P50/P95 比为 `1.0020/1.0244`，通过 5% 限制。所有基准配置的 root、频率摘要和完整频谱完全相同。冷启动首次含 JIT 求解为基线/当前 `20.668/21.060 s`；其余首次调用与预热后的逐对原始时间均存于 A/B JSON。

本轮没有修改生产 Fast 或 Reference 源码；这是对已发布统一稳定背景路径的新增 Fast 可求值区域验证。采样覆盖支持该矩形的区域改善，但不是连续域数学认证。

## 原始数据、图与复现

- 配对 lattice/Sobol：`screen_cr0_density_reheat_negative_tilt.py`、对应 `*_lattice.json` / `*_sobol.json`。
- Reference 与闭合：`validate_cr0_density_reheat_negative_tilt_z8.py/.json`、`close_reference_cr0_density_reheat_negative_tilt_z8.py/.json`。
- 完整 Fast A/B 与 1000 对控制复测：`benchmark_cr0_density_reheat_negative_tilt_ab.py`、`benchmark_cr0_density_reheat_negative_tilt_ab_500pairs.json`、`benchmark_cr0_reheat_trusted_control_1000pairs.py/.json`。
- guard 失败的宽框及无旧失败的 `kappa10` 探索数据保留在相邻 `screen_cr0_density_reheat_tilt_*`、`screen_cr0_tilt_temperature_*`、`screen_cr0_kappa_tilt_extension_*` 文件中，不纳入接受区。`n_t∈[.25,.5]` 的高倾斜温度探针 19/19 都是共享 N_eff guard 拒绝。
- 覆盖、失败分类、Reference 精度图：`cr0_density_reheat_negative_tilt_evidence.png`。
- 源码和原始数据哈希：`cr0_density_reheat_negative_tilt_manifest.json`。

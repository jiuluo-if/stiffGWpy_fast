# CR1 较高 `r` 相邻区域扩展

## 当前接受区域

在固定 `Omega_bh2=Omega_ch2=7e5`、`kappa10=.01`、`H0=67.32117`、`DN_eff=0`、`A_s=2.100549e-9`、`cr=1`、`DN_re=10` 下，当前 Fast 的连续参数片获得以下有限采样验证：

1. `T_re∈[1750,2000] GeV`，`r∈[.005,.0075]`。
2. `T_re∈[1875,2000] GeV`，`r∈[.005,.01]`。

两块与已验证的 `r∈[4e-43,.005]` 低幅度片相连。CR1 实际倾斜 `n_t=-r/8`；以上两个高 `r` 片 `n_t∈[-.00125,-.000625]`，为红倾斜，不触发蓝倾斜功率 cutoff。所有当前 screen 点 `N_inf` 有限、H² 正有限，且输出完整谱。

## 配对 screen 与旧结果分类

中等 `r` 矩形的 19 个输入中旧 Fast `1/19` 成功、当前 Fast `19/19` 成功；旧失败均为 `nonfinite_background_primitive`。

较高 `r` 矩形单轮 screen 为旧 Fast `0/19`、当前 Fast `19/19`，旧屏幕记录 18 个 `nonfinite_background_primitive` 和 1 个 15 秒 `resource_timeout`。同一参数的旧成功点 `(1968.75,.00875)` 在旧 Fast 独立重跑 20 次均成功；它在一轮 paired screen 中曾有一次不同的失败结果，因此不计为新增成功点。旧版资源超时点 `(1983.43,.006932)` 延长到 120 秒仍为 `resource_timeout`，且 H² 正有限；当前 Fast 求解成功。两类情况分别保留，未改成物理拒绝。

## 精度核验

- 中等 `r` 低温高 `r` 角 `(1750,.0075)` 对 Reference 的第三次 map 根差 `1.8944e-4`，二至三次 map 闭合 `1.29e-7`。
- 较高 `r` 片最紧独立 Sobol 留出第三次 map 根差 `1.9498e-4`，二至三次闭合 `7.51e-8`。
- 所有 Reference 样本最大谱差 `0.000903 dex`，低于 `1e-3 dex`；最大根差低于 `2e-4`。
- 更低温高 `r` 点 `(1750,.01)` 明确不接受：根差 `2.4147e-4` 超过门槛，尽管谱差 `0.000896 dex` 通过。

旧 Fast 的两个已成功点都经当前未修改 Reference 审计。高 `r` 点 `(1968.75,.00875)` 的旧根误差 `5.80%`、旧谱差 `4.64 dex`，当前 Fast 根误差 `1.754e-4`、谱差 `8.62e-4 dex`；中等 `r` Sobol 旧根误差 `5.91%`、旧谱差 `4.63 dex`，当前根误差 `1.483e-4`、谱差 `8.76e-4 dex`。两处旧/新频率网格均有 76 点但数组不同；旧与新输出分别用各自网格的 Reference 核验。

## 完整 Fast 性能 A/B

对照源 `d4ee1b9` 与当前生产 Fast 目录相同。同机、2 Numba 线程、TBB、BLAS 1 线程、独立常驻进程/cache，各侧预热 3 次并交替做 500 对完整求解。较高 `r` 片 P50/P95 最大比为 `1.0363/1.0286`，可信 CR0 控制点为 `0.9973/1.0030`。

中等 `r` 片低温高 `r` 角首轮 500 对的 P95 比是 `1.0568`，略过门槛。将该角扩到 1000 对并与首轮合并成 1500 对后，P50/P95 比为 `1.0109/1.0261`；其它 500 对点 P95 最高为 `1.0279`，CR0 控制点为 `1.0010`。稳定计时通过 5% 门槛。

## 证据与限制

这是格点、内部点、独立 Sobol 与代表 Reference 样本支持的有限采样结果，不是连续域认证。生产 Fast 和 Reference 源码均未改；扩展来自已发布 Fast 的实际求值能力。

- 区域 paired screen 与原始数据：`screen_cr1_r_midband_*`、`screen_cr1_r_upperband_candidate_*`。
- Reference、闭合与旧成功审计：`validate_cr1_r_midband_z8.*`、`validate_cr1_r_upperband_candidate_z8.*`、`close_reference_cr1_r_mid_upper_z8.*`、`close_reference_cr1_r_upperband_*`、`audit_cr1_r_upperband_old_success_z8.*`。
- 旧版超时与重复性原始记录：`audit_cr1_r_upperband_timeout_recheck.*`、`audit_cr1_r_upperband_old_repeatability_20.json`。
- 完整 Fast A/B：`benchmark_cr1_r_midband_*`、`benchmark_cr1_r_upperband_*`。
- 参数覆盖、失败分类与精度图、源码哈希：`cr1_r_upperband_evidence.png`、`cr1_r_upperband_manifest.json`。

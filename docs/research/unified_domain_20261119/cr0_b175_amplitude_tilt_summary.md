# CR0 正倾斜—小幅度二维 Fast 扩域核验（`c3c6ab3`）

## 参数区域与配对屏幕

固定 `Omega_bh2=1e6`、`Omega_ch2=0.1201075`、`T_re=2000`、`DN_re=10`、`kappa10=0.01`、`cr=0`。CR0 允许独立输入 `n_t`。目标矩形为 `log10(r)∈[-12,-8]` 与 `n_t∈[0.05,0.1]`；所有样本 `r>0`，`N_inf=47.7723412841` 有限，H² 全程正有限，`physical_guard=True`。

旧版 `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45` 与当前 `c3c6ab3fc34e705ebdea3d16739740a227d1acb7` 同输入配对：3×3 (`r=1e-12,1e-10,1e-8`; `n_t=.05,.075,.1`) 加两个内部点为 `0/11→11/11`；独立 scrambled Sobol 留出（seed `20261123`，8 点）为 `0/8→8/8`。旧版 19 点全部为 `nonfinite_background_primitive`，H² 正有限；候选全成功，75 个频率点完整且无警告。Fast 屏幕候选耗时中位数约 `77.3 ms`。

当前源 SHA-256：`exact_background.py=c55ff5755615a0eaa5a2f4c7ddbd938f84ffa1fe1cf722866c7b04225752e313`；`fast_sgwb.py=94936840905af9029504e369f312e0a6330ba03ef5fc2e104da5e0c9a9be9a42`；`_wkb_transfer.py=9fcb10055b2c9ac79b76bbbc2659c07d6eda0b82213218288eaa7bc087b7845e`。Reference SHA-256 `4e2eb80709da944923206a840c0864d454faad9b51dc3fb9ebfda90abe460536`，源码未改。

## 当前源码精度核验

早期 CR0 幅度条带的 Reference 工件哈希不同，未复用其误差。当前源码同 Fast 根和相同 76 点频率数组的 z8 Reference `rtol=1e-9` 结果如下。高倾斜端点及 Sobol holdout 的根值进一步迭代到第三次 Reference 映射；低倾斜两角报告首次映射相对误差，且其 map residual 很小。

| 样本 | 根相对误差 | Fast 根处最大谱差 |
|---|---:|---:|
| `(r=1e-12, n_t=.05)` | `9.2653e-5` | `8.0209e-4 dex` |
| `(r=1e-12, n_t=.1)` | `1.6119e-4` | `8.0206e-4 dex` |
| `(r=1e-8, n_t=.05)` | `9.2667e-5` | `8.0253e-4 dex` |
| `(r=1e-8, n_t=.1)` | `1.6122e-4` | `8.0212e-4 dex` |
| independent Sobol `(r≈6.3002e-9, n_t≈.09677)` | `1.5515e-4` | `8.0203e-4 dex` |

三点第三次映射的第二到第三次闭合差最大 `3.8958e-8`；三点 Fast 根相对第三次映射最大 `1.6122e-4`。全部采样通过根误差 `2e-4` 与谱误差 `1e-3 dex` 门槛。

## 性能 A/B

同一现源码下，以旧 WKB helper 对照候选 helper，在高幅高倾斜角点与独立 Sobol 留出各进行 100 组交替完整求解。资源为 CPU affinity `[0,1]`、Numba 2 线程、BLAS 1 线程、TBB、每侧 2 次预热、`goal/PCHIP`、`physical_guard=True`。候选/基线 side-wise P50/P95 比最大为 `1.0267/1.0145`，根、频率和谱逐位一致。成对时间比 P95 有明显尖峰（`1.156/1.182`），因此保留配对尾部抖动风险；边际 P50/P95 gate 通过。首点隔离 Numba cache 冷首调用时基线/候选耗时在原始 JSON 记录，包含 JIT。

## 数值表示与分类

旧失败点的 H² 为正有限，故该片扩展跨越的是背景 primitive 的表示溢出边界。沿用平滑支恒等式 `sigma=1-(d ln H²/dN)/3` 重建 `∫sigma dN`；再加热前仍按精确 `sigma=1`，物理转折不平滑。高动态范围时用对数/相对背景尺度，不改变 CR0 的倾斜定义。

分类定向回归另行复核：精确零 `r=0,cr=0` 返回精确零谱；`cr>0,r=0` 为 `invalid_r`；无法建立一致性 cutoff 为 `invalid_cutoff`；数学 H² 非正返回 `nonpositive_h2`；非零张量最终线性谱下溢返回 `nonrepresentable_output/final_spectrum_underflow`；中间乘积下溢但 log 谱可恢复时保持有限 log 输出。相关定向测试 `12 passed`。本轮未改生产或 Reference 代码。

额外回归断言覆盖输出表示边界：本参数点 14 个 `Ogw_today` 线性零格的 `log10OmegaGW` 仍有限，范围 `[-367.439,-325.375]`，均低于最小 subnormal 的 `log10=-323.306`；`r=0,cr=0` 的精确零谱仍为全 `-inf`。新增断言和 exact-zero/final-underflow 两项复测 `3 passed`，Ruff 通过。

这是一个 CR0 二维矩形的有限采样证据，不是连续域证明。原始数据、复现脚本、Reference arrays、map closure、计时、图和哈希清单见本目录 `screen_cr0_b175_amplitude_tilt*`、`validate_cr0_b175_amplitude_tilt_z8*`、`close_reference_cr0_b175_amplitude_tilt_z8*`、`benchmark_cr0_b175_amplitude_tilt_wkb_ab_100pairs*` 与 `cr0_b175_amplitude_tilt_manifest.json`。

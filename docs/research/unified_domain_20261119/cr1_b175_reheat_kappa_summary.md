# CR1 再加热温度 × kappa10 Fast 扩域核验（`8e71fb7`）

## 参数区域与配对屏幕

固定 `Omega_bh2=Omega_ch2=7e5`、`H0=67.32117`、`DN_eff=0`、`A_s=2.100549e-9`、`r=.01`、`DN_re=10`、`cr=1`。CR1 的实际一致性倾斜为 `n_t=-0.00125`；输入 `n_t=0` 不替代模型派生值。目标矩形为 `T_re∈[2000,2500]`、`kappa10∈[.01,.015]`。全部点的 `N_inf` 有限（约 `68.46–68.81`）、H² 正有限，`physical_guard=True`。

固定旧版基线 `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45` 与当前提交 `8e71fb74cadc50c1e3493a877cd19975a47e984b` 同输入配对：3×3 网格 + 两个内部点为 `0/11→11/11`；seed `20261126` 的 8 点独立 scrambled Sobol 为 `0/8→8/8`。旧版 19 点全部 `nonfinite_background_primitive`，且旧版 H² 正有限；候选全成功、76 个输出频率节点完整、无警告。候选 screen warm solver P50 约 `95 ms`。

候选 SHA-256：`exact_background.py=c55ff5755615a0eaa5a2f4c7ddbd938f84ffa1fe1cf722866c7b04225752e313`；`fast_sgwb.py=554d197b8c2271d85bf9ffb59f89ad51809af4da75a3e34cc6fc8d1e6a91534a`；`freq_adaptive.py=ee417a2f45c5c747fce36c95a0a124fc4a21824cc8d03431735a4f77f0cca28f`。Reference SHA `4e2eb80709da944923206a840c0864d454faad9b51dc3fb9ebfda90abe460536` 未变。

## 当前源码 Reference 精度

四个新角点和一个独立 Sobol holdout 使用当前 Fast 根、同一频率数组及 z8 Reference `rtol=1e-9`：

| 点 | 根误差代理 | Fast 根处最大谱差 |
|---|---:|---:|
| `(2000,.015)` | `1.5449e-4` | `7.8393e-4 dex` |
| `(2500,.01)` | `5.5334e-5` | `8.0727e-4 dex` |
| `(2500,.015)` | `7.1279e-5` | `7.4217e-4 dex` |
| Sobol `(2372.53,.01417)` | `8.2872e-5` | `8.0335e-4 dex` |

最大误差角 `(2000,.015)` 和 Sobol 留出继续到第三次 Reference map：Fast 根相对第三次 map 最大 `1.5449e-4`，第二至第三次 map 闭合差最大 `1.0502e-8`。共享角 `(2000,.01)` 使用当前-source Fast 复放校验后复用既有 z8 Reference：当前与保存根、76 点频率数组及 Fast log 谱逐位相同；该角根差 `1.0923e-4`，谱差 `8.6962e-4 dex`，map closure `2.9228e-7`。全部样本通过根 `2e-4` 和谱 `1e-3 dex` 门槛。

## 性能 A/B 与机制

旧 WKB helper 与当前 helper 在高角点和独立 Sobol 点各进行 100 组交替完整求解。CPU affinity `[0,1]`、Numba 2 线程、BLAS 1 线程、TBB、每侧 2 次预热、`goal/PCHIP`、`physical_guard=True`。候选/旧核最大 P50/P95 比为 `0.9997/0.9887`；根、频率和谱逐位相同。首点冷调用含 JIT，基线/候选为 `18.90/18.82 s`。成对时间比 P95 仍有调度尖峰（约 `1.07`），边际 P50/P95 gate 通过，配对尾部抖动保留为风险。

旧失败来自正 H² 区域上的背景 primitive 中间表示溢出。候选使用 H² 恒等式重建平滑支的 `∫sigma dN`，再加热前保留精确 `sigma=1`，不平滑物理转折；CR1 倾斜仍由模型约束派生。

## 证据边界

这是 CR1 `T_re×kappa10` 方域的有限样本证据，不是连续域证明；生产源码与 Reference 方程没有变化。本目录配对 lattice/Sobol、当前-source Reference、第三次 map closure、base-corner replay、100 对完整 A/B、图、复现脚本和 hash manifest 提供原始证据。

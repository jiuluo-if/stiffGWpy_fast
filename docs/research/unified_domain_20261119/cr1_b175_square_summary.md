# CR1 `[1e5, 7e5]²` Fast 区域证据（`b175f0b`）

## 范围与数学有效性

本轮核验固定 `cr=1`、`H0=67.32117`、`DN_eff=0`、`T_re=2000`、`DN_re=10`、`kappa10=0.01` 及其余屏幕配置；输入 `n_t=0` 仅作为 API 输入，模型按 CR1 一致性关系派生实际倾斜 `n_t=-0.00125`。所有抽样的 `N_inf` 有限、`physical_guard=True`，背景 `H²` 全程正且有限。再加热转折按原位置保留。

旧版基线为 `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45`，候选源码提交为 `b175f0b3ffceccd2d046808921cb2dceffda1611`。候选 `exact_background.py` / `fast_sgwb.py` SHA-256 分别为 `c55ff5755615a0eaa5a2f4c7ddbd938f84ffa1fe1cf722866c7b04225752e313` / `94936840905af9029504e369f312e0a6330ba03ef5fc2e104da5e0c9a9be9a42`；Reference SHA-256 为 `4e2eb80709da944923206a840c0864d454faad9b51dc3fb9ebfda90abe460536`，生产 Reference 未改。

## 配对区域筛查

- `[1e5,7e5] × [1e5,7e5]` 的 3×3 网格及 2 个内部点：旧版 `0/11`，候选 `11/11`。
- scrambled Sobol 独立留出（seed `20261120`，8 点）：旧版 `0/8`，候选 `8/8`。
- 旧版 19 个失败全部为 `nonfinite_background_primitive`；候选无失败。旧版与候选的 `H²` 均正且有限，说明这些样本不是 `H²≤0` 的数学无效点。
- 候选区域内 76 个频率节点完整输出，无警告；Fast 求解耗时中位数约 `94.6 ms`（按单次 solver 计时，含暖缓存）。

## 精度

当前源码哈希与早期 CR1 Reference 工件不匹配，且当前 Fast 复放未逐位重现旧根/频谱；旧 Reference 数值不作本轮证明。按相同当前 Fast 根、同一频率数组新做 `z_tail=8` Reference 映射，并继续迭代到第三次映射：

| 点 | Fast 根相对第三次 Reference 映射 | Fast 根处最大频谱误差 | 第二到第三次映射相对闭合差 |
|---|---:|---:|---:|
| 高角点 `(700k,700k)` | `1.0923e-4` | `8.6962e-4 dex` | `2.9228e-7` |
| 中角点 `(300k,300k)` | `2.3328e-5` | `8.7148e-4 dex` | `8.5404e-10` |
| 独立低区留出 `(200k,250k)` | `2.1809e-7` | `8.9352e-4 dex` | `9.6192e-8` |

三点均通过根误差 `2e-4`、频谱误差 `1e-3 dex` 门槛。Reference 映射闭合差另行列出，不把一次映射差误当作收敛证明。

## 性能 A/B 与机制

在相同现源码及同一 CR1 参数下，以旧 WKB 传递 helper 对照当前带奇异基底误差界的 helper；核心频率网格、PCHIP、物理 guard 相同。资源：CPU affinity `[0,1]`、Numba 2 线程、BLAS 1 线程、TBB、每侧 2 次预热、32 组交替完整 Fast 求解；首点使用隔离的新 Numba cache 记录冷启动。高角点与独立留出点的 P50 比率为 `1.003/1.022`，P95 比率为 `0.943/1.009`，均未超过 5%；根、频率数组及频谱均逐位相同。冷启动首点旧 helper/候选分别为 `21.44/20.08 s`，均包含 JIT 编译。

本片旧失败来自背景 primitive 的中间表示范围，而抽样 `H²` 本身为正有限。已有的背景实现利用平滑支上的恒等式 `sigma = 1 - (d ln H²/dN)/3`，从物理再加热分界积分得到 `F=∫sigma dN`；再加热前仍用精确 `sigma=1`，转折不做平滑。`Phi/S2/S2inv` 由该积分构造，指数接近 binary64 表示边界时按尺度转为相对背景表示。此结果验证这套统一表示在该 CR1 方域的有限样本上有效；本轮没有改生产代码。

## 证据边界与复现

这是有限样本的区域覆盖证据，不是连续域证明。独立精度核验覆盖 3 点；没有宣称逐点 Reference 精度或数学连续域认证。旧版耗时因失败不用于效率比较；CR1 A/B 只核验 WKB helper 的候选/旧核差异。生产源码无改动，Reference 方程和门槛无改动。

原始数据、脚本与图：

- `screen_cr1_b175_square.py`、`screen_h2_identity_region.py`
- `screen_cr1_b175_square_lattice.json`、`screen_cr1_b175_square_sobol.json`
- `validate_cr1_b175_square_z8.py`、`validate_cr1_b175_square_z8.json`
- `close_reference_cr1_b175_square_z8.py`、`close_reference_cr1_b175_square_z8.json`
- `benchmark_cr1_b175_wkb_guard_ab.py`、`benchmark_cr1_b175_wkb_guard_ab.json`
- `plot_cr1_b175_square.py`、`cr1_b175_square_evidence.png`、`cr1_b175_square_manifest.json`

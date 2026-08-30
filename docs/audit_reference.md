# 高精度参考管线与物理优先证据（2026-08-30）

## 0. 目的

`fast-vs-LSODA` 的对比有一个结构性缺陷：LSODA 本身不是真值，且 fast 与 LSODA 共用同一份固定步长
`sigma` 网格（瞬时再加热边界 σ=1→4/3 的样条离散），因此两者“一致”只证明引擎一致，不证明连续极限正确。
本阶段新建了一个独立高精度参考管线（`stiffgwpy/reference.py`），用「连续 σ(N) 求值 + 高阶自适应
DOP853」作为新的精度锚点，从而区分：

- 引擎误差（fast vs 同设置 LSODA）：~1e-5，已量化；
- 共享 σ 网格偏差（两者 vs 连续 σ 参考）：本阶段首次直接量化。

## 1. 参考管线的三个物理改进

1. **σ(N)/H(N) 精确求值**：在任意 N 上按物理分支分段求 σ（reheating 区 σ=1、辐射区 σ=4/3、
   stiff 区 σ→2、再加热边界为精确阶跃），不再用固定步长网格上的三次样条。σ(N) 穿过再加热边界的
   过冲被消除。
2. **高阶自适应 ODE**：用 `scipy.integrate.solve_ivp(method='DOP853')`，`rtol=1e-11`，停在可配置的
   `z_tail`，再进入解析深亚视界尾部。
3. **求积带误差界**：对 `Omega_GW(f)=Ogw-Oj` 用保形 PCHIP + 自适应 Gauss-Kronrod 求 `∫dlnf`，
   同时给出 quadrature error 与 interpolation error。

`run_reference` 的频率集现由 `freq_adaptive.grid_independent_freqs` 生成（~242 点，freq_res=1），
对 σ 网格分辨率不变（网格无关），自洽收敛后 `DN_gw=0.0022703`，与 construct_f 网格参考（0.0022708）
一致到 ~0.02%，确认参考锚点网格无关且准确。

解析极限自检（`tests/test_reference.py`）：

| 极限 | σ 理论值 | 参考实现 |
|---|---|---|
| reheating（MD, w=0） | 1 | 1.000000 |
| 纯辐射（w=1/3） | 4/3 | 1.333200 |
| 深 stiff（w=1） | 2 | 1.999773 |

σ(N) 在远离拐点处与 fast 网格一致（<5e-3），差异集中在拐点附近 → 印证误差源是网格样条而非模型。

## 2. default 点三方对比（full grid, 246 频）

复现：`python scripts/benchmark_reference.py --point default --out docs/reference --freq-full --with-lsoda`

| 引擎 | runtime (s) | DN_gw (SGWB 贡献) | 相对 reference | κ_r |
|---|---:|---:|---:|---:|
| reference（连续 σ, DOP853, z_tail=5） | 30.8 | 0.00227081 | 0（锚点） | 0.0019381 |
| lsoda（同网格, rtol=1e-6） | 18.1 | 0.00224921 | −0.95% | 0.0019454 |
| fast（production, h=0.01） | 0.20 | 0.00224081 | −1.32% | 0.0019381 |

关键结论：**fast 与 LSODA 都低估了积分 ΔN_eff 约 1%，且方向一致**，因为两者共用同一份带 σ-kink
样条偏置的 σ 网格。reference（连续 σ）给出的 0.00227081 更接近连续极限。

这一结论与既有 h 收敛扫描相互印证：LSODA/fast 在 h→0 时 DN_gw 单调爬向 ~0.002266–0.002269
（h=0.01→0.002249，h=0.00125→0.002266），reference 的 0.002271 正是去掉最后一点网格偏置后的位置。

## 3. 参考管线自身的误差估计（default 点）

| 误差源 | 估计 | 说明 |
|---|---:|---|
| ODE（rtol 1e-10 vs 1e-11） | dex_max 8.0e-8；ΔN_eff 相对 5.1e-10 | reference 已数值收敛 |
| 解析尾部（z_tail 5 vs 7） | dex_max 3.1e-3；ΔN_eff 相对 0.35% | 共享尾部公式；z_tail=5 仍有 ~0.3% 尾部误差 |
| 求积（quad abserr） | 7.1e-15 | 可忽略 |
| 插值（半密度 PCHIP） | 2.8e-11 | 可忽略 |

> 注意：reference 在上表使用 z_tail=5，因此 reference 自身带 ~0.3% 尾部误差。若以深尾
> （z_tail≥7）为真，reference 的 ΔN_eff 应再上调 ~0.35%，即 ~0.002279。fast 相对深尾参考的
> 低估约 ~1.7%，LSODA ~1.3%。

### 3.1 匹配误差估计（深亚视界渐近解，§2.1）

`solve_reference_mode` 现在在切换到解析尾部的点返回绝热参数
`eps_handoff = |1.5σ − 1| / e^{z_handoff}` 与 `matching_error_rel`（冻结振幅近似的首阶 WKB 缺陷，
取为该 ε）。实测（default 点，f=−8）ε 与冻结尾部误差的关系：

| z_tail | ε_handoff | 相对深参考(z=10)的 dex 误差 | ODE 阶步数(nfev) |
|---|---:|---:|---:|
| 3 | 4.9e-2 | 5.6e-3 | 4.1e3 |
| 5 | 6.7e-3 | 1.6e-3 | 1.5e4 |
| 7 | 9.0e-4 | 1.7e-4 | 9.3e4 |
| 10 | 4.5e-5 | 0（锚点） | 1.8e6 |

这直接印证 §2.1 的动机：**ODE 成本随 z_tail 爆炸（z=10 单频需 180 万次函数求值），而冻结尾部在
浅 z_tail 带 O(ε) 误差**。因此更便宜的深亚视界渐近解应在更小的 z 处利用 WKB 振幅演化（而非冻结）
来消除这部分振荡积分。`matching_error_rel` 现在可被上层读取，用于估计尾部贡献的误差。

### 3.2 WKB / 长波物理一致性测试

新增 `tests/test_physics_limits.py` 直接校验 §六 要求的物理极限并通过：

- 深亚视界：模式以 dθ/dN = e^z 振荡（WKB 频率），零交叉计数与 ∫e^z dN/(2π) 一致（<10%）。
- 长波：极低频模式到 today 仍超视界（`used_tail=False`），Ogw>0、Oj<0（超视界负 Omega_j 贡献）。
- 能量一致性：`DN_gw = Neff0*g2/Omega_nu` 与 `kappa_r` 定义式自洽（rel<1e-12）。
- σ 解析极限：MD=1、辐射=4/3、stiff=2（见 `tests/test_reference.py`）。
- no-stiff / stiff 增强方向：`kappa10→0` 时高频 Ω_GW 单调下降（stiff 增强消失，恢复标准 ΛCDM
  张量谱）；低频（stiff 影响带以下）Ω_GW 对 `kappa10` 几乎不敏感（`test_physics_limits.py`）。
- 标度律：纯辐射平台 `Ω_GW ∝ f^{n_t}`（n_t=−0.4 实测斜率 −0.41）；纯 stiff 增强带 `Ω_GW ∝ f^{1}`
  （实测斜率 +1.000），符合 stiff/kination（w=1）理论预期。

### 3.3 并行化

`spectrum_reference` 用线程池并行逐频 ODE（scipy DOP853 释放 GIL）；与串行逐位一致
（`max dex diff = 0`）。线程并发增益有限（GIL 在 Python 层仍受限），但保障了 reference 可直接在
满频网格上运行。Windows 下用 `mp.Pool` 会因 spawn/re-import 停滞，故弃用进程池。

## 4. 曲率自适应频率采样（§2.4）

新增 `stiffgwpy/freq_adaptive.py`：`adapt_refine_grid` 依据 `log10 Omega_GW` 的局部二阶导数
`y''` 估计插值误差 `|y''|h²/8`，只在超过 `target_dex` 的区间插入中点，增量求解（旧点不回算），
平坦区保持稀疏，尖锐特征（谱膝、stiff 峰、高频截止）自动加密。

default 点实测（reference 为求解器，z_tail=5）：

| 网格 | 点数 | g2（∫dlnf Ω_GW） | 相对细网格 |
|---|---:|---:|---:|
| 粗均匀 | 60 | 2.7577e-8 | −1.98% |
| 自适应（target 5e-3 dex） | 203 | 2.8172e-8 | +0.13% |
| 细均匀 | 220 | 2.8135e-8 | 0（锚点） |

自适应网格在 g2 积分上与细均匀网格一致到 ~0.13%，显著优于粗均匀网格的 ~2%，且仅用与细网格
相当的点数。单测 `tests/test_freq_adaptive.py` 用合成尖锐函数验证它能把插值误差压到目标以下、
并在特征附近聚类加密、平坦区稀疏。

> **重要物理澄清（低频超视界区）**：f 低于「今日地平线」(log10(aH/2π)≈−17.5) 的极低频模式到 today
> 仍超视界，其 `Ogw−Oj` 会因超视界负 `Omega_j` 贡献而符号可变、随 log f 急剧摆动，`log10(Omega_GW)`
> 在该区插值时病态。该区对 `Delta N_eff` 积分无实质贡献（审计已确认低频尾误差 ~1e-4 量级、不进入
> 积分量）。因此**曲率自适应采样应聚焦于已进入视界的谱区（低频再入拐点以上到高频截止）的物理特征**，
> 而非病态超视界区——后者不是物理可观测的静态 `Omega_GW`。`run_reference` 的 `log10OmegaGW` 已加
> `max(..., 1e-40)` 下限，避免负值产生 NaN 污染后续对比。

### 4.1 网格无关频率采样（消除 σ 网格依赖）

新增 `freq_adaptive.grid_independent_freqs`：用**连续背景**推导 `fmax`（膨胀起点视界）、`fmin`
（今日视界）、`fcmb`（CMB pivot）并用 `freq_res` 布点，不再读网格数组 `m.f_hor`。实测：在 gen_fast
（均匀 σ 网格）与 build_transition_grid（变量 σ 网格）两种模型上，生成的频率集**逐点一致**（967 点，
maxdiff=0）。这直接解决 §7.2 指出的“频率采样被 σ 网格分辨率污染”问题——`σ` 网格改动不再悄悄移动频点。
单测 `test_grid_independent_freqs_invariant_to_sigma_grid` 已固化。下一步：把该网格无关频点集 + §4
的曲率自适应细分接入 fast/RFT 生产路径，即可在启用 transition 变步长时不再污染 `ΔN_eff`。

## 6. Production accuracy modes 与错误预算（§9/§10/§11）

### 6.1 四档模式（§9）

`fast_sgwb.ACCURACY_MODES` 现提供 `debug / fast / production / reference` 四档，并保留
`ultra-fast` 作为 `fast` 的向后兼容别名：

| 模式 | h | col_step | z_tail | freq_res | outer tol | 用途 |
|---|---|---|---|---|---|---|
| debug | 0.005 | 1 | 10 | 2.0 | 1e-8 | 最高网格精度 + 诊断 |
| reference | 0.00125 | 1 | 10 | 2.0 | 1e-8 | 最紧固定网格档（最慢） |
| production | 0.01 | 4 | 7 | 1.0 | 1e-7 | 科研默认 |
| fast | 0.02 | 8 | 5 | 1.0 | 1e-6 | 快速探索 |

### 6.2 错误预算与自动升档（§10/§11）

`fast_sgwb.estimate_error(mode)` 返回按模式校准的分阶段相对误差：
`DN_gw_error` / `spectrum_error` / `quadrature_error` / `integration_error` /
`ODE_error` / `tail_error` / `model_bias_error`。其中引擎项（ODE/quadrature）由 h 收敛得到
（~1e-5）；`model_bias` 按实测 h 收敛曲线校准（h=0.02→3.9%、h=0.01→1.3%、h=0.005→0.55%、
h=0.00125→0.22%、h→0→0.10%）。fast 求解后 `m.error_estimates` 保存该预算。

`SGWB_iter(engine='fast', accuracy_mode=..., auto_escalate=True, error_tol=…)` 在
`DN_gw_error > error_tol` 时升到 `reference` 档并记录 `m.escalated_from` /
`m.reference_evals` / `m.escalations`；Cobaya `engine_stats` 已暴露这些字段与
`Delta_Neff_GW_error` derived 参数。

### 6.3 关键诚实声明：fast 的 `reference` 档 ≠ 连续 σ 真值

实测各档 `DN_eff`：fast 0.0021853、production 0.0022408、debug 0.0022464、
fast-reference(0.00125) 0.0022561，而独立 `reference.py`（连续 σ）为 0.0022708。
也就是说 **fast 即使跑到最紧的 `reference` 档仍带 ~0.2% 连续 σ vs 固定网格偏差**；只有
`reference.py` 的连续 σ 管线才消除它，而其成本 ~30 s/点，不适合 MCMC。因此：

- `auto_escalate` 到 fast-`reference` 只能把网格偏差压到 ~0.2%，达不到 §13 的 ΔN_eff<0.01% 目标；
- 要真正达到该目标，必须在 fast 内核里改用连续 σ(N) 求值（或拐点局部加密网格），这是后续阶段的
  主攻方向；
- 生产结论必须把 `model_bias` 计入误差，而不是把它藏在一个“fast 自洽收敛”的数字后面。

## 7. fast 内核接入连续 σ（`sigma_exact`，表 1 的直接回应）

本阶段把「连续 σ(N) 求值」注入了 fast 求解器（`stiffgwpy/exact_background.py`）：用分段精确的
σ（再加热拐点作为精确断点）重算 `F=∫σ dN`、`Phi`、`S2` 振幅，而不再用固定网格 σ 的三次样条，再把
这些精确表交给现有 numba 步进内核。`SGWB_iter(..., sigma_exact=True)` 开启。

default 点（z_tail=5，与 reference 基准同口径）：

| 路径 | DN_eff | 相对连续 σ 参考 (0.0022708) |
|---|---:|---:|
| fast 固定 σ 网格 | 0.0022494 | −0.94% |
| fast `sigma_exact` | 0.0022613 | −0.42% |

即 `sigma_exact` 把 model_bias 大致减半（0.94%→0.42%），直接印证表 1 的头号误差源就是 σ-kink 网格
偏差。`exact_background` 的 F 积分精度已核对到 ~2e-5（h/2 子网格 vs h/8 gold），不是残余误差来源。

### 7.1 残余误差：transition 区域步长相位振荡（§2.3）

`sigma_exact` 在不同 h 下的收敛呈**非单调**：

| h | sigma_exact DN_eff | 相对参考 |
|---|---:|---:|
| 0.01 | 0.0022613 | −0.42% |
| 0.005 | 0.0022660 | −0.21% |
| 0.0025 | 0.0022608 | −0.44% |
| 0.00125 | 0.0022651 | −0.25% |

这不随 h 单调消失，正是既有审计发现的共享网格路径在 σ 拐点处的相位振荡（σ=1→4/3 跳变导致的
样条/步长相位伪影），只不过本阶段把它从“样条过冲”移到了“Magnus 冻结中点跨越跳变段的 O(h) 误差”。
参考（DOP853）在该处靠自适应小步长消除；固定 h 的 Magnus 内核做不到。

**结论**：要消除最后的 ~0.2–0.4%，必须实现 §2.3 的 transition 感知自适应步长：在再加热拐点附近
（`N_re_abs` 邻域）对 ODE 自动细分步长，远处保持粗步长。这是下一步主攻，也是“物理自适应求解器”
区别于“LSODA 兼容近似器”的关键分水岭。

### 7.2 transition 网格步长实验（决定性结论：未接稳）

已搭好基座（`stiffgwpy/exact_background.build_transition_grid` 在 `N_re_abs` 插入网格节点、
`exact_phi_s2_grid` 在真中点算 `Phi`/`h_arr`，`fast_sgwb.solve_kernel` 支持可选 `h_arr` 变步长），
用**网格无关频率集**（§4.1，屏蔽频率网格）把残余歧义消掉后，得到决定性数字（default 点）：

| 路径 | DN_eff | 相对自洽连续 σ 参考 (0.0022708) |
|---|---:|---:|
| fast 固定 σ 网格 | 0.0022479 | −1.01% |
| fast `sigma_exact`（均匀网格，σ 拐点被步长跨越） | 0.0022617 | −0.40% |
| fast `transition`（kink 节点 + 变步长） | 0.0022954 | **+1.08%** |
| 自洽连续 σ 参考（DOP853 自适应，2 轮收敛） | 0.0022708 | 0（锚点） |

**决定性结论**：self-consistent 参考 = `0.0022708`（与单次值一致，因 DN_eff 小、反馈弱）。真相在
`sigma_exact` 与 `transition` 之间，**transition 的 kink 消解步长过校正**（+1.08%），而
`sigma_exact` 欠校正（−0.40%）。这说明**瞬时再加热的 σ 间断使二阶 Magnus（中点）方法无论网格相位如何
都带 ~0.4–1% 的系统误差，且随 h 非单调**——只有高阶自适应 ODE（DOP853）能正确解析该间断。

因此 `transition_refine` 显式抛 `NotImplementedError`，不作为生产路径；`sigma_exact` 是已验证的
可用改进（偏差 0.94%→0.40%）。**要真正突破到 <0.1%，需要不是“换个网格”而是“换个积分器/更高阶
Magnus”**，或接受该残差并在 production 的 `model_bias` 中如实报告（当前 `ERROR_BUDGET` 已把
production 的 model_bias 标为 1.3e-2，覆盖这一残差）。

补充负结果：把冻结中点 `z_mid` 换成端点平均 `z_avg=(z_k+z_{k+1})/2` 对结果几乎无影响
（sigma_exact：−0.4194% vs −0.4198%）。说明残差不是“步内 z 取值点”的选择，而是冻结-z Magnus
旋转矩阵在跨 σ 间断步处的**固有**近似——无论网格相位/步内 z 取法都带 ~0.4% 系统偏差。只有高阶自适应
ODE（reference 的 DOP853）或更高阶 Magnus 才能消除。

再补决定性负结果：实现**真正**的 transition 感知网格（`build_kink_refined_grid`：在 `N_re_abs`
±2 e-fold 内步长加密到 h/8、kink 为精确节点），配变步长 `h_arr` 重跑，得
`DN_eff=0.0022959`（+1.11%），与原 transition（+1.08%）几乎相同——**kink 网格细化并未减小残差**。
即 fast 的冻结-z Magnus 对进入 stiff 时代的高频模给错振幅，**与 kink 网格分辨率无关**；用 x 轴
分辨率无论怎么加，都停在 ~+1%。这判定了 §2.3 的“transition 感知网格”对当前 Magnus 方法不解决问题，
要达成 <0.1% 必须换成**高阶自适应 ODE**（不是加网格点）。

再补两个负结果（进一步收窄根因）：
- **深亚视界自适应子步**（`omega_crit`：ω>阈值时把步长细分到 ω≈阈值）对结果几乎无影响
  （sigma_exact −0.418% ~ −0.422%），说明残差并非深亚视界冻结-z 相位误差。
- **精确 z_tail 事件点切换尾部**：单频复刻把 f=−8 的偏差从 0.006 dex 降到 0.0005 dex，但**聚合
  DN_eff 反而变差**（sigma_exact −0.42% → −1.42%），单频与聚合行为矛盾——说明该修正在内核其它
  环节不一致，而非简单的手工修正。

综合所有负结果（z_avg、kink 节点、kink 细化、深亚视界子步、精确事件切换），fast 的冻结-z Magnus
在 σ 间断模型的聚合 ΔN_eff 上稳定停在 ~−0.4% 到 +1.1%，**无法通过任何网格/步长/切换点调节收敛到
连续 σ 真值**。结论不变且被充分验证：要达成 <0.1% 必须替换积分方法（高阶自适应 ODE），否则以
`sigma_exact`（−0.40%）为最优 fast 精度并在 `model_bias` 中如实报告。

### 7.3 真正的主导误差源：今日网格锚点量化（重大修复）

单频原型与生产内核的 ~10x 偏差最终指向**今日锚点**。模型网格用
`present = floor(N_inf/h)*h ≈ 68.81` 作为“今日”，而 `reference.py` 用连续 `N_inf = 68.81384`，
相差 `δN = N_inf − present ≈ 0.0038` e-fold。解析尾部幅值按 `exp(−δN)` 衰减，
`exp(−0.0038) ≈ 0.996` → 恰为 **~0.4% 系统偏差**。即：**真正的主导误差源是网格把今日锚点量化了**，
不是 σ-kink 样条，也不是 Magnus 阶数。

修复：`gen_fast` 与 `gen_expansion` 都把 `Nv` 末端锚定到连续 `N_inf`（真·今日），并相应更新
`index_re`。修复后（default 点，z_tail=5）：

| 配置 | DN_eff | 相对连续 σ 参考 (0.0022708) |
|---|---:|---:|
| fast grid 锚点前（h=0.01） | 0.0022494 / 0.0022408 | −0.94% / −1.32% |
| fast grid 锚点后（h=0.01，production） | 0.0022667 | **−0.18%** |
| fast grid 锚点后（h=0.00125，reference 档） | 0.0022722 | **+0.06%** |

**h 收敛的非单调幅度从 ~±1% 缩到 ~±0.2%**，且 reference 档（h=0.00125）达到 **+0.06% < 1e-3**。
这是不依赖 benchmark 调参的真实物理误差源修复（今日锚点量化）。剩余 ~±0.2% 的非单调性仍来自
σ-kink 的网格相位，需要 transition-aware 方法进一步压制。

### 7.4 transition-aware kink 细化（现已是 production 默认）

把 §7.3 的“今日锚点”与 transition-aware σ-kink（`build_kink_refined_grid` + `exact_phi_s2_grid`
变步长）结合（并将 `sigma_vec`/`H2_vec` 向量化以提速），得到（default 点，z_tail=5）：

| 路径 | DN_eff | 相对连续 σ 参考 (0.0022708) | 运行时间 |
|---|---:|---:|---:|
| 纯 grid（plain, h=0.01） | 0.0022667 | −0.18% | ~5 ms |
| transition-refine（production, h=0.01） | 0.0022698 | **−0.04%** | ~80 ms |

即 `transition-refine`（kink 节点 + 邻域细化 + 精确 Φ/S2 + 变步长）叠加今日锚点修复后达
**−0.04%（<1e-3）**，且 ~80 ms/点（相对 LSODA ~17 s 为 ~200×）。因此 `production` / `reference` /
`debug` 档已默认启用 `transition_refine`；`fast` 档保留纯 grid（最快、−0.18%）。

### 7.5 参数空间验证（§9，z_tail=5 同口径）

用 fast（`transition_refine`）与**同 z_tail=5 的连续 σ 自洽参考**对比：

| 参数点 | 参考 DN_gw | fast DN_eff | 相对参考 |
|---|---:|---:|---:|
| default（kappa10=1e-2） | 0.00227081 | 0.0022698 | **−0.04%** |
| stiff（kappa10=1.0） | 0.226372 | 0.226282 | **−0.04%** |

即求解器（`transition_refine`）在**同一 z_tail 口径**下对 moderate（default）与 strong（stiff）
stiff 都达 **−0.04%（<1e-3）**，说明不是 default 点特例。快速参数扫描（default/stiff/lowT/
radDominant/cr0_nt/tinyR/edge/norad/highT）均产出有限物理 ΔN_eff，`highT` 由共享 `DN_eff>5` 护栏拒绝。

> 注意：`production` 档配置 `z_tail=7`（更深的解析尾部，更接近连续极限），但与上述 `z_tail=5`
> 参考比较时会看到 ~0.35% 的差异——这是**尾部阈值（z_tail）**这一共享近似，不是求解器数值误差；
> 求解器在匹配的 z_tail 下为 −0.04%。

### 7.7 runtime vs physical-error Pareto（§13，default 点，z_tail=5 同口径）

`scripts/benchmark_pareto.py` 输出（参考为连续 σ 自洽锚点）：

| 引擎 | runtime | ΔN_eff 相对参考 |
|---|---:|---:|
| reference（连续 σ DOP853，锚点） | 171.7 s | 0 |
| fast transition-refine（production） | ~0.08–0.28 s | **−0.022%** |
| fast plain grid | 0.012 s | −0.161% |
| LSODA | 18.6 s | −0.168% |

关键结论：**fast transition-refine 在 runtime 与物理精度两个维度都优于 LSODA**（更快且更接近
连续 σ 参考），同时比 reference 快 ~600×。fast plain grid 最快（0.012 s）但 −0.16%。
数据存 `docs/reference/pareto_default.json`。

### 7.6 单频原型反例（提示根因在“内核 vs 单频”，而非数值方法）

用**独立单频原型**（连续精确 Φ(N)、2 阶 Magnus、精确 z_tail 尾部、深超视界起点）对照 reference：
单频误差仅 **~1.4e-4 到 8.3e-4 dex**（0.03%~0.19%），显著优于生产内核的 ~6e-3 dex (~1.4%)。且该
误差对**起点深度（z=-6.9 vs -12）、N_re 分段、4 阶 Magnus（mag2≈mag4）、细 Φ 子网格（h/8 vs
h/2）** 均几乎不敏感。即：**单频上 <1e-3 是可达的，但生产 numba 内核的聚合卡在 ~−0.4%，说明
根因不是数值方法/阶数/网格，而是“内核逐模实现 vs 单频原型”之间存在 ~10x 的偏差**（如 Φ/S2 表、
j0/z0 起始、尾部锚点的一致性），需定位修复。

## 4. fast 频谱的逐点误差（reference 为锚，频率子集）

`fast vs reference` 的 log10(Omega_GW) 绝对误差（dex，30 点频率子集）：

| 统计 | 值 |
|---|---:|
| dex_max | 0.156 |
| dex_p95 | 0.073 |
| dex_p50 | 0.0052 |
| 线性 Ω（信号区）相对 max | 0.432 |

最大误差集中在**最低频陡尾**（f≈10^-18 Hz），这与既有审计发现的 freq_res=1.0 低频欠采样一致：
reference 用连续 σ 逐点求解，暴露了 fast 在低频尾的 ~0.15 dex 偏差。良分辨区（中高频）fast 对
reference 仅 ~5e-3 dex。

## 5. 结论与后续

- **“fast≈LSODA”到 1e-5 是真的，但不足以作为精度声明**：两者共享 ~1% 的 σ-kink 网格偏差。
- **真实精度提升的第一优先项**是消除该共享偏差（连续 σ / 拐点局部加密），其次是提高 z_tail 以削减
  尾部 0.3% 误差、以及自适应频率采样解决低频陡尾 0.15 dex 欠采样。
- 参考管线自身 ODE 已收敛到 ~1e-7 dex，可作为 accuracy/reference 模式的基线；其成本约 30 s/点
  （246 频，无并行），适合 benchmark/病态点/收敛认证，不适合 MCMC 热路径。
- **尚未完成**：reference 的充分自洽收敛（本表用 fast 收敛的 DN_eff 做单次背景）、deep-tail
  （z_tail≥10）参考、以及低频陡尾的自适应频率采样。后续阶段将做 WKB 深亚视界渐近解以把 reference
  成本降到可做密集频谱，并用曲率感知采样消除低频尾误差。

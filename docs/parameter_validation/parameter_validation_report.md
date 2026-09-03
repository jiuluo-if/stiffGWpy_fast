# parameter_validation_report — fast vs 连续-σ reference 参数空间验证矩阵

> 生成日期：2026-09-03　git commit：`5125a5c`
>
> **本报告未重跑任何物理计算**：所有数字均回读自已提交的验证产物（见 §7 源文件），与认证运行完全一致，可仅凭仓库复现。配套机器可读文件：`docs/parameter_validation/validation_results.json`（逐点结构化）与`docs/parameter_validation/validation_results.csv`（逐点平表）。

## 0. 覆盖矩阵（引擎 / 层级 / oracle 对照）

| 层级 | tier | engine | 点数 | 参考对照 | 说明 |
| --- | --- | --- | --- | --- | --- |
| A_single_point | production | fast vs reference | 9 | 9/9 全谱 oracle（matched z8 grid_independent） | default/stiff/lowT/highT/rad_dominant/tiny_r/transition/cr0_blue/extreme |
| B_sobol240 | production | fast only | 240 | 无全谱 oracle（成本）；Layer C 提供 240 点 x 11 bin 参考对照 | 240 点 Sobol（r/n_t/cr/T_re/DN_re/kappa10），本地 a-posteriori 误差预算 |
| C_posterior_bulk | production | fast vs reference | 240 | 240/240 点 x 11 likelihood bin | likelihood bin 为原生求解节点；IS 后验 + e^ΔlogL 重加权 |
| A2_param_edges | production | fast-vs-reference | 16 | 13/16 全谱 oracle（matched grid_independent z8） | 参数轴边界（u~0.02/0.98）+ transition 敏感内部点 |
| P_plain9 | plain-grid | fast-plain-grid-vs-reference | 9 | 9/9 全谱 oracle（plain-grid 原生节点，matched z8） | plain-grid(fast) 引擎误差边界量化；1e-3 science gate 不满足 -> escalation 到 production/reference |
| anchor_default | plain-grid / production / reference | multi-engine | 3 | default 点同文件对照（commit 59563da9b527d80e26bd37f813b957f22ed2a168 时代设置） | plain-grid coarser 探索档 anchor；低于 1e-3 物理门槛，不作生产认证 |

## 1. 逐点分类统计与分类规则

状态编码（与需求 §6 一致）：`PASS` / `WARN` / `FAIL` / `PHYSICAL_INVALID` / `NUMERICAL_FAILURE`。分类规则：

- `PASS`：求解成功且该行门槛满足（Layer C 门槛 = 逐 bin dex < 1e-3 vs reference）。
- `WARN`：可运行但带明确警示（plain-grid 探索档 9-corner 边界：signal rel max 7.0e-2、DN rel abs med 9.1e-3，matched z8，远高于 1e-3 science gate）。
- `FAIL`：求解成功但认证门槛未达到（当前记录集中没有此类逐点行；集成 ΔNeff 门槛见 §5）。
- `PHYSICAL_INVALID`：显式物理/自洽拒绝 `shared_Neff_guard`（极端 r/DN_re/kappa10 角落），**不算 numerical failure**。
- `NUMERICAL_FAILURE`：异常 / 非有限 / 迭代失败（本批产物中 0 行）。

各层级计数：

| 层级 | PASS | WARN | FAIL | PHYSICAL_INVALID | NUMERICAL_FAILURE |
| --- | --- | --- | --- | --- | --- |
| A_single_point | 9 | 0 | 0 | 0 | 0 |
| A2_param_edges | 13 | 0 | 1 | 2 | 0 |
| B_sobol240 | 212 | 0 | 0 | 28 | 0 |
| C_posterior_bulk | 240 | 0 | 0 | 0 | 0 |
| P_plain9 | 0 | 9 | 0 | 0 | 0 |
| anchor_default | 2 | 1 | 0 | 0 | 0 |

## 2. Layer A — 9 个 matched z8 单点（production vs reference）

| label | r | cr | T_re | DN_re | kappa10 | signal rel max | transition rel max | DN_gw rel | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| default | 1.000e-02 | 1.000e+00 | 2.000e+03 | - | 1.000e-02 | 5.982e-04 | 5.250e-04 | 3.009e-04 | PASS |
| stiff | 1.000e-02 | 1.000e+00 | 2.000e+03 | - | 1.000e+00 | 6.749e-04 | 6.749e-04 | -4.323e-04 | PASS |
| lowT | 1.000e-02 | 1.000e+00 | 1.000e+01 | - | 1.000e-02 | 6.023e-04 | 5.897e-04 | 1.464e-03 | PASS |
| highT | 1.000e-02 | 1.000e+00 | 1.000e+04 | - | 1.000e-02 | 6.917e-04 | 6.917e-04 | 1.271e-04 | PASS |
| rad_dominant | 1.000e-02 | 1.000e+00 | 2.000e+03 | - | 1.000e-06 | 6.430e-04 | 6.430e-04 | 7.181e-04 | PASS |
| tiny_r | 1.000e-06 | 1.000e+00 | 2.000e+03 | - | 1.000e-02 | 7.093e-04 | 7.093e-04 | 3.819e-05 | PASS |
| transition | 5.000e-03 | 1.000e+00 | 5.000e+02 | - | 1.000e-01 | 6.636e-04 | 6.366e-04 | -7.559e-04 | PASS |
| cr0_blue | 1.000e-02 | 0.000e+00 | 1.000e+03 | 5.000e+00 | 1.000e-03 | 6.576e-04 | 6.576e-04 | 6.097e-04 | PASS |
| extreme | 3.000e-02 | 1.000e+00 | 1.000e+04 | - | 1.000e-01 | 6.921e-04 | 6.746e-04 | -4.340e-04 | PASS |

汇总：signal/transition 带 ΩGW 相对误差 max **7.093e-04**（门槛 <1e-3 → PASS）；集成 ΔNeff 相对误差 |DN| median **4.340e-04** / p95 **1.181e-03** / max **1.464e-03**（门槛 <1e-4 → FAIL，见 §5，诚实的架构极限）。

## 2b. Plain-grid tier — 9 个 matched z8 角落（fast plain-grid vs reference）

Plain-grid 引擎（`accuracy_mode='fast'`：h=0.02 / col_step=8 / 无 transition_refine / phase_max=0）在自身 construct 频率节点上与连续-sigma reference（z_tail=8, rtol=1e-9）逐点对照；reference 直接在 plain-grid 节点上求解，残差纯为引擎误差（无频率网格插值项）。

| label | signal rel max | transition rel max | DN_gw rel | classification |
| --- | --- | --- | --- | --- |
| stiff | 1.768e-02 | 1.768e-02 | 2.061e-03 | WARN |
| default | 1.867e-02 | 1.620e-02 | -7.345e-03 | WARN |
| lowT | 4.786e-02 | 3.069e-02 | 9.265e-04 | WARN |
| cr0_blue | 1.736e-02 | 1.736e-02 | -9.142e-03 | WARN |
| extreme | 1.718e-02 | 1.718e-02 | -1.017e-02 | WARN |
| highT | 1.660e-02 | 1.612e-02 | -7.681e-03 | WARN |
| tiny_r | 7.019e-02 | 2.609e-02 | -2.725e-02 | WARN |
| transition | 2.440e-02 | 2.440e-02 | -1.705e-02 | WARN |
| rad_dominant | 6.751e-02 | 6.751e-02 | -2.664e-02 | WARN |

明确精度包络（exploratory tier 边界）：signal 带 rel max **7.019e-02**（median 1.867e-02）、transition 带 rel max **6.751e-02**；集成 ΔNeff rel abs median **9.142e-03** / max **2.725e-02**；该验证 artifact 的 fast runtime median **0.76 s/点**（JIT 优化前的历史测量，不能作为当前性能口径；reference 中位 383 s/点）。1e-3 science gate 不满足 -> 该档仅用于探索；科学结论必须 escalation 到 production/reference（adapter 已实现 likelihood-aware auto_escalate，无 silent fallback）。

## 3. Layer B — 240 点 Sobol（production fast-only）

212/240 `ok`（PASS），28/240 显式 `shared_Neff_guard`（PHYSICAL_INVALID），0 numerical failure。
fast 遥测（212 ok 点；JIT 优化前的 validation artifact）：runtime median 5.34 s / p95 8.82 s / max 14.22 s；adaptive 频率网格节点 median 236；WKB handoff defect `handoff_eps` median 6.841e-04；本地估计 ΔNeff 相对误差 median 3.997e-04（DN→0 处饱和为绝对误差 ≤1e-5，物理不可观测）。当前 runtime 见 `docs/performance_comparison_20260903.md`。

## 4. Layer C — posterior-bulk（fast vs reference，240 点 x 11 bin）+ IS 后验

fast vs reference 逐 bin dex（11 个 likelihood bin 全为原生求解节点）：max **3.102e-04**，p95 **3.016e-04**，median **2.652e-04**（240/240 点 PASS <1e-3）。

- |ΔlogL| posterior bulk：max **7.303e-03**、p95 **4.695e-03**、mean **-4.622e-04**（n=240）
- IS 后验（9000 production draws, seed 20260903）：ESS **4167.4**（≥2000 PASS）。
- e^{ΔlogL} 重加权后验位移：log10 r **-0.0011 σ**、n_t **+0.0002 σ**（<0.1σ PASS；n_t 在 cr=1 下 prior-dominated，仅记录不作认证）。

## 5. 验收门槛（显式 PASS / FAIL / NOT YET VERIFIED）

> 门槛与数字**没有**为达标而调整；未达标的条目按实际极限如实报告。

| 门槛 | 状态 | 实测 | 证据 |
| --- | --- | --- | --- |
| signal-region Omega_GW rel err < 1e-3 (production vs reference, 9 matched z8 singles) | PASS | rel max 7.093e-04 | docs/paramsweep_z8/reference_points.jsonl |
| transition-region Omega_GW rel err < 1e-3 (same 9 points) | PASS | rel max 7.093e-04 | docs/paramsweep_z8/reference_points.jsonl |
| integrated Delta_Neff rel err < 1e-4 (production vs reference, matched z8) | FAIL | median 4.340e-04, p95 1.181e-03, max 1.464e-03 (lowT DN-of-DN; abs 5.5e-10); deep-oracle default -2.94e-4 | honest architecture limit ~3e-4..7.6e-4: frozen-z Magnus + grid envelope, not tuning-removable |
| posterior-bulk per-bin log10 Omega dex < 1e-3 (240 points x 11 bins, fast vs reference) | PASS | dex max 3.10e-4 | docs/mcmc_posterior/is_report.json + is_pointwise.json |
| posterior-bulk |Delta logL| < 0.1 | PASS | max 7.30e-3, mean -4.62e-4 (n=240) | docs/mcmc_posterior/is_report.json |
| posterior ESS >= 2000 (importance-sampled fast chain) | PASS | ESS 4167.4 (9000 production draws, seed 20260903) | docs/mcmc_posterior/is_report.json |
| posterior parameter shift < 0.1 sigma (fast vs reference-consistent posterior) | PASS | log10r -0.00110 sigma; n_t +0.00023 sigma (n_t prior-dominated under cr=1) | docs/mcmc_posterior/is_report.json |
| analytic limits + energy/scaling consistency | PASS | green | tests/test_physics_limits.py |
| production runtime >= 100x vs LSODA (matched-accuracy setting) | HISTORICAL | 旧 validation artifact 为 ~4.5x；JIT 后当前执行层 speedup 见性能报告，物理精度结论未因此改变 | docs/performance_comparison_20260903.md |
| fallback / escalation traceable; no silent fallback | PASS | FAST/FAST_ESCALATED/REFERENCE/LSODA_FALLBACK statuses + engine_stats; shared_Neff_guard explicit | tests/test_engine.py, tests/test_cobaya_adapter.py, docs/audit_acceptance.md |
| plain-grid tier: 9-corner accuracy boundary vs reference (matched z8, plain-grid own nodes) | FAIL | signal rel max 7.019e-02, transition rel max 6.751e-02; DN rel abs median 9.142e-03 / max 2.725e-02; runtime 0.76 s/pt 为 JIT 前历史值 | docs/paramsweep_plain/validation_summary.json (exploratory tier; science gate 1e-3 -> escalation) |
| plain-grid full parameter-space sweep vs oracle | NOT YET VERIFIED | 9 matched corners only; no 240-point plain-grid oracle sweep | docs/paramsweep_plain/plain_points.jsonl |
| extended oracle spots (16 axis-edge/interior, matched z8): signal/transition rel < 1e-3 | FAIL | 13 ok-PASS, 2 PHYSICAL_INVALID(guard), 1 gate-FAIL; signal rel max 1.641e-03; DN rel abs median 5.230e-04 | docs/paramsweep_z8b/reference_points.jsonl + validation_summary.json |
| production full-spectrum oracle coverage over the 240 Sobol points | NOT YET VERIFIED | oracle (360 s/pt) run for 9 singles + 240 posterior-bulk points at 11 bins only | docs/paramsweep_ref/fast_sweep.jsonl is fast-only |
| converged real-Cobaya 3-chain MCMC (plain / production / reference; KS/Wasserstein/KL/covariance, R-1) | NOT YET VERIFIED | bounded scaffold chains (~30 rows, plumbing only); certified Layer C = IS chain + exact e^dlogL reweighting | docs/mcmc_posterior/posterior_validation.md |
| per-parameter 1D scans (low/fid/mid/high/boundary/transition for every cosmology+inflation+reheating+stiff param) | NOT YET VERIFIED | 9 physics-corner singles + 240 Sobol cover regimes; no per-parameter 1D grid artifacts committed | docs/paramsweep_z8/reference_points.jsonl |

## 6. 覆盖边界 — 显式标记的 NOT YET VERIFIED / 诚实极限

以下条目按需求 §13 显式报告为 `NOT YET VERIFIED` / `FAIL` / `WARN`，均附原因与复现成本；没有为达标而移动任何门槛：

| 条目 | 状态 | 原因与成本 |
| --- | --- | --- |
| production 集成 Delta_Neff rel < 1e-4 | FAIL | frozen-z Magnus + z_tail/网格架构残差 ~3e-4..7.6e-4（deep-oracle default -2.94e-4 仍 >1e-4）；非步长或调参可消除，需换高阶/自适应 ODE 内核 |
| production matched-accuracy runtime >= 100x vs LSODA | HISTORICAL | 旧 z8 artifact 为 ~4.5x（4.1-5.3 s/点 vs LSODA 18.56 s）；JIT 后当前执行层对比见性能报告，matched-reference 精度边界仍为 signal rel max 7.0e-2 / DN rel abs med 9.1e-3 |
| plain-grid tier 全参数空间 vs oracle 扫描 | NOT YET VERIFIED | 9-corner matched z8 边界已量化（signal rel max 7.0e-2、DN rel abs med 9.1e-3，docs/paramsweep_plain/）；240 点 plain-grid 全谱 oracle 扫描未跑 |
| 240 Sobol 点全谱 oracle 对照 | NOT YET VERIFIED | reference ~360 s/点 -> 240 点约 24 CPU.h；当前 oracle 覆盖 9 singles 全谱 + 240 点 x 11 bin |
| 收敛的 real-Cobaya 三链 MCMC（plain/production/reference，KS/Wasserstein/KL/covariance、R-1） | NOT YET VERIFIED | reference 链 ~350-935 s/点；bounded scaffold 链仅验证 adapter plumbing（~30 行，未收敛）；已认证替代 = IS 后验 + e^{Delta logL} 精确重加权（ESS 4167） |
| 逐参数 1D 扫描网格（每个参数 low/fid/mid/high/boundary/transition） | NOT YET VERIFIED | 9 个物理 corner singles + 240 Sobol 覆盖 regime；逐参数 1D 网格产物未提交 |
| sub-horizon-today 极低频区（Ogw-Oj 变号附近） | WARN | 静态 Omega_GW 在该区物理定义受限且无 Delta_Neff 权重；文档化，非求解器数值误差 |

## 7. 复现与源文件

本矩阵由以下已提交产物生成（只读回放，无物理重算）：

- docs/paramsweep_z8/reference_points.jsonl
- docs/paramsweep_ref/fast_sweep.jsonl
- docs/mcmc_posterior/is_pointwise.json
- docs/mcmc_posterior/is_report.json
- docs/benchmark_current_fast_preset_20260903.json
- docs/paramsweep_plain/plain_points.jsonl
- docs/paramsweep_plain/validation_summary.json
- docs/paramsweep_z8b/reference_points.jsonl
- docs/paramsweep_z8b/validation_summary.json

重新生成：

    python scripts/build_validation_matrix.py

重跑底层物理验证的驱动：Layer A = `python scripts/validate_fast_vs_reference.py`；Layer C = `python scripts/importance_posterior.py --help`（draw/posterior/pointwise/report 阶段）；回归测试 = `python -m pytest tests/`（6 个 slow 门槛测试用 `-m slow` 加入）。更完整清单见 README `Reproduce the benchmark / validation`。

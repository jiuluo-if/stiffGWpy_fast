# 验收审计：逐条 VERIFIED / PARTIALLY VERIFIED / NOT VERIFIED

本文件把用户粘贴的 14 节审计/优化要求逐条对到当前仓库的**权威证据**，并按 §13 的要求如实标注
`VERIFIED` / `PARTIALLY VERIFIED` / `NOT VERIFIED`。不做“让 benchmark 数字更好看”的调参；凡声明都
指向可复现的文件/数字。

## 1. 重建 SGWB 误差预算 — VERIFIED

`docs/audit_error_budget.md`：把 inflationary tensor spectrum → … → ΔN_eff 全链路拆成 10 类误差
（truncation / ODE / quadrature / interpolation / tail / freq-grid / thermal / float / boundary /
model），逐项给出量级与缺口；并量化出**共享 σ 网格偏差是主导源**（fast 与 LSODA 同为约 −1%）。

## 2. 从物理结构优化算法 — PARTIALLY VERIFIED

* 深亚视界 WKB 渐近解：PARTIAL。`reference.solve_reference_mode` 已返回绝热参数
  `eps_handoff=|1.5σ−1|/e^z` 与 `matching_error_rel`；WKB 频率 `dθ/dN=e^z` 已由测试验证。但
  **WKB/渐近解尚未在 fast 内核落地**，停止振荡积分未实现。
* horizon-crossing 自适应步长：NOT DONE。当前 fast 用固定 `h` 步长；自适应未实现。
* stiff→radiation transition 感知网格：NOT DONE（见 §7.2：kink-resolved 变步长过校正，需更高阶
  积分器而非换网格）。
* 频率空间自适应：PARTIAL。`freq_adaptive.adapt_refine_grid`（曲率 `|y''|h²/8` 打分）与
  `grid_independent_freqs`（网格无关频点集、对 σ 网格不变）已实现并测试，但**未接入 fast 生产路径**。

## 3. 重设计求积 — PARTIALLY VERIFIED

`reference.integrate_spectrum` 用保形 PCHIP + 自适应 Gauss-Kronrod，并返回 `quadrature_error` 与
`interpolation_error`（default 点：quad ~7e-15、interp ~3e-11）。但 fast 路径仍用固定 Simpson
（`Wmat`）；fast 未返回 quadrature error estimate（reference 有，见 `estimate_error` 校准表）。

## 4. 物理收敛测试 — PARTIALLY VERIFIED

已有 h 收敛曲线（`docs/audit_phase2.md`）与 reference 的 ODE 收敛（rtol 1e-10 vs 1e-11 差 ~1e-7 dex）、
尾部灵敏度（z_tail 5→7 差 ~0.003 相对）、`freq_res` 灵敏度。但**相对连续 σ 真值**的收敛序
（h/freq/tail/quadrature tolerance 逐项 slope）尚未完整建立；且
`sigma_exact` 的模型偏差随 h 呈非单调（0.40%→0.21%→0.44%→0.25%），说明现阶段无干净二阶收敛。

## 5. 高精度参考（不把 LSODA 当真值）— VERIFIED

`stiffgwpy/reference.py` + `scripts/benchmark_reference.py`：连续 σ(N) + DOP853 + PCHIP/quad，
作为新锚点。default 点给出 `DN_eff=0.0022708`（self-consistent，2 轮收敛），并揭示 fast（−1.32%）
与 LSODA（−0.95%）都低估——**以参考为准**。`docs/audit_reference.md`、`docs/reference/`。

## 6. 物理守恒与解析极限 — VERIFIED

`tests/test_physics_limits.py` 已覆盖：WKB 深亚视界频率、长波超视界冻结、能量一致性
（`DN_gw=Neff0*g2/Ω_ν`、`κ_r`）、σ 解析极限（MD=1/辐射=4/3/stiff=2）、浮点稳健、stiff 增强/no-stiff
方向、`sigma_exact` 偏差下降、**纯辐射平台区 `Ω_GW ∝ f^{n_t}`**（实测斜率 −0.41 vs n_t=−0.4）、
**纯 stiff 增强带 `Ω_GW ∝ f^{1}`**（实测斜率 +1.000）。能量/标度一致性全部通过。

## 7. 浮点/取消误差 — VERIFIED

极端参数压力测试无 NaN/inf/下溢（Ω 最小 ~1e-30），`S2` 重标定避免大数相减与指数爆炸，
`log10OmegaGW` 加 `max(...,1e-40)` 下限，极端 stiff 由共享 `DN_eff>5` 护栏拒绝。回归测试
`test_float_robustness_extreme_params`、`docs/audit_error_budget.md` §5。

## 8. 速度优化建立在物理算法之上 — PARTIALLY VERIFIED

fast 在 MCMC 单点 ~5 ms（生产档），比 LSODA 快 ~1000×（warm）。但**精度上限被 σ 间断的二阶
Magnus 残差（~0.4–1%）限制**（§7.2），因此“accuracy↑”尚未在 fast 上完全达成——速度提升来自
重标定 + 固定步长 Magnus + Numba，而非物理结构（WKB+自适应）的落地。

> 更新（2026-08-30）：定位并修复**今日网格锚点量化**这一真正的主导误差源（§7.3）。网格把今日锚点
> 量化到 `floor(N_inf/h)*h`，与 `reference.py` 的连续 `N_inf` 差 ~0.0038 e-fold，恰造成 ~0.4%
> 系统偏差。锚定到连续 `N_inf` 后：production（h=0.01）误差 **−0.18%**、reference 档
> （h=0.00125）**+0.06%（<1e-3）**；h 收敛非单调幅度从 ~±1% 降到 ~±0.2%。这是真实物理误差源修复，
> 非 benchmark 调参。

## 9. 四档 accuracy modes — VERIFIED

`fast_sgwb.ACCURACY_MODES`：`debug/fast/production/reference`（`ultra-fast` 为兼容别名），
`apply_accuracy_mode` 校验并应用；单测 `tests/test_modes.py`。

## 10. production 返回 error estimates — VERIFIED

`fast_sgwb.estimate_error(mode)` 返回 `DN_gw_error/spectrum_error/quadrature_error/
integration_error/ODE_error/tail_error/model_bias_error`，`SGWB_iter(..., auto_escalate=True,
error_tol=…)` 超限升到 `reference` 并记录；`m.error_estimates` 挂模型。Cobaya derived 增
`Delta_Neff_GW_error`。注意：error budget 的 `model_bias` 是**按 h 收敛曲线校准**（覆盖 §7.2 的
~0.4–1% 残差），非逐点自适应估计。

## 11. Cobaya integration — PARTIALLY VERIFIED

支持 `engine: fast|lsoda`、`accuracy_mode`、`fallback`、`auto_escalate`、`error_tol`，`engine_stats`
暴露 fast/reference evals、escalations、failures、fallback fraction 与 error estimates。**未接入
`engine: reference`（连续 σ 的 reference.py）作为可直接调用的引擎**；生产→参考的错误驱动回退目前升到
fast-`reference` 档（网格仍有 ~0.2% 偏差），而非连续 σ 的 `reference.py`。

> 更新：`engine: reference` 已接入（`stiffgwpy/cobaya/stiffGW.py` 用连续 σ 的 `reference.py`
> 管线，`reference_rtol`/`reference_z_tail` 可配，`m.reference_evals` 计数，单测
> `test_reference_engine_sets_state` 通过）。生产→参考的错误驱动回退仍升到 fast-`reference` 档
> （网格仍带 ~0.2% 偏差）；Cobaya 直接以 `engine: reference` 可拿到连续 σ 高精度点，供 §13 的
> “后验 vs 高精度参考”对照。
> 再更新：`reference.apply_reference_to_model` 抽出为公共 helper，`LCDM_SG.SGWB_iter(engine=
> 'reference')` 与 Cobaya 适配器都复用它（去重），单测覆盖（`test_apply_reference_to_model`）。
> 因此连续 σ 高精度引擎现在既可在 Cobaya 用，也可在原生 `SGWB_iter` API 直接调用，`auto_escalate`
> 若要升到它只差一个配置开关。
> 再更新：`SGWB_iter(..., auto_escalate=True, escalate_to_reference=True)` 现可在误差超限时升到
> **真正的连续 σ `reference` 引擎**（`engine='reference'`，reference.py），而非 fast 网格档；升档后
> `DN_gw_error/spectrum_error` 置 0（自洽高精度），`m.escalations/escalated_from/reference_evals`
> 记录。单测 `test_auto_escalate_to_reference_engine` 通过。这补齐 §11 的
> “production→误差大→reference 计算→accept”完整回退链。

## 12. 四维 benchmark — VERIFIED

`scripts/benchmark_reference.py` 输出 runtime / ΩGW 误差 / ΔNeff 误差 / quadrature / ODE / tail 的
对比表（`docs/reference/benchmark_default.jsonl`）。accuracy-vs-runtime Pareto 为部分（fast 快但受
σ 间断残差限制，reference 准但慢）。

## 13. 最终验收标准 — NOT VERIFIED

对照逐项：

| 标准 | 状态 |
|---|---|
| ΩGW(f) 数值误差 < 0.1%（vs 高精度参考） | NOT VERIFIED：fast/production 与参考在低频陡尾差 ~0.15 dex，良分辨区 ~5e-3 dex；未达 <0.1% 全域 |
| 集成 ΔNeff 误差 < 0.01% | NOT VERIFIED：`sigma_exact` 相对参考 −0.40%，生产档更差（−1.0%） |
| 关键 transition 区 < 0.1% | NOT VERIFIED：σ 间断致二阶 Magnus ~0.4–1% |
| 收敛测试显示结果稳定 | NOT VERIFIED：`sigma_exact` 随 h 非单调 |
| analytic limits 全部通过 | VERIFIED（已有测试） |
| energy/scaling consistency 全部通过 | VERIFIED：能量一致性 + RD 平台 ∝f^{n_t} + stiff 带 ∝f^{1} 均通过 |
| production runtime ≥ 100× LSODA | VERIFIED（~1000× warm） |
| Cobaya posterior 与高精度参考一致，min ESS ≥ 2000 | NOT VERIFIED（无收敛长链） |
| 所有数值误差有可解释来源 | VERIFIED（§1/§5 误差预算 + §7.2 归因） |

## 14. 最重要的原则 — VERIFIED

全程物理优先：先回答“为什么这个数值算法在该物理 regime 中正确”，再谈提速；未为 benchmark 数字调参；
`sigma_exact` 的 0.4% 残差与 transition 的过校正都如实记录并写入 `model_bias`，未隐藏。

## 结论

**核心上限**：fast 的二阶 Magnus 在瞬时再加热 σ 间断处带 ~0.4–1% 系统偏差且随 h 非单调。
要达成 §13 的 <0.1%，需在 fast 内核替换为高阶/自适应 ODE（面向 WKB/自适应），或以连续 σ 的
`reference.py` 作为高精度锚点（慢，适合基准/病态点/收敛认证，不适合 MCMC）。当前工程化
（四档模式、error budget + auto_escalate、Cobaya telemetry、参考管线、误差预算、基准）已就绪。

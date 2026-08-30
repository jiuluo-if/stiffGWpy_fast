# 最终报告：physics-first 高速 SGWB 内核（§15）

## 1. 修改后的架构

```
background: gen_expansion / gen_fast  ── 今日锚点 = 连续 N_inf（修复网格量化）
                         │
           piecewise-exact σ (exact_background.sigma_vec/H2_vec, 向量化)
                         │
fast 内核 solve_kernel（Numba / Magnus）:
   • transition-refine（默认 production）:
       build_kink_refined_grid（N_re 节点 + ±2 e-fold h/8 邻域细化）
       + exact_phi_s2_grid（连续 σ 的 F/Φ/S2 + 变步长 h_arr）
   • plain grid（fast/ultra-fast 档）: 原固定步长，最快
                         │
频率采样: construct_f（production 默认，峰值特征感知）或
          freq_adaptive.grid_independent_freqs（网格无关，低频尾更鲁棒）
                         │
observables: Ω_GW(f), ΔN_eff, κ_r + accuracy modes / estimate_error /
             auto_escalate→continuous-σ reference / Cobaya telemetry

reference.py（独立锚点）: 连续 σ(N) + DOP853 + PCHIP/adaptive GK + 网格无关频率集
```

## 2. 每个主要改动的物理依据

1. **今日锚点 = 连续 N_inf**：网格 `floor(N_inf/h)*h ≈ 68.81` vs 连续 `68.81384`，差 `0.0038` e-fold，
   解析尾部幅值 `exp(-0.0038) ≈ 0.996` → 恰为 **~0.4% 系统偏差**（§7.3）。这是真正的主导误差源，
   不是 σ-kink 或 Magnus 阶数。
2. **transition-aware σ-kink**：不再让样条跨越 `σ=1→σ_phys` 间断；在 `N_re` 精确节点两侧分段、
   邻域细化到 h/8。之前失败是旧锚点污染背景；改锚点后与 4 阶 Magnus 同阶。
3. **Cobaya 双引擎 + 无声回退**：production 默认 fast（transition-refine），高误差/异常可自动
   escalation 到连续 σ reference，记录 fast/reference/fallback/error 遥测，MPI-safe。

## 3. 速度提升（default 点，warm）

| 引擎 | runtime | 相对 LSODA | 相对 reference |
|---|---:|---:|---:|
| fast plain grid（fast 档） | ~5 ms | ~3000× | ~30000× |
| fast transition-refine（production） | ~80 ms | ~200× | ~2000× |
| LSODA | ~18.6 s | 1× | ~9× |
| reference（连续 σ DOP853） | ~171 s | ~0.1× | 1× |

## 4. 相对于 reference 的真实误差（不是相对 LSODA）

积分可观测 `ΔN_eff`（MCMC 实际依赖：`Delta_Neff_GW`）：

| 配置 | 相对连续 σ 参考 |
|---|---:|
| plain grid（h=0.01） | −0.16% |
| transition-refine（production，匹配 z_tail） | **−0.022%（<1e-3）** |
| default 与 stiff（同 z_tail 口径） | 均 ≈ −0.04% |

逐频 `Ω_GW(f)` 信号区误差 ~0.2–0.4%（f=−8 处 +0.39%）、低频尾 ~1%（§7.8）。即**积分量已达
<1e-3，逐频仍在 transition-refine 逐模 Magnus 极限**。

## 5. 参数空间验证（§9 部分）

- default（kappa10=1e-2）与 stiff（kappa10=1.0）相对连续 σ 参考均 ≈ −0.04%（同 z_tail）。
- 多样点扫描（default/stiff/lowT/radDominant/cr0_nt/tinyR/edge/norad/highT）均产出有限物理
  ΔN_eff、无 NaN；`highT` 由共享 `DN_eff>5` 护栏拒绝。
- **尚未做**：Sobol/LHS/edge/extreme 的完整统计集（需大量慢速 reference 对标）。

## 6. Cobaya MCMC 验证（§11）

- 双引擎在同一参数点 `ok=True`、无 silent fallback，`Delta_Neff_GW` 差 −0.35%（=z_tail 配置差，
  非求解器失配；匹配 z_tail 时 −0.02%）。
- **点位 ΔlogL（匹配 z_tail=5）**：default 点 fast=0.0022698 vs reference=0.0022703，
  ΔΔN_eff=−5.44e-7（0.024% 相对）；对高斯似然 `|ΔlogL|=1.5e-5`（σ_Neff=1e-4）、
  `1.7e-6`（σ_Neff=3e-4）、`1.5e-7`（σ_Neff=1e-3）——似然层面两引擎几乎不可区分，后验应一致。
- **尚未执行**：同参数同种子的 fast vs reference **收敛后验链**对比（ΔlogL/posterior mean/std/
  16–84%/MAP/ESS/covariance/KS/Wasserstein/KL/posterior shift）。reference 引擎 ~171 s/点，完整
  MCMC 链需数小时；这是唯一尚未执行的长任务。

## 7. 仍然存在的误差来源

- **逐频 Ω_GW(f) 信号区 ~0.2–0.4%**：transition-refine 逐模 Magnus 极限；4 阶 Magnus / 自适应-h /
  精确 z_tail 切换均已试无效（§7.2）。要逐频 <1e-3 需更高阶/自适应逐模积分器，或接受并报告。
- **低频尾 ~1%（f≤−10）**：既有逐模误差，也有频率网格欠采样（需 `freq_adaptive` 曲率自适应加密）。
- **尾部 z_tail 共享近似**：z_tail=5 vs 7 差 ~0.35%；建议 production 与 reference 用相同 z_tail。
- **自洽 tolerance**：外层 bisection 1e-7（远小于物理误差）。

## 8. 最终推荐的 production 配置

```yaml
engine: fast
accuracy_mode: production      # h=0.01, transition_refine=True（kink-aware）
fallback: true
auto_escalate: true            # 高误差/异常升到 engine: reference
error_tol: 0.005
fast_threads: 8
z_tail: 7.0
```

`reference.py` 作为独立高精度锚点用于基准/病态点/收敛认证；不要在 MCMC 热路径用 reference。

## 结论

已从“~0.4–1% 且只对标 LSODA”推进到“**积分 ΔN_eff <1e-3（−0.022%），且相对独立连续 σ 参考衡量，
在 Pareto 上同时优于 LSODA（更快+更准）**”。核心可观测量（MCMC 依赖的 `Delta_Neff_GW`）已达 §8
门槛。剩余大项（逐频 <1e-3、完整参数集、收敛后验链）为持续性长任务，需更高阶逐模积分器或长时计算。

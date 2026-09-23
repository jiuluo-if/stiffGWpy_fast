# Task Plan: stiffgwpy_fast 单一 fast 生产求解器

## Goal

在当前 `codex/fast_v0.2` 分支，将 `fast` 演进为唯一用户可见的科学生产路径，并以独立 reference、参数扫描和性能证据证明精度、速度与稳定性达到目标。

## Current Session: 2026-09-23 — fast_v0.2 tensor propagation optimization

### Goal

在 `fast_v0.2` 目标分支上，以用户指定的 `465196c82b1c938af7bfce125933a4dfda22f4a2` 为基线，先完成同资源配置的 fresh attribution，再逐项验证 P0 repeated-exp elimination、P1 modulo/division assembly elimination 与 P2 LLVM/ASM 热点证据。任何 production candidate 必须先通过 kernel 输出 bitwise gate，再通过 full-outer、edge/Sobol、确定性和正式线程环境门禁；失败候选只保留独立 artifact，不进入生产路径。

### Acceptance criteria

- 不修改 fast profile 的数学公式、精度参数、输出字段、failure/guard 语义或确定性。
- baseline/candidate 只改变一个可验证因素；`f`、`log10OmegaGW`、`DN_gw`、`g2`、`w2` 以及 `handoff_eps` 通过 bitwise gate 后才能继续。
- full-outer 至少覆盖 default、lowT、highT、stiff、high_kappa，并记录 total、tensor kernel、expansion、Phi-S2、frequency preparation、integration 的 median/p95、outer iterations、kernel calls、每通道传播步数。
- 正式资源环境明确记录 Numba threads、threading layer、BLAS budget、affinity、warmup/repeats；不跨资源配置比较。
- 只有稳定收益且所有门禁通过的最小 production diff 才能提交；无收益或不等价候选不得强行合入。
- 有效修改完成后使用 `2966684515@qq.com` 提交，并推送到远端 `fast/fast_v0.2`。

### Phases

- [x] P0: fresh baseline/profile and production call-shape audit
- [x] P1: standalone repeated-exp elimination spike with bitwise gate (REJECTED_FOR_PRODUCTION; no production diff)
- [x] P2: standalone counted assembly-state spike with digest/guard gate (REJECTED_FOR_PRODUCTION; no production diff)
- [x] P3: LLVM/ASM hotspot audit and candidate decision (recompile-backed; specialized spike rejected for formal high-kappa regression)
- [x] P4: final artifact verification and review complete; commit/push is the remaining delivery action (production source unchanged)

## Autonomous continuation: Round 29+ — grouped SoA/AoSoA propagation feasibility

### Acceptance criteria

- Reuse the fresh canonical profile contract: `kink_split=True`, goal grid, PCHIP, fixed affinity/resources, alternating A/B measurements, median and p95.
- Keep the production Cartesian transfer map, phase subdivision, tail handoff, assembly schedule, guards, failure reasons, output ordering, and deterministic semantics unchanged inside the standalone candidate.
- Candidate must account for bucket-induced extra work. A bucket width is admissible only if its lockstep work overhead is measured and the full outer result is not slower after the expected SIMD benefit.
- Kernel outputs must be bitwise equal on all five named cases before outer timing. Full outer must preserve spectra/DN/g2/w2/failure/convergence and deterministic replay; no production dispatch without formal 2/16/20-thread stability and edge/Sobol evidence.

### Phases

- [x] P5a: fresh Round 29 canonical profile, Amdahl analysis, literature scan, and j0/tail divergence screen
- [x] P5b: TDD red/green standalone grouped SoA/AoSoA kernel with explicit scatter and work-overhead telemetry
- [x] P5c: named-case bitwise gate and 2-thread full-outer gate; formal 16/20-thread expansion skipped after decisive regression
- [x] P5d: REJECT grouped SoA/AoSoA; LLVM showed no vector double lanes; select a materially different next experiment

### Next Step

Select the next standalone experiment from the surviving candidate pool; do not retune bucket width or reopen the rejected grouped-lockstep implementation.

## Autonomous continuation: Round 30 — coarse fixed-point predictor feasibility

### Hypothesis

The first outer map at the original `DN_eff` is a stable contraction across named
regimes. A sparse 32-frequency Cartesian pre-solve may estimate the fixed-point
`DN_eff` cheaply enough that one full solve at the extrapolated point replaces the
current two full solves on hard cases.

### Acceptance criteria

- Standalone only; no production API or solver source changes.
- Use the exact production Cartesian transfer, `kink_split=True`, goal-grid ordering,
  PCHIP integration and existing physical guards. The only changed factor is the
  coarse frequency pre-solve plus fixed-point extrapolation.
- First screen: five named cases, candidate full outputs against the converged
  production baseline with spectrum max difference `<=1e-3 dex`, `DN_gw` relative
  difference `<=2e-4`, equal failure/convergence status, and deterministic replay.
- If the first screen passes, expand to edge/Sobol plus independent Cartesian/Prüfer/WKB
  checks and 2/16/20-thread median/p95 end-to-end accounting. If it fails, reject
  immediately and retain the artifact; do not tune the factor on the same sample.
- Report coarse-prepass and full-solve costs separately, then compare total candidate
  end-to-end time against the two-solve production baseline.

### Next Step

Write the failing predictor contract test, then implement a standalone one-pass evaluator and coarse prepass; production remains untouched.

### Round 30 outcome

- [x] TDD caught and fixed the prototype's missing `ln(10)` conversion in the
  PCHIP `d ln(f)` measure.
- [x] Re-ran the corrected gain-1 first-iterate predictor across all five
  named cases. Status/determinism gates passed, but spectrum max error was
  `3.342e-3`–`5.905e-3 dex`, above the `1e-3 dex` gate in every case.
- [x] Audited the residual mismatch to production's one-iteration output
  contract: final scalar `DN_eff` can be updated without regenerating the
  returned spectrum/background. A predictor full solve at the new scalar
  cannot be contract-equivalent while saving the first propagation.
- [x] Reject the predictor family for this output contract; do not tune gain or
  reopen a secant variant on the same path.

### Next experiment: residual-controlled local phase/adiabatic transfer

- [ ] Build a standalone residual-controlled local phase/adiabatic transfer
  twin around the existing Cartesian mode equation. It must estimate local
  phase/commutator defect, use larger blocks only when the defect gate passes,
  and fall back exactly to the production transfer otherwise.
- [ ] First gate: named/edge/Sobol mode-level comparison to the independent
  reference/oracle, physical guards and failure semantics; no production edit.
- [ ] Only if it reduces actual propagation steps and passes the correctness
  gate, run alternating full end-to-end median/p95 A/B at fixed resources.
- [ ] Record ACCEPT/REJECT and select the following concrete candidate in this
  ledger immediately.

### Round 31 outcome

- [x] TDD fallback contract passed after reproducing the production midpoint,
  kink, tail and resource-lifecycle semantics; threshold `0` was bitwise equal.
- [x] Independent reference screen passed on five named and eight edge/Sobol
  points; `edge_tre_hi` was recorded as a pre-existing `shared_Neff_guard`.
- [x] Full outer alternating A/B completed at 2/16/20 threads with warmup and
  median/p95. Low-thread kernel/outer gains did not survive formal thread
  resources; high-kappa regressed at 20 threads.
- [x] Reject residual-controlled two-step midpoint transfer for production;
  retain standalone artifacts and do not retune its threshold.

### Next experiment: fourth-order commutator-free Magnus block

- [ ] Implement a standalone two-exponential fourth-order commutator-free
  Magnus composition over four native intervals, using Gaussian-node operator
  combinations for the existing 2x2 Cartesian transfer system.
- [ ] Use a local two-vs-four native-step defect/commutator estimate as the
  acceptance gate; if the gate fails, fall back exactly to production transfer.
- [ ] First gate: threshold-zero bitwise fallback, named/edge/Sobol independent
  reference/oracle checks, guard/failure/determinism semantics, and measured
  reduction in transfer evaluations. No production edit.
- [ ] Only if the first gate passes, run alternating full end-to-end median/p95
  at 2/16/20 threads. Reject on any formal-regime regression or lack of stable
  >5% total improvement.

### Baseline correction note

The first fresh profiler invocation omitted `--kink-split`; its output is explicitly non-canonical and excluded from evidence. Re-run the same five cases with `kink_split=True` before selecting a candidate.

### Verification blocker note

The new candidate contracts pass. The existing fast test scope has three baseline failures unrelated to this round's unchanged production source: two stale monkeypatch signatures and one pre-existing shared-N_eff guard case. Do not claim a green full suite; proceed with scoped review/push while reporting these failures.

### Next Step

读取本节并运行 fresh baseline；保留既有未跟踪实验产物，不纳入本轮提交。

## Next Step

当前 fresh HEAD `b8e71bf` 的 2-thread stage profile 显示 tensor kernel 仍是目标区
最大阶段（high-T/stiff/high-kappa warm median `4.131/4.682/4.375 ms`），随后是
`fast_phi_s2_split` `1.350/1.535/1.475 ms` 与 background `0.647/0.717/0.710 ms`。
本轮 analytic branch、local radiation jump、单独一阶 adiabaticity handoff 均已
standalone rejected；下一单一动作是设计并验证带高阶 curvature/phase safety guard
的 hybrid handoff prototype，先过 focused accuracy，再考虑 runtime。

二阶 adiabaticity trigger focused A/B 已完成：精度改善但目标 runtime 不稳定，已
rejected；kink post-transition sigma probe cache 也已 bitwise 安全但目标 runtime
无收益，已 rejected。下一单一动作改为 nonoscillatory phase/Riccati 的结构性
propagation prototype。

fresh observable-aware outer reuse 也已完成：强制复用可降低目标工况总耗时约
`23–27%`，但 spectrum 最大偏差为 `1.05e-3–1.55e-2 dex`，DN-only proxy 不足以
放行；不再继续单一阈值调参。

Oracle C 一阶 tail correction 的 z=4 early-handoff standalone 也已拒绝：kernel
约快 29%，但 default spectrum p95 约 `1.83e-2`、DN relative `7.4e-3`；下一步
只研究更完整的 nonoscillatory carrier/phase prototype，不调同一 correction。

以下为历史阶段记录：

已完成 Phase A 测试去重审计与 Phase B/C 的本地实现：compatibility 五版本只跑轻量 smoke，3.11 承担一次 canonical regression，static/package/Cobaya 分离；默认 Numba=2、BLAS=1、reference/oracle workers=1。benchmark 脚本已统一低压力默认值并记录资源 telemetry，主要参数/参考扫描入口已增加 `--threads` 与外层并行时的内层单线程保护。oracle tail 脚本新增 commit/schema/reference-version 绑定的 cache key、逐 z_tail 原子 checkpoint 与 `--resume`；default 完整 76 点 Stage B 已完成，Stage C 已以 low-T/high-T/stiff 各 8 点最小网格完成并记录到 `docs/oracle_tail_convergence_stageC_min8.json`。default full-grid tail systematic `3.6295e-3` 且非单调；Stage C 三点仍为 oracle-sensitive，正式 tail correction 仍未放行。新增有限 phase-window Oracle B 原型后，四点最大 z=5→7 observable 变化为 `3.927e-3–6.897e-3`，证明 handoff 仍敏感但未形成独立 truth anchor，候选不晋升。

Prüfer full native-grid certification 已于 2026-09-11 完成（完整 76 频率、
4 named + 10 edge + 5 Sobol、z=5/7、频谱 + outer 重放 + default 确定性
重放）：full-grid `DN_gw` 相对差 median `6.44e-10`、max `2.44e-9`；
26 accepted outer 全部迭代一致、4 个显式物理 guard 双实现一致；PASS /
VERIFIED，Prüfer 保持 reference-only oracle，不切换正式 kernel。下一步按
当前优先队列转向 **higher-order / 独立 tail oracle（Oracle C）** 或
**nested native-frequency quadrature（真实嵌套求积）**。

Oracle C 解析高阶 WKB tail 也已完成（2026-09-11）：推导
`transfer² × (1 + sin(2θ_f)/ω_f)`，在四点全 76 频率网格上把
frozen-vs-deep 的 `1.34e-3..3.65e-3` tail systematic 压到
`5.46e-6..1.25e-5`（216-665 倍），并解释了 Stage B 的非单调来源。
PASS / VERIFIED。下一轮建议：**nested native-frequency quadrature**
（第九原则，验证 E_nested 覆盖，服务 release gate DN<2e-4）或把 Oracle C
修正晋升为 reference tail 的独立第二锚点（需更严格验证）。

Fast 真实误差标定也已完成（2026-09-11）：`scripts/benchmark_fast_true_error.py`
以 Oracle C 的 WKB 锚点测得 fast 真实 `DN_gw` 误差为 default `4.31e-4`、
highT `7.50e-4`、stiff `1.29e-3`、lowT `1.24e-2`（同约定 PCHIP `2.94e-4`
低估了真实误差）；default/highT/stiff 的剩余误差由残余 fast-vs-WKB 主导，
lowT 是非 tail 独立源。ACCEPTED as calibration finding。下一实验建议：
**按频率分解 fast-vs-WKB 残差**（区分 tail / deep-subhorizon stepping /
frequency quadrature 三源），acceptance criteria：`spectrum max < 1e-3`、
`DN rel < 5e-4`、no new failure、determinism pass；若残差 tail-dominated，
则在 fast tail assembly 内升格 Oracle C 的 `1 + sin(2θ_f)/ω_f` 修正。

频率残差分解已完成（2026-09-11）：
`scripts/benchmark_fast_residual_decomposition.py` 用 `build_Wmat` 权重复现
fast 的 Simpson DN（逐位一致），并在同一积分器下比较 fast 与 Oracle C WKB
锚点的 per-node 值：default `9.49e-06`、highT `9.32e-06`、stiff `1.10e-05`、
lowT `1.41e-04`；而 Simpson-vs-PCHIP 积分器差为
`4.20e-04/7.40e-04/1.28e-03/1.26e-02`。结论：真实 DN 误差由默认 Simpson 频率
积分器主导，fast 传播 kernel 已 ~1e-5。下一实验（先写 acceptance criteria）：
默认 `frequency_quadrature` 换 PCHIP（必要时 Numba 化），要求 DN rel `<2e-4`、
spectrum max 不退化、warm median 增幅 `<10%`、no new failure、
determinism pass。

频率积分器 A/B 已完成（2026-09-11）：`scripts/benchmark_fast_quadrature_ab.py`
显示 PCHIP 把四个命名点的真实 DN 误差降到 `1.07e-05..1.66e-04`（全部进入
`2e-4` gate），但 scipy PCHIP 路径的 warm runtime 比值 `1.204..1.307` 超出
`<10%` 预算；成本来自 `PchipInterpolator.integrate`（`0.141 ms/次`）与局部
estimator（`1.39 ms/次`）。原假设“PCHIP 积分是节点值的线性泛函、可预计算固定
权重向量”已被否证：Fritsch–Carlson 斜率是节点值的非线性函数，用单位矩阵探针
得到的“权重”与逐列积分一致，但 `weights @ y` 与 `PchipInterpolator(x, y)
.integrate()` 相差约一个数量级。因此改为**单次拟合共享**：
`pchip_integral_breakdown` 构造一次 `PchipInterpolator` 并取其原函数在节点上的
差值，同时返回全程积分与逐区间积分，供 DN 全局积分、`g2c[-1]` 与局部
estimator 复用同一次拟合。

单次拟合共享已完成（2026-09-11）（实验前已写定 acceptance criteria：
四点 DN 与现 scipy PCHIP `rel < 1e-12`、DN vs WKB `< 2e-4`、PCHIP/Simpson
warm median 比值 `< 1.10`、spectrum 不退化、no new failure、determinism
pass，否则拒绝切换默认）。实现：`pchip_integral_breakdown` 一次构造
`PchipInterpolator` + `antiderivative()`，同时返回全程与逐区间积分；
`estimate_frequency_quadrature_local(..., candidate_intervals=...)` 复用它；
非重叠三点 Simpson 基线改为等价向量化（实测逐位一致）。同会话配对 A/B
（`docs/fast_quadrature_reuse_ab.json` 与 `..._before.json`，各 50 repeats、
NUMBA=2、BLAS=1）：PCHIP 专属开销 `2.0-2.3 ms -> 0.6-1.0 ms`，比值
`1.23-1.37 -> 1.07-1.17`；四点 DN 与现 scipy PCHIP 逐位一致（rel `0.0`）、
vs WKB 不变。ACCEPTED：opt-in PCHIP 路径加速（同等观测值、约 60% 更少
PCHIP 专属开销）。**REJECTED：本阶段不切换默认 `frequency_quadrature`**，
因为 `<1.10` 未稳健满足（default 同会话 `1.165`，另两次会话 `1.093/1.123`；
16 线程生产刻度探针 `1.101-1.154`）。详见
`docs/fast_quadrature_reuse_assessment.md`。

向量化 PCHIP 积分核也已完成（2026-09-11）（实验前已写定 acceptance：
与 scipy PCHIP `<1e-12` rel、2 线程与 16 线程比值均 `<1.10`、DN 不变、
no new failure、determinism pass）。实现：`_pchip_integrals_vectorized` 复刻
scipy 的 Fritsch-Carlson 斜率与形状保持 guard，并用 Hermite 分段积分闭式
`h/2*(y_i+y_{i+1}) + h^2*(m_i-m_{i+1})/12`；`pchip_integral_breakdown` 保持
scipy 参考实现；estimator 的逐面板分摊改为等价向量化（三种 allocation 逐位
一致）。配对 A/B（`--repeats 50`，2 线程与 16 线程各一次前后配对）：
2 线程比值 `1.1622/1.1261/1.0495/1.0901 -> 1.0105/1.0551/1.0397/1.0411`，
16 线程 `1.1987/1.1137/1.1270/1.1203 -> 1.0428/1.0342/1.0180/1.0413`；
同会话第二次配对复现同一模式（before `1.050-1.184`、after `1.002-1.056`）。
PCHIP 专属开销降到 `0.06-0.42 ms`（2 线程）/`0.15-0.33 ms`（16 线程）；
实测 `DN_gw` 相对 scipy PCHIP 路径变化 `3.83e-16/4.98e-16/0/0`（1-2 ulp）。
ACCEPTED：预算达标且默认 Simpson 路径逐位不变。

下一实验（acceptance criteria 已写定，见
`docs/fast_quadrature_reuse_assessment.md` 的 “Next experiment (pre-registered)”）：
把 `SGWB_iter_fast` 的默认 `frequency_quadrature` 切到 `pchip`，要求四点
`DN_gw` vs Oracle C WKB `< 2e-4`、Sobol/edge 屏幕在预算内、spectrum max 不退化、
no new failure、determinism pass、2/16 线程比值 `< 1.10`，并同步刷新
`validation_manifest.json`、README、`ERROR_BUDGET` 与 estimator coverage
artifact，重新满足 “validation artifact == release HEAD”。

**本轮执行（2026-09-11，默认切换已落地）**：`_SGWB_iter_fast_impl` 与
`SGWB_iter_fast` 的默认 `frequency_quadrature` 改为 `pchip`（`simpson` 仍可
显式选择）；PCHIP 热路径遇到非有限 integrand 时改为把 NaN 传给统一的
`math.isfinite(DN_gw_new)` guard，保持既有 `fast_failure_reason='nonfinite'`
的 abort/restore 语义，`_pchip_integrals_vectorized` 继续 fail-loud。
四点 `DN_gw` vs Oracle C WKB 为 `1.07e-5/1.66e-4/9.88e-6/1.15e-5`（全
`< 2e-4`）；同会话 50-repeat 配对的 PCHIP/Simpson warm 比值 2 线程
`1.0105/1.0551/1.0397/1.0411`、16 线程 `1.0428/1.0342/1.0180/1.0413`（全
`< 1.10`）；canonical full pytest `159 passed, 6 deselected`。评测与刷新
清单见 `docs/fast_quadrature_default_switch_assessment.md`。

追加（同一轮，证据与文档刷新后定稿）：默认切换的权威 50-repeat 配对为
`docs/fast_quadrature_default_ab_t{2,16}.json`，2 线程比值
`1.019/1.035/1.034/1.024`、16 线程 `1.069/1.028/0.987/1.017`
（default/lowT/highT/stiff）；20 线程 25-repeat 矩阵 default warm median/p95
`4.932/5.467 ms`、highT/stiff/high-kappa `7.79/7.79/7.27 ms`；stability
`guard_count=3`、`failure_count=0`；LHS 400 点 `254 success/146 guard/0
numerical failure`。`ERROR_BUDGET` 频率积分项 `1.0e-3 -> 2.0e-4`，
README/README_zh/CHANGELOG/`docs/fast_v02_audit_report.md` 与 manifest 同步
刷新；结论 ACCEPT。下一方向（按预期科学价值排序）：high-T/stiff/high-kappa
runtime 分解（`7.8/7.8/7.3 -> <=4 ms`）、nested native-frequency 真实嵌套
求积（区分 same-grid 残差来源）、analytic stiff/RD branch。

## Pre-registered experiment: fast Python 层准备冗余消除（2026-09-12，改动前登记）

**Hypothesis**：20 线程下 fast 的 runtime 由 Python/背景准备层主导
（`tensor_solve_kernel` 仅约 15%），其中存在纯冗余计算与冗余内存写入；消除它们可在
**不改变任何数值输出**的前提下降低 warm runtime。

**改动前证据**（六点、25 repeats、`docs/fast_pyoverhead_before_t20.json` 与
`docs/fast_pyoverhead_before_t2.json`）：20 线程 warm median default `5.290 ms`、
highT `6.228`、stiff `6.965`、high_kappa `7.067`、lowT `4.298`、low_r `4.076`；
2 线程 default `6.885 ms`、highT `9.987`、high_kappa `9.764`、lowT `6.360`、
low_r `5.820`、stiff `11.308`。default 20 线程阶段分解：
`expansion_background` 1.035、`tensor_solve_kernel` 0.813、`fast_phi_s2_split` 0.627、
`goal_frequency_construction` 0.482、`column_integration` 0.170、`pchip_fine` 0.170、
`frequency_preparation` 0.109、`correct_kink_background` 0.214 ms。

**候选改动**（全部为等价冗余消除）：(a) `Ogw/Oj/Opgw` 在 `np.zeros` 之后立即
`fill(0.0)` 的重复写入；(b) `grid_independent_freqs.f_hor_cont` 每次重复求值的
`N`-无关 `H2_vec`；(c) `gen_fast` 的 `np.unique(np.sort(...))` 整网格重排；
(d) `_sigma_node_limits`/`fast_phi_s2_split` 内重复的 `m.derived_param` 属性求值。

**Acceptance（全部满足才接受）**：

1. 六点 `f`/`log10OmegaGW`/`DN_gw`/`g2`/`w2` 的 SHA256 digest 与改动前**逐位一致**；
2. 2 线程与 20 线程 warm median，default 改善 ≥ 5%，六点均不退化超过 2%；
3. `converged`/`fast_failure_reason`/guard 计数不变，无新失败；
4. 同参数重复调用 digest 稳定（determinism pass）；
5. full pytest、Cobaya、Ruff、mypy、compileall、manifest、
   `git diff --check`、wheel + installed-wheel smoke 全部通过。

未满足则 REJECT，并把负面结论写入 `docs/`。

**结果（2026-09-12，改动后）**：ACCEPT。

1. PASS：六点 `f`/`log10OmegaGW`/`DN_gw`/`g2`/`w2` SHA256 digest 逐位一致
   （20 线程 3 cycle x 4 次运行共 360 项 + A/B 交叉、2 线程 3 cycle x 4 次
   运行共 288 项 + A/B 交叉，0 mismatch）；
2. PASS：A,B,B,A 对称序配对 default ratio `0.9456`（改善 `5.44%`，20 线程）、
   `0.9177`（`8.23%`，2 线程）；六点中位 ratio `0.816-1.011`，无点退化超 `2%`；
3. PASS：`converged`/`n_freq`/`fast_failure_reason` 全部不变，无新 guard；
4. PASS：同一字段在所有 12 次（20 线程）与 6 次（2 线程）运行中 digest 相同；
5. PASS：full pytest `159 passed, 6 deselected`、`pytest -m cobaya` `1 passed`、
   Ruff、mypy、compileall、manifest、`git diff --check`、
   wheel + installed-wheel smoke 全部通过。

证据：`docs/fast_pyoverhead_ab_symmetric.json`、
`docs/fast_pyoverhead_before_t20.json`、`docs/fast_pyoverhead_after_t20.json`、
`docs/fast_pyoverhead_before_t2.json`，以及 `findings.md` 的
「Fast Python preparation-layer de-duplication (2026-09-12)」。

## Current Phase

Phase 3: Implementation and evidence-driven optimization

## Phases

### Phase 1: Requirements & Discovery
- [x] 读取目标文件并提取可验证要求
- [x] 确认当前分支、工作树、远端和推送约束
- [x] 完成 fast/production 真实代码路径审计
- [x] 建立 Phase A exact kink split 的可复现实验基线
- **Status:** complete

### Phase 2: Design & Planning
- [x] 比较 uniform-grid split、variable-grid refine 与 exact primitive 方案
- [x] 形成 reference-only oracle、分阶段 benchmark、兼容 guard 边界
- [x] 编写并自审实现计划
- **Status:** complete

### Phase 3: Implementation
- [x] 以 TDD 实现 Phase A exact kink split 和 one-sided primitive
- [x] 实现 Phase B split-segment phase cap
- [x] 实现 Phase C goal frequency grid 和 native eval nodes
- [x] 合并为唯一正式 fast preset；高层旧 production/transition-refine 别名弃用并映射到 fast，底层 validation 入口保留（`c9110c9`）
- [ ] 继续优化 DN/速度；本轮必须先重建远端 HEAD fresh profiling 与 DN_gw 误差分解，不得在 profiler 证据前 micro-optimize；已接受物理尾部匹配修正（gamma=1）、背景节点缓存、formal kink 的 frequency-only preparation、受门限保护的 exact primitive 复用、smooth-node `fast_phi_s2_split` 和受背景稳定性门控的 outer full-solve 复用，四阶 Magnus/曲率子步、low-T 局部加密、无门控删除第二次 full solve、批量 sigma 采样、phase_max 加密及 z_tail 加深已拒绝
- [x] 继续优化 DN/速度；已完成固定环境 fresh profiling、六点 25-repeat runtime/积分差矩阵、同网格 oracle 与 PCHIP/插值法探针；默认 telemetry 的 Simpson-trapezoid estimator 已以 50-repeat profiler 证明约 8.7% 局部收益，DN-driven midpoint 排序原型因默认点收敛不单调暂不接受；有限 phase-window Oracle B 原型已完成但不晋升；Prüfer standalone 已完成固定 8 频率完整 `Ogw/Oj/Opgw/DN_gw` 与 outer self-consistency 复核，并完成 10 个参数轴 edge + 5 个固定 Sobol 点的 guard/outer 复核（26 个 accepted comparison、4 个 physical guard、accepted 最大 outer DN 差 `6.44e-9`）
- [x] Prüfer full native-grid certification：`scripts/benchmark_prufer_fullgrid.py` 在完整 76 频率网格完成 19 点（4 named + 10 edge + 5 Sobol）x `z_tail=5/7` 的频谱与 outer 自洽重放；full-grid `DN_gw` 相对差 median `6.44e-10`、max `2.44e-9`；26 accepted outer 迭代一致、4 个显式 physical guard 双实现一致；default 重放 bitwise 一致；PASS / VERIFIED，Prüfer 保持 reference-only oracle，不切换正式 kernel
- [x] Python 准备层冗余消除（2026-09-12）：`gen_fast` 目标网格去重降到单次归并、`grid_independent_freqs` 的 `N`-无关量提到闭包外、新分配缓冲跳过冗余 zero-fill、重复 `derived_param` 求值单次绑定；六点 × 五字段 digest 逐位一致（648 项比较 + A/B 交叉 0 mismatch），A,B,B,A 对称序配对 default 改善 20 线程 `5.4%` / 2 线程 `8.2%`，六点无超 `2%` 退化；新增 `scripts/benchmark_fast_runtime_ab.py` 计时工具，ACCEPT
- [x] Nested native-frequency quadrature（2026-09-12）：新增 `scripts/benchmark_nested_frequency_quadrature.py`，把局部估计器 top-N 区间的真实中点并入**积分** support grid；PCHIP 六点 `E_nested` 绝对 `1.7e-14`-`6.3e-10`、相对 `7.4e-12`-`2.8e-9`，比同点 `actual_rel`（`1.83e-4`-`2.97e-4`）小 5-7 个量级；simpson 阳性对照 `1.362e-9` 证明测量有效；Oracle A 同网格 `z_tail` 8->10 使 default fast-vs-oracle DN 差 `2.944e-4 -> 4.956e-6`（与 oracle 自报 frozen-handoff 缺陷同阶），fast 与两个独立 oracle 在 z=10 上一致到 `5e-6`；证据性 ACCEPT，未改 kernel
- **Status:** in_progress

### Phase 4: Testing & Verification
- [x] 运行单元、API、完整频率 reference 和 native likelihood-node 验证
- [x] 运行 cold/warm、p95、线程和 profiler 对照
- [ ] 完成参数空间、likelihood bins、guards、NaN/fallback、Cobaya derived 全矩阵，并形成 accuracy-runtime Pareto 数据
- [ ] 完成 76/80/90/110 节点的 DN convergence、spectrum convergence、runtime 与 estimator coverage；验证 eval_freqs 不改变 DN 的科研生产 invariant
- **Status:** pending

### Phase 5: Delivery
- [x] 更新 README、README_zh、Cobaya YAML 和 tests
- [ ] 更新与 release HEAD 完全一致的 manifest、benchmark docs、README 表格、README_zh、CHANGELOG 并生成完整优化报告
- [x] 记录 goal refinement、局部传播和 Richardson 探针的独立 A/B 与回退结论
- [x] Phase A commit `07895e6` 已 fresh verification 后推送 `fast/fast_v0.2`
- [x] Phase B/C 和正式 preset 变更 commit `78cc5b9` 已 fresh verification 后推送；CI 修复 commit `4436b9d` 已通过远端全矩阵
- **Status:** pending

## Decisions Made

| Decision | Rationale |
|----------|-----------|
| 先做 Phase A exact kink split，再决定后续算法阶段 | 目标明确要求逐项隔离实验，不能一次性合并多个未经证明的优化 |
| `stiffgwpy_fast.reference` 是唯一 precision oracle | 目标明确禁止用 LSODA 作为精度真值 |
| 保留 float64、physical guards、fallback 语义和输出契约 | API 合并不得改变物理参数和下游 Cobaya 行为 |
| 重大 Phase 单独提交并推送到 `fast` 远端同名分支 | 用户明确要求每次重大修改后同步 GitHub |

## Errors Encountered

| Error | Resolution |
|-------|------------|
| 首次 `create_goal` 失败：当前线程已有未完成 goal | 继续使用用户更新后的 active goal，不重复创建 |
| 初次计划补丁匹配失败 | 读取 init-session 实际模板后按当前内容重新生成计划文件 |
| 本轮对修改文件执行全量 ruff 暴露 84 个既有 E701/E702/F841 | 未扩大范围格式化；compileall 与 diff check 通过，保留该既有门禁问题并单独报告 |
| Q2 dense-reference 首次运行暴露局部 estimator 的 Chebyshev spline 未初始化 | 补充局部 Chebyshev 回归测试并初始化 PCHIP spline 后重跑通过 |
| 全局 Simpson 权重分摊造成局部 estimator 假保守 | 改为非重叠 local composite-Simpson panel；default 仍约 15.3% 过保守，保留为诊断候选，不进入 production |

## Current Session Notes: 2026-09-17

- 已从正确的 `fast` 远端 fetch `fast_v0.2`，本地 HEAD 与远端同为 `f411bbb252e20b7626b33d007ea73d40bf34f46c`，工作树初始干净。
- 已重新读取 `progress.md`、`findings.md`、`docs/validation/validation_manifest.json` 与当前 round benchmark/rejection artifacts；manifest 仍绑定旧 commit `e1d2dde28404827f2c7c43107fab15420edf5ffb`，不能直接作为当前 HEAD 证据。
- fresh profile 首次尝试失败：`run_fixed_profile.py` 外层参数只允许 `A/B`，无法接受内层已定义的 `default/highT/stiff/high_kappa`；未据此形成性能结论，改用支持六点的 runtime 矩阵与 `A/B` 阶段 profiler。
- outer reuse headroom 首次包含 `edge_tre_hi` 后被既有 `total N_eff too large` physical guard 中止，脚本未将 guard 转为记录而直接访问缺失的 `DN_gw`；该点不计 numerical failure，后续先排除显式 guard 点完成 accepted-point 诊断。

## Next Step

已完成当前 HEAD fresh runtime/阶段 profile、outer-reuse baseline、kernel 输入局部响应原型，以及 phase-substep exponential recurrence 的 standalone A/B、outer solve、Oracle C focused 对照和 14 点 named/edge/Sobol 矩阵；当前下一步为扩大 recurrence 的 full-grid Oracle A/Prüfer/WKB 与 total-runtime candidate 认证，仍未进入 production。

## Current session update: 2026-09-20

- [x] 完成全 named/edge/Sobol 的解析 S2 端点项 + horizon-start 项边界诊断；`14/24`
  eligible，安全因子 1.1 下 coverage `4/14`、false-safe `0/4`，目标
  high-T/stiff/high-kappa `0/3`，因此不进入 production。
- [ ] 补齐 `j0/z0/tail/phase-path` 响应项，先以 standalone 形式和 Cartesian/Prüfer/WKB
  oracle 对照；未获得预算内 target-regime coverage 前不做 runtime A/B。

- [x] standalone no-Psi preparation A/B：五个 regime 各 50 次逐位一致，但 target
  best runtime gain 仅 `4.44%` 且其他目标点变慢，REJECTED，不改 production。
- [ ] 下一候选继续从 fresh profile 的 allocation/copy 或 tensor-kernel 结构中寻找可证明
  的收益；不得重复已拒绝的 Psi、derived-param、goal-grid 或固定 z 解析分支实验。
- [x] 三对角恒定 `aa/cc` buffer standalone production-candidate A/B：digest 全一致，
  high-T/stiff 收益 `3.95%/2.67%`，high-kappa 退化 `3.85%`，REJECTED 并回退；测试与
  benchmark 线程环境改为独立进程启动时固定，避免复现 Numba 环境冲突。

- [x] 完成一阶绝热性触发 carrier prototype；`eps=3e-4` 精度可接受但无稳定 >5%
  runtime 收益，`eps>=1e-3` 误差或速度不合格，故 rejected。
- [x] 定位并修复本轮 CI 根因：Ruff `I001` import order；Ruff、mypy 已复核通过。
- [ ] 下一轮若继续 hybrid，加入二阶 adiabaticity/error bound 与 full-grid Oracle A/
  Prüfer/WKB；仍不修改 production path。

- [x] 完成 outer goal-grid invariant 与 25-repeat named/Sobol/edge A/B；`10/13`
  digest 不一致且 speed ratio `0.960..1.003`，候选 rejected。
- [ ] 研究只缓存不影响最终 `f` digest 的 background/primitive 中间量；先建立逐位
  输出契约，再做目标 regime runtime A/B。

- [x] 完成 `derived_param` 单次 expansion 局部缓存 A/B；四个 regime digest 一致但
  收益最高约 1.9%，按 >5% 门槛拒绝并撤回正式实现。
- [x] 接受 H2 endpoint 函数内缓存：13 点 digest mismatch 为 0，high-kappa 50-repeat
  runtime 改善 6.54%，Cartesian/Prüfer oracle 与 error budget 通过。
- [x] 完成 analytic branch eligibility boundary；目标 regime 可覆盖比例仅
  `14.31%–17.66%`，低于 >30% 优先门槛，closed-form prototype 暂缓。
- [x] 清理无关中文注释/编码门禁及同类 workflow 检查；将 Numba 资源环境测试隔离到
  导入前子进程，避免同一环境污染原因重复触发 CI 失败。
- [x] 完成 HEAD `759d5cc` fresh round-5 分层 profile；exact primitive allocation/fusion
  与 phase exponential reuse 均按 digest/稳定 runtime gate rejected，production 未改。
- [ ] 下一候选聚焦 tensor phase kernel 的结构性工作量削减；先做 standalone exact
  prototype，再进行 25–50 repeat A/B 与独立 Oracle A/Prüfer/WKB 检验。
- [x] 完成 round-6 fresh profile；static-grid reuse、no-assembly probe 及组合均因目标
  区间未稳定超过 5% 而 rejected，未修改 production。
- [ ] 下轮转向 phase transfer 的跨 channel 结构复用或独立 adiabaticity estimator；
  不重复本轮两个 specialization。
- [x] Phase-substep exponential recurrence standalone prototype：2/16-thread 30-repeat
  kernel A/B、五工况 outer solve、14 点 named/edge/Sobol 矩阵与当前 HEAD fresh
  Oracle C 完成；速度通过，digest 非逐位一致，保留为 prototype，不进入 production。
- [ ] 对 recurrence candidate 完成 full-grid Oracle A/Prüfer/WKB、reheating/kink
  edge、扩大 Sobol coverage 和 total runtime gate；未完成前不得接入正式 fast。

### Round continuation: full native-grid audit

- [x] Fetch latest `fast_v0.2` and reread required state/artifacts/rejections.
- [x] Add and run standalone full native-grid recurrence versus Prüfer/DOP853 audit.
- [x] Record five-regime spectrum/DN results; keep production path unchanged.
- [ ] Run edge/kink and expanded Sobol coverage.
- [ ] Run formal 16/20-thread total-runtime candidate and independent Oracle A/WKB gate.

### Recurrence decision

- [x] Full 24-point edge/Sobol native-grid audit with explicit physical-guard handling.
- [x] Full outer matrix: status match and false-safe gate passed (`0` mismatches / `0` false-safe).
- [x] Formal 16/20-thread total-runtime A/B completed; target-wide stable `>5%` gate failed.
- [x] Reject recurrence for production; preserve artifacts and precision/guard evidence.
- [ ] Select the next candidate only after a new fresh stage profile; do not repeat this path.

### Round 14 decision

- [x] Fresh HEAD `cf183c4` 上完成 WKB carrier `z_match=4 -> z_tail=5` standalone
  25-repeat A/B；default/high-T/high-kappa 有局部收益，但 stiff 退化约 5.1%。
- [x] Reject carrier preintegration for production；保留 artifact 和 twin accuracy
  evidence，production source 未修改。
- [ ] 下一轮重新读取 fresh profile 与全部 rejection evidence，再选择新的
  propagation 工作量候选；不得继续微调同一 handoff。

### Round 15 decision

- [x] Fresh stage profile 重新确认 tensor kernel 为主要热点；完成 sigma-node
  no-copy standalone spike 与 formal/exact-kink TDD contract。
- [x] 五工况 25-repeat digest/guard/收敛一致，但最高 runtime 改善仅 1.81%，拒绝
  production；不再重复 primitive allocation 微调。
- [ ] 下一轮重新 fetch/profile，优先寻找尚未验证的严格等价 kernel 外围重复工作，
  或在 standalone 中提出新的数学依据后再做结构性 prototype。

### Round 16 decision

- [x] Fresh cProfile 定位并验证 FD interpolator lookup cache 假设；50-repeat
  determinism 与逐位输出契约通过。
- [x] 25-repeat 局部收益无法在 50 repeats 复现，拒绝 production；不再重复同一
  Python lookup/cache 微优化。
- [ ] 下一轮重新 fetch/profile，转向尚未验证的 outer allocation/assembly 或新的
  standalone 数学路径，并保持独立 precision gate。

### Round 17 decision

- [x] Fresh stage profile 与 channel-overlap diagnostic 完成；`(j0, tail_index)`
  pair 无重复，简单 cross-channel exact transfer reuse 排除。
- [ ] 下一轮重新 fetch/profile，重点拆解 outer allocation/assembly；若尝试结构性
  propagation，必须先处理 mode-local phase subdivision 的 exactness。
## Round 18: outer assembly attribution (2026-09-17)

- [x] Fresh stage attribution at current HEAD；确认 high-T/stiff/high-kappa
  每次运行执行两次 `assemble=1` full solve，default/low-T 执行一次。
- [x] 保留 artifact `docs/outer_attribution_round18_20260917.json` 与脚本
  `scripts/benchmark_outer_attribution.py`；production path 未修改。
- [x] 排除把 outer assembly shortcut 当作新候选：既有独立 profile 显示总收益
  仅约 2--4%，低于速度门槛。
- [ ] 下一阶段只研究能减少真实 propagation arithmetic/steps 的 standalone
  phase-function/Riccati/WKB 方向；先做数学残差和独立 Cartesian/Prüfer/WKB
  对照，再考虑 25--50 repeat runtime gate。

## Round 19: tail-factor cache spike (2026-09-17)

- [x] 先完成 TDD identity contract，再做只读 tail-factor cache twin。
- [x] 五工况 25-repeat kernel A/B 完成；candidate/base 为
  `0.9633/0.9763/0.8911/0.9610/0.8892`。
- [x] 发现五工况 digest 均 mismatch；拒绝 production，不进入 Oracle 或 outer
  gate，保留失败 artifact 与精度边界。
- [ ] 下一候选必须保持原始浮点运算顺序，或提出新的 standalone 数学依据；不再
  直接推广预计算乘积到正式 kernel。

## Round 20: phase loop-state reuse (2026-09-17)

- [x] fresh stage profile 与 high-T/stiff/high-kappa outer iteration 归因完成。
- [x] standalone exact twin 通过 kernel/full-outer bitwise、spectrum/DN、failure
  gate；2-thread 50-repeat default/high-kappa full-outer 改善 `7.4%/8.6%`。
- [x] 16/20-thread formal A/B 无明显退化，确认正式资源下收益受 preparation
  层限制；不把候选宣称为 10% breakthrough。
- [x] 接入 production 的最小 loop-state reuse patch；待新 HEAD post-commit
  benchmark、全门禁和 artifact provenance 复核。
- [x] 新 HEAD `94144c1` post-commit full pytest 通过（`165 passed, 6 deselected`），
  五点 fresh profile artifacts 已绑定该 SHA；production 优化保留。

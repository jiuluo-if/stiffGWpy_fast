# Task Plan: stiffgwpy_fast 单一 fast 生产求解器

## Goal

在当前 `codex/fast_v0.2` 分支，将 `fast` 演进为唯一用户可见的科学生产路径，并以独立 reference、参数扫描和性能证据证明精度、速度与稳定性达到目标。

## Next Step

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
5. full pytest、Cobaya、Ruff、mypy、compileall、中文注释门禁、manifest、
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
   Ruff、mypy、compileall、中文注释门禁、manifest、`git diff --check`、
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

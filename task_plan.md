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

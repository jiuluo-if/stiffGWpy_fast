# Task Plan: stiffgwpy_fast 单一 fast 生产求解器

## Goal

在当前 `codex/fast_v0.2` 分支，将 `fast` 演进为唯一用户可见的科学生产路径，并以独立 reference、参数扫描和性能证据证明精度、速度与稳定性达到目标。

## Next Step

已在 `e122fc7` 完成 Phase 1 CI/package audit；当前为 oracle tail convergence 诊断准备阶段。新增 `z_tail=5/6/7/8/10` 汇总器与可复现实验脚本，但本轮完整 native-grid 扫描因 reference 计算过慢停止，尚无新的数值结论；不得把旧 tail 片段当作本轮证据。

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
- [ ] 继续优化 DN/速度；已完成固定环境 fresh profiling、六点 25-repeat runtime/积分差矩阵、同网格 oracle 与 PCHIP/插值法探针；默认 telemetry 的 Simpson-trapezoid estimator 已以 50-repeat profiler 证明约 8.7% 局部收益，DN-driven midpoint 排序原型因默认点收敛不单调暂不接受；下一步是参数空间 oracle/estimator coverage，不得切换 PCHIP 默认
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

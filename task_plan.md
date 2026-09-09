# Task Plan: stiffgwpy_fast 单一 fast 生产求解器

## Goal

在当前 `codex/fast_v0.2` 分支，将 `fast` 演进为唯一用户可见的科学生产路径，并以独立 reference、参数扫描和性能证据证明精度、速度与稳定性达到目标。

## Next Step

继续以 profiler 和独立 reference 优化正式 fast 的 DN/速度：已接受背景节点缓存，outer probe 审计确认不能直接删掉第二次 full solve；下一步优先研究能减少显式振荡步数或频率通道的算法，再做跨参数验证；任何有证据的重大阶段单独提交并推送。

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
- [ ] 继续优化 DN/速度并决定后续阶段是否合并；已接受物理尾部匹配修正（gamma=1）和背景节点缓存，四阶 Magnus/曲率子步、low-T 局部加密、直接删除第二次 full solve 及批量 sigma 采样已拒绝
- **Status:** in_progress

### Phase 4: Testing & Verification
- [x] 运行单元、API、完整频率 reference 和 native likelihood-node 验证
- [x] 运行 cold/warm、p95、线程和 profiler 对照
- [ ] 完成参数空间、likelihood bins、guards、NaN/fallback、Cobaya derived 全矩阵
- **Status:** pending

### Phase 5: Delivery
- [x] 更新 README、README_zh、Cobaya YAML 和 tests
- [ ] 更新 manifest、benchmark docs、CHANGELOG 并生成完整优化报告
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

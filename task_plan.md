# Task Plan: stiffgwpy_fast 单一 fast 生产求解器

## Goal

在当前 `codex/fast_v0.2` 分支，将 `fast` 演进为唯一用户可见的科学生产路径，并以独立 reference、参数扫描和性能证据证明精度、速度与稳定性达到目标。

## Next Step

基于已完成的代码路径审计和 fresh baseline，确认 Phase A 设计取舍并等待用户批准后进入实现。

## Current Phase

Phase 1: Requirements & Discovery

## Phases

### Phase 1: Requirements & Discovery
- [x] 读取目标文件并提取可验证要求
- [x] 确认当前分支、工作树、远端和推送约束
- [x] 完成 fast/production 真实代码路径审计
- [x] 建立 Phase A exact kink split 的可复现实验基线
- **Status:** complete

### Phase 2: Design & Planning
- [ ] 比较 2–3 个 kink split / hybrid solver 方案
- [ ] 形成经用户确认的设计与验收边界
- [ ] 编写并自审实现计划
- **Status:** in_progress

### Phase 3: Implementation
- [ ] 以 TDD 实现 Phase A exact kink split
- [ ] 单独验证并决定接受或回退 Phase A
- [ ] 依次推进 Phase B–E，仅合并有证据收益的阶段
- **Status:** pending

### Phase 4: Testing & Verification
- [ ] 运行单元、API、reference、参数空间和 likelihood bins 验证
- [ ] 运行 cold/warm、p95、线程 scaling 和 profiler before/after
- [ ] 审计 guards、NaN/fallback、Cobaya derived quantities 与兼容别名
- **Status:** pending

### Phase 5: Delivery
- [ ] 更新 README、Cobaya YAML、schema、tests、manifest、benchmark docs、CHANGELOG
- [ ] 生成完整优化报告
- [ ] fresh verification 后按阶段提交并推送远端同名分支
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

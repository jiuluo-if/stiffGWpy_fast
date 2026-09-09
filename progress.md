# Progress Log: stiffgwpy_fast 单一 fast 生产求解器

## Session: 2026-09-09

### Current Status

- **Phase:** 3 - Implementation and evidence-driven optimization
- **Started:** 2026-09-09

### Actions Taken

- 读取目标文件，确认最终目标为单一 fast 生产 solver，并保留 reference oracle。
- 读取并启用 `using-superpowers`、`self-evolution`、`planning-with-files`、`brainstorming` 与 `writing-plans` 的适用流程。
- 确认分支 `codex/fast_v0.2`、远端目标 `fast/fast_v0.2`；重大修改按阶段提交并推送。
- 完成首轮关键词和入口检索，发现 README、Cobaya、benchmarks、tests 仍广泛暴露双档位。
- 核对远端：`fast` 指向目标仓库；用户要求重大 Phase 提交并推送远端同名分支。
- Phase A 已加入 uniform-grid 单区间 transfer split 与 one-sided reheating primitive；目标测试覆盖断点属性和单侧积分约定。
- Phase A 的 sparse subset 结果已降级为诊断；benchmark 现默认在候选完整 246 点频率网格上运行独立 reference，避免把 subset 频率积分误报成总 `DN_gw` 误差。
- Phase B 已完成 split transfer segment 的 phase-cap 实现；`phase_max=.25` target test 通过，但单点 oracle 收益有限。
- Phase C 已接入 `freq_grid='goal'`：seed + reheating feature reserve + native `eval_freqs`；测试和 solver native-node 验证通过。
- 正式 fast preset 已切换为单一 `goal-kink-hybrid` 组合；旧 production 保留作 validation/compatibility mode。
- CI #53 的根因是新增英文代码注释未满足中文注释门禁；`4436b9d` 已修复并推送，CI run `34304423555` 的 3.9–3.13 与 Cobaya jobs 全部通过。
- seed64 频率预算已完成独立 reference A/B：76 点比 seed80 的 90 点更省节点，默认点 spectrum p95 `2.719e-3`、DN rel `1.368e-3`；保留 seed64 作为当前实验基线。
- 四阶 Magnus 与背景曲率局部子步实验均回退：没有把 `h=.01` 的 spectrum p95 压到 `<3e-3`，且引入额外热路径成本。
- 新增安全热路径优化：外层首轮 `solve_kernel` 只组装最终列，后续收敛轮再组装完整列；新增测试通过，like-for-like plain warm median 约 `6.727 -> 5.051 ms`，正式 kink 路径尚未达到 `<=4 ms/point`。
- 节点 sigma 复用实验虽数值等价，但没有稳定耗时收益，已回退，不纳入主线。
- exact split 的节点 sigma 缓存已通过逐项等价测试并接入 formal kink 路径；profile steady warm median 约 `10.62 -> 8.90 ms`。外层审计显示第二轮 DN 更新很小，但一次 full solve 会让最终 `g2/w2` 在 high-T 点偏离约 `3%`，所以只保留首轮 probe 优化，不跳过第二次 full solve。

### Test Results

| Test | Expected | Actual | Status |
|------|----------|--------|--------|
| `git branch --show-current` | `codex/fast_v0.2` | `codex/fast_v0.2` | PASS |
| `git diff --stat` | 初始无代码 diff | 无输出 | PASS |
| `git remote -v` | 目标 `fast` 远端存在 | `https://github.com/jiuluo-if/stiffGWpy_fast.git` | PASS |
| `git ls-remote --heads fast` | 远端同名目标存在 | `fast_v0.2`=`4436b9d`，`main`=`b83aa89` | PASS |
| `python -m pytest -q` | 现有回归通过 | `110 passed, 6 deselected` | PASS |
| `python scripts/bench_fast.py --reps 3 --cases 0` | fresh fast runtime 基线 | warm median `2.461 ms`，fallback `0` | PASS |
| `python scripts/benchmark_reference.py --point default --z-tail 5 --no-ode-error --no-tail-error` | fresh oracle 对照 | production fast vs reference 已输出，误差仍超目标 | PASS |
| `python -m pytest -q tests/test_freq_adaptive.py::test_split_primitive_uses_left_limit_at_reheating_kink tests/test_freq_adaptive.py::test_breakpoint_phi_s2_is_accurate_without_dense_subgrid tests/test_fast_sgwb.py::test_kink_split_variant_inserts_only_reheating_breakpoint` | Phase A target tests | 2 passed, 1 passed after endpoint correction; current rerun pending | PASS |
| `python scripts/benchmark_phase_a.py --reps 3 --subset-reference` | Phase A sparse diagnostic | plain `7.389 ms`, kink `8.305 ms`; transition max `1.281% -> 0.776%`; total DN comparison invalid by design | DIAGNOSTIC |
| `python scripts/benchmark_phase_a.py --reps 3` | Phase A full-grid oracle | plain `7.713 ms`, kink `10.020 ms`; DN rel `1.069% -> 0.411%`; spectrum dex p95 `1.004e-2 -> 4.745e-3`; transition max `3.066% -> 1.317%` | PARTIAL |
| `python scripts/profile_fast_breakdown.py --case A --reps 3` and `--kink-split` | Phase A profiler | steady median total `3.701 ms` plain / `7.441 ms` kink at 16 threads; kink exact primitive overhead is measurable | PARTIAL |
| `python scripts/benchmark_phase_a.py --reps 3 --freq-grid goal --phase-max 0.25` | Phase B/C at `h=.02` | 76-point grid; kink spectrum dex p95 `4.792e-3`, DN rel `4.344e-3` | PARTIAL |
| `python scripts/benchmark_phase_a.py --reps 2 --h 0.005 --freq-grid goal --phase-max 0.25` | Phase B/C precision probe | seed64 76-point grid; kink spectrum p95 `2.719e-3`, max `3.040e-3`; DN rel `1.368e-3`; thread=1 median `15.030 ms` | PARTIAL |
| `FAST_THREADS=16 python scripts/bench_fast.py --reps 3 --cases 0` | Formal fast runtime probe | warm median `6.727 ms`, cold `1.919 s`, fallback `0`; speed target not met | PARTIAL |
| `python -m pytest -q` after CI fix and experiment rollback | full regression | `111 passed, 6 deselected in 25.15s` | PASS |
| `python scripts/bench_fast.py --reps 3 --cases 0` with `FAST_THREADS=16` | outer-probe hot-path comparison | plain warm median `5.051 ms`; formal kink probe `8.872 ms`; fallback `0` | PARTIAL |
| scoped ruff + mypy + Chinese comment gate + manifest | CI-maintained checks | all PASS; comments base `93645ec` | PASS |

### Errors

| Error | Resolution |
|-------|------------|
| 计划文件首次补丁未匹配 init-session 生成的简化模板 | 重新读取实际文件后以完整替换方式更新 |

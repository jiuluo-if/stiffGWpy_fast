# Progress Log: stiffgwpy_fast 单一 fast 生产求解器

## Session: 2026-09-09

### Current Status

- **Phase:** 1 - Requirements & Discovery
- **Started:** 2026-09-09

### Actions Taken

- 读取目标文件，确认最终目标为单一 fast 生产 solver，并保留 reference oracle。
- 读取并启用 `using-superpowers`、`self-evolution`、`planning-with-files`、`brainstorming` 与 `writing-plans` 的适用流程。
- 确认分支 `codex/fast_v0.2`、HEAD `b83aa89`，初始工作树无代码 diff。
- 完成首轮关键词和入口检索，发现 README、Cobaya、benchmarks、tests 仍广泛暴露双档位。
- 核对远端：`fast` 指向目标仓库；用户要求重大 Phase 提交并推送远端同名分支。
- Phase A 已加入 uniform-grid 单区间 transfer split 与 one-sided reheating primitive；目标测试覆盖断点属性和单侧积分约定。
- Phase A 的 sparse subset 结果已降级为诊断；benchmark 现默认在候选完整 246 点频率网格上运行独立 reference，避免把 subset 频率积分误报成总 `DN_gw` 误差。
- Phase B 已完成 split transfer segment 的 phase-cap 实现；`phase_max=.25` target test 通过，但单点 oracle 收益有限。
- Phase C 已接入 `freq_grid='goal'`：seed + reheating feature reserve + native `eval_freqs`；测试和 solver native-node 验证通过。
- 正式 fast preset 已切换为单一 `goal-kink-hybrid` 组合；旧 production 保留作 validation/compatibility mode，尚未完成 README/Cobaya 全面清理。

### Test Results

| Test | Expected | Actual | Status |
|------|----------|--------|--------|
| `git branch --show-current` | `codex/fast_v0.2` | `codex/fast_v0.2` | PASS |
| `git diff --stat` | 初始无代码 diff | 无输出 | PASS |
| `git remote -v` | 目标 `fast` 远端存在 | `https://github.com/jiuluo-if/stiffGWpy_fast.git` | PASS |
| `git ls-remote --heads fast` | 远端同名目标存在 | `fast_v0.2` 与 `main` 均在 `b83aa89` | PASS |
| `python -m pytest -q` | 现有回归通过 | `103 passed, 6 deselected` | PASS |
| `python scripts/bench_fast.py --reps 3 --cases 0` | fresh fast runtime 基线 | warm median `2.461 ms`，fallback `0` | PASS |
| `python scripts/benchmark_reference.py --point default --z-tail 5 --no-ode-error --no-tail-error` | fresh oracle 对照 | production fast vs reference 已输出，误差仍超目标 | PASS |
| `python -m pytest -q tests/test_freq_adaptive.py::test_split_primitive_uses_left_limit_at_reheating_kink tests/test_freq_adaptive.py::test_breakpoint_phi_s2_is_accurate_without_dense_subgrid tests/test_fast_sgwb.py::test_kink_split_variant_inserts_only_reheating_breakpoint` | Phase A target tests | 2 passed, 1 passed after endpoint correction; current rerun pending | PASS |
| `python scripts/benchmark_phase_a.py --reps 3 --subset-reference` | Phase A sparse diagnostic | plain `7.389 ms`, kink `8.305 ms`; transition max `1.281% -> 0.776%`; total DN comparison invalid by design | DIAGNOSTIC |
| `python scripts/benchmark_phase_a.py --reps 3` | Phase A full-grid oracle | plain `7.713 ms`, kink `10.020 ms`; DN rel `1.069% -> 0.411%`; spectrum dex p95 `1.004e-2 -> 4.745e-3`; transition max `3.066% -> 1.317%` | PARTIAL |
| `python scripts/profile_fast_breakdown.py --case A --reps 3` and `--kink-split` | Phase A profiler | steady median total `3.701 ms` plain / `7.441 ms` kink at 16 threads; kink exact primitive overhead is measurable | PARTIAL |
| `python scripts/benchmark_phase_a.py --reps 3 --freq-grid goal --phase-max 0.25` | Phase B/C at `h=.02` | 76-point grid; kink spectrum dex p95 `4.792e-3`, DN rel `4.344e-3` | PARTIAL |
| `python scripts/benchmark_phase_a.py --reps 2 --h 0.005 --freq-grid goal --phase-max 0.25` | Phase B/C precision probe | 90-point grid; kink spectrum p95 `2.598e-3`, max `3.040e-3`; DN rel `1.677e-3`; thread=1 median `16.128 ms` | PARTIAL |
| `python scripts/bench_fast.py --reps 3 --cases 0` after fast preset merge | Formal fast runtime probe | median `6.799 ms`, fallback `0`; precision target partly met, speed target not met | PARTIAL |
| `python -m pytest -q` after Phase B/C and preset merge | full regression | `109 passed, 6 deselected in 15.93s` | PASS |

### Errors

| Error | Resolution |
|-------|------------|
| 计划文件首次补丁未匹配 init-session 生成的简化模板 | 重新读取实际文件后以完整替换方式更新 |

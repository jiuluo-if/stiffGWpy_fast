# stiffGWpy fast_v0.2 Tensor Propagation Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task with verification checkpoints.

**Goal:** 在用户指定的 `fast_v0.2` 基线 `465196c82b1c938af7bfce125933a4dfda22f4a2` 上，以严格逐位等价为前提，验证并在确有稳定收益时优化 `SGWB_iter_fast()` 的 tensor propagation 热路径，最后推送到远端 `fast/fast_v0.2`。

**Architecture:** 先用固定资源环境建立五工况 fresh stage profile，再分别在 standalone kernel twin 中验证 P0 repeated-exp elimination 和 P1 counted assembly state；任何 candidate 只有在 bitwise output、guard/failure、determinism 和 full-outer gate 全部通过后才允许最小化接入 production。P2 只根据 Numba LLVM/ASM 证据选择下一候选，不预先改数学公式或启用 `fastmath`。

**Tech Stack:** Python 3, NumPy, Numba/LLVM, pytest, PowerShell, Git.

**Spec:** `C:\Users\联想\.codex\attachments\d413dab6-fcb0-4898-ab67-4e2070a031ae\pasted-text-1.txt`

## Global Constraints

- 工作基线必须是 `465196c82b1c938af7bfce125933a4dfda22f4a2`，目标远端分支为 `fast/fast_v0.2`。
- 不改变 fast profile 的公式、精度参数、输出要求、failure/guard 行为、传播步序或确定性。
- 禁止以 `fastmath=True` 作为本轮生产优化手段；禁止把 relaxed experimental candidate 接入 production。
- baseline/candidate 必须在相同 Numba thread、CPU affinity、threading layer、BLAS budget、warmup/repeat 配置下交替测量。
- P0/P1 先做 standalone kernel bitwise gate；full-outer 依次覆盖 default、lowT、highT、stiff、high_kappa、edge/Sobol 与正式线程环境。
- Git 提交邮箱必须为 `2966684515@qq.com`，提交说明使用“英文：中文内容”格式；有效修改完成后推送目标分支。
- 既有未跟踪实验产物不自动纳入本轮提交；只添加本轮明确产生且需要交付的脚本、artifact、测试或生产改动。

### Task 1: Establish fresh baseline and call-shape audit

**Files:**
- Read: `stiffgwpy_fast/fast_sgwb.py`, `stiffgwpy_fast/exact_background.py`, `stiffgwpy_fast/stiff_SGWB.py`
- Read: `scripts/profile_fast_breakdown.py`, `scripts/_resource_budget.py`, `docs/benchmarks.md`
- Create: `docs/profile_fast_breakdown_round28_fresh_20260923_<case>.json`
- Modify: `task_plan.md`, `findings.md`, `progress.md`

- [ ] Confirm local/remote SHA, branch, clean tracked diff, configured email, and preserve unrelated untracked artifacts.
- [ ] Record the actual formal resources: Numba threads, affinity, threading layer, BLAS environment, warmup policy, repeats, Python/NumPy/Numba versions.
- [ ] Run fresh five-regime stage profiles with the production fast configuration and capture total, tensor kernel, expansion, Phi-S2, frequency preparation, integration, outer iterations, kernel calls, and steps per channel.
- [ ] Read existing 2026-09-17 profiles and rejected artifacts only as historical context; do not compare across resource configurations.
- [ ] Decide P0/P1 priority from the fresh tensor attribution and save the decision before code changes.

### Task 2: P0 exact repeated-exp standalone spike

**Files:**
- Read: `scripts/benchmark_phase_exp_hoist_spike.py`
- Create: `scripts/benchmark_phase_exp_hoist_strict_spike.py`
- Create: `docs/phase_exp_hoist_strict_round28_20260923.json`
- Test: standalone kernel digest/failure contract in the new script

- [ ] Build the baseline by calling the current production `FS.solve_kernel` directly.
- [ ] Build the candidate with the production `z_node`, `z_mid_step`, `z_end`, kink endpoint, phase subdivision, tail matching, assembly, and no `fastmath` changes.
- [ ] Compute the production `z_mid` expression once, reuse one `w_mid = exp(z_mid)` for `n_sub` and only the `n_sub == 1` scaled step, and retain the original substep loop for `n_sub > 1`.
- [ ] Preserve the exact floating-point operation order inside `_scaled_step_with_w` relative to `scaled_step`.
- [ ] Run kernel A/B first; if any required digest differs, report the first mode/k/variable divergence and stop this candidate before outer timing.
- [ ] Only after bitwise equality, run alternating full-outer A/B on the five named regimes, then edge/Sobol and 16/20-thread formal environments with 50-repeat or specified stable counts.
- [ ] Accept only a stable end-to-end gain; otherwise mark REJECTED and keep the candidate standalone.

### Task 3: P1 counted assembly-state standalone spike

**Files:**
- Read: current production `solve_kernel` and `assemble_main` call sites
- Create: `scripts/benchmark_counted_assembly_spike.py`
- Create: `docs/count_assembly_round28_20260923.json`

- [ ] Replace only `k % col_step`/`k // col_step` assembly scheduling with `next_output_k` and `output_slot` state in the twin.
- [ ] Keep tensor state arithmetic, propagation step order, kink interval, z-tail crossing, tail matching, and assembly node set unchanged.
- [ ] Compare per-mode assembly node indices and output digests before measuring speed.
- [ ] Run the same named, edge/Sobol, and formal-thread gates; keep rejected artifacts if there is no stable total-runtime improvement.

### Task 4: P2 LLVM/ASM audit and decision

**Files:**
- Create: `scripts/audit_fast_kernel_llvm.py`
- Create: `docs/fast_kernel_llvm_audit_round28_20260923.json`

- [ ] Compile the canonical Numba signatures and inspect LLVM/ASM for `scaled_step`, `_phase_segment`, and `solve_kernel`.
- [ ] Count ordinary-step `exp`, separate `sin`/`cos` calls, `sqrt`/`ceil`, integer division/modulo, invariant branches, spills/duplicate loads, and whether inlining occurred.
- [ ] Record specializations for `h_arr is None`, `Sv is not None`, and canonical configuration branches.
- [ ] Choose at most one new candidate only when the generated code proves a remaining cost; otherwise stop optimization and report the next hypothesis.

### Task 5: Production gate, review, commit and push

**Files:**
- Modify only the accepted minimal production file(s), if any
- Update: artifact provenance and relevant docs

- [ ] Run focused contract tests, full relevant pytest, compile/import checks, Ruff/mypy scope checks, and `git diff --check`.
- [ ] Use fresh verification after all changes; confirm artifact `commit` equals final HEAD and no unrelated untracked files are staged.
- [ ] Run `awesome-code-review` on the final diff and resolve any correctness/security findings.
- [ ] Configure `user.email=2966684515@qq.com`, commit with English prefix plus Chinese content, push explicitly to `fast fast_v0.2`, and verify exact remote SHA.
- [ ] Report ACCEPT/REJECT for each hypothesis, evidence paths, resource configurations, residual risks, and whether GitHub CI was independently observed.

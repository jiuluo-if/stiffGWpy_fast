# Findings & Decisions: fast 单档位重构

## Requirements

- 当前分支 `codex/fast_v0.2`，远端 `fast` 指向 `https://github.com/jiuluo-if/stiffGWpy_fast.git`；重大 Phase 完成并 fresh verification 后提交、推送到远端同名分支。
- 目标是彻底合并 `fast` / `production`，最终只保留唯一正式 `fast`；`reference` 仅作 precision oracle，validation/debug 仅内部使用。
- Phase A 必须只加入 exact `N_re` breakpoint；保持 `h=0.02`、`z_tail=5`、`phase_max=0`、原频率网格、outer tolerance、`col_step` 不变，并测量 runtime、spectrum/DN_gw 误差、transition 最大误差、400 点 guard/failure。
- 后续优先级：短振荡区间的 phase/WKB envelope、goal-oriented sparse frequency grid、outer self-consistency 加速、解析 stiff/RD branch，最后才做 profiler 证明的底层优化。
- 默认 float64；不得静默 NaN、放松 guard、改变 fallback、物理参数含义、输出字段或 Cobaya derived quantities。
- 精度目标：spectrum median <1e-3、p95 <3e-3、尽量 max <1e-2；DN_gw median <5e-4、p95 尽量 <1e-3；速度不明显慢于当前约 4.4 ms/point，第一阶段目标 ≤4 ms/point。
- 每个 Phase 单独提供 commit、runtime、oracle accuracy、parameter sweep、profiler、A/B 结果及接受/回退结论；最终更新 README、Cobaya YAML、schema、tests、manifest、benchmark docs、CHANGELOG，并生成优化报告。

## Research Findings

- `README.md` 当前仍明确写着“两种 user-facing fast profiles”：`fast` plain-grid 与 `production` transition-refine；用户 API、Cobaya YAML、benchmark 脚本和 tests 仍暴露 `production`/`transition_refine`。
- 当前关键代码入口：`fast_sgwb.py` 的 `normalize_accuracy_mode`/`apply_accuracy_mode`/`resolve_config`/`_SGWB_iter_fast_impl`/`SGWB_iter_fast`；`stiff_SGWB.py` 的 `SGWB_iter`；`reference.py` 的 `run_reference`。
- 当前仓库 HEAD 为 `b83aa89`，最近提交为“统一项目名称为 stiffgwpy_fast”；工作树在审计开始时无未提交代码改动，之后仅生成了计划文件。
- 当前仓库已有 `exact_background.py`，名称与注释表明它支持 continuous-sigma expansion integrals，但是否已接入轻量 breakpoint 仍需沿调用链核实。
- 当前测试已有 `test_modes.py`、`test_engine.py`、`test_fast_sgwb.py`、`test_reference.py`、`test_freq_adaptive.py`、`test_cobaya_adapter.py`，可作为兼容性与 oracle 骨架。
- 文档基线记录 plain-grid warm median `4.442 ms`、transition-refine `21.772 ms`；这些是历史快照，必须重新运行后才能作为本轮证据。
- 本机 fresh baseline（2026-09-09，Python 3.11.9 / NumPy 2.4.4 / Numba 0.67.0 / SciPy 1.17.1，32 CPUs）：`bench_fast.py --reps 3 --cases 0` 的 plain fast warm median `2.461 ms`、min `2.281 ms`、cold `0.200 s`，fallback `0`；这是当前机器快照，不替代历史基线。
- 本机 fresh reference 对照（default、reference `z_tail=5`、30 个 subset bins）：当前 production fast cold `0.313 s`，fast `DN_gw=0.002270462664`；independent reference `DN_gw=0.002796635640`；spectrum linear relative median `3.44e-3`、max `4.06e-1`，dex p50 `1.50e-3`、p95 `3.39e-2`、max `1.48e-1`。该结果再次证明当前 fast/production 仍未达到目标，且低幅度/transition 区域需按 signal mask 与全量指标分别解释。
- 远端 `fast` 实际 refs 为 `refs/heads/fast_v0.2` 和 `refs/heads/main`，均指向 `b83aa89`；本地分支名为 `codex/fast_v0.2`。因此推送目标应明确使用 `fast HEAD:fast_v0.2`，不能假定本地 upstream 已配置。
- `scripts/benchmark_phase_a.py` 原先固定使用 31 点 reference subset；该 subset 会改变总 `DN_gw`，不能作为总量 oracle。脚本现默认以候选完整频率网格调用独立 reference，稀疏 subset 仅通过显式参数保留作诊断。
- `exact_background` 的 reheating primitive 已修正为保留 `N_re` 的左右单侧 sigma 极限：结束于断点的 Simpson panel 用 pre-transition `sigma=1`，开始于断点的 panel 用 post-transition 值。
- 完整 246 点 reference（固定候选 `DN_eff`，`z_tail=5`, `rtol=1e-11`）显示 Phase A plain -> kink split：`DN_gw` 相对误差 `1.069% -> 0.411%`，spectrum dex median `2.653e-3 -> 1.690e-3`，p95 `1.004e-2 -> 4.745e-3`，max `1.311e-2 -> 5.925e-3`，transition max `3.066% -> 1.317%`。这是有效的完整频率网格结果，但尚未满足最终 spectrum p95 `<3e-3`。
- Phase A thread=1 median timing（fixed `h=.02,z_tail=5,phase_max=0,col_step=8`）为 plain `7.713 ms`、kink `10.020 ms`；16-thread profiler 的 JIT 后 steady sample 为 plain `3.701 ms`、kink `7.441 ms`，主要额外成本落在 exact background primitive，而非 tensor kernel。

## Technical Decisions

| Decision | Rationale |
|----------|-----------|
| 先定位真实数据流再设计 split | 目标要求精确计算 `N_re` 且不得让 transfer step 跨 kink；不能凭文件名猜测已有实现 |
| 设计阶段必须区分“接口合并”和“数值算法变更” | 这样可以保持物理契约，且每一项算法收益可独立归因 |

## Issues Encountered

| Issue | Resolution |
|-------|------------|
| `create_goal` 返回已有未完成 goal | 使用当前 active goal；不创建第二个 goal |
| 远端 `fast` 尚未在本地分支跟踪信息中显示 | 已确认 remote 名称和 URL；下一步核对远端同名 ref |

## Resources

- `C:\Users\联想\.codex\attachments\79e36cae-8136-418c-9325-fde68b9fc4ff\goal-objective.md`
- `F:\codex\stiffGWpy\stiffgwpy_fast\fast_sgwb.py`
- `F:\codex\stiffGWpy\stiffgwpy_fast\stiff_SGWB.py`
- `F:\codex\stiffGWpy\stiffgwpy_fast\reference.py`

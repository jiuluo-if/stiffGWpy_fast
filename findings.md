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
- 远端 `fast` 实际 refs 为 `refs/heads/fast_v0.2` 和 `refs/heads/main`；审计初始快照均在 `b83aa89`，当前 `fast_v0.2` 已推进到 `4436b9d`。本地分支名为 `codex/fast_v0.2`，推送目标必须明确使用 `fast HEAD:fast_v0.2`。
- `scripts/benchmark_phase_a.py` 原先固定使用 31 点 reference subset；该 subset 会改变总 `DN_gw`，不能作为总量 oracle。脚本现默认以候选完整频率网格调用独立 reference，稀疏 subset 仅通过显式参数保留作诊断。
- `exact_background` 的 reheating primitive 已修正为保留 `N_re` 的左右单侧 sigma 极限：结束于断点的 Simpson panel 用 pre-transition `sigma=1`，开始于断点的 panel 用 post-transition 值。
- 完整 246 点 reference（固定候选 `DN_eff`，`z_tail=5`, `rtol=1e-11`）显示 Phase A plain -> kink split：`DN_gw` 相对误差 `1.069% -> 0.411%`，spectrum dex median `2.653e-3 -> 1.690e-3`，p95 `1.004e-2 -> 4.745e-3`，max `1.311e-2 -> 5.925e-3`，transition max `3.066% -> 1.317%`。这是有效的完整频率网格结果，但尚未满足最终 spectrum p95 `<3e-3`。
- Phase A thread=1 median timing（fixed `h=.02,z_tail=5,phase_max=0,col_step=8`）为 plain `7.713 ms`、kink `10.020 ms`；16-thread profiler 的 JIT 后 steady sample 为 plain `3.701 ms`、kink `7.441 ms`，主要额外成本落在 exact background primitive，而非 tensor kernel。
- Phase B `phase_max=.25` 在 `h=.02` 下对完整 oracle 几乎不改变 spectrum p95（kink `4.7629e-3`），但已修复 split segment 未应用 phase cap 的逻辑缺口。
- Phase C goal grid 的 seed80 历史结果为 90 点、spectrum dex p95 `2.598e-3`、DN rel `1.677e-3`；当前 seed64 结果与节点/误差取舍见 Latest evidence。
- 当前正式 `fast` preset（`h=.005,col_step=8,z_tail=5,phase_max=.25,freq_grid=goal,kink_split=True`）在本机 16 threads 为 warm median `6.727 ms`、cold `1.919 s`；高于原先 plain `2.654 ms` 与第一阶段 `<=4 ms` 目标，必须继续 profiler/outer iteration 优化，不能宣称速度达标。

## Latest evidence (2026-09-09)

- 外层首轮只需计算最终频率列来得到 `DN_gw`，不需要写入所有中间 e-fold 列；已改为首轮 `assemble=0`、后续完整轮 `assemble=1`。新增回归测试确认首轮/最终轮调用顺序，最终数组摘要保持一致。
- 该热路径改动的 like-for-like plain benchmark（16 threads，A，3 repeats）warm median 从旧记录 `6.727 ms` 降至 `5.051 ms`；仍未达到 `<=4 ms/point`，不能视为速度目标完成。正式 kink 路径仍需继续优化。
- 本轮新增代码经完整回归 `111 passed, 6 deselected`、ruff、mypy、中文注释门禁和 manifest 校验通过；一次未形成稳定收益的 sigma 节点复用实验已回退。

- CI run `34303864035` 的五个 test matrix job 均在“中文注释门禁”失败；根因是 `78cc5b9` 新增 5 行英文 `#` 注释。仅补充中文语义后提交 `4436b9d` 并推送 `fast_v0.2`，CI run `34304423555` 的 Python 3.9–3.13、ruff、mypy、pytest、manifest、wheel、distribution、smoke 和 Cobaya 全部 PASS。
- 正式 fast 的 seed64 A/B（`h=.005,z_tail=5,phase_max=.25,kink_split,goal`，76 点，full candidate-grid reference）为 spectrum dex median `6.528e-4`、p95 `2.719e-3`、max `3.040e-3`，DN rel `1.368e-3`；seed80 的 90 点结果为 p95 `2.598e-3`、DN rel `1.677e-3`。seed64 保留为较低节点且 DN 更好的当前基线，但最终 DN gate 仍未满足。
- 当前 seed64 formal warm benchmark（16 threads，3 repeats，default A）为 median `6.727 ms`、cold `1.919 s`、fallback `0`；仍高于 `<=4 ms/point` 目标。
- 尾部隔离（reference 固定 `z_tail=5`）显示 fast `z_tail=3/4/5` 的 spectrum p95 分别约 `2.10e-2/5.07e-3/2.72e-3`，DN rel 约 `1.48e-2/3.31e-3/1.37e-3`；不能以浅尾部换速。
- 四阶两点 Magnus A/B：`h=.01` 下 DN rel 可到 `1.90e-4`，但 spectrum p95 `3.998e-3`；`h=.005` 热路径约增至 32.3 ms（thread=1）而精度无实质收益，已回退。基于 `Phi` 曲率的局部子步也未改善 h=.01 谱 p95，已回退。
- 高层 API 合并提交 `c9110c9`：`accuracy_mode='fast'` 成为唯一正式用户档；旧 `production`/`ultra-fast`/transition-refine 别名在高层发出 `DeprecationWarning` 并映射到 fast，底层 validation preset 仍可供内部脚本使用。完整测试 `111 passed, 6 deselected`，CI run `34313611135` 全部 PASS。
- rejected numerical A/B：从视界前 3 个十进制数量级提前到 4 个只改变 DN 约 `3e-9`；标准四阶 RK4、两节点四阶 Magnus、视界局部 2/4 子步和 sigma 变化触发的局部半步均未缩小 `h=.005` 与 `h=.0025` 的 DN 差异，且部分方案变慢，均已回退。
- rejected sparse goal refinement：76 点增加到 87 点后，spectrum dex p95 从 `2.719e-3` 改善到 `2.675e-3`，但 integrated DN rel 从 `1.368e-3` 恶化到 `1.755e-3`，已回退。32 个 h/2 探针通道可以估计默认点的 h 偏差约 `7.35e-4`（完整 h/2 对照 `7.44e-4`），但单次探针仍约 `10 ms` 且会造成 endpoint-only 修正，尚未形成可接受的输出一致性方案。
- 尾部匹配 A/B：当前 `gamma=(3-1.5*sigma)/2` 是阻尼去除变量的系数；由物理张量解 `T~a^-1` 推导，深亚视界振幅匹配的首阶组合应为 `x+y/omega`，即 `gamma=1`。只改这一项，在默认点 16 个代表 native 节点上，z_tail=5 的 spectrum dex p95 由约 `2.63e-3` 降至 `3.20e-4`，max 由约 `2.66e-3` 降至 `4.64e-4`，warm fast 约 `9-10 ms`，无新增计算分支；保留该修改。
- 尾部修正跨参数探针：默认/high-T_re/stiff 的 30 个 native 节点 spectrum dex p95 分别约 `9.59e-4/9.14e-4/1.02e-3`；low-T_re 约 `3.35e-3`，误差集中在 `log10(f/Hz)=-18..-16.5`，这些模式未进入 analytic tail。low-T_re 的 h=.005 -> .0025 将该局部 p95 约 `3.42e-3 -> 1.55e-3`，但局部 phase 二分、起点提前、终点直接取 Phi_grid[k+1] 均未带来稳定收益，均回退。
- 尾部修正改变了旧 validation 测试中的固定数值口径：production transition-refine 的 continuous-sigma 分支为 `DN_eff=0.00227259`，普通 grid 分支为 `0.00225706`；测试现改为要求 continuous-sigma 分支接近 reference 且优于普通 grid，避免沿用旧 tail match 的过时绝对阈值。
- 尾部匹配修正 commit `da4fd02` 已按 `2966684515@qq.com` 提交并推送到 `fast/fast_v0.2`；远端 CI run `34318085109` 的 Python 3.9–3.13、ruff、mypy、manifest、wheel、distribution、smoke 与 Cobaya 全部 PASS。
- 背景节点缓存 A/B：`exact_phi_s2_split` 现在可直接复用 `m.sigma` 的节点值，同时保留原有形状检查和 reheating 左右单侧约定；与原始 `sigma_vec` 路径逐项比较到 `2e-13` 以内。formal kink profile 的 steady warm median 约由 `10.62 ms` 降至 `8.90 ms`，但仍高于 `<=4 ms` 目标。
- outer self-consistency 审计：default、low-T、high-T、stiff 和 low-r 五个代表点的首轮 probe 到第二轮 full solve 的 `DN_gw` 绝对变化分别约为 `7.7e-14`、`4.3e-9`、`7.7e-12`、`6.2e-13`、`7.8e-16`；说明 DN 本身对外层更新很不敏感。但临时强制“一次 full solve 即结束”时，最终 `g2/w2` 与当前两轮流程仍出现差异，high-T 点 `g2` 总和约差 `3.0%`，因此不能仅凭 DN 收敛就删除第二次 full solve，当前 outer probe 方案保留。
- rejected background batching：将 exact split 的 midpoint 与 quarter-point sigma 合并为一次更大的向量调用，数值结果保持等价，但 formal kink profile 的稳定 warm median 约由 `8.4 ms` 上升至 `10.6 ms`；已回退，不纳入主线。线程扫描 `1/2/4/8/16/32` 的 warm median 约为 `13.4/9.96/9.90/7.08/8.76/9.91 ms`，稳定 profiler 的 8/16 线程均约 `8.4 ms`，暂不改变默认线程数。
- formal kink preparation 优化：正式 kink 路径不再先计算随后被 exact primitive 覆盖的 uniform-grid `Phi/S2` spline，只保留频率起点、`f_hor` 和 tail 所需准备；新增测试确认两轮 outer solve 均使用 frequency-only preparation。与临时恢复 full preparation 的同环境 10 次 profiler 对照，steady warm median 约 `9.68 -> 8.98 ms`，所有 spectrum/DN/g2/w2 digest 保持一致；保留该修改。

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

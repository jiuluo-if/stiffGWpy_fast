# stiffgwpy fast/plain-grid 性能工程报告

日期：2026-09-03  
基线：`80af01b`（solver 与当前 HEAD 相同）  
环境：Windows，Python 3.11.9，NumPy 2.4.4，Numba 0.67.0，SciPy 1.17.1；基准命令固定 `FAST_THREADS=4`。

## 结论

本轮发现并修复了一个实现级热点：`solve_kernel` 外层缺少 Numba JIT 装饰器，导致 `prange` 在 Python 层执行。补上 `@njit(parallel=True, cache=True)` 后，fast/plain-grid 获得数量级加速；公式、步长、频率网格、outer criterion 和公开行为均未改变。

## 基线与 breakdown

使用 7 次 warm 重复，首个调用单独记录 JIT/cold。baseline exact（`80af01b`）默认点、4 threads：

| 项目 | median | p95 | 说明 |
|---|---:|---:|---|
| warm total | 398.9 ms | 2386.6 ms | p95 受 Windows/JIT/调度噪声影响 |
| expansion/background | 0.375 ms | 470.0 ms | 每次 outer iteration 2 次调用 |
| frequency grid | 0.573 ms | 0.694 ms | 每次 outer iteration 2 次调用 |
| kernel prepare | 0.238 ms | 481.3 ms | 常态低于 1 ms |
| tensor solve kernel | 396.7 ms | 566.8 ms | 主要热点，约 99% |
| column integration | 0.203 ms | 312.3 ms | 最终一次调用 |

首个 cold call 总计约 3.24 s，其中包含 Numba JIT；未混入 warm runtime。

修复后的 public benchmark（`FAST_THREADS=4`，15 次 warm 重复）：A 点 cold `0.253 s`、warm median `6.879 ms`、p95 `7.642 ms`；B 点 warm median `4.085 ms`、p95 `4.684 ms`。相对 baseline warm median A `398.9 ms`、B `200.6 ms`，speedup 约 `58.0x` / `49.1x`。修复后 cold 仍单独计量，未混入 warm。

线程探索的 3 次 warm 小样本受并发测量噪声影响，不能作为正式 speedup：1/2/4/8/16 threads 的 median total 分别约 528/496/521/389/522 ms。它没有显示单调扩展，正式结论以固定线程、串行复测为准；不能声称更多线程必然更快。

## 候选 AB 结果

### 候选 1：复用 `exp(Phi_mid)`，减少热循环 `exp`

- hotspot：`solve_kernel` 每个 channel/step 的 `exp(z_mid)`。
- 结果：失败并回退；baseline median `398.9 ms`，candidate median `543.3 ms`，慢约 36%。新增 Numba 函数层未带来收益。
- 数值：因性能门已失败，未推进到完整 oracle AB；候选已删除，最终 solver 无此差异。
- 风险：中等（改变浮点计算路径）；不采用。

### 候选 2：tail coarse-slot while 改为等价整数索引

- hotspot：每个 mode 的 tail handoff bookkeeping。
- 结果：失败并回退；一次 7-repetition 实测 candidate median `437.1 ms`，高于 baseline exact `398.9 ms`，理论可省的工作占比过小。
- 数值：代数上等价，未形成可接受性能收益；候选已删除。
- 风险：低；不采用。

## 采用的优化：JIT 并行化外层 tensor kernel

- change：为 `solve_kernel` 添加 `@njit(parallel=True, cache=True)`。
- hotspot：Python 外层 channel loop；原有 `prange` 未进入 Numba 编译。
- baseline runtime：A `398.9 ms`、B `200.6 ms` warm median。
- candidate runtime：A `6.879 ms`、B `4.085 ms` warm median。
- speedup：A 约 `58.0x`、B 约 `49.1x`。
- numerical max diff：A/B 的 `DN_eff`、频率、完整 spectrum、`DN_gw`、g2、w2 digest 均逐位相同；max diff 为 `0`（逐位比较）。
- tests：定向 fast/engine、A/B baseline digest、线程 scaling 已通过；全量门禁见下节。
- risk level：低到中；改变执行层和并行调度，但每个 channel 独立写入，未改变 channel 内浮点顺序。

线程 scaling（A 点，warm median，15 次；每个线程数独立进程）：1 `6.123 ms`，2 `4.124 ms`，4 `3.568 ms`，8 `3.386 ms`，16 `3.233 ms`。本机在允许范围内没有反向变慢。

### 调度 chunk-size 实验

为验证 `prange` 负载均衡，额外测量了 `FAST_CHUNKSIZE=1/2/4/8/16/32/64`。短样本曾显示 A 点 chunk16 约有 8% 优势，但长样本串行复测不支持全局采用：A 点 baseline 默认 chunk 的 median `401.3 ms`、chunk16 `405.5 ms`；B 点分别为 `198.2 ms`、`192.9 ms`。方向不一致且接近机器调度噪声，因此没有改变 solver 默认调度；chunk 仅作为 profiler 的实验参数保留。

### coarse-column assembly 上限与 history 重建实验

临时关闭非最终 coarse-column assembly 后，A 点 kernel 约从 `408 ms` 降至 `307 ms`，B 点约从 `187 ms` 降至 `149 ms`，说明 assembly 数学本身约占 20–25%。随后实现了“保存每个 coarse slot 的 x/y，收敛后一次性重建”的候选；A 点实际 warm median `437.9 ms`、B 点 `201.8 ms`，由于 history 写入、清零和最终重建开销而失败，已完整回退。该实验没有作为数值结果使用。

随后将该方案改为专用的无 optional-参数 deferred kernel，避免通用 kernel 的 optional 分支。结果仍未改善：A 点 warm median `400.2 ms`、B 点 `200.7 ms`，与 baseline 基本持平；完整输出的 history 写入成本约等价于原 assembly。候选已回退，说明下一步不能继续在相同 coarse-state 存储策略上叠加复杂度。

## 数值与兼容性

本轮唯一 solver 变化是执行层的 Numba 装饰器；A/B fast 与 production AB 的频率、完整 spectrum、`DN_gw`、g2、w2 均逐位一致。已有 plain-grid matched-reference artifact（9 点）记录的 envelope 为：signal relative median `1.867e-2`、max `7.019e-2`；integrated `DN_gw` relative median `9.142e-3`、max `2.725e-2`。本轮没有修改 reference、validation artifact、preset、alias、Cobaya adapter 或 telemetry。

## 测试与构建

- `python -m pytest`：JIT 候选最终门禁 `99 passed, 6 deselected`。
- fast/engine 定向回归：`29 passed, 3 deselected`。
- `python -m pytest -m cobaya`：JIT 候选最终门禁 `1 passed, 104 deselected`。
- `ruff check .`：未通过；仓库既有约 399 条 lint 问题，非本轮 profiler 文件引入。新增 profiler 单文件需单独检查。
- `wheel build`：JIT 候选最终门禁成功。
- plain-grid validation：已有 artifact envelope 如上；本轮未重新消耗完整连续-sigma oracle 运行，不能声称本轮 validation PASS。

## 修改文件

- `scripts/profile_fast_breakdown.py`：新增阶段级 profiler/benchmark 工具。
- `docs/profile_*`：保存 baseline、候选和 chunk-size 实验记录。
- `docs/performance_fast_plain_20260903.md`：本报告。
- `stiffgwpy/fast_sgwb.py`：为原有 `solve_kernel` 添加 `@njit(parallel=True, cache=True)`。
- solver 公式、公开 API、preset、alias、Cobaya 配置均未修改。

## 后续方向（本轮未采用）

下一轮应在独立分支中针对 `tensor_solve_kernel` 做编译级实验：检查 Numba IR 是否内联 `assemble_*`，评估按 work 量重排 channel 对 `prange` 负载均衡的影响，并用独立进程串行测量 1/2/4/8/16 threads。任何改变浮点 evaluation order 的 kernel fusion、exp 分解或 loop 重排，都必须先完成完整 spectrum、`DN_gw`、g2/w2、kappa_r、迭代次数和 failure telemetry AB，再决定是否采用。

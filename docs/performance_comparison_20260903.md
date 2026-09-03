# stiffgwpy 优化前后全方法对比报告

日期：2026-09-03  
基线：`80af01b`  
候选：`fb54190`  
环境：Windows，Python 3.11.9，NumPy 2.4.4，Numba 0.67.0，SciPy 1.17.1；fast benchmark 固定 `FAST_THREADS=4`。

## 1. 结论

本轮唯一采用的性能修改是在现有 `solve_kernel` 外层增加
`@njit(parallel=True, cache=True)`。原有 channel-level `prange` 因缺少外层
JIT 时实际在 Python 层执行；修复后进入 Numba 并行 kernel。公式、preset、
频率节点生成、outer `Delta N_eff` 收敛条件、返回对象、fallback 和 telemetry
均未改变。

具体代码变化只有一行：

```python
@njit(parallel=True, cache=True)
def solve_kernel(...):
    ...  # 原有 channel-level prange 与计算顺序保持不变
```

也就是说，本轮没有放大 `h`、减少频率节点、降低 `tol`、提前 tail、使用
`float32`、删除物理项或改变 outer convergence criterion；优化目标是让已有
`prange` 真正进入 Numba 编译并行执行。

默认 fast/plain-grid 点的 warm 单点 median 从 `398.9 ms` 降至 `4.442 ms`
（约 `89.8x`）；第二个 validation point 从 `200.6 ms` 降至 `2.501 ms`
（约 `80.2x`）。candidate 的首个 A 点 cold call 单独计量为 `0.325 s`，不混入 warm
runtime。production 同一 kernel 的当前测量为 cold `0.226 s`、warm median
`21.772 ms`、p95 `22.149 ms`。

## 2. 方法和兼容性范围

| 用户入口 | 实际路径 | 本轮处理 |
|---|---|---|
| `engine='fast', accuracy_mode='fast'` | fast plain-grid | 优化，preset 不变 |
| `plain_grid` / `plain-grid` / `plain` | `fast` 别名 | 保持兼容，输出指纹相同 |
| `ultra-fast` | `fast` 历史别名 | 保持兼容，输出指纹相同 |
| `accuracy_mode='production'` | transition-refine fast | 同受 kernel JIT 影响，行为不改 |
| `transition_refine` / `transition-refine` / `tr` | `production` 别名 | 保持兼容 |
| `engine='lsoda'` | 原始 LSODA 路径 | 未修改 |
| `engine='reference'` | 独立 continuous-sigma + DOP853 oracle | 未修改 |
| `debug` / `deep` / fast 内部 `reference` | 验证档位 | 未修改，非新增生产档位 |

`fast` preset 仍为 `h=0.02`、`col_step=8`、`z_tail=5.0`、`freq_res=1.0`、
`tol=1e-6`、`transition_refine=False`、`phase_max=0.0`、`freq_grid=construct`。

## 3. Runtime 对比

### 3.1 fast/plain-grid：默认点和已有 plain-grid 点

| 点位 | baseline warm median | candidate cold | candidate warm median | candidate p95 | speedup |
|---|---:|---:|---:|---:|---:|
| A：`r=1e-2, T_re=2e3, kappa10=1e-2` | 398.9 ms | 0.325 s | 4.442 ms | 5.105 ms | 89.8x |
| B：`r=1e-3` | 200.6 ms | 0.003 s* | 2.501 ms | 3.060 ms | 80.2x |

candidate 使用当前提交 `b0909d0` 重新运行的 15 次 warm 重复，原始记录为
`docs/benchmark_current_fast_preset_20260903.json`；speedup 采用 median，不采用最佳一次。
`*` B 点的 cold 字段发生在同一进程 A 点之后，Numba cache 已可用，不能解释为
独立的全冷启动；全冷启动以独立首个 A 点记录为准。
baseline 的 breakdown profiler 使用 7 次重复；其较大的 p95 受到独立进程
首次 JIT/Windows 调度噪声影响，因此 public benchmark 的 15 次 p95 作为
单点尾延迟口径。

### 3.2 warm batch

在 12 个既有 benchmark 参数点上先完成一次 warm-up，再连续执行一轮：

| batch | 总耗时 | 每点平均 | 对照 |
|---|---:|---:|---|
| candidate fast/plain-grid，12 点 | 28.205 ms | 2.350 ms/点 | 同一进程、4 threads |
| candidate fast/plain-grid，单点重复 15 次 | — | median 3.113 ms/点，p95 3.603 ms | 点位 A |

batch 平均低于单点重复是参数点负载不同和连续调度的结果，不作为额外
算法 speedup；它说明 JIT kernel 适合扫描场景。

### 3.3 全部引擎/档位的现行量级

| 方法 | 当前 runtime 证据 | 数值定位 |
|---|---:|---|
| fast/plain-grid | warm median `4.442 ms`（A） | 探索档；plain-grid oracle envelope 未改变 |
| fast/production | warm median `21.772 ms`（A），p95 `22.149 ms` | 科学生产档；仍需按现有 oracle/误差预算解释 |
| LSODA | 近期 A 点 `22.137 s` | 回归/fallback/历史 runtime anchor，未修改 |
| independent reference | 历史 full-oracle 约 `360–383 s/点` | 精度锚点，未修改、不进 MCMC 热路径 |

### 3.4 速度对比：只与 LSODA 比

以下是同一当前提交、同一机器、同一参数点的 `bench_fast.py` 结果；LSODA
没有被优化，fast 的 cold/JIT 不计入 warm speedup：

| 点位 | LSODA median/单次 | fast warm median | LSODA / fast |
|---|---:|---:|---:|
| A | `21.117182 s` | `4.442 ms` | `4753.87x` |
| B | `8.621666 s` | `2.501 ms` | `3446.87x` |

这里的速度结论不使用 independent reference 的耗时，也不使用 fast 与
reference 的精度差异推导 speedup。旧的 `0.37 s`、`3.7–4.1 s` 和 `~1000x` 只代表 JIT 修复前的历史口径，
已从活动文档移除并保留在 archive 的历史材料中。

### 3.5 线程 scaling

candidate A 点，独立进程、15 次 warm 重复的 median：

| FAST_THREADS | 1 | 2 | 4 | 8 | 16 |
|---:|---:|---:|---:|---:|---:|
| runtime | 6.123 ms | 4.124 ms | 3.568 ms | 3.386 ms | 3.233 ms |

在本机允许范围内未观察到反向变慢。`FAST_CHUNKSIZE` 长样本实验方向不一致，
未改变默认调度。

## 4. Runtime breakdown 和热点排序

baseline A 点 7 次 warm profiler 的 median：

| 阶段 | median | 占比/结论 |
|---|---:|---|
| tensor solve kernel | 396.7 ms | 主要热点，约 99% |
| frequency-grid construction | 0.573 ms | 非热点 |
| expansion/background | 0.375 ms | 非热点 |
| kernel prepare | 0.238 ms | 非热点 |
| column/integration | 0.203 ms | 非热点 |

因此本轮没有猜测式地改 allocation、背景缓存、积分公式或 outer iteration；
先由 profiler 锁定 Python 外层 channel loop，再只改变执行层。P1/P2 审计中
尝试的 `exp` 复用、tail index 化、helper inline、step normalization、
channel reorder、chunk 调度和 deferred assembly/history 均未形成稳定收益，
已回退。

## 5. 精度对比：只与 independent reference 比

精度基准是 `stiffgwpy.reference` 的 continuous-sigma + DOP853 独立实现，
不是 LSODA。已有 matched-reference 结果：

| fast 档位 | 对 independent reference 的结果 | 结论 |
|---|---|---|
| plain-grid | signal relative median `1.867e-2`、max `7.019e-2`；integrated `DN_gw` median `9.142e-3`、max `2.725e-2` | 仅探索档，未通过 `1e-3` science gate |
| production/transition-refine | signal relative max `7.09e-4`（dex max `3.1e-4`）；integrated `DN_gw` median `4.3e-4`、p95 `1.2e-3` | signal gate `<1e-3` 通过；integrated `<1e-4` 未通过 |

这些是物理方法之间的精度结果，不是 JIT 优化造成的误差；本轮 candidate
与 baseline 的 AB 结果如下。

## 6. 数值和行为 AB

candidate 与 `80af01b` 的 A/B fast 输出逐位一致：

- `DN_eff`、returned frequency nodes、完整 `Omega_GW(f)`、`DN_gw`、`g2`、`w2`：digest 全部相同，max diff 为 `0`；
- production/transition-refine：`DN_eff=0.002261731150563835`、`n_freq=247`，`f`、`log10OmegaGW`、`DN_gw`、`g2`、`w2` 的 SHA-256 全部相同；
- fast failure、physical guard、LSODA fallback telemetry：公开 smoke test 中 `fast_evals=1`、`fast_failures=0`、`lsoda_fallbacks=0`，各 alias 与 canonical fast 输出相同；
- 未修改 `LCDM_SG.SGWB_iter(engine='fast', accuracy_mode='fast')`、Cobaya 参数/YAML、array shape、frequency ordering 或 derived quantities。

已有 plain-grid matched-reference oracle envelope 仍为 signal relative median
`1.867e-2`、max `7.019e-2`，integrated `DN_gw` relative median `9.142e-3`、
max `2.725e-2`。本轮没有修改 reference 或 validation artifact，也没有重新
消耗完整 continuous-sigma oracle，因此这里报告“既有 envelope 未改变”，不把
它表述为本轮重新跑出的 oracle PASS。

## 7. 验证状态

| 检查 | 结果 |
|---|---|
| `python -m pytest` | `99 passed, 6 deselected` |
| `python -m pytest -m cobaya` | `1 passed, 104 deselected` |
| wheel build | `python -m build --wheel` 成功 |
| `ruff check .` | 未通过；仓库既有约 399 条 lint 问题，本轮未扩大范围 |
| plain-grid oracle | 既有 artifact envelope 未改变；本轮未重跑完整 oracle |

## 8. 本轮修改和风险

采用修改只有 `stiffgwpy/fast_sgwb.py` 的一行 JIT 装饰器；本报告和活动文档
同步更新了现行 runtime 口径。风险为低到中：执行顺序只在独立 channel 之间
并行，channel 内浮点 evaluation order 未改变；Numba cache/cold startup
会影响首次调用延迟，因此生产部署和 benchmark 必须分开报告 cold 与 warm。

## 9. 后续方向

本轮不采用：继续做 kernel fusion、按 work 量重排 channel、改变 chunk policy、
减少中间 array 或调整 outer iteration。这些方向若改变浮点 evaluation order，
必须建立独立 commit，完成全输出、迭代次数、guard/failure telemetry 和
plain-grid matched-reference AB 后再决定。

复现入口：`scripts/bench_fast.py`、`scripts/profile_fast_breakdown.py`，以及
`docs/performance_fast_plain_20260903.md` 的实验记录。

# 测试重复审计

本审计以当前 `fast_v0.2` 工作树为准，目标是减少重复求解，而不是减少独立风险覆盖。

## 分层结果

| 类别 | 主要文件或入口 | reference/LSODA | 责任边界 |
| --- | --- | --- | --- |
| 兼容性 | `test_compatibility_smoke.py`, `test_resource_budget.py` | 否 | 导入、配置、quadrature、资源默认和 tiny fast smoke |
| fast canonical | `test_fast_sgwb.py`, `test_engine.py`, `test_freq_adaptive.py`, `test_physics_limits.py` | 仅显式小交互 | solver invariant、状态恢复、grid、guard、determinism |
| reference unit | `test_reference.py` | 仅轻量解析检查 | background、积分、tail summary；长 reference 标为 slow |
| Cobaya | `test_cobaya_adapter.py` | 独立 job | adapter API、derived names、fallback telemetry、一次真实导入 smoke |
| science certification | `tests` 中的 `slow`、slow workflow、oracle scripts | 是 | LSODA/reference/尾部收敛，不进入版本矩阵重复执行 |

## 重复性结论

- 5 个 Python 版本只运行兼容性 smoke，不重复 canonical numerical、reference、build 或 Cobaya。
- canonical regression 只在 Python 3.11 执行一次；Cobaya 单独执行一次，避免把集成依赖重复装入其他版本。
- `eval_freqs` 的 native-node 约束与 DN 不变性分别由一个 canonical 测试和一个 interaction 测试负责。
- thread-state restore、A→B→A 参数隔离、physical guard 分类和 repeated-call determinism 已有独立断言；没有再添加同一 full solve 的参数笛卡尔积。
- 资源契约测试只检查环境与 telemetry，不启动 reference 或额外 ODE solve。

## 资源审计

默认 correctness/reference 路径使用 Numba/fast 线程 2、BLAS 类线程 1、reference workers 1；benchmark 脚本的显式 `--threads` 或 `--workers` 才能提高预算。Cobaya 真实集成 smoke 是有固定 180 秒上限的独立例外，CI 显式使用 `SGWB_POOL_SIZE=2`，不改变库默认值。benchmark artifact 记录 logical CPU、affinity、线程层、BLAS caps、workers、版本、commit 和并发进程数。

本文件是执行结构审计，不把 benchmark 数值结果当作 correctness 证据；数值认证仍必须读取带 commit/schema 的独立 artifact。

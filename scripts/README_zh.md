# 仓库脚本

[English](README.md) | 中文说明

本目录包含验证、基准、性能剖析、打包和历史诊断脚本。请从仓库根目录运行，确保程序包导入和产物路径按说明解析。

## 脚本分类

| 前缀或文件 | 用途 |
|---|---|
| `validate_*.py` | 重放 fast/reference 对比、参数边缘、参数扫描和档位检查 |
| `build_*.py` | 根据已有产物汇总验证矩阵和 manifest |
| `bench_*.py`、`benchmark_*.py` | 测量运行时间、数值精度或指定诊断假设 |
| `profile_*.py` | 分解求解器阶段耗时，或测量指定路径的属性/调用开销 |
| `smoke_*.py`、`verify_distribution.py` | 检查安装后 wheel 的资源和打包边界 |
| `plot_*.py` | 从已保存的验证或基准数据绘图 |
| `_resource_budget.py`、`_smoke_probe.py` | 线程资源限制和隔离冒烟检查的共享工具 |

`validate_fast.py`、`validate_random.py` 和 `convergence_study.py` 均按参照谱
（`log10(Omega_GW) >= -30`）选择信号区频点；线性相对误差以候选谱为分母。

`benchmark_phase_a.py` 针对一个内置参数点，对比普通固定网格与精确再加热转折点切分。
频谱相对误差统一调用共享函数，并以候选谱作分母。默认参考计算使用候选的完整频率网格；
`--subset-reference` 是可能改变参考解 `DN_gw` 的历史稀疏网格诊断。该脚本只覆盖指定比较，
不能作为全参数空间精度认证。

`check_random_freq.py` 在一个固定参数点比较两个引擎，标出 10 个带种子的频率通道并同时报告全网格残差。
这 10 个通道仅用于展示，不是 10 个独立参数案例。全网格相对误差包括近零频点，应结合 dex 差解读；
PNG 默认写入 `scripts/`，也可通过输出路径参数指定其他位置。

`bench_fast.py` 对比一次 LSODA 耗时和 fast 冷/热调用时间，只用于运行时间与回归观察，不比较求解精度。
计时调用使用 fast 档位的网格和求积方法，容差与频率分辨率可通过 `--tol`、`--freq-res` 覆盖。
默认计时不切分精确再加热转折点；`--kink-split` 会加入该转折点。单独的安全包装器诊断始终走正式 fast 的转折点切分路径，
并在 JSON 中记录其配置。P95 使用线性分位数插值，默认仅重复 3 次时只能粗略观察耗时尾部。

`benchmark_transition_frontier_ab.py` 在两个 `DN_re` 边界上交替运行正式 fast 与匹配的转折细化配置，
记录每次求解耗时及各配置内部的输出摘要。它本身不认证精度；精度需查看匹配的 Oracle C 记录。

许多 `benchmark_*_spike.py` 是历史独立候选方案测试；它们不会选择或修改生产求解路径。运行前先看脚本的 `--help` 和输出路径。保留带日期的实验产物；新实验应使用新文件名，除非确实要替换原文件。

SageNet 对比脚本还需要单独的 SageNet 代码仓库和运行环境；本项目的 extras 不会安装 SageNet。

## 可复现实验

按研究问题使用对应的环境、线程预算、随机种子、预热方式和重复次数。冷启动/JIT 开销与预热耗时分开记录；数值精度比较与速度比较也分别报告。记录代码提交。详见[`docs/reproducibility.md`](../docs/reproducibility.md)、[`docs/benchmarks.md`](../docs/benchmarks.md)和[`docs/experiment_catalog_zh.md`](../docs/experiment_catalog_zh.md)。

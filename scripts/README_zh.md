# 仓库脚本

[English](README.md) | 中文说明

本目录包含验证、基准、性能剖析、打包和历史诊断脚本。请从仓库根目录运行，确保程序包导入和产物路径按说明解析。

## 脚本分类

| 前缀或文件 | 用途 |
|---|---|
| `validate_*.py` | 重放 fast/reference 对比、参数边缘、参数扫描和档位检查 |
| `build_*.py` | 根据已有产物汇总验证矩阵和 manifest |
| `bench_*.py`、`benchmark_*.py` | 测量运行时间、数值精度或指定诊断假设 |
| `profile_*.py` | 将全流程耗时分解到求解器各阶段 |
| `smoke_*.py`、`verify_distribution.py` | 检查安装后 wheel 的资源和打包边界 |
| `plot_*.py` | 从已保存的验证或基准数据绘图 |
| `_resource_budget.py`、`_smoke_probe.py` | 线程资源限制和隔离冒烟检查的共享工具 |

许多 `benchmark_*_spike.py` 是历史独立候选方案测试；它们不会选择或修改生产求解路径。运行前先看脚本的 `--help` 和输出路径。保留带日期的实验产物；新实验应使用新文件名，除非确实要替换原文件。

SageNet 对比脚本还需要单独的 SageNet 代码仓库和运行环境；本项目的 extras 不会安装 SageNet。

## 可复现实验

按研究问题使用对应的环境、线程预算、随机种子、预热方式和重复次数。冷启动/JIT 开销与预热耗时分开记录；数值精度比较与速度比较也分别报告。记录代码提交。详见[`docs/reproducibility.md`](../docs/reproducibility.md)、[`docs/benchmarks.md`](../docs/benchmarks.md)和[`docs/experiment_catalog_zh.md`](../docs/experiment_catalog_zh.md)。

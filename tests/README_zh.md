# 测试目录

[English](README.md) | 中文说明

测试用于保护公共 API、数值不变量、求解器档位、打包资源以及选定的科学对照门槛。请从仓库根目录运行命令。

## 主要测试范围

- `test_fast_sgwb.py`、`test_engine.py`、`test_freq_adaptive.py`：fast 求解结果、档位解析、物理保护、网格和原生评估节点。
- `test_reference.py`、`test_numerics.py`、`test_physics_limits.py`：reference 辅助函数和物理极限。
- `test_cobaya_adapter.py`：适配器选项、派生参数、回退和遥测。
- `test_compatibility_smoke.py`、`test_resource_budget.py`：支持的 Python 版本导入和运行资源默认值。
- `test_*_spike.py`、`test_*_round*.py`：隔离候选方案或历史诊断契约；spike 测试通过不表示候选方案已进入生产路径。
- `test_validation_manifest.py`、`test_posterior_tools.py`：产物结构和分析辅助函数。

## 运行

```bash
python -m pytest -q                 # 默认测试集；慢速测试会被排除
python -m pytest -m slow -q         # 耗时较长的 LSODA/reference 检查
python -m pytest -m cobaya -q       # Cobaya 集成检查
```

测试风险与耗时分层见[`docs/test_coverage_matrix.md`](../docs/test_coverage_matrix.md)；CI 如何避免重复昂贵求解见[`docs/test_duplication_audit.md`](../docs/test_duplication_audit.md)。测试结果绑定实际验证的提交和环境，不能替代有原始产物支撑的精度或参数空间结论。

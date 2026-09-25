# `stiffgwpy_fast` 程序包

[English](README.md) | 中文说明

本目录包含安装到 Python 环境中的程序包。项目概览和用户快速开始见[根目录 README](../README_zh.md)。

## 引擎职责

- `fast` 是唯一正式的用户可见快速档位，使用 goal 频率网格、精确再加热 kink 拆分、Numba 张量传播和 PCHIP 玻尔积分。
- `reference` 是独立的连续 `sigma(N)` + DOP853 精度参照，用于有明确范围的对比和认证运行。
- `lsoda` 保留原始自适应积分路径，用于回归检查和显式故障回退。

`production`、`ultra-fast` 等历史名称仅为验证兼容保留。高层 API 会将弃用的快速档位名称映射为 `fast`；它们不是额外的生产档位。

## 模块索引

| 模块 | 职责 |
|---|---|
| `stiff_SGWB.py` | 公共入口 `LCDM_SG.SGWB_iter` 和引擎分派 |
| `LCDM_stiff_Neff.py` | 宇宙学模型参数与派生量 |
| `fast_sgwb.py` | fast 配置、Numba 内核、迭代、频率积分和局部误差预算 |
| `exact_background.py` | 连续背景求值及 fast 准备阶段使用的 kink 感知积分原语 |
| `freq_adaptive.py` | 自适应和面向目标的频率网格 |
| `reference.py` | 独立连续背景 DOP853 参照流程 |
| `config.py` | 不可变且经过验证的逐次调用配置 |
| `_resources.py` | 以 context manager 访问包内数据文件 |
| `_metrics.py` | 共享谱误差和相对误差计算 |
| `global_param.py`、`functions.py` | 物理常数和旧版数值辅助函数 |
| `cobaya/` | Cobaya theory 适配器和似然模块 |

## 配置边界

`LCDM_SG.SGWB_iter()` 默认选择正式 `fast` 档位。具名配置按次解析，`FastSolverConfig` 不可变。直接构造它会使用旧手动设置默认值，并非正式预设；需要具名档位时使用 `fast_sgwb.resolve_config('fast')`。只有明确要使用旧模块设置时才传入 `accuracy_mode=None`。显式参数覆盖只作用于当前调用。

```python
from stiffgwpy_fast import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
model.SGWB_iter()  # 使用正式 fast 档位
```

使用数值结果前，请先阅读[数值方法](../docs/numerical_method.md)、[精度边界](../docs/accuracy.md)和[复现说明](../docs/reproducibility.md)，并核对每项带日期结果的代码提交与实验范围。

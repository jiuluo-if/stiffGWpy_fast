# stiffGWpy_fast

[English](README.md) | 中文说明

`stiffGWpy` 是一个计算 LCDM、刚性物质和原初随机引力波背景（SGWB）的宇宙学程序包。
项目提供快速求解器、独立连续 sigma 高精度参考流程，以及仅用于回归和运行时间对比的
LSODA 流程。

## 功能概览

对于给定的宇宙学模型 `LCDM_SG`，程序会：

1. 构造包含辐射、中微子、刚性物质和 Lambda 的背景；
2. 按频率通道求解张量模扰动方程；
3. 迭代 `Delta N_eff`，使 SGWB 对额外辐射的贡献达到自洽。

主要输出包括今天的 `Omega_GW(f)`、积分量 `DN_gw`、辐射物种数 `kappa_r`、膨胀曲线，
以及 Cobaya 理论模块使用的派生量。

## 快速求解档位

项目对外提供两个快速档位：

| 档位 | 用途 | 特点 |
|---|---|---|
| `fast` | 快速探索、信号形状筛查 | plain-grid，速度优先，精度包络仍需按验证结果解释 |
| `production` | 正式计算和 MCMC 热路径 | transition-refine，处理 kink、相位和自适应频率网格 |

独立的 `stiffgwpy.reference` 连续 sigma 流程是精度锚点。LSODA 只用于回归、故障回退和
运行时间比较，不是第三个生产档位。

## 安装

在项目目录执行：

```bash
pip install .
```

如需 Cobaya 集成：

```bash
pip install .[cobaya]
```

## 基本用法

```python
from stiffgwpy import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
print(model.derived_param)
```

高层接口可通过 `engine="fast"` 使用快速求解；需要精度锚点时使用
`engine="reference"`。具体参数、异常处理和输出契约请查看英文 README 以及 `docs/` 中的
专题文档。

## 验证

```bash
python -m pytest -q
python -m pytest -q -m slow
python -m pytest -m cobaya -q
python scripts/validate_manifest.py
python -m build --wheel
python scripts/smoke_installed_wheel.py dist/stiffgwpy-*.whl
```

当前验证结果和精度边界以 [`docs/validation/validation_manifest.json`](docs/validation/validation_manifest.json)
为准；可复现命令见 [`docs/reproducibility.md`](docs/reproducibility.md)。

## 重要限制

- 快速求解器的速度优势不等于普适精度认证；请根据验证产物和局部误差预算解释结果。
- `reference` 使用冻结尾部近似，深尾部认证的计算成本很高。
- 完整的独立 reference-engine MCMC 尚未运行，现有后验结论基于 importance reweighting。

中英文 README 的入口和关键信息应保持同步；新增功能、参数或验证结论时，请同时更新本页
和 [`README.md`](README.md)，详细主题文档仍以 `docs/` 为准。

## 目录

- `stiffgwpy/`：Python 程序包和 Cobaya 适配器；
- `tests/`：单元测试、慢速数值门禁和 Cobaya 测试；
- `scripts/`：验证、构建 manifest 和 wheel 冒烟脚本；
- `docs/`：验证产物、复现说明和专题文档。

## 许可证

GPL-3.0，详见 [`LICENSE.md`](LICENSE.md)。

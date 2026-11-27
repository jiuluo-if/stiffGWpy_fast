# stiffgwpy_fast

`stiffgwpy_fast` 是用于计算含刚性物质宇宙学中随机引力波背景的 Python 软件包。默认使用 `fast` 引擎；另提供独立的 `reference` 和 `lsoda` 引擎，用于对照和回退。

## 安装

```bash
python -m pip install .
```

使用 Cobaya 适配器时安装 `.[cobaya]`；开发和运行测试时安装 `.[dev]`。

## 快速开始

```python
from stiffgwpy_fast import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
model.SGWB_iter()  # 默认 fast 引擎
```

模型对象提供计算得到的频谱和派生宇宙学量。可选的 Cobaya 集成参见 [Cobaya 适配说明](docs/cobaya.md)。

## 数值方法与验证

正式 `fast` 档位使用 Numba 传播核、考虑再加热拐点的背景预处理和自适应频率积分。`reference` 引擎使用独立的连续背景 DOP853 求解；`lsoda` 用于回归对照和回退。物理模型、[数值方法](docs/numerical_method.md)和[精度边界](docs/accuracy.md)分别说明模型假设、求解流程及验证范围。

核心回归测试使用 `tests/data/` 中固定的参考样本。这些样本仅作为测试夹具；批量扫描、实验脚本、报告和生成图表保存在本地，不纳入正式发布仓库。

## 开发

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
ruff check stiffgwpy_fast scripts tests
python -m mypy
python -m build --sdist --wheel
```

打包要求见 [PACKAGING.md](PACKAGING.md)，贡献要求见 [CONTRIBUTING.md](CONTRIBUTING.md)，复现规范见 [docs/reproducibility.md](docs/reproducibility.md)。

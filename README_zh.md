# stiffGWpy_fast

[English](README.md) | 中文说明

`stiffGWpy_fast` 是一个计算 LCDM、刚性物质和原初随机引力波背景（SGWB）的宇宙学程序包。
项目提供快速求解器、独立连续 sigma 高精度参考流程，以及仅用于回归和运行时间对比的
LSODA 流程。

PyPI 安装名为 `stiffgwpy_fast`，Python 导入名为 `stiffgwpy_fast`，两者都与旧的
`stiffgwpy` 项目区分开来。

## 功能概览

对于给定的宇宙学模型 `LCDM_SG`，程序会：

1. 构造包含辐射、中微子、刚性物质和 Lambda 的背景；
2. 按频率通道求解张量模扰动方程；
3. 迭代 `Delta N_eff`，使 SGWB 对额外辐射的贡献达到自洽。

主要输出包括今天的 `Omega_GW(f)`、积分量 `DN_gw`、辐射物种数 `kappa_r`、膨胀曲线，
以及 Cobaya 理论模块使用的派生量。

## 快速求解档位

项目对外只提供一个正式快速档位 `fast`，组合了：

- `h=0.005`、`col_step=8`、`z_tail=5`、`phase_max=0.25`；
- 精确 `N_re` 断点拆分，transfer step 不跨越 reheating kink；
- 通常约 70–120 点的 goal-oriented 频率网格，并保留 `eval_freqs` 原生节点。

`production`/`transition_refine` 仅作为旧调用的弃用兼容别名，高层接口会将它们映射到
`fast`；它们不再是第二个用户档位，底层验证脚本仍可保留内部预设。
独立的 `stiffgwpy_fast.reference` 连续 sigma 流程是精度锚点。LSODA 只用于回归、故障回退和
运行时间比较。

## 安装

在项目目录执行：

```bash
pip install stiffgwpy-fast
```

如需 Cobaya 集成：

```bash
pip install .[cobaya]
```

打包和发布到 PyPI 的完整说明见 [`PACKAGING.md`](PACKAGING.md)。发布包会主动排除
`docs/`、测试文件、验证脚本、CI 配置和仅用于研究的配置文件。

## 本地隐私与仓库卫生

仓库会忽略本地凭据和自动生成文件，包括 `.pypirc`、`.env*`、密钥/证书、缓存、覆盖率
报告、构建目录、MCMC 链、检查点和临时参数扫描结果。PyPI Token 只能保存在本地，不能
写入 README、Issue、命令日志或 Git 提交。提交前请运行 `git status --short`；对可能含有
凭据或私有实验结果的文件，可用 `git check-ignore -v <文件路径>` 确认其已被忽略。

Git 中的 `docs/` 只保留选定的验证说明和较小的复现摘要；可再生成的大型或未完成研究结果
不提交，也不会进入 PyPI 发布包。

## 基本用法

```python
from stiffgwpy_fast import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
model.SGWB_iter()
print(model.derived_param)
```

高层接口 `SGWB_iter()` 默认使用唯一的 `fast` goal-kink-hybrid 档位。旧的
`accuracy_mode="production"` 仅用于兼容性验证；需要精度锚点时使用
`engine="reference"`；需要原始回归路径时使用 `engine="lsoda"`。

## 验证

```bash
python -m pytest -q
python -m pytest -q -m slow
python -m pytest -m cobaya -q
python scripts/validate_manifest.py
python -m build --wheel
python scripts/smoke_installed_wheel.py dist/stiffgwpy_fast-*.whl
```

当前验证结果和精度边界以 [`docs/validation/validation_manifest.json`](docs/validation/validation_manifest.json)
为准；可复现命令见 [`docs/reproducibility.md`](docs/reproducibility.md)。

当前 fresh HEAD 使用正式 fast 路径（`h=.005`、`col_step=8`、`z_tail=5`、
`phase_max=.25`、exact kink split、goal grid）。固定 20 threads/workqueue 的六点
25 次矩阵中，default warm median/p95 为 `5.36/5.73 ms/point`，尚未达到 `<4 ms`。
同一 76-node 网格的独立 reference 对照中，Simpson/PCHIP 的 `DN_gw` 相对误差为
`7.14e-4/2.94e-4`，所以 PCHIP 仍只作为 opt-in。`eval_freqs` 已与积分 support
nodes 解耦，六点 DN invariant 实测变化为 0。详细证据见 `findings.md`、
`progress.md` 和对应 JSON artifact。

## 重要限制

- 快速求解器的速度优势不等于普适精度认证；请根据验证产物和局部误差预算解释结果。
- `reference` 使用冻结尾部近似，深尾部认证的计算成本很高。
- 完整的独立 reference-engine MCMC 尚未运行，现有后验结论基于 importance reweighting。

中英文 README 的入口和关键信息应保持同步；新增功能、参数或验证结论时，请同时更新本页
和 [`README.md`](README.md)，详细主题文档仍以 `docs/` 为准。

## 目录

- `stiffgwpy_fast/`：Python 程序包和 Cobaya 适配器；
- `tests/`：单元测试、慢速数值门禁和 Cobaya 测试；
- `scripts/`：验证、构建 manifest 和 wheel 冒烟脚本；
- `docs/`：验证产物、复现说明和专题文档。

## 许可证

GPL-3.0，详见 [`LICENSE.md`](LICENSE.md)。

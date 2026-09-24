# stiffGWpy_fast

[English](README.md) | 中文说明

`stiffGWpy_fast` 是计算 LCDM、刚性物质与原初随机引力波背景（SGWB）的 Python 程序包。它提供一个用户可见的 `fast` 快速求解档位、独立的连续 `sigma(N)` 高精度参考流程，以及仅用于回归、显式故障回退和耗时比较的原始 LSODA 路径。

PyPI 安装名为 `stiffgwpy_fast`，Python 导入名也为 `stiffgwpy_fast`；两者都与旧的 `stiffgwpy` 项目区分开来。

## 功能概览

对于模型 `LCDM_SG`，程序构造含辐射、中微子、刚性物质和 Lambda 的背景，逐频率求解张量模扰动，并迭代 `Delta N_eff` 直到 SGWB 对额外辐射的贡献达到自洽。主要输出包括今天的 `Omega_GW(f)`、积分量 `DN_gw`、`kappa_r`、膨胀曲线和 Cobaya 理论模块使用的派生量。

## 当前快速档位

对外只有一个正式快速档位 `fast`（goal-kink-hybrid）：`h=0.005`、`col_step=8`、`z_tail=5`、`phase_max=0.25`、goal 频率网格、精确 reheating kink 拆分和默认 PCHIP 频率积分。旧的 `production`、`transition_refine` 等名称只保留兼容或验证用途，不代表第二个用户档位。具体配置以 `stiffgwpy_fast.fast_sgwb.ACCURACY_MODES['fast']` 为准。

## 安装与运行

```bash
pip install stiffgwpy-fast
```

在项目目录安装源码和依赖，或添加可选 Cobaya 支持：

```bash
pip install .
pip install ".[cobaya]"
```

```python
from stiffgwpy_fast import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
model.SGWB_iter()  # 默认使用唯一的 fast goal-kink-hybrid 档位
print(model.DN_gw[-1])
```

需要独立精度参照时使用 `model.SGWB_iter(engine="reference")`；要运行原始 LSODA 路径时使用 `model.SGWB_iter(engine="lsoda")`。打包和 PyPI 发布说明见 [`PACKAGING.md`](PACKAGING.md)。发布包不包含 `docs/`、测试、CI 和研究专用配置。

## 最新性能与精度证据

最新全流程性能数据来自 2026-09-23 的 solver commit `c7d4766`，固定 2 个 Numba `workqueue` 线程、1 个 BLAS 线程、CPU affinity `[0, 1]`，重复 25 次。预热后 median/p95：default `7.189/7.916 ms`、low-T `7.116/7.782 ms`、high-T `10.475/10.890 ms`、stiff `10.260/10.936 ms`、high-kappa `9.929/11.062 ms`。每个情景的首个样本单独记录；只有 default 首样本是进程冷启动并包含 JIT，后续情景是在同一进程内核已初始化后测得。当前仓库 HEAD 后续补充了 MCMC 证据，求解器源码与该 profile commit 相同。完整协议见 [`docs/benchmarks.md`](docs/benchmarks.md)。

最新有明确精度归属的审计记录日期为 2026-09-12：在四个指定点，PCHIP 对 Oracle C WKB anchor 的误差为 `1.07e-5`–`1.66e-4`；六点 `z_tail=8` same-grid reference 对比的 `DN_gw` 差值 median 为 `2.93e-4`。这些都是有限点位和指定协议下的记录，不能视为全参数空间认证。更多边界见 [`docs/accuracy.md`](docs/accuracy.md) 与 [`docs/experiment_catalog_zh.md`](docs/experiment_catalog_zh.md)。

2026-09-24 的 LVK MCMC 对比报告使用三种方法、三个情景和 4 条链；当前 fast 的采样步骤约比 SageNet+ 快 `2.9–4.9` 倍。SageNet+ 有两种情景未达到报告采用的链诊断参考值；低重加热温度情景的观测频段覆盖为 0%，不能据此比较数据拟合表现。参见 [`docs/mcmc_sagenet_compare/report.md`](docs/mcmc_sagenet_compare/report.md)。

参数验证产物仍是带日期的抽样：400 点 LHS 屏查为 254 个成功、146 个物理保护、0 个数值失败；另一个 240 点 Sobol 运行记录 212 个成功和 28 个物理保护。它们不等于当前完整参数空间认证。原始 manifest 的 commit 是 `083fdf5`，查看时应同时核对文件中记录的提交和实验范围。

## 验证与仓库说明

```bash
python -m pytest -q
python -m pytest -q -m slow
python -m pytest -m cobaya -q
python scripts/validate_manifest.py
python -m build --wheel
python scripts/smoke_installed_wheel.py dist/stiffgwpy_fast-*.whl
```

测试数量会随提交和可选标记变化；新鲜结果以 CI 和 [`docs/reproducibility.md`](docs/reproducibility.md) 的验证命令为准。README 中的固定读数均注明对应日期、提交或实验范围，不能直接外推为其他机器或整个参数空间的保证。

仓库 `docs/` 保存主题说明、带日期的研究报告和验证数据；本机链文件及检查点按 `.gitignore` 管理。提交前检查 `git status --short`，不要提交凭据、Token 或未审核的本地实验数据。`.pypirc`、`.env*`、密钥/证书、缓存、覆盖率和构建输出等本机文件应保持忽略；PyPI Token 只保存在本机。README 的中英文入口和关键结论应同步维护。

## 目录

- `stiffgwpy_fast/`：Python 包、求解器与 Cobaya 适配器；
- `tests/`：默认回归、慢速数值、Cobaya 和兼容性测试；
- `scripts/`：验证、基准、manifest 与打包冒烟脚本；
- `docs/`：当前指南、历史实验报告、验证产物和原始实验数据。

## 许可证

GPL-3.0，详见 [`LICENSE.md`](LICENSE.md)。

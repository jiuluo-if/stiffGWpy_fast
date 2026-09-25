# stiffGWpy_fast

[English](README.md) | 中文说明

`stiffGWpy_fast` 是计算 LCDM、刚性物质与原初随机引力波背景（SGWB）的 Python 程序包。它提供唯一正式的用户可见 `fast` 档位、独立连续 `sigma(N)` 高精度参照流程，以及仅用于回归、显式回退和耗时比较的原始 LSODA 路径。PyPI 发行包名和 Python 导入名均为 `stiffgwpy_fast`，与旧的 `stiffgwpy` 项目区分。

> **文档说明：** 带日期的精度和验证结论须结合原始数据与代码提交阅读。manifest 可能汇总较早且范围有限的运行，不能视为当前 HEAD 的认证结果。`fast` 配置以 `stiffgwpy_fast.fast_sgwb.ACCURACY_MODES['fast']` 为准；最新性能数据见 [`docs/benchmarks.md`](docs/benchmarks.md)，各实验的证据状态见 [`docs/experiment_catalog_zh.md`](docs/experiment_catalog_zh.md)。

## 计算内容和物理参数

给定 `LCDM_SG`，程序构造包含辐射、中微子、刚性物质和 Lambda 的背景，逐频率求解张量模扰动，并迭代 `Delta N_eff`，使 SGWB 对额外辐射的贡献达到自洽。主要输出包括今天的 `Omega_GW(f)`、积分量 `DN_gw`、辐射组分数 `kappa_r`、膨胀曲线 `g2`/`w2` 和 Cobaya 所需派生量。

| 参数 | 含义 |
|---|---|
| `Omega_bh2`、`Omega_ch2`、`H0` | 重子/冷暗物质密度和哈勃常数 |
| `DN_eff` | 常量额外辐射（`Delta N_eff`） |
| `A_s`、`r`、`n_t` | 标量振幅、张量标量比、张量倾斜 |
| `cr` | 大于 0 时强制单场一致性关系；否则使用 `n_t`/`DN_re` |
| `T_re`、`DN_re` | 再加热温度（GeV）和类物质再加热 e-fold 数 |
| `kappa10` | 10 MeV 时刚性物质与光子的能量密度比 |

瞬时再加热转换在 `sigma(N)` 中表现为 kink。正式 `fast` 档位包含精确 `N_re` 传播拆分、相位限制子步进和面向目标区域的稀疏频率网格。旧的 `production`、`transition_refine` 等名称只为兼容验证保留，不是其他用户档位。精度锚点是独立的 `stiffgwpy_fast.reference` 流程。

## 当前 fast 档位

| 配置项 | `fast` |
|---|---|
| `h`、`col_step` | `0.005`、`8` |
| `z_tail`、`phase_max` | `5.0`、`0.25` |
| `freq_grid`、`kink_split` | goal（通常 70–120 点）、开启 |
| `frequency_quadrature` | 默认 `pchip`；可显式使用 `simpson` |
| 外层容差 | `1e-6` |

goal 网格在再加热特征附近保留节点，并将 `eval_freqs` 作为原生求解节点；`z_tail` 之后由解析 WKB 包络处理。局部误差预算见 `stiffgwpy_fast.fast_sgwb.estimate_local_error`。

最新全流程 profile 于 2026-09-23 在 solver commit `c7d4766` 生成：25 次重复、2 个 Numba `workqueue` 线程、1 个 BLAS 线程、CPU affinity `[0, 1]`。预热 median/p95（毫秒）为 default `7.189/7.916`、low-T `7.116/7.782`、high-T `10.475/10.890`、stiff `10.260/10.936`、low-r `5.923/6.313`、high-kappa `9.929/11.062`；六个情景均收敛。只有 default 首样本是进程冷启动并包含 JIT。当前 HEAD 后续增加了 MCMC 证据，没有改动该 profile 提交之后的求解器源码。完整协议见 [`docs/benchmarks.md`](docs/benchmarks.md)。

## 精度参照与已知边界

`reference.py` 独立实现连续 `sigma(N)`（kink 为精确断点）、逐频率自适应 `DOP853`，以及带误差估计的 PCHIP 和自适应 Gauss–Kronrod 求积。LSODA 不作为精度真值。参照流程在 `z_tail` 使用冻结解析尾部：默认点从 `z_tail=7` 到 8 的 `DN_gw` 相对变化为 `4.2e-4`，从 8 到 10 为 `3.0e-4`；`z_tail=14` 的无截断求解因深次视界刚性而不可行。因此参照本身有约 `3e-4` 的尾部敏感性，正式路径约 `4e-4` 的残差与之同量级。`reference.oracle_variants()` 运行 Oracle A/B/C 并报告 `CONSISTENT` 或 `ORACLE-SENSITIVE`。

2026-09-12 的指定点审计中，四点 PCHIP 对 Oracle C WKB anchor 的误差为 `1.07e-5`–`1.66e-4`；六点 `z_tail=8` 同网格对照的 `DN_gw` 差值 median 为 `2.93e-4`。原生 76 点默认网格上，fast 对连续 `sigma` 参照的相对误差为 `2.94e-4`，Simpson 为 `7.14e-4`。76/80/90/110 点扫描并非单调；这些读数均不是全参数空间认证。完整范围见 [`docs/accuracy.md`](docs/accuracy.md)。

## 安装和 Python 用法

```bash
pip install stiffgwpy-fast
pip install .
pip install '.[cobaya]'
pip install '.[cobaya,mpi]'
pip install '.[dev]'
```

```python
from stiffgwpy_fast import LCDM_SG

model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
model.SGWB_iter()  # 默认使用唯一的 fast goal-kink-hybrid 档位
print(model.DN_gw[-1])
model.SGWB_iter(engine="reference")  # 独立精度参照
```

使用 `engine="lsoda"` 可运行原始回归路径。旧的 `production`、`transition_refine` 和 `ultra-fast` 名称在高层 API 中映射到正式 `fast` 档位。显式传入 `h`、`col_step`、`z_tail`、`freq_res` 或 `tol` 可覆盖对应值；仅当有意读取旧模块手动设置快照时，才显式传入 `accuracy_mode=None`。底层调用也可按次传入不可变的 `fast_sgwb.FastSolverConfig`。

打包和 PyPI 发布说明见 [`PACKAGING.md`](PACKAGING.md)。发布包不包含 `docs/`、测试、验证脚本、CI 文件和研究专用配置。本机链文件、凭据、缓存、覆盖率和构建输出由 `.gitignore` 管理；PyPI Token 只保存在本机。提交前检查 `git status --short`，不要提交凭据或未审核的本地数据。

## Cobaya 用法

```yaml
theory:
  stiffgwpy_fast.cobaya.stiffGW.stiffGW:
    engine: fast
    fallback: True
    accuracy_mode: fast
    fast_threads: 8
```

`stiffGW.yaml` 中 fast 配置的 `0` 是“采用所选 accuracy mode”的哨兵值，解析顺序为 `accuracy_mode -> 预设默认值 -> 用户显式覆盖`。因此默认 `engine: fast` 会使用组合后的 fast 设置。

在理论 YAML 中指定 `eval_freqs: [数值列表或文件路径]` 可把 `log10(f/Hz)` 频率作为原生求解节点，从而避免陡峭谱特征处的插值误差。Layer C 在 11 个 PTA 类频点上测得逐 bin dex 插值误差不超过 `3.1e-4`。默认 `eval_freqs: null` 时，似然对返回频谱插值。原生评估节点不进入玻尔积分支撑网格，因此不会改变自洽的 `DN_gw`。

## 参数空间与后验验证

已有参数扫描是带日期的样本，不代表当前全参数空间认证：

- 历史 plain-grid LHS 400 点：254 成功、146 个共享 `Delta_Neff` 物理保护、0 个数值失败。36% 的保护比例来自总 `N_eff > 5` 的物理拒绝；记录生成于 `f87e969`。
- 历史 fast Sobol 240 点：212 成功、28 个物理保护，见 `docs/paramsweep_ref/fast_sweep.jsonl`。
- 16 个参数轴边缘及转换区域点见 `docs/paramsweep_z8b/`。

拒绝状态分为 `PHYSICAL_INVALID`、`PHYSICAL_GUARD`、`NUMERICAL_FAILURE`、`FAST_ERROR` 和 `ORACLE_ERROR`；物理保护不是数值失败。

mock-data 重要性采样验证（Layer C）使用固定种子的 9,000 个 fast 样本：ESS **4167**（门槛 2000，通过），`log10 r` 后验移动 **-0.0011 sigma**（门槛小于 0.1 sigma，通过），逐 bin dex 最大值 **3.1e-4**，`|Delta logL|` 最大值 **7.3e-3**（门槛 0.1，通过）。这是重要性重加权结果，不是独立 reference MCMC 链。详见 [`docs/mcmc_posterior/posterior_validation.md`](docs/mcmc_posterior/posterior_validation.md)。

另一个 2026-09-24 LVK 报告在三个情景下比较 plain-grid、fast 和 SageNet+，每种方法每个情景运行四条链。fast 采样单步约快 `2.9–4.9` 倍；SageNet+ 在两个情景中未达到报告采用的链诊断参考值。低再加热温度情景对所有方法的 LVK 频段覆盖均为 0%，不能据此比较数据拟合。详见 [`docs/mcmc_sagenet_compare/report.md`](docs/mcmc_sagenet_compare/report.md)。

## 性能边界与复现

2026-09-23 profile 没有同资源 LSODA 对照，因此不报告当前速度比。default warm runtime 为 `7.189 ms`，高于 `4 ms` 目标。reference 在有日期的深尾研究中约需每点 `360–579 s`，具体取决于尾部设置；深尾无截断认证因方程刚性而不可行。

贡献和发布门槛见 [`CONTRIBUTING.md`](CONTRIBUTING.md)，分阶段工程审计见 [`docs/engineering_audit.md`](docs/engineering_audit.md)。常用验证命令如下：

```bash
python scripts/validate_two_modes.py --phase convergence
python scripts/validate_two_modes.py --phase param_sweep --n 400
python scripts/build_two_mode_manifest.py
python -m pytest
python -m pytest -m cobaya
python scripts/validate_manifest.py
python -m build --wheel
python scripts/smoke_installed_wheel.py dist/stiffgwpy_fast-*.whl
```

驱动脚本会记录 Git 提交和环境信息。测试数量依赖提交和可选标记；核验当前状态时，以 CI 和 [`docs/reproducibility.md`](docs/reproducibility.md) 为准。历史测试数量只代表对应提交。

## 目录与许可证

`stiffgwpy_fast/` 提供求解器、引擎分派和 Cobaya 适配器；`tests/` 包含默认、慢速、Cobaya 和兼容性测试；`scripts/` 提供验证、基准、manifest 和安装冒烟脚本；`docs/` 保存当前指南、带日期的报告、验证数据与实验记录。

各目录的开发说明见[程序包模块索引](stiffgwpy_fast/README_zh.md)、[脚本目录](scripts/README_zh.md)、[测试目录](tests/README_zh.md)和[Cobaya 适配器](stiffgwpy_fast/cobaya/README_zh.md)。

GPL-3.0，详见 [`LICENSE.md`](LICENSE.md)。

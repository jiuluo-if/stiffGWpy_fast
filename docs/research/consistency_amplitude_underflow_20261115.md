# `cr>0` 张量振幅下溢与一致性背景分类

日期：2026-11-15

## 根因

一致性分支 `cr>0` 用张量振幅决定暴胀尺度和暴胀终点。输入 `r>0`、`A_s>0` 时数学乘积 `A_t=r A_s` 非零；但当两者乘积低于 float64 最小次正规数时，派生参数中 `A_t` 舍入为 0。原背景代码仍执行 `math.log(45/2*A_t/(rho_rad_re+rho_stiff_re))`，在 `derived_param` 阶段抛出 `ValueError: math domain error`。这发生在 Fast 模态传播、`S2` 原函数构造和终端组装之前，不能归类成 `S2` 溢出。

`cr>0` 下不能直接采用 `cr=0` 的 Fast 单位振幅回退：这里 `A_t` 同时决定背景 cutoff。正确处理是仅在存储乘积为零、且两个输入仍为正时，以

```text
log_A_t = log(r) + log(A_s)
```

计算背景公式中的正值表达式。原有公式

```text
V_inf = (3 A_t / 2)^(1/4) * sqrt(pi) * M_Pl
Delta_N = 4/3 (N_10 - N_re) + 4/3 log(M_Pl)
          + 1/3 log[(45/2) A_t / (rho_rad_re + rho_stiff_re)]
```

在对数域等价为

```text
log(V_inf) = 1/4 [log(3/2) + log_A_t] + 1/2 log(pi) + log(M_Pl)
Delta_N = 4/3 (N_10 - N_re) + 4/3 log(M_Pl)
          + 1/3 [log(45/2) + log_A_t - log(rho_rad_re + rho_stiff_re)]
```

这只修复了派生背景公式的浮点表示；若由此得到 `Delta_N < 0`，`N_inf` 仍不存在，Fast 必须拒绝 cutoff。不会把依赖 `A_t` 的一致性背景替换成自由倾斜背景，也不会在此分支传播单位功率。

## 四分类试验

分类脚本在相同 cosmology、`T_re=2000`、`kappa10=.01` 下检查：

- `cr=0, r=1e-316, A_s>0`：存储 `A_t=0`，但 `log(r)+log(A_s)` 有限；Fast 继续现有 log-power 路径并完成。
- `cr=0, r=0`：由输入定义的精确零振幅；Fast 使用精确零路径并完成。
- `cr>0, r=1e-316, A_s>0`：`log_A_t=-747.598`；对数背景计算得到有限 `V_inf=3.23e-63`，但 `Delta_N=-209.962`，故 `N_inf=None`，Fast 按 `invalid_cutoff` 拒绝。
- `cr>0, r=0`：一致性分支不支持零 tensor-to-scalar ratio，Fast 按 `invalid_r` 拒绝。

这些拒绝分别位于背景 cutoff 和输入一致性边界。此处没有进入 `S2` 原函数检查；`S2=exp(Psi)` 的指数溢出仍是独立失败类别。该改动没有把 `cr>0` 下溢样本纳入适用域。

## 回归证据

新增测试先按 TDD 在旧代码上因 `ValueError: math domain error` 失败；实现对数域背景表达式后，三项 `cr>0` 回归通过，覆盖非零乘积下溢、精确 `r=0` 和精确 `A_s=0`。在同一修订上，完整测试为 `406 passed, 9 deselected, 2 warnings`；`ruff check .`、项目配置的 Mypy、sdist 和 wheel 构建均通过。

另外重跑了一点 `cr=1,r=.01` 的常规振幅尾深验证；`A_t` 在 Fast 和 Reference 中完全相等，z5/z8 精度结果与改动前一致。非零 `cr>0` 下溢案例的原始输入、推导量、最终 `N_inf` 和 Fast failure reason 保存在原始 JSON。

## 判定和边界

接受这项背景计算稳定性修复：对数域求值与原数学表达式等价，且现有物理 cutoff 约束仍生效。它修复的是由下溢触发的未处理数学异常，不是对 `cr>0` 下溢振幅进行 Fast 频谱恢复。只有当对数推导出的 `Delta_N` 非负、背景其余约束有效，而且一致性分支的功率传播也具有独立的稳定表达式及误差证据时，才能讨论适用域扩展。

可复查产物：[分类脚本](../../scripts/screen_consistency_amplitude_underflow_20261115.py)、[原始结果](results/consistency_amplitude_underflow_20261115.json)、[回归测试](../../tests/test_fast_log_amplitude.py)。

## `T_re`/`kappa10` 边界和 `S2` 对照

又对 `T_re={10,2000,10000}` 与 `kappa10={1e-3,1e-2,1}` 的 3×3 离散组合重复 `r=1e-316, A_s=2.100549e-9, cr=1`。九点的 log 派生 `Delta_N` 均为负，范围 `−198.43` 至 `−214.72`，都在 Fast 背景原函数构造前以 `invalid_cutoff` 返回。该扫描说明在这些离散温度和刚性密度点，log 域背景公式可算且边界拒绝一致；不能据此声称连续范围认证，也没有证明其他参数可恢复 `N_inf`。

为确认 `S2` 是不同阶段，另在当前代码重跑仓库已有的普通振幅压力点 `Omega_bh2=1e6, cr=1, r=.01, n_t=.5, T_re=2000, kappa10=.01`。它的 `N_inf=68.81384` 有效，Fast 实际进入背景原函数并返回 `nonfinite_background_primitive`，detail=`S2`、index=`13763`、记录的对数尺度=`12678.0768`。对照下溢点的 `invalid_cutoff` 且没有 primitive failure detail，证实二者是独立分类：前者是 cutoff 物理前置拒绝，后者是 `S2=exp(Psi)` 的有限值表示溢出。这个控制点只验证分类，不是对该压力参数的精度或适用性认可。

矩阵、控制点输入、逐点输出、环境及当前相关源码哈希见[边界扫描脚本](../../scripts/scan_consistency_underflow_boundary_20261115.py)和[原始 JSON](results/consistency_underflow_boundary_20261115.json)。

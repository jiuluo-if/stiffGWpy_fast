# 第一阶段审计：物理等价性与实现正确性（2026-08-29）

审计对象：`stiffgwpy.fast_sgwb`（fast）vs `stiffgwpy.stiff_SGWB` + `functions.py`（参考 LSODA）。
方法：逐项对照物理方程、初始条件、视界穿越条件、尾部解析近似、频谱积分；并推导 fast 步进算子的数学结构。

## A. 物理方程等价性（连续层面逐一验证）

参考系统（`functions.py:tensor`，N = ln a，prime = d/dN）：

```
z' = 1.5*sigma - 1
x' = -3*x + 1.5*sigma*x - exp(z)*y
y' = -y + 1.5*sigma*y + exp(z)*x
```

fast 使用重标定变量（`assemble_main` 逆推）：x = sqrt(S2)*xh, y = sqrt(S2)*yh，
其中 S2 = exp(Psi), Psi = 3*F - 4*N, F' = sigma, Phi = 1.5*F - N + N0, Phi' = 1.5*sigma - 1 = z'。

代入可得 xh/yh 的精确 ODE（sigma 项被解析吸收，完全消失）：

```
xh' = -xh - exp(z)*yh
yh' = +yh + exp(z)*xh
```

步进算子 `scaled_step` 的矩阵为 [[c-si, -w*si],[w*si, c+si]]，w = exp(z_mid)，
Om = sqrt(w^2-1)，c = cos(Om*h)，si = sin(Om*h)/Om（w>=1）；
这正是常系数系统 A = [[-1,-w],[w,1]]（A^2 = (1-w^2)I）的矩阵指数 e^(A h)。
w<1 时用 cosh/sinh 级数截断（c = 1+x^2/2, si = h(1+x^2/6), x = Om*h）。

结论：**fast 步进是该重标定系统在“z 冻结于区间中点”下的精确解**（Magnus 型中点法），
离散误差仅来自 z 的中点求积（理论上全局二阶）+ w<1 时级数截断（量级小）。

逐项对照结果（均一致）：

| 项目 | 参考实现 | fast | 结论 |
|---|---|---|---|
| 膨胀历史公式 | `gen_expansion` | `gen_kernel` | 逐项相同；仅 FD 积分 quad vs 三次样条查表（测试差 ~1.7e-11） |
| 热历史样条 | `gp.spl_rho/rhop` | 同一 spline | 完全一致 |
| 初始条件 | (z0, 0, e^{z0}) | xh=0, yh=e^{z0}*S2inv, z=z0 | x=0, y=e^{z0} ✓ |
| 起始视界条件 | `find_index_hc(freq+3)`（k/aH<=1e-3） | prep_kernel 二分找 f_hor>=freq+3 | 语义一致 ✓ |
| 深亚视界阈值 | `subhorizon` 事件 z=5 | 循环条件 zz<5.0 | 一致 ✓ |
| 解析尾部 | `Th_hf = coeff*exp(-zf+N_now-N_hf)` | `assemble_tail` 同式 | 一致 ✓（离散停止点略差，见 Major 1） |
| 组装公式 | x^2+y^2、x*Th/3、(-5x^2+7y^2)/72 | 重标定后同式 | 一致 ✓ |
| Simpson 积分 | scipy simpson 逐列 | 预计算 Wmat | 测试差 ~1e-12 ✓ |
| 频率网格 | `construct_f`（共用） | 同一函数 | 一致 ✓ |
| DN_gw/kappa_r 输出 | 同式 | 同式 | 一致 ✓ |

## B. 四级问题清单

### Critical（必须修，阻塞认证）

- 无。连续层面未发现物理方程/初始条件/尾部/积分口径的实质性错误；
  现有 12 点 + 随机点交叉验证与组件级测试（1e-11~1e-15）与本文推导自洽。

### Major（影响认证结论，需数值研究量化）

1. **h 与 z_tail 硬编码不可配置**：`solve_kernel` 内 `scaled_step(..., 0.01)` 与
   `zz < 5.0`、`gen_fast` 的 0.01 网格间距均写死；无法做收敛研究，也无法暴露
   离散停止点（last_z<5 与事件 z=5）差异。→ 第二阶段参数化。
2. **收敛阶与误差来源未量化**：中点法理论二阶、PCHIP 插值、COL_STEP 粗列、
   FD 查表、尾部离散停止——各自误差贡献没有分解。→ 第二、五阶段。
3. **参数空间覆盖不足**：仅 12 手选点 + 少量随机点；无边界/失败集/极端参数。
   → 第三、四阶段。
4. **无严格参考**：默认 LSODA（rtol=1e-6）本身有误差，且外层 1e-4 只是停止条件；
   需收紧容差 LSODA / Radau 独立交叉验证。→ 第二、三阶段。
5. **likelihood 级未验证**：Cobaya 已接入 engine 参数，但未做 ΔlogL / posterior 对比。
   → 第六阶段。

### Minor（小缺陷，顺手修）

1. `scaled_step` 的 w<1 分支是 cosh/sinh 级数截断而非解析延拓；仅超视界初期出现，
   量级 ~O((Om*h)^4)，预计 <1e-10，但应在收敛研究中确认。
2. fast 尾部初始振幅用“最后一个 zz<5 的点”而参考用事件 z=5 精确点；差异
   ~O(h)，进入误差预算（已并入 Major 1 的研究）。
3. `m.DN_gw` 列表在 fast 路径不保留迭代历史（局部 DN_gw_list），与参考行为不同，
   不构成物理问题，但影响下游调试。
4. `_build_fd_table` 用户缓存目录 `~/.cache/stiffgwpy` 无并发写保护（多进程首启可能
   竞争写同一 npz）；建议写临时文件再原子 rename。
5. `global_param.py`/`stiff_SGWB.py` 的 star import 与单行多语句风格遗留（ruff 历史项），
   不影响运行。

### Nice-to-have

1. fast 路径未暴露 engine 元数据（h/col_step/z_tail/threads/版本）→ 第七阶段补
   `m.fast_metadata` / `m.engine_used` / `m.approximate`。
2. `set_threads` 运行中修改只对下次 kernel 启动生效，文档未说明。
3. 收敛研究、随机验证、benchmark 的 JSONL 输出可加 schema 版本号。
4. `MAX_ITER=60` 与“>5 截断”路径的打印信息可包含更丰富的诊断。

## C. 与原始审计报告的衔接

本阶段确认：fast 是同一物理方程的精确重标定 + 中点 Magnus 离散 + 粗列 PCHIP +
解析尾部，**不是另一套物理模型**。因此上一轮“同物理方程”的定位成立；当前差距
集中在“离散误差未量化、参数空间未覆盖、无严格参考、likelihood 未验证”，
全部列入后续阶段。

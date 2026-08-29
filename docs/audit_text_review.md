# 粘贴审计文本逐项复核（2026-08-30）

复核对象：用户粘贴的独立审计报告文本（落款 2026-08-29，审计快照 `6b01ca8`）。
复核基准：当前 fast 分支 `444988d`（相对上游 `origin/main = c8a8733` 领先 18 个提交（写本文档时 `HEAD=444988d`））。
结论：**文本中大部分“当前仓库没有/做不到”的指控已被后续提交修复，仅“全参数认证”与
“posterior shift 收敛链”两项仍为待办（因用户要求停止长跑而未做）。**

## 1. 已过时/已被修复的指控（文本纰漏）

| 文本主张 | 复核判定 | 当前证据 |
|---|---|---|
| “fast 仓库比上游只多一个提交，物理/Cobaya 代码未修改” | 过时 | `fast/main=444988d` 已领先 18 提交，新增 Cobaya 适配器、测试、CI、打包 |
| “脚本把 `F:\codex\stiffGWpy` 硬编码进 sys.path[0]”（bench_fast.py:17 / validate_fast.py:3） | 已修复 | 脚本改用仓库相对路径 |
| “Cobaya 仍调用原求解器 / 未接入” | 已修复 | `stiffgwpy/cobaya/stiffGW.py`：engine=lsoda\|fast、fallback、全部旋钮、MPI-safe；`tests/test_mcmc_compare.py` |
| “pytest 0 tests / coverage 0%” | 已修复 | `tests/` 6 个文件共 38 项测试全绿 |
| “无 pyproject.toml / pip install . 失败” | 已修复 | `pyproject.toml` + wheel 构建冒烟通过 |
| “README 只列 NumPy/SciPy” | 已修复 | README 列出 numpy/scipy/astropy/pyyaml/numba |
| “Ruff 133 项，fast_sgwb.py 26 项含重复导入/重复定义” | 过时 | 维护集 ruff 全绿；`fast_sgwb.py` 的 F/I 类问题为 0（仅 E701/E702 风格债，CI 有意排除并注释说明） |
| “没有 JIT 失败回退到 LSODA 的机制” | 已修复 | `m.SGWB_iter(engine='fast', fallback=True)` + Cobaya `fallback` |
| “while True 无最大迭代” | 已修复 | `MAX_ITER = 60` |
| “DN_gw_new 无 isfinite 检查” | 已修复 | 非有限值直接中止 + `DN_eff>5` 护栏 |
| “修改传入 DN_eff 无 try/finally” | 已修复 | `finally` 中恢复 `DN_eff` |
| “set_col_step(0) 触发 ZeroDivisionError” | 已修复 | setter 校验 1..8 |
| “_COL_STEP 全局可变导致并发布局错位” | 已修复 | 调用入口快照 `col_step` |
| “import 时强改宿主线程数” | 已修复 | 默认 = numba 默认；仅显式 `FAST_THREADS` 时设置并校验 |
| “无随包 FD 表” | 已修复 | `stiffgwpy/fd_table.npz` 随包发布 |
| “报告数字无 commit/环境链路” | 已修复 | bench/validate 等脚本记录 env + git commit 元数据 |
| “报告 10 次取最小值 vs 脚本 3 次中位数” | 已声明 | `docs/benchmark_report.md` 顶部修订声明承认口径差异；`bench_fast.py` 记录 min/median/p95 |
| “报告无正式测试/打包/CI” | 已修复 | `.github/workflows/ci.yml`：pytest + ruff + wheel build |

## 2. 仍成立（非纰漏，属诚实声明）

- **全参数精度认证未完成**：`docs/audit_phase3.md` 标注 PARTIALLY CERTIFIED，扫描 918/1030（91%），
  extreme 角点未跑；此为用户要求停止长跑所致，非能力缺失。
- **posterior shift 无收敛链证据**：`docs/audit_mcmc.md` 明确 NOT CERTIFIED，仅 N=20 短链，
  max|ΔlogL|=0.09、failure rate=0；≥2000 样本收敛链待跑。
- **fast 非 drop-in**：不产出 `N_hc/Th/Oj/Ogw/Opgw`，README 已如实说明。
- **1e-4 只是外层停止条件**，不是统一误差界；README/审计文档均已说明。

## 3. 处置

文本中的指控已作为“部分修改建议”逐条落实（见 §1）；剩余两项（§2 前两条）为长跑任务，
按用户“不用继续跑下去”的要求挂起，续跑命令见 `docs/audit_summary.md` 末尾。

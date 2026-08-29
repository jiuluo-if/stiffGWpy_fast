# 第三阶段审计：大规模参数空间验证（2026-08-30）

审计对象：`stiffgwpy.fast_sgwb`（fast）vs `stiffgwpy.stiff_SGWB`（LSODA 参考）。
目标：建立 ≥1000 个确定性参数点的 LSODA-vs-fast 交叉验证，输出 worst cases / 误差分布 / 参数误差图。
工具：`scripts/param_sweep.py`、`scripts/plot_param_sweep.py`（均已通过冒烟验证与 ruff）。

**结论：第三阶段 NOT CERTIFIED。** 采样设计与工具脚本已入库，但全量扫描在本机两次尝试均因
**系统提交内存墙**（commit-charge 上限）失败，仅 1 点完成全指标对比，不能代表参数空间。
恢复路径见文末。项目整体维持 **PARTIALLY CERTIFIED**（阶段一物理等价性、阶段二 default 点收敛已完成）。

## 1. 设计与协议

- 采样（确定性，seed=20260830，共 1030 点）：Sobol 400 + LHS 350 + 边界过采样 LHS 200
  （unit 坐标向 0/1 推 `|u-0.5|^0.6`）+ 手选极端/角点 80（单参数上下界 22、联合极端 22、
  cr=1 的 (r,T_re) 角格 24、cr=0 的 (n_t,DN_re) 角格 12）。
- 参数覆盖：Omega_bh2 [0.018,0.026]、Omega_ch2 [0.09,0.15]、H0 [60,76]、DN_eff [0,2]、
  A_s [1e-9,4e-9]（log）、r [1e-4,1e-1]（log）、n_t [-0.5,0.5]、cr {0,1}、T_re [10,1e6]（log）、
  DN_re [0,30]、kappa10 [1e-4,1]（log）。对数参数 log 采样，线性参数 linear 采样。
- 每点：同一 σ 网格上 LSODA（rtol=1e-8, atol=[1e-12,1e-22,1e-22], outer tol=1e-7）vs fast
  （h=0.01, col_step=4, z_tail=5, freq_res=1.0, outer tol=1e-7）；记录 status、耗时、ΔN_eff_final、
  DN_gw[-1]、κ_r、频谱 dex 误差（max/p50/p95/p99）、信号区线性 Ω 相对误差（max/median）、
  DN_gw 曲线相对误差（max/median）。
- 绘图脚本产出：`summary.json`、`worst_cases.json`（worst 20）、`error_distribution.png`、
  `parameter_error_map/`（11 张，含 Spearman 相关）。

## 2. 全量扫描受阻：环境内存墙（非代码缺陷）

- 现象：`--workers 2` 下前 7 点全部失败（`MaybeEncodingError` / `MemoryError: Unable to allocate ~1 MB`
  / `OSError(22, 系统资源不足)`），随后系统连新建进程都报 `页面文件太小 (os error 1455)`；
  同一批点在单点冒烟中均正常。
- 根因量化：每个 Python 子进程导入本栈（numpy+scipy+numba+`stiffgwpy`，其中 `global_param`
  顶层 import astropy）私有提交 ~1.66 GB；参考求解器 `run_SGWB` 每次求解开 `mp.Pool(4)`，
  2 个扫描 worker → 峰值 11 进程 ≈ 18 GB，超出本机提交上限（RAM + 自动管理页面文件，
  实测 ~29.6–40.7 GB，基线占用 ~24 GB）。

## 3. 目前唯一完整数据点（sobol-0000）

LSODA 9.2 s + fast 0.07 s：ΔN_eff 相对差 3.2e-13；DN_gw[-1] 相对差 1.06e-5；κ_r 相对差 3.2e-13；
频谱 dex_max 2.97e-3（p50 3.0e-6）；线性 Ω 相对误差 max 6.8e-3 / median 6.9e-6；
DN_gw 曲线相对误差 median 1.06e-5 / max 0.239（早期近零区，与阶段二一致）。
与阶段二 default 点量级一致，但不能外推。

## 4. NOT CERTIFIED 清单与恢复路径

- 未完成：全参数空间 ≥1000 点交叉验证、误差分布与参数误差图、worst 20 cases、各参数误差敏感性、
  极端/边界区域失败率（阶段四输入）。
- 恢复：释放内存后 `python scripts/param_sweep.py --out docs/paramsweep --workers 1 --retry-failed --no-warmup`
  （checkpoint 续跑，`--retry-failed` 重跑失败点），完成后 `python scripts/plot_param_sweep.py`。
- 低内存主机建议：将 `stiff_SGWB.run_SGWB` 的硬编码 `mp.Pool(4)` 改为可配置池大小
  （env var，默认 4；纯进程管理改动，不改变任何数值结果），或延迟化 `global_param` 的
  astropy 顶层 import（每子进程私有提交可降至 ~0.8 GB）。

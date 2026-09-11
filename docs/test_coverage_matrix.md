# Test coverage matrix

本矩阵按风险和成本组织测试，避免在多个 Python 版本重复昂贵 reference 求解。

| Layer | 主要职责 | 运行范围 | 代表入口 |
| --- | --- | --- | --- |
| A | 配置、导入、quadrature helper、资源与 API 兼容性 | Python 3.9–3.13 | `tests/test_compatibility_smoke.py` |
| B | fast 数值、kink、tail、outer reuse、eval_freqs、guard 与 Cobaya adapter 的 canonical 回归 | Python 3.11 一次 | `python -m pytest -q`、`python -m pytest -m cobaya -q` |
| C | LSODA/reference、oracle tail 与科学精度认证 | 手动或 scheduled，默认串行 | `.github/workflows/slow.yml` 与 `scripts/benchmark_oracle_tail_convergence.py` |
| D | sdist/wheel、归档边界和安装后 smoke | Python 3.11 一次 | `python -m build --sdist --wheel` |
| Static | 中文注释、Ruff、mypy、compileall、diff | Python 3.11 一次 | `.github/workflows/ci.yml` 的 `static` job |

## 正交覆盖

| 风险 | canonical test | interaction case |
| --- | --- | --- |
| reheating kink | `tests/test_freq_adaptive.py` | low-T numerical regression |
| tail matching | `tests/test_fast_sgwb.py` | `tests/test_reference.py` 与 oracle 子集 artifact |
| outer reuse | `tests/test_fast_sgwb.py` | high-T/stiff fallback checks |
| frequency grid / eval nodes | `tests/test_fast_sgwb.py` | `tests/test_engine.py` |
| quadrature | `tests/test_fast_sgwb.py` | default/Sobol diagnostic scripts |
| physical guards | `tests/test_engine.py` | guard-edge stability artifact |
| thread state / determinism | `tests/test_fast_sgwb.py`, `tests/test_engine.py` | sequential configuration calls |

## 资源预算

普通 CI 和 slow workflow 默认设置 `NUMBA_NUM_THREADS=2`、`FAST_THREADS=2`，并将
OMP、BLAS、MKL、VECLIB、NUMEXPR 限制为 1。oracle benchmark 默认 `workers=1`；
任何更高并行度必须由调用者显式指定，并在 artifact 中记录。

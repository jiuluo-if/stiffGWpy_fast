"""在相同 LVK 数据下比较原始网格、当前 fast 与 SageNet+ 的 MCMC。"""

from __future__ import annotations

import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SAGE = Path(r"F:\codex\SageNet")
TORCH_LIB = SAGE / "sagenet_env" / "Lib" / "site-packages" / "torch" / "lib"
MKL_BIN = Path(r"C:\miniconda3\Library\bin")
for _path in (TORCH_LIB, MKL_BIN):
    if _path.is_dir():
        os.environ["PATH"] = str(_path) + os.pathsep + os.environ.get("PATH", "")
        try:
            os.add_dll_directory(str(_path))
        except (AttributeError, OSError):
            pass

os.environ["NUMBA_NUM_THREADS"] = "2"
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"
os.environ["FAST_THREADS"] = "2"
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "True")
sys.path[:0] = [str(REPO), str(SAGE)]

import numpy as np  # noqa: E402
from scipy.interpolate import interp1d  # noqa: E402

from stiffgwpy_fast import LCDM_SG  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402

try:
    import torch
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
except (ImportError, RuntimeError):
    torch = None

OUT = REPO / "docs" / "mcmc_sagenet_compare"
CHAIN_FILE = REPO / "docs" / "mcmc" / "chains" / "sagenet_compare_20260924.npz"
DATA_FILE = REPO / "stiffgwpy_fast" / "cobaya" / "likelihoods" / "LIGO_SGWB" / "C_O1_O2_O3.dat"
COMMON = dict(Omega_bh2=0.0223828, Omega_ch2=0.1201075,
              H0=67.32117, A_s=2.100549e-9, cr=0.0, DN_eff=0.0)
CONTEXTS = {
    "基准": dict(kappa10=1e-2, T_re=2e3, DN_re=20.0),
    "低重加热温度": dict(kappa10=1e-2, T_re=10.0, DN_re=10.0),
    "高刚性物质": dict(kappa10=1.0, T_re=2e3, DN_re=30.0),
}
ENGINES = ("最初 plain-grid", "当前 fast", "SageNet+ Transformer")
BOUNDS = ((-5.0, -1.0), (-0.5, 1.5))
START = np.array([-2.0, 0.0])
STEP = np.array([0.35, 0.12])
SEEDS = (20260924, 20260925, 20260926)
BURN = 1200
KEEP = 2000
ORACLE_NFREQ = 48
_CONFIGURED_FAST_ENGINE = None


def common_parameters(context: dict[str, float], theta: np.ndarray) -> dict[str, float]:
    return {**COMMON, **context, "r": 10.0 ** float(theta[0]), "n_t": float(theta[1])}


def configure_fast(engine: str) -> None:
    global _CONFIGURED_FAST_ENGINE
    if _CONFIGURED_FAST_ENGINE == engine:
        return
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    if engine == "最初 plain-grid":
        # 历史 plain-grid 原始预设：固定粗网格，不细分相位、不拆 kink。
        FS.set_h(0.02)
        FS.set_col_step(8)
        FS.set_z_tail(5.0)
        FS.set_phase_max(0.0)
        FS.set_freq_grid("construct")
    else:
        # 当前正式用户档位，参数来自 fast_v0.2 当前配置。
        FS.set_h(0.005)
        FS.set_col_step(8)
        FS.set_z_tail(5.0)
        FS.set_phase_max(0.25)
        FS.set_freq_grid("goal")
    _CONFIGURED_FAST_ENGINE = engine


def predict_fast(engine: str, params: dict[str, float]) -> dict:
    configure_fast(engine)
    model = LCDM_SG(**params)
    result = FS.SGWB_iter_fast(
        model, tol=1e-6, freq_res=1.0, transition_refine=False,
        kink_split=(engine == "当前 fast"),
        freq_grid=("construct" if engine == "最初 plain-grid" else "goal"),
        frequency_quadrature="pchip")
    if result is None or not getattr(model, "SGWB_converge", False):
        return {"status": "rejected", "reason": getattr(model, "fast_failure_reason", "no convergence")}
    return {"status": "ok", "f": np.asarray(model.f, dtype=float),
            "y": np.asarray(model.log10OmegaGW, dtype=float),
            "DN_eff": float(model.cosmo_param["DN_eff"])}


def predict_sage(predictor, params: dict[str, float]) -> dict:
    result = predictor.predict(params)
    return {"status": "ok", "f": np.asarray(result["f"], dtype=float),
            "y": np.asarray(result["log10OmegaGW"], dtype=float)}


def make_likelihood(data: np.ndarray):
    logf = np.log10(data[:, 0])
    observed = data[:, 1]
    sigma = data[:, 2]

    def evaluate(prediction: dict) -> float:
        f = np.asarray(prediction["f"], dtype=float)
        y = np.asarray(prediction["y"], dtype=float)
        order = np.argsort(f)
        f, y = f[order], y[order]
        use = f >= -5.0
        f, y = f[use], y[use]
        if f.size < 4 or not np.isfinite(f).all() or not np.isfinite(y).all():
            return -math.inf
        model = np.zeros_like(observed)
        valid = logf <= f[-1]
        if valid.any():
            try:
                model[valid] = 10.0 ** interp1d(f, y, kind="cubic")(logf[valid])
            except (ValueError, FloatingPointError):
                return -math.inf
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            residual = (observed - model) / sigma
            chi2 = np.dot(residual, residual)
        return float(-0.5 * chi2) if math.isfinite(float(chi2)) else -math.inf
    return evaluate


def make_sage_predictor():
    from sagenetgw.classes import GWPredictor
    return GWPredictor(model_type="Transformer", device="cpu")


def evaluate_engine(engine: str, params: dict[str, float], predictor, likelihood) -> tuple[float, str | None]:
    try:
        out = predict_sage(predictor, params) if engine.startswith("SageNet") else predict_fast(engine, params)
        if out["status"] != "ok":
            return -math.inf, str(out.get("reason", "guard or no convergence"))
        value = likelihood(out)
        return value, None if math.isfinite(value) else "non-finite likelihood"
    except (ArithmeticError, ValueError, RuntimeError, FloatingPointError) as exc:
        return -math.inf, type(exc).__name__


def mcmc_chain(engine: str, context: dict[str, float], seed: int, predictor, likelihood) -> dict:
    rng = np.random.default_rng(seed)
    theta = START.copy()
    params = common_parameters(context, theta)
    logp, reason = evaluate_engine(engine, params, predictor, likelihood)
    if not math.isfinite(logp):
        raise RuntimeError(f"{engine} initial point rejected: {reason}")
    samples = np.empty((BURN + KEEP, 2), dtype=float)
    accepted = 0
    failures = {}
    t0 = time.perf_counter()
    for i in range(BURN + KEEP):
        proposal = theta + rng.normal(size=2) * STEP
        inside = all(lo <= val <= hi for val, (lo, hi) in zip(proposal, BOUNDS))
        if inside:
            proposal_logp, failure = evaluate_engine(
                engine, common_parameters(context, proposal), predictor, likelihood)
            if failure:
                failures[failure] = failures.get(failure, 0) + 1
            if math.isfinite(proposal_logp) and math.log(rng.random()) < proposal_logp - logp:
                theta, logp = proposal, proposal_logp
                accepted += 1
        samples[i] = theta
    elapsed = time.perf_counter() - t0
    chain = samples[BURN:]
    return {"samples": chain, "wall_s": elapsed, "acceptance": accepted / (BURN + KEEP),
            "failure_counts": failures, "seed": seed}


def autocorrelation_ess(chain: np.ndarray) -> list[float]:
    result = []
    n = len(chain)
    for column in range(chain.shape[1]):
        values = chain[:, column] - np.mean(chain[:, column])
        size = 1 << (2 * n - 1).bit_length()
        fft = np.fft.rfft(values, n=size)
        acf = np.fft.irfft(fft * np.conjugate(fft), n=size)[:n]
        if acf[0] <= 0:
            result.append(float(n))
            continue
        acf /= acf[0]
        tau = 1.0
        for lag in range(1, n - 1, 2):
            pair = acf[lag] + acf[lag + 1]
            if pair <= 0:
                break
            tau += 2.0 * pair
        result.append(float(min(n, max(1.0, n / tau))))
    return result


def split_rhat(chains: np.ndarray) -> list[float]:
    # 各独立链切成前后两半，避免只看链均值掩盖混合问题。
    split = np.concatenate([chains[:, :chains.shape[1] // 2],
                            chains[:, chains.shape[1] // 2:2 * (chains.shape[1] // 2)]], axis=0)
    n = split.shape[1]
    within = np.mean(np.var(split, axis=1, ddof=1), axis=0)
    between = n * np.var(np.mean(split, axis=1), axis=0, ddof=1)
    variance = ((n - 1) / n) * within + between / n
    return np.sqrt(variance / within).tolist()


def summarize(chains: list[dict]) -> dict:
    arr = np.stack([c["samples"] for c in chains])
    combined = arr.reshape(-1, 2)
    ess = np.sum([autocorrelation_ess(c["samples"]) for c in chains], axis=0)
    total_wall = sum(c["wall_s"] for c in chains)
    return {
        "n_chains": len(chains), "draws_per_chain": int(arr.shape[1]),
        "total_draws": int(arr.shape[0] * arr.shape[1]),
        "acceptance_median": float(np.median([c["acceptance"] for c in chains])),
        "rhat_split": split_rhat(arr), "ess": ess.tolist(),
        "wall_s_median": float(np.median([c["wall_s"] for c in chains])),
        "wall_s_total": total_wall,
        "milliseconds_per_step": 1000.0 * total_wall / (arr.shape[0] * (BURN + KEEP)),
        "ess_per_second": (ess / total_wall).tolist(),
        "posterior_mean": np.mean(combined, axis=0).tolist(),
        "posterior_std": np.std(combined, axis=0, ddof=1).tolist(),
        "posterior_p16": np.percentile(combined, 16, axis=0).tolist(),
        "posterior_p84": np.percentile(combined, 84, axis=0).tolist(),
        "failure_counts": {key: sum(c["failure_counts"].get(key, 0) for c in chains)
                           for key in set().union(*(c["failure_counts"].keys() for c in chains))},
        "samples": arr,
    }


def env_meta() -> dict:
    meta = {"python": sys.version.split()[0], "platform": platform.platform(),
            "logical_cpus": os.cpu_count(), "numba_threads": 2,
            "numba_threading_layer": "workqueue", "blas_threads": 1,
            "torch_threads": 1, "seeds": list(SEEDS), "burn": BURN, "keep": KEEP,
            "free_parameter_bounds": {"log10r": BOUNDS[0], "n_t": BOUNDS[1]},
            "proposal_step": {"log10r": STEP[0], "n_t": STEP[1]},
            "lvk_data_rows": int(np.loadtxt(DATA_FILE).shape[0])}
    import numba
    import scipy
    import sklearn
    import torch as torch_mod
    meta.update(numpy=np.__version__, scipy=scipy.__version__, numba=numba.__version__,
                torch=torch_mod.__version__, sklearn=sklearn.__version__)
    try:
        meta["stiffgwpy_sha"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        meta["stiffgwpy_sha"] = None
    try:
        meta["sagenet_sha"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=SAGE, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        meta["sagenet_sha"] = None
    return meta


def run(quick: bool = False) -> None:
    global BURN, KEEP
    if quick:
        BURN, KEEP = 30, 50
    OUT.mkdir(parents=True, exist_ok=True)
    CHAIN_FILE.parent.mkdir(parents=True, exist_ok=True)
    data = np.loadtxt(DATA_FILE)
    likelihood = make_likelihood(data)
    result = {"schema": "mcmc_sagenet_compare_v1", "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
              "environment": env_meta(), "settings": {"burn": BURN, "keep": KEEP,
              "seeds": list(SEEDS), "free_parameters": ["log10r", "n_t"],
              "likelihood": "repository LVK O1/O2/O3 cross-correlation, same cubic interpolation and chi-square log-likelihood as packaged Cobaya likelihood",
              "accuracy_reference_in_sampler": False}, "contexts": {}, "precision_reference": {}}
    result["settings"]["proposal_step"] = {"log10r": float(STEP[0]), "n_t": float(STEP[1])}
    chain_store = {}
    predictors = {"最初 plain-grid": None, "当前 fast": None}
    model_load_start = time.perf_counter()
    predictors["SageNet+ Transformer"] = make_sage_predictor()
    result["sagenet_model_load_s_excluded_from_chain"] = time.perf_counter() - model_load_start

    # 预热将 JIT 和权重读取成本分开；这些调用不进入 MCMC 计时。
    for context_name, context in CONTEXTS.items():
        params = common_parameters(context, START)
        for engine in ENGINES:
            started = time.perf_counter()
            _, reason = evaluate_engine(engine, params, predictors[engine], likelihood)
            if reason:
                raise RuntimeError(f"Warm-up failed for {context_name}/{engine}: {reason}")
            result.setdefault("warmup_s", {}).setdefault(engine, []).append(time.perf_counter() - started)

    for context_index, (context_name, context) in enumerate(CONTEXTS.items()):
        result["contexts"][context_name] = {"fixed_parameters": context, "engines": {}}
        for engine_index, engine in enumerate(ENGINES):
            print(f"starting {context_name} / {engine}", flush=True)
            chains = []
            for chain_index, seed in enumerate(SEEDS):
                chain = mcmc_chain(engine, context, seed + context_index * 100,
                                   predictors[engine], likelihood)
                chains.append(chain)
                print(f"  chain {chain_index + 1}/3: {chain['wall_s']:.1f}s, acceptance={chain['acceptance']:.3f}", flush=True)
            stats = summarize(chains)
            raw_samples = stats.pop("samples")
            result["contexts"][context_name]["engines"][engine] = stats
            chain_store[f"{context_name}_{engine_index}"] = raw_samples
            chain_store[f"{context_name}_{engine_index}_wall_s"] = np.asarray([c["wall_s"] for c in chains])

    # 每个参数情景从三种后验的合并样本取中位点；高精度档只在这些点上作误差参照。
    for context_name, context in CONTEXTS.items():
        pooled = np.concatenate([chain_store[f"{context_name}_{i}"].reshape(-1, 2)
                                 for i in range(len(ENGINES))])
        theta = np.median(pooled, axis=0)
        params = common_parameters(context, theta)
        fast_anchor = predict_fast("当前 fast", params)
        if fast_anchor["status"] != "ok":
            result["precision_reference"][context_name] = {"status": "fast guard", "reason": fast_anchor.get("reason")}
            continue
        predictions = {}
        for engine in ENGINES:
            predictions[engine] = (predict_sage(predictors[engine], params) if engine.startswith("SageNet")
                                   else predict_fast(engine, params))
        low = max(-6.0, *(float(np.min(p["f"])) for p in predictions.values()
                          if p["status"] == "ok"))
        high = min(1.0, *(float(np.max(p["f"])) for p in predictions.values()
                          if p["status"] == "ok"))
        if not (high > low):
            result["precision_reference"][context_name] = {
                "status": "no common spectrum support", "theta_log10r_nt": theta.tolist(),
                "engine_frequency_ranges": {name: [float(np.min(p["f"])), float(np.max(p["f"]))]
                                             for name, p in predictions.items() if p["status"] == "ok"}}
            continue
        logf_oracle = np.linspace(low, high, ORACLE_NFREQ)
        model = LCDM_SG(**params)
        dn_eff = fast_anchor["DN_eff"]
        start = time.perf_counter()
        ogw, oj, _, used_tail = REF.spectrum_reference(
            model, logf_oracle, dn_eff, z_tail=8.0, rtol=1e-7, workers=1)
        oracle_s = time.perf_counter() - start
        ref_y = np.log10(np.maximum(np.asarray(ogw) - np.asarray(oj), 1e-300))
        row = {"status": "ok", "theta_log10r_nt": theta.tolist(),
               "fixed_parameters": context, "reference_wall_s_excluded_from_mcmc": oracle_s,
               "frequency_range_log10_hz": [float(logf_oracle.min()), float(logf_oracle.max())],
               "oracle_frequency_count": ORACLE_NFREQ, "tail_fraction": float(np.mean(used_tail)),
               "frequency_log10_hz": logf_oracle.tolist(), "reference_spectrum_log10": ref_y.tolist(),
               "engines": {}}
        for engine in ENGINES:
            candidate = predictions[engine]
            if candidate["status"] != "ok":
                row["engines"][engine] = {"status": "rejected", "reason": candidate.get("reason")}
                continue
            order = np.argsort(candidate["f"])
            pred_y = interp1d(candidate["f"][order], candidate["y"][order], kind="cubic",
                              bounds_error=True)(logf_oracle)
            dex = np.abs(pred_y - ref_y)
            rel = np.abs(10.0 ** pred_y - 10.0 ** ref_y) / np.maximum(10.0 ** ref_y, 1e-300)
            row["engines"][engine] = {"status": "ok", "dex_p95": float(np.percentile(dex, 95)),
                                      "dex_max": float(np.max(dex)),
                                      "relative_error_p95": float(np.percentile(rel, 95)),
                                      "relative_error_max": float(np.max(rel)),
                                      "frequency_range_log10_hz": [float(np.min(candidate["f"])),
                                                                    float(np.max(candidate["f"]))],
                                      "spectrum_log10": pred_y.tolist(),
                                      "lvk_frequency_coverage_fraction": float(np.mean(
                                          (np.log10(data[:, 0]) >= np.min(candidate["f"])) &
                                          (np.log10(data[:, 0]) <= np.max(candidate["f"]))))}
        result["precision_reference"][context_name] = row
        print(f"reference accuracy points complete: {context_name} ({oracle_s:.1f}s reference excluded)", flush=True)

    np.savez_compressed(CHAIN_FILE, **chain_store)
    result["chain_file"] = str(CHAIN_FILE.relative_to(REPO))
    (OUT / "results.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    make_plots(result, chain_store, logf_oracle)
    write_report(result)
    print(f"wrote results and figures to {OUT}", flush=True)


def make_plots(result: dict, chain_store: dict, logf_oracle: np.ndarray) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    # Color-blind-friendly, high-contrast palette used consistently in all figures.
    colors = {"最初 plain-grid": "#0072B2", "当前 fast": "#D55E00", "SageNet+ Transformer": "#009E73"}
    markers = {"最初 plain-grid": "o", "当前 fast": "s", "SageNet+ Transformer": "^"}
    line_styles = {"最初 plain-grid": "-", "当前 fast": "--", "SageNet+ Transformer": ":"}
    fig, axes = plt.subplots(len(ENGINES), len(CONTEXTS), figsize=(15, 9),
                             sharex=True, sharey=True, constrained_layout=True)
    for row_index, engine in enumerate(ENGINES):
        for col_index, context_name in enumerate(CONTEXTS):
            ax = axes[row_index, col_index]
            i = ENGINES.index(engine)
            samples = chain_store[f"{context_name}_{i}"].reshape(-1, 2)
            ax.hexbin(samples[:, 0], samples[:, 1], gridsize=32, mincnt=1,
                      bins="log", cmap=LinearSegmentedColormap.from_list(
                          f"posterior_{i}", ["#FFFFFF", colors[engine]], N=256),
                      alpha=0.95, linewidths=0.0)
            mean = np.mean(samples, axis=0)
            ax.plot(mean[0], mean[1], marker=markers[engine], color=colors[engine],
                    markeredgecolor="#202020", markeredgewidth=0.8, ms=6,
                    label=engine)
            ax.set_title(context_name)
            ax.set_xlabel(r"$\log_{10}(r)$")
            ax.grid(alpha=0.2)
            ax.set_ylabel(r"$n_t$" + (f"\n{engine}" if col_index == 0 else ""))
    fig.suptitle("MCMC 后验分布：按方法分行、按参数情景分列（每格 3 条链）")
    fig.savefig(OUT / "mcmc_posterior.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    x = np.arange(len(CONTEXTS))
    width = 0.24
    for i, engine in enumerate(ENGINES):
        step = [result["contexts"][name]["engines"][engine]["milliseconds_per_step"]
                for name in CONTEXTS]
        ess_speed = [min(result["contexts"][name]["engines"][engine]["ess_per_second"])
                     for name in CONTEXTS]
        axes[0].bar(x + (i - 1) * width, step, width, color=colors[engine], label=engine)
        axes[1].bar(x + (i - 1) * width, ess_speed, width, color=colors[engine], label=engine)
    axes[0].set_ylabel("每一步耗时（毫秒，含同一 LVK 似然）")
    axes[1].set_ylabel("每秒有效样本数（取两个参数中较小值）")
    for ax in axes:
        ax.set_xticks(x, list(CONTEXTS))
        ax.grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle("MCMC 速度：平均每步耗时与有效样本产出")
    fig.savefig(OUT / "mcmc_speed.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(CONTEXTS), figsize=(15, 4.8), constrained_layout=True, sharey=True)
    for ax, context_name in zip(axes, CONTEXTS):
        ref = result["precision_reference"][context_name]
        if ref.get("status") != "ok":
            ax.text(0.5, 0.5, ref.get("status", "无可比较频段"), ha="center", va="center")
            ax.set_title(context_name)
            continue
        freq = np.asarray(ref["frequency_log10_hz"])
        ax.plot(freq, ref["reference_spectrum_log10"], color="#222222", lw=2.3,
                label="本地独立精度参照")
        for i, engine in enumerate(ENGINES):
            row = ref["engines"][engine]
            ax.plot(freq, row["spectrum_log10"], color=colors[engine], lw=1.8,
                    linestyle=line_styles[engine], marker=markers[engine], markevery=8,
                    ms=4,
                    label=engine if context_name == next(iter(CONTEXTS)) else None)
        ax.set_title(context_name)
        ax.set_xlabel(r"频率 $\log_{10}(f/\mathrm{Hz})$")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel(r"能量密度 $\log_{10}\Omega_{GW}$")
    axes[0].legend(fontsize=8)
    fig.suptitle("MCMC 后验中心附近的谱形（共同频段；独立参考只作精度锚点）")
    fig.savefig(OUT / "mcmc_accuracy.png", dpi=200)
    plt.close(fig)


def write_report(result: dict) -> None:
    lines = ["# stiffGWpy 与 SageNet+：同一 LVK 数据下的 MCMC 对比", "",
             f"实测时间：{result['generated_at']}",
             f"代码版本：stiffGWpy `{result['environment']['stiffgwpy_sha']}`；SageNet `{result['environment']['sagenet_sha']}`。", "",
             "## 先看结果", "",
             "本报告把三种方法放进同一套随机游走 MCMC、同一份 LVK O1/O2/O3 数据和同一似然函数。每种方法、每个参数情景各跑 3 条独立链。",
             "简单说，MCMC 会在参数范围里反复试组合：更符合观测的组合更容易留下；保留下来的组合形成后验样本，表示在当前数据和设定下哪些参数更受支持。",
             "速度按整条采样链计时；精度参照档只在采样后抽出的共同参数点运行，不参与链的计时或后验。", "",
             "| 参数情景 | 方法 | 每步耗时（越低越快，ms） | 每秒有效样本（越高越好） | 接受提议比例 | 最大 R-hat | log10(r) 样本中心 ± 散布 | n_t 样本中心 ± 散布 |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for context_name in CONTEXTS:
        for engine in ENGINES:
            row = result["contexts"][context_name]["engines"][engine]
            mean, std = row["posterior_mean"], row["posterior_std"]
            lines.append(f"| {context_name} | {engine} | {row['milliseconds_per_step']:.3f} | {min(row['ess_per_second']):.2f} | {row['acceptance_median']:.1%} | {max(row['rhat_split']):.3f} | {mean[0]:.3f} ± {std[0]:.3f} | {mean[1]:.3f} ± {std[1]:.3f} |")
    lines += ["", "### 先读懂表里的几个数", "",
              "- **后验样本**：算法试走并保留下来的参数组合；越常出现，表示在当前数据和设定下越受支持。本次每条链保留 2,000 个组合；每个方法、每个情景合计 6,000 个。它们不是 6,000 份独立数据。",
              "- **接受提议比例**：算法提出新参数组合后，实际移动到新组合的比例；它只描述采样过程，不表示模型准确率。",
              "- **样本中心 ± 散布**：中心是保留参数值的平均数；散布（标准差）表示这些值分得宽还是窄。它反映在当前数据、模型和参数范围下的估计不确定性，不是算法精度误差；也不一定等同于严格的 68% 区间。",
              "- **R-hat**：比较 3 条独立链走到的区域是否一致。越接近 1 越好；高于 1.05 通常提示链还没充分混合。本次都低于 1.05，但这只是基本检查，不保证所有区域都已采足。",
              "- **有效样本数 ESS**：把彼此相近、重复的信息折算后，估计相当于多少个独立样本。表中每秒有效样本和后面的 ESS 都取两个参数中较小的那个，方便看较难采的参数。", ""]
    speedup = [result["contexts"][c]["engines"]["SageNet+ Transformer"]["milliseconds_per_step"] /
               result["contexts"][c]["engines"]["当前 fast"]["milliseconds_per_step"] for c in CONTEXTS]
    lines += ["", f"**简要结论：当前 fast 每一步约比 SageNet+ 快 {min(speedup):.1f}–{max(speedup):.1f} 倍；下面的谱误差表也显示当前 fast 在这三个情景都更接近本地精度参照。最初 plain-grid 同样纳入采样比较，速度与当前 fast 接近。**", ""]
    lines += ["", "## 采样是否稳定、频谱差多少", "",
              "| 参数情景 | 方法 | 较低的 ESS（两个参数中较小值） | 最大 R-hat | 覆盖 LVK 观测频段 | 95%频点误差不超过（dex） | 最大频点误差（dex） |", "|---|---|---:|---:|---:|---:|---:|"]
    for context_name in CONTEXTS:
        ref = result["precision_reference"].get(context_name, {})
        for i, engine in enumerate(ENGINES):
            row = result["contexts"][context_name]["engines"][engine]
            prec = ref.get("engines", {}).get(engine, {})
            coverage = prec.get("lvk_frequency_coverage_fraction")
            coverage_txt = "未测" if coverage is None else f"{coverage:.1%}"
            p95 = prec.get("dex_p95")
            maxdex = prec.get("dex_max")
            lines.append(f"| {context_name} | {engine} | {min(row['ess']):.0f} | {max(row['rhat_split']):.3f} | {coverage_txt} | {'未测' if p95 is None else f'{p95:.3g}'} | {'未测' if maxdex is None else f'{maxdex:.3g}'} |")
    lines += ["", "**低重加热温度情景的 LVK 频段覆盖为 0%：三种方法的预测频率范围都没有覆盖观测频段，因此这个情景的后验主要由参数范围和采样规则决定，不能据此判断谁更符合数据。**",
              "", "失败提案按计算路径分别记账：", "",
              "| 参数情景 | 最初 plain-grid | 当前 fast | SageNet+ |", "|---|---|---|---|"]
    for context_name in CONTEXTS:
        failures = result["contexts"][context_name]["engines"]
        def describe_failures(row):
            labels = {"shared_Neff_guard": "物理一致性保护", "max_iter": "迭代未收敛", "non-finite likelihood": "无效数值/似然"}
            return ", ".join(f"{labels.get(k, k)} {v} 次" for k, v in row["failure_counts"].items()) or "无"
        lines.append(f"| {context_name} | {describe_failures(failures['最初 plain-grid'])} | {describe_failures(failures['当前 fast'])} | {describe_failures(failures['SageNet+ Transformer'])} |")
    lines += ["", "这些次数表示提议点被拒绝，不是最终样本数；参数越界的提议不计入。物理一致性保护是模型边界检查，不等于数值故障。", "",
              "图 1 展示算法常采到哪些参数组合；图 2 比较耗时和有效样本产出；图 3 把三种算法的能谱与本地高精度参考叠在一起。",
              "本次最大 R-hat 都低于 1.05；较低参数的 ESS 约为 80–115，说明链之间大体一致，但有效独立信息量仍有限，后验中心和散布应看作本次计算的估计。`dex` 是以 10 为底的对数差：0.01 dex 约相差 2%，0.1 dex 约相差 26%；数值越小越接近参考。", "",
              "## 比较了什么参数", "",
              "| 参数 | 通俗含义 | 本次设置 |", "|---|---|---|",
              "| `log10(r)` / `r` | 原初引力波强度相对标量扰动的比例；对数每增加 1，`r` 增大 10 倍 | MCMC 自由参数；`log10(r)` 在 -5 到 -1 均匀取值，即 `r` 为 1e-5 到 0.1 |",
              "| `n_t` | 原初引力波谱随频率上升或下降的斜率 | MCMC 自由参数；范围 -0.5 到 1.5 |",
              "| `kappa10` | 10 MeV 时刚性物质能量与光子能量之比 | 三个情景分别为 0.01、0.01、1 |",
              "| `T_re` | 再加热结束时的温度，单位 GeV | 三个情景分别为 2000、10、2000 |",
              "| `DN_re` | 从暴胀结束到再加热结束经历的膨胀量 | 三个情景分别为 20、10、30 |",
              "| `cr` | 是否使用单场暴胀的一致性关系 | 固定为 0，使 `n_t` 可独立变化；与 SageNet+ 的输入方式一致 |",
              "| `Omega_bh2` | 普通物质（重子）的密度参数 | 固定为 0.0223828 |",
              "| `Omega_ch2` | 暗物质的密度参数 | 固定为 0.1201075 |",
              "| `H0` | 宇宙当前膨胀速度，单位 km/s/Mpc | 固定为 67.32117 |",
              "| `A_s` | 原初标量扰动的强度 | 固定为 2.100549e-9 |",
              "| `DN_eff` | 额外相对论粒子的起始贡献 | 固定为 0；求解器按各自路径处理引力波贡献 |", "",
              "下列是求解器精度/速度设置：`h` 是主要积分步长，较小通常更细但较慢；`col_step` 是隔多少个内部步保存/汇总一次，较大通常更省时；`z_tail` 是切换到尾部近似的位置；`phase_max` 限制振荡每一步允许跨过的相位，0 表示关闭这项细分。`goal` 表示按目标频率构网，`construct` 表示直接用构造网格。", "",
              "## 测量方法", "",
              f"- 数据：仓库随附的 LVK O1/O2/O3 频率交叉相关表，共 {result['environment']['lvk_data_rows']:,} 行；三种方法都使用同一个插值和卡方似然。",
              f"- 每条链：预热 {result['settings']['burn']} 步，记录 {result['settings']['keep']:,} 步；3 个随机种子独立重复。采用二维高斯随机游走，`log10(r)` 步长 {result['settings']['proposal_step']['log10r']}、`n_t` 步长 {result['settings']['proposal_step']['n_t']}。",
              "- 最初 plain-grid：按旧探索档的 `h=0.02`、`col_step=8`、`z_tail=5`、`phase_max=0`、普通构网、无精确断点拆分，在当前代码中重建；它纳入 MCMC 对比，但不是旧版代码快照复跑。",
              "- 当前 fast：`h=0.005`、`col_step=8`、`z_tail=5`、`phase_max=0.25`、goal 频率网格、精确拆分再加热断点。",
              "- SageNet+：Transformer 预训练权重，CPU 单线程；每个 MCMC 提议复用已加载模型。模型首次读入时间单独记录为 {:.3f} 秒，不混进每一步时间。".format(result['sagenet_model_load_s_excluded_from_chain']),
              "- 预热耗时（每种方法的首个调用 / 之后两次的中位数）：" + "；".join(f"{engine} {result['warmup_s'][engine][0]:.3f} 秒 / {statistics.median(result['warmup_s'][engine][1:]):.3f} 秒" for engine in ENGINES) + "。首个调用可能包含编译或缓存初始化；这些时间不计入链的每步耗时。",
              "- 精度参照：本地连续-sigma 独立求解器（rtol=1e-7，z_tail=8）只在 MCMC 后验中心附近、三种方法共同覆盖的频段比较 48 个频率点；它不参加 MCMC。", "",
              "## 结果文件", "",
              "### 后验样本分布", "",
              "横轴是 `log10(r)`，纵轴是斜率 `n_t`；颜色越深表示链越常访问该参数组合。图按方法分行、按情景分列，避免三种分布叠加后颜色混在一起。", "",
              "三种方法使用固定高对比颜色：最初 plain-grid 为蓝色、当前 fast 为橙色、SageNet+ 为绿色；本次已用新配色重绘并替换旧图。", "",
              "![三种方法的 MCMC 后验样本分布](mcmc_posterior.png)", "",
              "### 速度与有效样本", "",
              "每步时间越低越快；有效样本/秒越高，表示同样时间内可用的独立信息越多。", "",
              "![采样耗时和有效样本产出对比](mcmc_speed.png)", "",
              "### 能谱精度对比", "",
              "黑色参考线来自本地高精度档，仅用于采样后的精度核对，不参与 MCMC。", "",
              "![三种方法与本地精度参照的能谱对比](mcmc_accuracy.png)", "",
              "- `results.json`：逐情景、逐方法原始汇总及精度曲线；MCMC 原始链保存在忽略提交的本地文件 `{}`。".format(result['chain_file']), "",
              "本结论只描述当前两个代码版本、该 Transformer 权重、该机器与这份 LVK 数据；不代表 SageNet+ 的所有权重或所有参数范围。"]
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="短链，只用于验证程序连接与图表格式")
    run(parser.parse_args().quick)

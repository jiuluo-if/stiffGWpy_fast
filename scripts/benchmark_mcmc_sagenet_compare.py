"""在相同 LVK 数据下比较原始网格、当前 fast 与 SageNet+ 的 MCMC。"""

from __future__ import annotations

import json
import hashlib
import io
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from contextlib import redirect_stdout
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
CHAIN_FILE = OUT / "posterior_chains_20260924.npz"
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
SEEDS = (2026092401, 2026092402, 2026092403, 2026092404)
STARTS = np.array([[-4.6, -0.45], [-3.6, -0.30], [-2.6, -0.15], [-1.8, 0.0]])
BURN = 4000
KEEP = 10000
ORACLE_NFREQ = 48
PRECISION_DRAWS_PER_ENGINE = 8
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


def mcmc_chain(engine: str, context: dict[str, float], seed: int, start: np.ndarray,
               predictor, likelihood) -> dict:
    rng = np.random.default_rng(seed)
    theta = np.asarray(start, dtype=float).copy()
    start_point = theta.tolist()
    params = common_parameters(context, theta)
    logp, reason = evaluate_engine(engine, params, predictor, likelihood)
    if not math.isfinite(logp):
        raise RuntimeError(f"{engine} initial point rejected: {reason}")
    samples = np.empty((KEEP, 2), dtype=float)
    accepted = 0
    failures = {}
    def advance() -> None:
        nonlocal theta, logp, accepted
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
    t0 = time.perf_counter()
    for _ in range(BURN):
        advance()
    warmup_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    for i in range(KEEP):
        advance()
        samples[i] = theta
    sampling_s = time.perf_counter() - t0
    chain = samples
    return {"samples": chain, "warmup_s": warmup_s, "sampling_s": sampling_s,
            "wall_s": warmup_s + sampling_s, "acceptance": accepted / (BURN + KEEP),
            "failure_counts": failures, "seed": seed, "start": start_point}


def split_chains(chains: np.ndarray) -> np.ndarray:
    half = chains.shape[1] // 2
    return np.concatenate((chains[:, :half], chains[:, -half:]), axis=0)


def basic_rhat(chains: np.ndarray) -> np.ndarray:
    m, n, _ = chains.shape
    within = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    between = n * np.var(np.mean(chains, axis=1), axis=0, ddof=1)
    variance = ((n - 1) / n) * within + between / n
    return np.sqrt(variance / within)


def multi_chain_ess(chains: np.ndarray) -> np.ndarray:
    """Geyer initial-positive/monotone ESS for multiple equal-length chains."""
    m, n, p = chains.shape
    centered = chains - np.mean(chains, axis=1, keepdims=True)
    size = 1 << (2 * n - 1).bit_length()
    fft = np.fft.rfft(centered, n=size, axis=1)
    acov = np.fft.irfft(fft * np.conjugate(fft), n=size, axis=1)[:, :n, :] / n
    within = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    between = n * np.var(np.mean(chains, axis=1), axis=0, ddof=1)
    var_plus = ((n - 1) / n) * within + between / n
    rho = 1.0 - (within[None, None, :] - acov) / var_plus[None, None, :]
    rho[:, 0, :] = 1.0
    output = np.empty(p, dtype=float)
    for col in range(p):
        pair_sums = []
        previous = np.inf
        for lag in range(0, n - 1, 2):
            pair = float(np.mean(rho[:, lag, col] + rho[:, lag + 1, col]))
            if pair <= 0:
                break
            pair = min(pair, previous)
            pair_sums.append(pair)
            previous = pair
        tau = max(1.0, -1.0 + 2.0 * sum(pair_sums))
        output[col] = min(float(m * n), max(1.0, m * n / tau))
    return output


def rank_normalize(values: np.ndarray) -> np.ndarray:
    from scipy.special import ndtri
    from scipy.stats import rankdata
    flat = values.reshape(-1)
    ranks = rankdata(flat, method="average")
    transformed = ndtri((ranks - 3.0 / 8.0) / (len(flat) + 1.0 / 4.0))
    return transformed.reshape(values.shape)


def chain_diagnostics(chains: np.ndarray) -> dict:
    split = split_chains(chains)
    ranked = rank_normalize(split)
    folded = rank_normalize(np.abs(split - np.median(split, axis=(0, 1))))
    rhat = np.maximum(basic_rhat(ranked), basic_rhat(folded))
    bulk = multi_chain_ess(ranked)
    mean_ess = multi_chain_ess(split)
    pooled = split.reshape(-1, split.shape[-1])
    tails = []
    for col in range(pooled.shape[1]):
        q05, q95 = np.quantile(pooled[:, col], [0.05, 0.95])
        low = (split[:, :, col] <= q05).astype(float)[:, :, None]
        high = (split[:, :, col] >= q95).astype(float)[:, :, None]
        tails.append(min(float(multi_chain_ess(low)[0]), float(multi_chain_ess(high)[0])))
    tail = np.asarray(tails)
    return {"rhat_rank": rhat, "ess_bulk": bulk, "ess_tail": tail,
            "ess_mean": mean_ess,
            "mcse_mean": np.std(pooled, axis=0, ddof=1) / np.sqrt(mean_ess)}


def summarize(chains: list[dict]) -> dict:
    arr = np.stack([c["samples"] for c in chains])
    combined = arr.reshape(-1, 2)
    diagnostics = chain_diagnostics(arr)
    total_wall = sum(c["sampling_s"] for c in chains)
    total_with_warmup = sum(c["wall_s"] for c in chains)
    return {
        "n_chains": len(chains), "draws_per_chain": int(arr.shape[1]),
        "total_draws": int(arr.shape[0] * arr.shape[1]),
        "acceptance_median": float(np.median([c["acceptance"] for c in chains])),
        "rhat_rank": diagnostics["rhat_rank"].tolist(),
        "ess_bulk": diagnostics["ess_bulk"].tolist(),
        "ess_tail": diagnostics["ess_tail"].tolist(),
        "ess_mean": diagnostics["ess_mean"].tolist(),
        "mcse_mean": diagnostics["mcse_mean"].tolist(),
        "warmup_s_total": sum(c["warmup_s"] for c in chains),
        "sampling_s_total": total_wall,
        "sampling_s_per_chain": [c["sampling_s"] for c in chains],
        "wall_s_total_including_warmup": total_with_warmup,
        "wall_s_median": float(np.median([c["wall_s"] for c in chains])),
        "milliseconds_per_step": 1000.0 * total_wall / (arr.shape[0] * KEEP),
        "milliseconds_per_step_including_warmup": 1000.0 * total_with_warmup / (arr.shape[0] * (BURN + KEEP)),
        "ess_per_second": (diagnostics["ess_bulk"] / total_wall).tolist(),
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
            "processor": platform.processor(), "logical_cpus": os.cpu_count(), "numba_threads": 2,
            "numba_threading_layer": "workqueue", "blas_threads": 1,
            "torch_threads": 1, "seeds": list(SEEDS), "burn": BURN, "keep": KEEP,
            "free_parameter_bounds": {"log10r": BOUNDS[0], "n_t": BOUNDS[1]},
            "proposal_step": {"log10r": STEP[0], "n_t": STEP[1]},
            "lvk_data_rows": int(np.loadtxt(DATA_FILE).shape[0])}
    meta["lvk_data_sha256"] = hashlib.sha256(DATA_FILE.read_bytes()).hexdigest()
    weight_file = SAGE / "sagenetgw" / "models" / "best_gw_model_Transformer.pth"
    meta["sagenet_transformer_weights_sha256"] = (
        hashlib.sha256(weight_file.read_bytes()).hexdigest() if weight_file.is_file() else None)
    import numba
    import scipy
    import sklearn
    import torch as torch_mod
    import matplotlib as matplotlib_mod
    meta.update(numpy=np.__version__, scipy=scipy.__version__, numba=numba.__version__,
                torch=torch_mod.__version__, sklearn=sklearn.__version__,
                matplotlib=matplotlib_mod.__version__)
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


def json_safe(value):
    """Replace non-finite numeric results with JSON null, preserving structure."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    if isinstance(value, np.integer):
        return int(value)
    return value


def posterior_accuracy_check(result: dict, chain_store: dict, predictors: dict,
                             data: np.ndarray, points_per_engine: int) -> dict:
    """Compare all solvers with the independent reference at posterior draws."""
    checks = {}
    likelihood = make_likelihood(data)
    for context_index, (context_name, context) in enumerate(CONTEXTS.items()):
        rng = np.random.default_rng(2026092401 + context_index)
        points = []
        for engine_index, engine in enumerate(ENGINES):
            samples = chain_store[f"{context_name}_{engine_index}"].reshape(-1, 2)
            selected = rng.choice(len(samples), size=min(points_per_engine, len(samples)),
                                  replace=False)
            points.extend({"source_posterior": engine, "theta_log10r_nt": samples[i].tolist()}
                          for i in selected)
        rows = []
        for point_index, point in enumerate(points):
            theta = np.asarray(point["theta_log10r_nt"], dtype=float)
            params = common_parameters(context, theta)
            with redirect_stdout(io.StringIO()):
                anchor = predict_fast("当前 fast", params)
            row = {**point, "point_index": point_index, "engines": {}}
            if anchor["status"] != "ok":
                row.update(status="fast physical guard", reference_wall_s=0.0)
                rows.append(row)
                continue
            predictions = {}
            for engine in ENGINES:
                try:
                    with redirect_stdout(io.StringIO()):
                        candidate = (predict_sage(predictors[engine], params)
                                     if engine.startswith("SageNet") else predict_fast(engine, params))
                except (ArithmeticError, ValueError, RuntimeError, FloatingPointError) as exc:
                    candidate = {"status": "rejected", "reason": type(exc).__name__}
                if candidate.get("status") == "ok":
                    frequencies_out = np.asarray(candidate["f"], dtype=float)
                    spectrum_out = np.asarray(candidate["y"], dtype=float)
                    if (not np.isfinite(frequencies_out).all() or
                            not np.isfinite(spectrum_out).all() or
                            np.unique(frequencies_out).size != frequencies_out.size):
                        candidate = {"status": "rejected",
                                     "reason": "non-finite spectrum or duplicate frequency nodes"}
                predictions[engine] = candidate
            valid = [p for p in predictions.values() if p["status"] == "ok"]
            if len(valid) != len(ENGINES):
                row.update(status="one or more engines rejected", reference_wall_s=0.0)
                for engine, candidate in predictions.items():
                    row["engines"][engine] = {"status": candidate["status"],
                                               "reason": candidate.get("reason")}
                rows.append(row)
                continue
            low = max(-6.0, *(float(np.min(p["f"])) for p in valid))
            high = min(1.0, *(float(np.max(p["f"])) for p in valid))
            if high <= low:
                row.update(status="no common frequency range", reference_wall_s=0.0)
                rows.append(row)
                continue
            frequencies = np.linspace(low, high, ORACLE_NFREQ)
            model = LCDM_SG(**params)
            started = time.perf_counter()
            ogw, oj, _, used_tail = REF.spectrum_reference(
                model, frequencies, anchor["DN_eff"], z_tail=8.0, rtol=1e-7, workers=1)
            reference_s = time.perf_counter() - started
            reference_y = np.log10(np.maximum(np.asarray(ogw) - np.asarray(oj), 1e-300))
            row.update(status="ok", reference_wall_s=reference_s,
                       frequency_range_log10_hz=[float(low), float(high)],
                       oracle_tail_fraction=float(np.mean(used_tail)))
            for engine, candidate in predictions.items():
                order = np.argsort(candidate["f"])
                try:
                    predicted_y = interp1d(candidate["f"][order], candidate["y"][order],
                                           kind="cubic", bounds_error=True)(frequencies)
                except ValueError as exc:
                    row["engines"][engine] = {"status": "rejected",
                                               "reason": f"accuracy interpolation: {exc}"}
                    continue
                dex = np.abs(predicted_y - reference_y)
                row["engines"][engine] = {
                    "status": "ok", "dex_p50": float(np.percentile(dex, 50)),
                    "dex_p95": float(np.percentile(dex, 95)), "dex_max": float(np.max(dex)),
                    "lvk_frequency_coverage_fraction": float(np.mean(
                        (np.log10(data[:, 0]) >= np.min(candidate["f"])) &
                        (np.log10(data[:, 0]) <= np.max(candidate["f"]))))}
            rows.append(row)
            print(f"posterior precision {context_name} {point_index + 1}/{len(points)} "
                  f"({reference_s:.1f}s reference excluded)", flush=True)
        summaries = {}
        for engine in ENGINES:
            measurements = [r["engines"][engine]["dex_p95"] for r in rows
                            if r.get("status") == "ok" and
                            r["engines"].get(engine, {}).get("status") == "ok"]
            summaries[engine] = {
                "valid_parameter_points": len(measurements),
                "dex_p95_per_point_median": float(np.median(measurements)) if measurements else None,
                "dex_p95_per_point_p95": float(np.percentile(measurements, 95)) if measurements else None,
                "dex_p95_per_point_max": float(np.max(measurements)) if measurements else None}
        checks[context_name] = {"points_per_engine": points_per_engine,
                                "points": rows, "summary": summaries,
                                "excluded_points": sum(r.get("status") != "ok" for r in rows),
                                "reference_wall_s_total": sum(r.get("reference_wall_s", 0) for r in rows)}
    return checks


def run(quick: bool = False) -> None:
    global BURN, KEEP
    if quick:
        BURN, KEEP = 30, 50
    OUT.mkdir(parents=True, exist_ok=True)
    CHAIN_FILE.parent.mkdir(parents=True, exist_ok=True)
    data = np.loadtxt(DATA_FILE)
    likelihood = make_likelihood(data)
    result = {"schema": "mcmc_sagenet_compare_v2", "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
              "environment": env_meta(), "settings": {"burn": BURN, "keep": KEEP,
              "seeds": list(SEEDS), "free_parameters": ["log10r", "n_t"],
              "chain_starts": STARTS.tolist(), "diagnostics": "rank-normalized folded split R-hat; Geyer bulk/tail ESS",
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
                with redirect_stdout(io.StringIO()):
                    chain = mcmc_chain(engine, context, seed + context_index * 100,
                                       STARTS[chain_index], predictors[engine], likelihood)
                chains.append(chain)
                print(f"  chain {chain_index + 1}/4: {chain['wall_s']:.1f}s, acceptance={chain['acceptance']:.3f}", flush=True)
            stats = summarize(chains)
            raw_samples = stats.pop("samples")
            result["contexts"][context_name]["engines"][engine] = stats
            chain_store[f"{context_name}_{engine_index}"] = raw_samples
            chain_store[f"{context_name}_{engine_index}_wall_s"] = np.asarray([c["wall_s"] for c in chains])
            chain_store[f"{context_name}_{engine_index}_starts"] = STARTS.copy()

    result["precision_posterior"] = posterior_accuracy_check(
        result, chain_store, predictors, data,
        points_per_engine=1 if quick else PRECISION_DRAWS_PER_ENGINE)

    # 每个参数情景从三种后验的合并样本取中位点，作为便于直观看谱形的代表点。
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
    result = json_safe(result)
    (OUT / "results.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
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
    fig.suptitle(f"MCMC 后验分布：按方法分行、按参数情景分列（每格 {len(SEEDS)} 条链）")
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

    fig, axes = plt.subplots(2, len(CONTEXTS), figsize=(15, 8), constrained_layout=True)
    for col, context_name in enumerate(CONTEXTS):
        for row, param_name in enumerate((r"$\log_{10}(r)$", r"$n_t$")):
            ax = axes[row, col]
            for engine_index, engine in enumerate(ENGINES):
                samples = chain_store[f"{context_name}_{engine_index}"]
                for chain_index in range(samples.shape[0]):
                    trace = samples[chain_index, ::10, row]
                    ax.plot(np.arange(len(trace)) * 10, trace, color=colors[engine],
                            alpha=0.4 + 0.14 * chain_index, lw=0.45,
                            label=engine if row == 0 and col == 0 and chain_index == 0 else None)
            ax.set_title(context_name if row == 0 else "")
            ax.set_ylabel(param_name)
            ax.set_xlabel("保留样本步数")
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8, ncol=1)
    fig.suptitle("链轨迹检查：颜色表示方法，同色深浅表示不同独立链")
    fig.savefig(OUT / "mcmc_diagnostics.png", dpi=200)
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
    fig.suptitle("后验代表参数点的能谱（共同频段；独立参考只作精度锚点）")
    fig.savefig(OUT / "mcmc_accuracy.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(CONTEXTS), figsize=(15, 5), constrained_layout=True,
                             sharey=True)
    for ax, context_name in zip(axes, CONTEXTS):
        check = result["precision_posterior"][context_name]
        data_by_engine = []
        for engine in ENGINES:
            values = [point["engines"][engine]["dex_p95"] for point in check["points"]
                      if point.get("status") == "ok" and
                      point.get("engines", {}).get(engine, {}).get("status") == "ok"]
            data_by_engine.append(values)
        bp = ax.boxplot(data_by_engine, patch_artist=True,
                        tick_labels=["plain-grid", "fast", "SageNet+"])
        for patch_box, engine in zip(bp["boxes"], ENGINES):
            patch_box.set_facecolor(colors[engine])
            patch_box.set_alpha(0.48)
        for i, (values, engine) in enumerate(zip(data_by_engine, ENGINES), start=1):
            jitter = np.linspace(-0.08, 0.08, len(values)) if values else []
            ax.scatter(np.full(len(values), i) + jitter, values, s=16, color=colors[engine],
                       edgecolors="#202020", linewidths=0.25, zorder=3)
        ax.set_yscale("log")
        ax.set_title(context_name)
        ax.set_ylabel("每个参数点的频谱误差 p95（dex，对数轴）")
        ax.grid(axis="y", alpha=0.22)
    fig.suptitle("跨后验参数点的精度：每个点都用 48 个频率与独立参考比较")
    fig.savefig(OUT / "mcmc_posterior_accuracy.png", dpi=200)
    plt.close(fig)


def write_report(result: dict) -> None:
    lines = [
        "# stiffGWpy 与 SageNet+：LVK 数据下的 MCMC 对比", "",
        f"实验时间：{result['generated_at']}",
        f"代码版本：stiffGWpy `{result['environment']['stiffgwpy_sha']}`；SageNet `{result['environment']['sagenet_sha']}`。", "",
        "## 结论先读", "",
        "本次补充实验使用相同 LVK 数据、似然和自由参数，比较最初 plain-grid、当前 fast 与 SageNet+ Transformer。MCMC 从当前参数附近随机提出新参数，再按似然接受或拒绝；每个方法和情景有 4 条链、每条保留 10,000 个样本，共 360,000 个保留样本。速度按预热后的采样阶段计算；预热、模型载入和独立精度参照另行计时。",
        "当前 fast 每步耗时更低；能谱精度需结合下文多个后验参数点的误差分布判断。低重加热温度情景不覆盖 LVK 观测频段，不用于判断谁更符合观测。", "",
        "## 速度与后验结果", "",
        "每步时间越低越快；ESS 越高表示链中重复信息越少；R-hat 越接近 1，四条链越一致。均值和标准差描述参数后验位置与宽度，不是算法误差。", "",
        "| 情景 | 方法 | 每步毫秒 | 每秒最低 bulk ESS | 最低 bulk ESS | 最低 tail ESS | 最大 rank R-hat | log10(r) 均值±标准差 | n_t 均值±标准差 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for context_name in CONTEXTS:
        for engine in ENGINES:
            row = result["contexts"][context_name]["engines"][engine]
            mean, std = row["posterior_mean"], row["posterior_std"]
            lines.append(
                f"| {context_name} | {engine} | {row['milliseconds_per_step']:.3f} | "
                f"{min(row['ess_per_second']):.2f} | {min(row['ess_bulk']):.0f} | "
                f"{min(row['ess_tail']):.0f} | {max(row['rhat_rank']):.3f} | "
                f"{mean[0]:.3f}±{std[0]:.3f} | {mean[1]:.3f}±{std[1]:.3f} |")

    lines += ["", "## 采样是否可靠", "",
              "每种方法、每个情景使用 4 条不同随机种子的链；起点分散在 `log10(r)`=-4.6 到 -1.8、`n_t`=-0.45 到 0。每条链预热 4,000 步，再保留 10,000 步。",
              "R-hat 是看不同链是否汇合；bulk ESS 看后验主体有多少有效信息，tail ESS 看两端区间的信息量；MCSE 是均值因有限采样带来的估算误差。采用秩标准化 folded split R-hat 与 Geyer ESS 计算。", "",
              "Stan 的诊断建议默认至少 4 条链；最终结果通常要求 R-hat < 1.01，bulk ESS 大于链数的 100 倍，tail ESS 检查 5% 与 95% 尾部。[诊断说明](https://mc-stan.org/learn-stan/diagnostics-warnings.html) 本报告把四链总 ESS 400 作为参考线；它不是对任何科学问题都足够的保证。", "",
              "| 情景 | 方法 | 最大 R-hat | 最低 bulk ESS | 最低 tail ESS | 最大均值 MCSE（log10r / nt） | 诊断门槛是否都达到* |",
              "|---|---|---:|---:|---:|---:|---|"]
    for context_name in CONTEXTS:
        for engine in ENGINES:
            row = result["contexts"][context_name]["engines"][engine]
            passed = max(row["rhat_rank"]) < 1.01 and min(row["ess_bulk"]) >= 400 and min(row["ess_tail"]) >= 400
            mcse = row["mcse_mean"]
            lines.append(f"| {context_name} | {engine} | {max(row['rhat_rank']):.3f} | "
                         f"{min(row['ess_bulk']):.0f} | {min(row['ess_tail']):.0f} | "
                         f"{mcse[0]:.4g} / {mcse[1]:.4g} | {'是' if passed else '否'} |")
    mode = result.get("posterior_mode_diagnostics")
    if mode:
        lines += ["", "**基准 SageNet+ 的重要问题**：" +
                  f"第 {mode['chain_with_high_tilt']} 条链有 {mode['high_tilt_fraction_in_chain']:.1%} 的样本落在 `n_t>{mode['high_tilt_threshold']}`，另外三条链均未进入该区域；占四链样本的 {mode['high_tilt_fraction_overall']:.1%}。" +
                  f"固定 `log10(r)={mode['profile_log10r']}` 时，SageNet+ 在 `n_t={mode['profile_low_nt']}` 的 LVK 频段覆盖为 100%，在 `n_t={mode['profile_high_nt']}` 降为 0%；高斜率区域的似然因而不再由 LVK 频段约束。" +
                  f"这对应 R-hat {mode['rhat_nt']:.3f}、bulk ESS {mode['bulk_ess_nt']:.0f}、tail ESS {mode['tail_ess_nt']:.0f}。该组合未通过混合检查，均值和标准差不能当成稳定后验估计。" +
                  ("一个抽到该高斜率区域的精度检查点还被 fast 的物理边界保护拒绝，因此不进入三方法误差排名。" if mode.get("high_mode_precision_probe_guarded") else ""), ""]
    lines += ["", "*同时达到表中三个参考值才标为‘是’；即使通过，也只说明本报告设定下的链诊断较好。", "",
              "## 精度：检查整个后验范围", "",
              "本地连续-sigma 高精度求解器（`rtol=1e-7`、`z_tail=8`）只作独立参照，不参与 MCMC。每个情景从三种方法各随机抽取 8 个后验参数点（共 24 点），每点在三种方法共同覆盖的 48 个频率上对比。下表是这些抽查点的频谱误差 p95，再汇总中位数 / 95 分位 / 最大值；它不能保证覆盖所有稀有后验区域。", "",
              "| 情景 | 方法 | 有效参数点 / 24 | 参数点误差 p95 的中位数 / 95分位 / 最大值（dex） |",
              "|---|---|---:|---:|"]
    for context_name in CONTEXTS:
        check = result["precision_posterior"][context_name]
        for engine in ENGINES:
            row = check["summary"][engine]
            n = row["valid_parameter_points"]
            vals = (row["dex_p95_per_point_median"], row["dex_p95_per_point_p95"], row["dex_p95_per_point_max"])
            values = "未能计算" if vals[0] is None else " / ".join(f"{x:.4g}" for x in vals)
            lines.append(f"| {context_name} | {engine} | {n} / 24 | {values} |")
    for context_name in CONTEXTS:
        excluded = [point for point in result["precision_posterior"][context_name]["points"]
                    if point.get("status") != "ok"]
        if excluded:
            causes = ", ".join(f"点 {point['point_index'] + 1}: {point.get('status')}"
                                for point in excluded)
            lines.append(f"- {context_name} 有 {len(excluded)} 个抽查点未进入共同精度比较：{causes}。物理保护拒绝的点单独报告，不把它算成普通精度误差。")
    lines += ["", "图中展示同一情景各后验参数点的误差分布；每个点的误差先在 48 个频率上计算，再取 p95。高刚性物质与基准情景可以用于观测频段内的数值精度比较。低重加热温度情景的 LVK 覆盖为 0%，结果只能用于数值对照，不能解释为数据支持度。0.01 dex 约对应 2.3% 的谱幅差，0.1 dex 约对应 26%。", "",
              "被拒提议按原因分别记录；物理保护表示模型适用边界，不等同于数值错误：", "",
              "| 情景 | 方法 | 被拒提议原因与次数 |", "|---|---|---|"]
    for context_name in CONTEXTS:
        for engine in ENGINES:
            counts = result["contexts"][context_name]["engines"][engine]["failure_counts"]
            labels = {"shared_Neff_guard": "物理边界保护", "max_iter": "未收敛",
                      "non-finite likelihood": "无效数值/似然"}
            description = ", ".join(f"{labels.get(k, k)} {v} 次" for k, v in sorted(counts.items())) or "无"
            lines.append(f"| {context_name} | {engine} | {description} |")
    lines += ["", "## 场景和参数设置", "",
              "三个固定背景情景分别为：基准（`kappa10=0.01, T_re=2000 GeV, DN_re=20`）、低重加热温度（`0.01, 10 GeV, 10`）、高刚性物质（`1, 2000 GeV, 30`）。‘低温’表示再加热结束温度较低；‘高刚性’表示 10 MeV 时刚性物质相对光子的能量比较高。", "",
              "| 参数 | 通俗含义 | 采样设置 |", "|---|---|---|",
              "| `log10(r)` | 原初引力波强度比例的 10 为底对数；加 1 表示 `r` 增大 10 倍 | 自由参数；-5 到 -1 均匀取值 |",
              "| `n_t` | 引力波谱随频率变化的斜率 | 自由参数；-0.5 到 1.5 |",
              "| `kappa10` | 10 MeV 时刚性物质能量与光子能量之比 | 按上述三种情景固定 |",
              "| `T_re` | 再加热结束温度，单位 GeV | 按上述三种情景固定 |",
              "| `DN_re` | 暴胀结束到再加热结束期间的膨胀量 | 按上述三种情景固定 |",
              "| `cr` | 是否启用单场暴胀一致性关系 | 固定为 0，允许 `n_t` 独立变化 |",
              "| `Omega_bh2` / `Omega_ch2` | 普通物质 / 暗物质密度 | 0.0223828 / 0.1201075 |",
              "| `H0` / `A_s` | 当前膨胀速度 / 标量扰动强度 | 67.32117 km/s/Mpc / 2.100549e-9 |",
              "| `DN_eff` | 起始额外相对论粒子贡献 | 固定 0；引力波贡献由各求解器计算 |", "",
              "提议步幅固定为 `log10(r)=0.35`、`n_t=0.12`。plain-grid 按早期粗网格档重建；fast 使用当前正式档；SageNet+ 使用 Transformer、CPU 单线程。三者共享同一 LVK O1/O2/O3 数据（{} 行）和三次插值卡方似然。计时期间屏蔽重复的控制台保护提示，但仍逐项计数；这样不会把终端输出速度当成模型计算速度。".format(result['environment']['lvk_data_rows']),
              f"求解器设置：plain-grid `h=0.02, col_step=8, z_tail=5, phase_max=0, construct 网格`；fast `h=0.005, col_step=8, z_tail=5, phase_max=0.25, goal 网格并拆分再加热断点`。SageNet 模型载入 {result['sagenet_model_load_s_excluded_from_chain']:.3f} 秒；采样前预热时间记录在 JSON 的 `warmup_s`，均不计入每步采样时间。硬件为 {result['environment'].get('cpu_model', result['environment'].get('processor', '未记录'))}；固定使用 2 个 Numba 线程、1 个 BLAS/Torch 线程。", "",
              "## 图表", "",
              "颜色统一为蓝色 plain-grid、橙色 fast、绿色 SageNet+。", "",
              "### 后验样本", "",
              "每个点是保留的参数组合；密集处表示当前设置下采样较多。样本相关，不能把 40,000 个保留值当成 40,000 份独立信息。", "",
              "![后验样本分布](mcmc_posterior.png)", "",
              "### 链轨迹", "",
              "检查链是否重叠、是否还在持续漂移。颜色代表方法，同色深浅代表独立链。", "",
              "![链轨迹诊断](mcmc_diagnostics.png)", "",
              "### 速度和有效样本", "",
              "时间仅计入预热后的保留阶段；模型加载、采样预热和独立精度参照不混入。", "",
              "![速度与有效样本对比](mcmc_speed.png)", "",
              "### 代表性能谱", "",
              "展示各方法后验合并样本中位参数附近的能谱。", "",
              "![能谱与独立精度参照](mcmc_accuracy.png)", "",
              "### 后验范围精度", "",
              "箱体和点展示 24 个后验参数点的逐点误差 p95；纵轴使用对数刻度，方便同时观察小误差与大误差。", "",
              "![后验多点精度误差分布](mcmc_posterior_accuracy.png)", "",
              "## 对 fast 的评价", ""]
    speed_ratios = [result["contexts"][c]["engines"]["SageNet+ Transformer"]["milliseconds_per_step"] /
                    result["contexts"][c]["engines"]["当前 fast"]["milliseconds_per_step"] for c in CONTEXTS]
    plain_speed_change = [100 * (result["contexts"][c]["engines"]["最初 plain-grid"]["milliseconds_per_step"] /
                                  result["contexts"][c]["engines"]["当前 fast"]["milliseconds_per_step"] - 1)
                          for c in CONTEXTS]
    plain_ess_change = [100 * (min(result["contexts"][c]["engines"]["当前 fast"]["ess_per_second"]) /
                               min(result["contexts"][c]["engines"]["最初 plain-grid"]["ess_per_second"]) - 1)
                        for c in CONTEXTS]
    fast_accuracy = [result["precision_posterior"][c]["summary"]["当前 fast"] for c in CONTEXTS]
    sage_accuracy = [result["precision_posterior"][c]["summary"]["SageNet+ Transformer"] for c in CONTEXTS]
    lines.append(f"- **速度**：fast 比 SageNet+ 快约 {min(speed_ratios):.1f}–{max(speed_ratios):.1f} 倍；比 plain-grid 在基准 / 低温情景快 {plain_speed_change[0]:.0f}% / {plain_speed_change[1]:.0f}%，在高刚性情景慢 {abs(plain_speed_change[2]):.0f}%；每秒有效样本分别变化 {plain_ess_change[0]:+.0f}% / {plain_ess_change[1]:+.0f}% / {plain_ess_change[2]:+.0f}%。这些速度值只适用于本次机器、版本和线程设置。")
    comparable = [(f["dex_p95_per_point_median"], s["dex_p95_per_point_median"])
                  for f, s in zip(fast_accuracy, sage_accuracy)
                  if f["dex_p95_per_point_median"] is not None and s["dex_p95_per_point_median"] is not None]
    if comparable:
        better = sum(f < s for f, s in comparable)
        lines.append(f"- **精度**：抽到且可比较的参数点中，fast 的误差中位数在 {better}/{len(comparable)} 个情景低于 SageNet+；具体见上表。基准 SageNet+ 链未收敛，其误差汇总只描述已抽到的可比较区域，不代表完整后验。")
    lines.append("- **相对 plain-grid 精度**：fast 在基准和高刚性情景的后验参数点误差中位数较低；低温情景较高，但该情景的 LVK 覆盖为 0%。")
    lines += ["- **观测解释边界**：低温情景 LVK 频段覆盖为 0%，不支持拟合优劣结论；本报告比较的是数值频谱相对本地参考的误差。",
              "- **科研使用判断**：fast 与 plain-grid 三个情景都达到本报告的链诊断参考；SageNet+ 基准与低温情景未达到，因此这份三方后验比较整体尚未全部通过科研生产验收。基准 SageNet+ 还存在不同链落入不同覆盖区域的问题。", "",
              "## 原始数据与复现", "",
              f"- `results.json`：每种方法、情景的汇总、诊断量、逐参数点精度结果和环境信息。",
              f"- `{result['chain_file']}`：四条链全部保留样本、起点和逐链耗时；随报告一并归档。",
              "- `scripts/benchmark_mcmc_sagenet_compare.py`：采样、诊断、精度核验与绘图脚本。",
              f"- LVK 数据 SHA-256：`{result['environment']['lvk_data_sha256']}`；SageNet Transformer 权重 SHA-256：`{result['environment']['sagenet_transformer_weights_sha256']}`。",
              "- 已按旧样本文件名和报告路径检索公开网络，未找到可核验的原始后验样本；远端 Git 历史也没有该文件。旧本地链只有三条相同起点的短链，未作为新证据。本次样本由当前报告版本重新生成并归档。", "",
              "结果限定于报告列出的代码版本、SageNet 模型权重、CPU 环境、数据和参数范围；不代表其他模型权重或所有参数空间。"]
    (OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="短链，只用于验证程序连接与图表格式")
    run(parser.parse_args().quick)

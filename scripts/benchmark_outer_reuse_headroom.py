"""测量正式 fast 路径的 outer full-solve reuse 余量（只读 A/B，不改默认行为）。

`_OUTER_FULL_REUSE_SIGMA_TOL`/`_OUTER_FULL_REUSE_FHOR_TOL = 1e-4` 目前只让低
DN_eff 更新点跳过第二次完整 tensor solve；high-T/stiff/high-kappa 因 DN_eff
更新更大而必须跑满两次，warm runtime 约为前者的两倍（见
`docs/benchmark_head_matrix.json` 与固定 profiler 的 kernel 调用次数）。本脚本对
同一代表点**交替**测量「现状判据」与「强制 reuse」，记录 DN/频谱偏差、runtime
比值、每轮背景最大变化与 kernel 调用次数，用来判断阈值是否过紧、能否结构化省掉
第二次求解。

所有被替换的函数与阈值都在 finally 中恢复，solver 默认行为不变。
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

REUSE_TOL = 1.0e-4
FORCED_TOL = 1.0e9


class Trace:
    """只读插桩：记录每轮背景与一次求解中的 kernel 调用次数。"""

    def __init__(self):
        self.sigma = []
        self.f_hor = []
        self.kernel_calls = 0
        self._saved = None

    def reset(self):
        self.sigma = []
        self.f_hor = []
        self.kernel_calls = 0

    def install(self):
        saved = (FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel)
        self._saved = saved

        def gen_fast(*args, **kwargs):
            out = saved[0](*args, **kwargs)
            sigma = getattr(args[0], 'sigma', None)
            if sigma is not None:
                self.sigma.append(np.asarray(sigma, dtype=float).copy())
            return out

        def prep(*args, **kwargs):
            out = saved[1](*args, **kwargs)
            self.f_hor.append(np.asarray(out[1], dtype=float).copy())
            return out

        def solve(*args, **kwargs):
            self.kernel_calls += 1
            return saved[2](*args, **kwargs)

        FS.gen_fast = gen_fast
        FS.prep_frequency_only = prep
        FS.solve_kernel = solve

    def restore(self):
        if self._saved is not None:
            FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel = self._saved
            self._saved = None


def _step_sizes(traces):
    """逐轮最大绝对变化；形状不一致的轮次跳过（视为不可比）。"""
    out = []
    for i in range(1, len(traces)):
        a, b = traces[i], traces[i - 1]
        if a.shape == b.shape:
            out.append(float(np.max(np.abs(a - b))))
    return out


def _relative(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def _solve(point, forced, trace):
    """按模式设置 reuse 阈值后跑一次正式 fast 路径。"""
    if forced:
        FS._OUTER_FULL_REUSE_SIGMA_TOL = FORCED_TOL
        FS._OUTER_FULL_REUSE_FHOR_TOL = FORCED_TOL
    else:
        FS._OUTER_FULL_REUSE_SIGMA_TOL = REUSE_TOL
        FS._OUTER_FULL_REUSE_FHOR_TOL = REUSE_TOL
    model = LCDM_SG(**CASES[point])
    trace.reset()
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
    dn = float(np.asarray(model.DN_gw, dtype=float).reshape(-1)[-1])
    spec = np.asarray(model.log10OmegaGW, dtype=float).copy()
    return model, dn, spec


def _trace_row(trace):
    return {'kernel_calls': trace.kernel_calls,
            'sigma_steps': _step_sizes(trace.sigma),
            'f_hor_steps': _step_sizes(trace.f_hor)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default', 'lowT', 'highT', 'stiff', 'low_r',
                                 'high_kappa'])
    parser.add_argument('--warmup', type=int, default=2)
    parser.add_argument('--repeats', type=int, default=15)
    parser.add_argument('--out',
                        default='docs/benchmark_outer_reuse_headroom_head.json')
    args = parser.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    trace = Trace()
    trace.install()
    rows = []
    try:
        for point in args.points:
            _, dn_current, spec_current = _solve(point, False, trace)
            current_trace = _trace_row(trace)
            _, dn_forced, spec_forced = _solve(point, True, trace)
            forced_trace = _trace_row(trace)
            for mode in (False, True):
                for _ in range(args.warmup):
                    _solve(point, mode, trace)
            times = {'current': [], 'forced': []}
            for repeat in range(args.repeats):
                order = (False, True) if repeat % 2 == 0 else (True, False)
                for mode in order:
                    start = time.perf_counter()
                    _solve(point, mode, trace)
                    key = 'forced' if mode else 'current'
                    times[key].append((time.perf_counter() - start) * 1e3)
            stats = {}
            for key, values in times.items():
                arr = np.asarray(values)
                stats[key] = {'warm_median_ms': float(np.median(arr)),
                              'warm_p95_ms': float(np.percentile(arr, 95)),
                              'warm_min_ms': float(arr.min()),
                              'n_repeats': int(arr.size)}
            rows.append({
                'point': point,
                'DN_gw_current': dn_current,
                'DN_gw_forced': dn_forced,
                'DN_rel_forced_vs_current': _relative(dn_forced, dn_current),
                'spectrum_max_abs_dex_diff': float(
                    np.max(np.abs(spec_forced - spec_current))),
                'current': dict(stats['current'], **current_trace),
                'forced': dict(stats['forced'], **forced_trace),
                'runtime_ratio_forced_over_current': (
                    stats['forced']['warm_median_ms']
                    / stats['current']['warm_median_ms']),
                'times_ms': times,
            })
    finally:
        trace.restore()
        FS._OUTER_FULL_REUSE_SIGMA_TOL = REUSE_TOL
        FS._OUTER_FULL_REUSE_FHOR_TOL = REUSE_TOL
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(workers=1, threads=2),
        'warmup': args.warmup,
        'repeats': args.repeats,
        'reuse_tol_current': REUSE_TOL,
        'reuse_tol_forced': FORCED_TOL,
        'rows': rows,
        'semantics': ('fast SGWB_iter_fast wall time, self-consistent DN_gw, '
                      'spectrum and per-iteration background change for the '
                      'current reuse gate versus a forced-reuse gate on the '
                      'same native grid'),
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()

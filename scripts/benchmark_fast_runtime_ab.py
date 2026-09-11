"""固定资源的 fast warm runtime 与逐位输出 digest A/B 计时。"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time

import numpy as np

ROOT = os.environ.get('STIFFGWPY_ROOT',
                      os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}

DIGEST_FIELDS = ('f', 'log10OmegaGW', 'DN_gw', 'g2', 'w2')


def p95(values):
    return float(np.percentile(np.asarray(values), 95))


def digest(model, name):
    arr = np.ascontiguousarray(np.asarray(getattr(model, name), dtype=np.float64))
    return hashlib.sha256(arr.tobytes()).hexdigest()


def run_once(case):
    model = LCDM_SG(**case)
    start = time.perf_counter()
    FS.SGWB_iter_fast(model, kink_split=True)
    elapsed = time.perf_counter() - start
    row = {'seconds': elapsed,
           'DN_eff': float(model.cosmo_param['DN_eff']),
           'DN_gw_last': float(np.asarray(model.DN_gw)[-1]),
           'n_freq': int(len(model.f)),
           'reason': getattr(model, 'fast_failure_reason', None),
           'converged': bool(getattr(model, 'SGWB_converge', False))}
    for name in DIGEST_FIELDS:
        row['digest_' + name] = digest(model, name)
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--cases', default='all')
    parser.add_argument('--json', default=None)
    parser.add_argument('--affinity', type=int, default=None,
                        help='限制进程到前 N 个逻辑 CPU，降低移动 CPU 的频率噪声')
    args = parser.parse_args()
    if args.affinity:
        import psutil
        psutil.Process().cpu_affinity(list(range(args.affinity)))
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    names = sorted(CASES) if args.cases == 'all' else args.cases.split(',')
    summary = {'threads': args.threads, 'repeats': args.repeats,
               'affinity': args.affinity, 'points': {}}
    for name in names:
        # 首次调用含 JIT 编译与首次分配，不计入 warm 统计。
        run_once(CASES[name])
        rows = [run_once(CASES[name]) for _ in range(args.repeats)]
        times = [row['seconds'] for row in rows]
        point = {'median_s': statistics.median(times),
                 'p95_s': p95(times),
                 'min_s': min(times)}
        point.update({k: v for k, v in rows[-1].items() if k != 'seconds'})
        summary['points'][name] = point
        print('%-11s median %7.3f ms  p95 %7.3f ms  DN_gw %.6e' % (
            name, point['median_s'] * 1e3, point['p95_s'] * 1e3, point['DN_gw_last']))
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as handle:
            json.dump(summary, handle, ensure_ascii=False, indent=2)
    digests = {name: {k: v for k, v in point.items() if k.startswith('digest')}
               for name, point in summary['points'].items()}
    print(json.dumps(digests, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

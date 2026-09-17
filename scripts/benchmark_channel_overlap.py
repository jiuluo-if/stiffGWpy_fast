"""Read-only diagnostic for exact cross-channel propagation reuse potential."""
from __future__ import annotations

import json
import os
import sys

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def run_case(name):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    model = LCDM_SG(**CASES[name])
    FS.SGWB_iter_fast(model, kink_split=True)
    nv = model.Nv.astype(np.float64)
    _, _, j0s, z0s, _ = FS.prep_frequency_only(model, nv, model.f)
    phi_grid = EB.fast_phi_s2_split(
        model, nv, model.cosmo_param['DN_eff'], sigma_nodes=model.sigma)[0]
    tail_idx = np.asarray([
        np.searchsorted(phi_grid, phi_grid[j0] + 5.0 - z0,
                        side='left')
        for j0, z0 in zip(j0s, z0s)
    ], dtype=np.int64)
    unique, counts = np.unique(j0s, return_counts=True)
    unique_tail, tail_counts = np.unique(tail_idx, return_counts=True)
    pairs = np.stack((j0s, tail_idx), axis=1)
    unique_pairs, pair_counts = np.unique(pairs, axis=0, return_counts=True)
    return {
        'case': name,
        'n_modes': int(j0s.size),
        'unique_j0': int(unique.size),
        'j0_reuse_fraction': float(1.0 - unique.size / max(j0s.size, 1)),
        'max_j0_multiplicity': int(counts.max()),
        'unique_tail_index': int(unique_tail.size),
        'tail_reuse_fraction': float(1.0 - unique_tail.size / max(j0s.size, 1)),
        'max_tail_multiplicity': int(tail_counts.max()),
        'unique_j0_tail_pairs': int(unique_pairs.shape[0]),
        'pair_reuse_fraction': float(1.0 - unique_pairs.shape[0] /
                                     max(j0s.size, 1)),
        'max_pair_multiplicity': int(pair_counts.max()),
        'mean_steps': float(np.mean(tail_idx - j0s)),
        'median_steps': float(np.median(tail_idx - j0s)),
        'j0_digest': int(np.sum(j0s.astype(np.int64))),
    }


def main():
    rows = [run_case(name) for name in ('default', 'highT', 'stiff', 'high_kappa', 'lowT')]
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'channel_overlap_read_only_diagnostic',
        'production_path_changed': False,
        'resources': telemetry(workers=1, threads=2),
        'rows': rows,
    }
    out = 'docs/channel_overlap_round17_20260917.json'
    with open(out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()

from __future__ import annotations

import json
import math
import statistics
import sys
import timeit
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from stiffgwpy_fast.fast_sgwb import (  # noqa: E402
    _phase_segment_base,
    _phase_segment_wkb3,
    _phase_substeps,
)

H = 0.00125
Q = 1.0
PHASE_MAX = 0.25
TOLERANCE = 1.0e-5
LOOPS = 200
REPEATS = 7
rows = []

for target in (16, 32, 64, 128, 256, 512, 1024):
    z_mid = math.log(target * PHASE_MAX / H)
    z_start = z_mid - 0.5 * Q * H
    z_end = z_mid + 0.5 * Q * H
    actual = _phase_substeps(H, z_mid, PHASE_MAX)
    accepted = _phase_segment_wkb3(
        0.0, 1.0, z_start, z_end, H, PHASE_MAX, TOLERANCE, 0
    )[2]
    if not accepted:
        raise RuntimeError(f"WKB rejected the break-even probe at {actual} steps")

    wkb_samples = timeit.repeat(
        lambda: _phase_segment_wkb3(
            0.0, 1.0, z_start, z_end, H, PHASE_MAX, TOLERANCE, 0
        ),
        number=LOOPS,
        repeat=REPEATS,
    )
    phase_samples = timeit.repeat(
        lambda: _phase_segment_base(0.0, 1.0, z_start, z_end, H, PHASE_MAX),
        number=LOOPS,
        repeat=REPEATS,
    )
    wkb_s = statistics.median(wkb_samples) / LOOPS
    phase_s = statistics.median(phase_samples) / LOOPS
    rows.append(
        {
            "target_substeps": target,
            "actual_substeps": int(actual),
            "accepted": bool(accepted),
            "wkb_us": wkb_s * 1.0e6,
            "phase_us": phase_s * 1.0e6,
            "phase_over_wkb": phase_s / wkb_s,
        }
    )
    print(rows[-1])

output = {
    "status": "complete",
    "source": "constant-q positive-frequency segment, q=100, h=0.00125",
    "phase_max": PHASE_MAX,
    "wkb_defect_tolerance": TOLERANCE,
    "resources": {"loops": LOOPS, "repeats": REPEATS},
    "rows": rows,
}
path = ROOT / "docs/research/unified_domain_20261119/benchmark_wkb3_break_even.json"
path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
print(path)

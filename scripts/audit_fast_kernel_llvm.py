"""Inspect canonical Numba LLVM/ASM for the fast propagation hot spots."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, limit_affinity, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, limit_affinity, telemetry

apply_environment()

import psutil  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import _make_args, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402

TOKENS = {
    "exp": (r"\bexp\b", r"llvm\.exp", r"__svml_exp"),
    "sin": (r"\bsin\b", r"llvm\.sin", r"__svml_sin"),
    "cos": (r"\bcos\b", r"llvm\.cos", r"__svml_cos"),
    "sqrt": (r"\bsqrt\b", r"llvm\.sqrt", r"__svml_sqrt"),
    "ceil": (r"\bceil\b", r"llvm\.ceil", r"__svml_ceil"),
    "integer_div": (r"\b[su]?div\b", r"\bdiv\b"),
    "integer_mod": (r"\b[su]?rem\b", r"\brem\b"),
    "alloca": (r"\balloca\b",),
    "stack_spill": (r"\bspill\b", r"\brsp\b", r"\bymm\d+\b"),
}


def _count(text: str, patterns: tuple[str, ...]) -> int:
    return sum(len(re.findall(pattern, text, flags=re.IGNORECASE)) for pattern in patterns)


def _excerpt(text: str) -> list[str]:
    lines = text.splitlines()
    interesting = ("exp", "sin", "cos", "sqrt", "ceil", "div", "rem", "call", "alloca", "spill")
    return [line[:240] for line in lines if any(token in line.lower() for token in interesting)][:80]


def _inspect(fn, name: str) -> dict:
    if not fn.signatures:
        return {"name": name, "error": "no compiled signature"}
    signature = fn.signatures[-1]
    # Numba's on-disk cache deliberately disables inspection for cache-loaded
    # code. Recompile the already-specialized signature in this process so the
    # returned LLVM/ASM is real rather than an empty placeholder.
    recompiled = False
    recompile_error = None
    try:
        fn.recompile()
        recompiled = True
        signature = fn.signatures[-1]
    except Exception as exc:  # pragma: no cover - version dependent
        recompile_error = repr(exc)
    try:
        llvm = fn.inspect_llvm(signature)
        llvm_error = None
    except Exception as exc:  # pragma: no cover - version dependent
        llvm = ""
        llvm_error = repr(exc)
    try:
        asm = fn.inspect_asm(signature)
        asm_error = None
    except Exception as exc:  # pragma: no cover - version dependent
        asm = ""
        asm_error = repr(exc)
    return {
        "name": name,
        "signature": repr(signature),
        "recompiled": recompiled,
        "recompile_error": recompile_error,
        "llvm_error": llvm_error,
        "asm_error": asm_error,
        "llvm_sha256": hashlib.sha256(llvm.encode()).hexdigest(),
        "asm_sha256": hashlib.sha256(asm.encode()).hexdigest(),
        "llvm_bytes": len(llvm),
        "asm_bytes": len(asm),
        "llvm_counts": {key: _count(llvm, pats) for key, pats in TOKENS.items()},
        "asm_counts": {key: _count(asm, pats) for key, pats in TOKENS.items()},
        "llvm_has_fast_math": bool(
            re.search(r"\bfast(?:\s+(?:nnan|ninf|nsz|arcp|contract|afn|reassoc))+", llvm, flags=re.IGNORECASE)
        ),
        "asm_excerpt": _excerpt(asm),
        "llvm_excerpt": _excerpt(llvm),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--affinity-count", type=int, default=None)
    parser.add_argument("--output", default="docs/fast_kernel_llvm_audit_round28_20260923.json")
    args = parser.parse_args()
    process = psutil.Process()
    affinity_count = args.affinity_count or args.threads
    selected = limit_affinity(process, affinity_count)
    os.environ["FAST_THREADS"] = str(args.threads)
    FS.apply_accuracy_mode("fast")
    FS.set_threads(args.threads)

    model, common = _prepared("default", args.threads)
    solve_args = _make_args(common)
    FS.solve_kernel(*solve_args)
    FS.scaled_step(0.0, 1.0, 0.1, 0.005)
    FS._phase_substeps(0.005, 0.1, 0.25)
    FS._phase_segment(0.0, 1.0, -1.0, 0.0, 0.005, 0.25)
    payload = {
        "commit": os.popen("git rev-parse HEAD").read().strip(),
        "profile": "fast",
        "threads": args.threads,
        "affinity_selected": selected,
        "resources": telemetry(process, threads=args.threads),
        "n_freq": int(model.f.size),
        "n_background": int(model.Nv.size),
        "kernels": [
            _inspect(FS.scaled_step, "scaled_step"),
            _inspect(FS._phase_substeps, "_phase_substeps"),
            _inspect(FS._phase_segment, "_phase_segment"),
            _inspect(FS.solve_kernel, "solve_kernel"),
        ],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

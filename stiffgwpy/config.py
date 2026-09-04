"""Immutable configuration records for the fast SGWB solver.

The numerical kernels still expose legacy module setters for compatibility,
but a solver invocation should consume a snapshot rather than a mutable module
namespace. Keeping this record dependency-light avoids importing the solver
module from configuration code and makes it safe to pass through adapters.
"""

from dataclasses import dataclass
from math import isfinite
from typing import Optional


@dataclass(frozen=True)
class FastSolverConfig:
    """Validated per-call settings for the fast numerical kernels."""

    h: float = 0.01
    col_step: int = 4
    z_tail: float = 8.0
    phase_max: float = 0.5
    freq_grid: str = "construct"
    threads: Optional[int] = None

    def __post_init__(self) -> None:
        raw_col_step = self.col_step
        raw_threads = self.threads
        h = float(self.h)
        z_tail = float(self.z_tail)
        phase_max = float(self.phase_max)
        col_step = int(self.col_step)
        threads = None if self.threads is None else int(self.threads)
        if not isfinite(h) or not 1e-4 <= h <= 0.1:
            raise ValueError("h must be finite and in [1e-4, 0.1], got %r" % self.h)
        if isinstance(raw_col_step, bool) or float(raw_col_step) != col_step:
            raise ValueError("col_step must be an integer in [1, 8], got %r" % raw_col_step)
        if not 1 <= col_step <= 8:
            raise ValueError("col_step must be an integer in [1, 8], got %r" % self.col_step)
        if not isfinite(z_tail) or not 2.0 <= z_tail <= 15.0:
            raise ValueError("z_tail must be finite and in [2.0, 15.0], got %r" % self.z_tail)
        if not isfinite(phase_max) or not 0.0 <= phase_max <= 10.0:
            raise ValueError("phase_max must be finite and in [0, 10], got %r" % self.phase_max)
        if self.freq_grid not in ("construct", "grid_independent", "adaptive"):
            raise ValueError("freq_grid must be construct/grid_independent/adaptive, got %r" % self.freq_grid)
        if (raw_threads is not None and
                (isinstance(raw_threads, bool) or float(raw_threads) != threads)):
            raise ValueError("threads must be a positive integer, got %r" % raw_threads)
        if threads is not None and threads < 1:
            raise ValueError("threads must be a positive integer, got %r" % self.threads)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "z_tail", z_tail)
        object.__setattr__(self, "phase_max", phase_max)
        object.__setattr__(self, "col_step", col_step)
        object.__setattr__(self, "threads", threads)

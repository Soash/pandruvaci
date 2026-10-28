#!/usr/bin/env python3
"""
Hardware Detection & Multi-Core Optimization Module
===================================================
Provides robust hardware inspection, auto-detection of available physical
and logical CPU cores (respecting OS process affinity / cgroups), and dynamic
thread allocation for high-performance bioinformatics pipelines.
"""

import os
import logging
from typing import Any, Dict, Optional, Tuple

logger = logging.getLogger("hardware_detector")


def detect_available_cores() -> int:
    """
    Detects the maximum number of CPU cores available to the current process.
    Prioritizes OS scheduler affinity (os.sched_getaffinity) to respect container
    or taskset boundaries, falling back to os.cpu_count() or 4.
    """
    try:
        if hasattr(os, "sched_getaffinity"):
            cores = len(os.sched_getaffinity(0))
            if cores > 0:
                return cores
    except Exception:
        pass

    try:
        cores = os.cpu_count()
        if cores and cores > 0:
            return cores
    except Exception:
        pass

    return 4


def parse_percentage_ratio(val: Any) -> Optional[float]:
    """
    Parses a percentage string (e.g. "50%", "75.5%"), a float ratio (e.g. 0.5 for 50%),
    or string float into a normalized fraction (0.0 < ratio <= 1.0).
    Returns None if the value is not a percentage or ratio representation.
    """
    if val is None:
        return None
    if isinstance(val, float) and 0.0 < val <= 1.0:
        return val
    if isinstance(val, str):
        v = val.strip()
        if v.endswith("%"):
            try:
                pct = float(v[:-1].strip())
                if pct > 0:
                    return min(1.0, pct / 100.0)
            except ValueError:
                pass
        else:
            try:
                f_val = float(v)
                if 0.0 < f_val <= 1.0:
                    return f_val
            except ValueError:
                pass
    return None


def resolve_threads(
    step_threads: Any = "auto",
    config: Optional[Dict[str, Any]] = None,
    default: int = 4,
    task_name: str = "task",
) -> int:
    """
    Resolves the exact number of threads to utilize for a pipeline step based
    on system capability and user configuration.

    Resolution Hierarchy:
    1. Detect available physical/logical cores (e.g. 32).
    2. Check global 'hardware' configuration in config.yaml:
       - 'auto_detect_cores': If true, enables auto-tuning.
       - 'reserve_cores': Number of cores reserved for OS responsiveness (default 0).
       - 'max_threads': Global upper bound cap (if set).
       - 'threads_percent' / 'cpu_percent': Global percentage of CPU cores (e.g. 50, 75, "50%").
    3. If step_threads is a percentage (e.g. "50%", "75%", 0.5):
       Allocates round(usable_cores * percentage).
    4. If step_threads is "auto", "max", 0, or None:
       Inherits global threads or threads_percent (or allocates 100% of usable cores).
    5. If step_threads is an explicit integer:
       Respects the user's manual override (clamped to max_threads if configured).
    """
    total_cores = detect_available_cores()
    hw_cfg = config.get("hardware", {}) if config else {}

    auto_detect = hw_cfg.get("auto_detect_cores", True)
    reserve_cores = int(hw_cfg.get("reserve_cores", 0))
    max_cap = (config.get("max_threads") if config else None) or hw_cfg.get("max_threads")
    # Support both config['hardware']['threads'] and top-level CLI flag: --config threads="50%"
    global_threads = (config.get("threads") if config else None) or hw_cfg.get("threads", "auto")

    # If step_threads is not specified, None, or "auto", inherit global threads
    if step_threads is None or (isinstance(step_threads, str) and step_threads.strip().lower() == "auto"):
        step_threads = global_threads

    usable_cores = max(1, total_cores - reserve_cores)

    # 1. Check if step_threads specifies a percentage directly (e.g. "50%", "75%", 0.5)
    ratio = parse_percentage_ratio(step_threads)

    # 2. If not in step_threads, check if step_threads is auto/None and global has a percentage
    if ratio is None and (
        isinstance(step_threads, str) and step_threads.strip().lower() in ("auto", "max", "all")
        or step_threads == 0
        or step_threads == global_threads
    ):
        ratio = parse_percentage_ratio(global_threads)
        if ratio is None:
            # Check both root-level CLI flag and hardware config block
            raw_pct = (
                (config.get("threads_percent") or config.get("cpu_percent") or config.get("threads_ratio"))
                if config else None
            ) or hw_cfg.get("threads_percent") or hw_cfg.get("cpu_percent") or hw_cfg.get("threads_ratio")
            if raw_pct is not None:
                if isinstance(raw_pct, (int, float)) and raw_pct > 1.0:
                    ratio = min(1.0, float(raw_pct) / 100.0)
                else:
                    ratio = parse_percentage_ratio(raw_pct)

    if ratio is not None:
        allocated = max(1, round(usable_cores * ratio))
    else:
        is_auto = False
        if isinstance(step_threads, str):
            if step_threads.strip().lower() in ("auto", "max", "all"):
                is_auto = True
        elif step_threads == 0:
            is_auto = True

        if is_auto or (auto_detect and step_threads in (None, "auto")):
            allocated = usable_cores
        else:
            try:
                allocated = int(step_threads)
                if allocated <= 0:
                    allocated = usable_cores
            except (ValueError, TypeError):
                allocated = default

    # Apply global ceiling if defined
    if max_cap is not None:
        try:
            cap = int(max_cap)
            if cap > 0:
                allocated = min(allocated, cap)
        except (ValueError, TypeError):
            pass

    allocated = max(1, allocated)
    logger.debug(
        "Resolved %d threads for '%s' (System cores: %d, Configured: %s)",
        allocated,
        task_name,
        total_cores,
        step_threads,
    )
    return allocated


def optimize_parallel_workers_and_threads(
    total_items: int,
    step_threads: Any = "auto",
    config: Optional[Dict[str, Any]] = None,
    preferred_threads_per_worker: int = 4,
) -> Tuple[int, int]:
    """
    Optimizes the trade-off between parallel workers (processes/threads)
    and internal tool threads (e.g. MMseqs2 threads per worker) such that:
        workers * threads_per_worker <= total_budget

    Returns:
        (workers, threads_per_worker)
    """
    budget = resolve_threads(step_threads, config=config, default=4, task_name="parallel_pool")
    
    if total_items <= 1:
        return 1, budget

    # Choose threads per worker
    threads_per_worker = min(budget, preferred_threads_per_worker)
    workers = max(1, budget // max(1, threads_per_worker))
    
    # Don't spawn more workers than total items
    workers = min(workers, total_items)
    
    # Give remaining budget to threads per worker if items are fewer than workers
    threads_per_worker = max(1, budget // workers)

    return workers, threads_per_worker

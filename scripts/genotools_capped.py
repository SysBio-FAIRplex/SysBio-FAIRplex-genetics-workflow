#!/usr/bin/env python3
"""Run GenoTools with its worker pool sized to the SLURM allocation instead of the node.

WHY THIS EXISTS
---------------
genotools/ancestry.py:513-518 sizes the GridSearchCV worker pool from the machine:

    available_ram_gb = psutil.virtual_memory().total / (1024**3)
    max_workers_by_ram = ...
    n_jobs = min(os.cpu_count(), max_workers_by_ram)

Both inputs report NODE totals and ignore the cgroup. That is correct on a dedicated VM or
container — how GenoTools is normally run — and wrong on a shared SLURM node, where the job
owns 64 of 192 CPUs and 256 of 755 GB. n_jobs is then passed explicitly to GridSearchCV, so
LOKY_MAX_CPU_COUNT does not apply, and os.cpu_count() ignores CPU affinity on Python 3.11,
so taskset does not either. Patching the two readings is the only lever from outside.

WHAT ACTUALLY BROKE
-------------------
biowulf sets RLIMIT_NPROC (`ulimit -u`) to 1024, per-user and node-wide, shared across every
job that user has on the node. 192 loky workers — each a full interpreter with numba and
llvmlite loaded — do not fit in that budget. pthread_create returns EAGAIN, the C++ layer
throws an uncaught std::runtime_error, the worker aborts on SIGABRT, joblib raises
TerminatedWorkerError, and the parent exits 1. SLURM records a bare `ExitCode 1:0` with no
OUT_OF_MEMORY and no MaxRSS, which is why this read as anything but a resource limit.

Note what is NOT the fix: joblib already caps threads INSIDE each worker. Its
MAX_NUM_THREADS_VARS list includes NUMBA_NUM_THREADS, and the per-worker default is
max(cpu_count() // n_jobs, 1) — which is 1 when n_jobs equals the core count. Setting those
env vars by hand changes nothing joblib was not already doing. The worker COUNT is the lever.

UPSTREAM
--------
The real fix belongs in ancestry.py: use len(os.sched_getaffinity(0)), or honour
SLURM_CPUS_PER_TASK, instead of os.cpu_count(). This wrapper is the same fix applied from
outside so it works without waiting for a release. Delete it once upstream lands.

USAGE
-----
Drop-in for the `genotools` console script — identical arguments:

    python3 scripts/genotools_capped.py --bfile ... --out ... --ancestry ...

Override the cap explicitly if needed:

    GENOTOOLS_MAX_WORKERS=32 python3 scripts/genotools_capped.py ...
"""

import os
import sys


def _allocated_cpus():
    """CPUs this job actually owns, most authoritative source first."""
    override = os.environ.get("GENOTOOLS_MAX_WORKERS")
    if override and override.isdigit() and int(override) > 0:
        return int(override)

    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm and slurm.isdigit() and int(slurm) > 0:
        return int(slurm)

    # Affinity is right when something has pinned us; os.cpu_count() never is on a
    # shared node. os.sched_getaffinity is Linux-only, hence the fallback.
    try:
        return max(len(os.sched_getaffinity(0)), 1)
    except (AttributeError, OSError):
        return max(os.cpu_count() or 1, 1)


def _allocated_bytes():
    """Memory this job actually owns. SLURM reports MB; 0 means 'all of the node'."""
    for var in ("SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU"):
        raw = os.environ.get(var)
        if not (raw and raw.isdigit()):
            continue
        mb = int(raw)
        if mb <= 0:                      # --mem=0 asks for the whole node
            continue
        if var == "SLURM_MEM_PER_CPU":
            mb *= _allocated_cpus()
        return mb * 1024 * 1024
    return None                          # not under SLURM — leave psutil alone


def main():
    _node_cpus = os.cpu_count()          # read the truth before shadowing it
    cpus = _allocated_cpus()
    mem = _allocated_bytes()

    # Patch BEFORE importing genotools so ancestry.py reads the constrained values.
    os.cpu_count = lambda: cpus

    if mem is not None:
        import psutil

        _real_virtual_memory = psutil.virtual_memory

        def _capped_virtual_memory():
            # svmem is a namedtuple, so _replace keeps every other field intact.
            vm = _real_virtual_memory()
            return vm._replace(total=min(vm.total, mem),
                               available=min(vm.available, mem))

        psutil.virtual_memory = _capped_virtual_memory

    gb = "unconstrained" if mem is None else f"{mem / 1024**3:.1f}GB"
    print(f"[genotools_capped] reporting {cpus} CPUs / {gb} to GenoTools", flush=True)
    print(f"[genotools_capped] node actually has {_node_cpus} CPUs — the difference is the "
          f"point; ulimit -u is the binding constraint here.", flush=True)

    from genotools.__main__ import handle_main
    return handle_main()


if __name__ == "__main__":
    sys.exit(main())

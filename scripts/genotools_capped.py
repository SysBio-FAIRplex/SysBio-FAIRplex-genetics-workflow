#!/usr/bin/env python3
"""Run GenoTools with its worker pool sized to the SLURM allocation instead of the node.

GenoTools (ancestry.py:513-518) sizes its GridSearchCV pool from os.cpu_count() and
psutil.virtual_memory(), which report NODE totals. On a shared biowulf node that pool exceeds the
per-user `ulimit -u` of 1024: workers abort and SLURM logs a bare `ExitCode 1:0`. The worker count
is the lever — thread env vars are not (joblib already sets them per worker). This patches both
readings before GenoTools imports. Delete it once upstream uses os.sched_getaffinity or
SLURM_CPUS_PER_TASK.

Drop-in for the `genotools` console script, identical arguments. GENOTOOLS_MAX_WORKERS overrides
the cap.
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

    # Affinity reflects any pinning; os.cpu_count() ignores it.
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

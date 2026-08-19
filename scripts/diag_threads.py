#!/usr/bin/env python3
"""Measure what a loky worker ACTUALLY sees and how many OS threads it creates.

Written because three successive fixes to the genotools thread-exhaustion bug were reasoned
about rather than measured, and two of the three were wrong. This reproduces the exact shape
of genotools' GridSearchCV call — joblib Parallel over N loky workers, each touching numba
through UMAP — and reports the numbers instead of inferring them.

Run it under the SAME allocation the real job uses:

    srun --cpus-per-task=64 --mem=256g -t 15 python3 scripts/diag_threads.py
    srun --cpus-per-task=64 --mem=256g -t 15 python3 scripts/diag_threads.py 16   # n_jobs=16

Questions it answers:
  1. Does the parent's os.cpu_count() monkeypatch survive into a spawned worker? (expected: NO)
  2. Does joblib actually set NUMBA_NUM_THREADS in worker environments? (docs say no, source says
     yes — MAX_NUM_THREADS_VARS includes it; this settles which is true for THIS joblib version)
  3. How many OS-level threads does one worker really hold once numba has launched? Python's
     threading.active_count() cannot see numba's native threads, so read /proc/self/task.
  4. Peak total tasks for this user on this node, against `ulimit -u`.

The product of (3) and n_jobs versus `ulimit -u` is the whole bug. Everything else is commentary.
"""

import os
import subprocess
import sys


def probe(i):
    """Runs INSIDE a loky worker. Reports what that process can see about itself."""
    import os

    before = len(os.listdir("/proc/self/task"))

    # Force numba's lazy thread-pool launch — the exact call UMAP makes at umap_.py:2295,
    # and the call that appears at the abort site in every failed run.
    import numba
    numba.get_num_threads()

    after = len(os.listdir("/proc/self/task"))

    return {
        "worker": i,
        "os_cpu_count": os.cpu_count(),
        "sched_affinity": len(os.sched_getaffinity(0)),
        "env_NUMBA": os.environ.get("NUMBA_NUM_THREADS"),
        "env_OMP": os.environ.get("OMP_NUM_THREADS"),
        "numba_threads": numba.get_num_threads(),
        "os_threads_before": before,
        "os_threads_after": after,
    }


def user_tasks():
    """Total tasks (processes + threads) for this user on this node — what RLIMIT_NPROC counts."""
    try:
        out = subprocess.run(["ps", "-eLf", "-u", str(os.getuid())],
                             capture_output=True, text=True, check=True)
        return max(len(out.stdout.strip().splitlines()) - 1, 0)
    except Exception as exc:                                  # noqa: BLE001 - diagnostic only
        return f"unavailable ({exc})"


def main():
    n_jobs = int(sys.argv[1]) if len(sys.argv) > 1 else 64

    import joblib
    from joblib import Parallel, delayed

    print("=" * 78)
    print("PARENT")
    print(f"  joblib               {joblib.__version__}")
    print(f"  os.cpu_count()       {os.cpu_count()}")
    print(f"  sched_getaffinity    {len(os.sched_getaffinity(0))}")
    print(f"  SLURM_CPUS_PER_TASK  {os.environ.get('SLURM_CPUS_PER_TASK')}")
    print(f"  ulimit -u (soft)     {__import__('resource').getrlimit(__import__('resource').RLIMIT_NPROC)[0]}")
    print(f"  user tasks at start  {user_tasks()}")
    print(f"  n_jobs for this run  {n_jobs}")
    print("=" * 78)

    results = Parallel(n_jobs=n_jobs)(delayed(probe)(i) for i in range(n_jobs))
    peak = user_tasks()   # workers are torn down by now, but catches anything left behind

    ok = [r for r in results if isinstance(r, dict)]
    sample = ok[0] if ok else {}

    print("\nWORKER (first one; all should agree)")
    for key in ("os_cpu_count", "sched_affinity", "env_NUMBA", "env_OMP",
                "numba_threads", "os_threads_before", "os_threads_after"):
        print(f"  {key:<20} {sample.get(key)}")

    per_worker = sample.get("os_threads_after")
    print("\nARITHMETIC")
    if isinstance(per_worker, int):
        print(f"  threads per worker   {per_worker}")
        print(f"  x n_jobs ({n_jobs})".ljust(23) + f"{per_worker * n_jobs}")
        print(f"  vs ulimit -u         {__import__('resource').getrlimit(__import__('resource').RLIMIT_NPROC)[0]}")
    print(f"  user tasks after     {peak}")

    print("\nREAD IT LIKE THIS")
    print("  env_NUMBA is None      -> joblib is NOT capping numba; set NUMBA_NUM_THREADS yourself")
    print("  os_cpu_count == node   -> the parent's monkeypatch did NOT cross the spawn boundary")
    print("  threads*n_jobs >= ulim -> that product is the bug; lower n_jobs or threads per worker")


if __name__ == "__main__":
    sys.exit(main())

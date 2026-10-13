"""Phase-3 follow-up ROUND 2, preflight (docs/phase3_followup.md).  Gates the
invalidation and the run.  ~5-10 min on >= 3 cores, cosmo_env active, from the
repo root:

    python GCE/scripts/phase3_followup_preflight.py          # all checks
    python GCE/scripts/phase3_followup_preflight.py -n 3     # cap the worker pool

Every check runs node_joint of the GENERATED round-2 drivers, so it tests what
the job will actually run: the 1e-3 first try, the PHASE3_REQUIRE retry ladder,
the early-stopping alpha sweep and the local crossing, end to end.

  G1  three round-1 glitch nodes come back repaired.  In round 1 the 1e-3 trigger
      left one or more of their alphas on the legacy path (Omega 2-7x low) and
      the crossing fell by 3.4-5.7x.  Each must now give the crossing measured
      on 2026-10-13 (all alphas on Phase 3, crossing from the four points around
      it) to 0.5 %, with:
        - no alpha legacy or failed;
        - at least one alpha engaged by the retry ladder (Sdot < 1e-3);
        - the alphas past the crossing skipped, not solved.

Exit code 0 iff G1 passes on all three.  Prints a GATE line last.
"""

import argparse
import importlib.util
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
POOL_TIMEOUT = 3600          # s per solve
G1_TOL = 0.005

# (driver, i, j, round-1 crossing, repaired crossing, alphas solved), measured
# 2026-10-13 through solver 9098886f73 (scratch fu_glitch_axis.py, test_sweep_live.py)
G1_NODES = [
    ("d1_e09", 44, 55, 1.7140e-03, 9.7117e-03, 7),    # mX=95.87 rv=0.582
    ("d1_e09", 12, 52, 6.0581e-04, 1.7966e-03, 5),    # mX=24.88 rv=0.509
    ("d1_e10",  1, 62, 6.9504e-04, 3.3703e-03, 6),    # mX=15.65 rv=0.795
]


def _load_driver(key, stem):
    """Import a generated driver under its own module name, so mp.Pool can
    pickle its functions."""
    spec = importlib.util.spec_from_file_location(key, HERE / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[key] = mod
    spec.loader.exec_module(mod)
    mod.SOLVE_TIMEOUT = POOL_TIMEOUT
    mod._set_alphax(mod.EPS_JOINT[0])
    return mod


DRIVERS = {
    "d1_e09": _load_driver("d1_e09", "make_relic_grid_vP_fu2_e09"),
    "d1_e10": _load_driver("d1_e10", "make_relic_grid_vP_fu2_e10"),
}


def _node(task):
    """Pool worker: (driver key, i, j) -> node_joint's (alpha, path, sdot), seconds."""
    key, i, j = task
    mod = DRIVERS[key]
    t0 = time.perf_counter()
    _, _, alpha, _, path, sdot = mod.node_joint(
        (i, j, float(mod.mX_grid[i, j]), float(mod.mY_grid[i, j]), mod.EPS_JOINT[0]))
    return alpha, path, sdot, time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser(description="Phase-3 follow-up round-2 preflight")
    ap.add_argument("-n", "--ncores", type=int, default=min(3, mp.cpu_count()))
    args = ap.parse_args()

    for key, mod in DRIVERS.items():
        kw = mod.PHASE3_KW
        assert kw.get("enable_phase3") and kw.get("Sdot_switch_thresh") == 1e-3, \
            f"{key}: PHASE3_KW is not the follow-up's forced setting: {kw}"
        assert mod.PHASE3_REQUIRE, f"{key}: PHASE3_REQUIRE is off"
        assert mod.ALPHAX_SWEEP_STOP, f"{key}: ALPHAX_SWEEP_STOP is off"
        assert mod.ALPHAX.max() == 1.0, f"{key}: alpha_X axis is not (2e-5, 1.0)"

    tasks = [(k, i, j) for k, i, j, *_ in G1_NODES]
    print(f"preflight: {len(tasks)} nodes on {args.ncores} workers "
          f"(cap {POOL_TIMEOUT} s/solve)", flush=True)
    with mp.Pool(args.ncores) as pool:
        res = list(pool.imap(_node, tasks, chunksize=1))

    go = True
    print("\nG1 round-1 glitch nodes, repaired:")
    for (key, i, j, old, want, n_solved), (alpha, path, sdot, secs) in zip(G1_NODES, res):
        mod = DRIVERS[key]
        solved = path != mod.PATH_SKIPPED
        checks = {
            "crossing": bool(np.isfinite(alpha) and abs(alpha / want - 1) <= G1_TOL),
            "no legacy/failed": not np.any(np.isin(path, (mod.PATH_LEGACY, mod.PATH_FAILED))),
            "retry used": bool(np.any(sdot[solved & (path == mod.PATH_P3)] < 1e-3)),
            "tail skipped": int(solved.sum()) == n_solved and bool(np.all(solved[:n_solved])),
        }
        ok = all(checks.values())
        go &= ok
        print(f"   eps={mod.EPS_JOINT[0]:.0e} mX={mod.mX_grid[i, j]:.2f} "
              f"rv={10 ** mod.log_rv[j]:.3f}: alpha_relic={alpha:.4e} (repaired {want:.4e}, "
              f"round 1 {old:.4e}) [{secs:.0f} s] {'PASS' if ok else 'FAIL'}")
        print(f"      paths {path.tolist()}  Sdot {sdot.tolist()}")
        for name, c in checks.items():
            if not c:
                print(f"      FAILED: {name}")

    print(f"\nGATE  invalidate + run: {'GO' if go else 'NO-GO'}")
    sys.exit(0 if go else 1)


if __name__ == "__main__":
    main()

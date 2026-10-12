"""Regression test for relic_coupling against the campaign's validated anchor.

The 2026-10 Phase-3 campaign measured alpha_relic = 8.3862e-03 at
mX = 34.86, rv = 0.950, xi_ini = 2.5, epsX = 1e-9 with Phase 3 on
(docs/phase3_rerun_campaign.md, probe P1).  phase3_preflight.py checks that
through the legacy make_relic_grid_vP*.py drivers, which carry their own copy
of the crossing code; this checks it through make_relic_grid.py, the shared
one, so a change to relic_coupling is actually exercised.

    python GCE/scripts/anchor_check.py -n 8

Exit code 0 iff the anchor comes back within ANCHOR_TOL.
"""

import argparse
import importlib.util
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

HERE = Path(__file__).resolve().parent

ANCHOR_MX = 34.856        # GeV
ANCHOR_RV = 0.950         # mY / mX
ANCHOR_XI = 2.5           # T_D/T_SM at Tinf
ANCHOR_EPS = 1e-9
ANCHOR_VALUE = 8.3862e-03  # the campaign's measured crossing
ANCHOR_TOL = 0.02          # |ratio - 1|


def _load():
    spec = importlib.util.spec_from_file_location("drv", HERE / "make_relic_grid.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["drv"] = mod
    spec.loader.exec_module(mod)
    return mod


DRV = _load()


def _task(arg):
    """Pool worker: Omega h^2 at one coupling, through the shared driver."""
    mX, mY, coupling = arg
    DRV._set_portal("vP", ANCHOR_XI)
    try:
        return DRV.omega_h2_joint(mX, mY, coupling, ANCHOR_EPS)
    except Exception:
        return np.nan


def main():
    ap = argparse.ArgumentParser(description="anchor regression for relic_coupling")
    ap.add_argument("-n", "--ncores", type=int, default=min(8, mp.cpu_count()))
    args = ap.parse_args()

    DRV._set_portal("vP", ANCHOR_XI)
    DRV._set_coupling(ANCHOR_EPS)

    i = int(np.argmin(np.abs(10.0 ** DRV.log_mX - ANCHOR_MX)))
    j = int(np.argmin(np.abs(10.0 ** DRV.log_rv - ANCHOR_RV)))
    mX, mY = float(DRV.mX_grid[i, j]), float(DRV.mY_grid[i, j])

    print(f"anchor node mX={mX:.3f} mY={mY:.3f} rv={mY / mX:.4f}  "
          f"xi_ini={ANCHOR_XI} eps={ANCHOR_EPS:g}  Phase 3 "
          f"{'ON' if DRV.PHASE3_KW else 'OFF'}", flush=True)

    tasks = [(mX, mY, float(c)) for c in DRV.COUPLING]
    with mp.Pool(args.ncores) as pool:
        omega = np.array(pool.map(_task, tasks))

    for c, w in zip(DRV.COUPLING, omega):
        print(f"    {c:10.4g}  {w:12.4g}" + ("" if np.isfinite(w) else "   <- failed"))

    crossing = DRV.relic_coupling(omega)
    ratio = crossing / ANCHOR_VALUE
    ok = bool(np.isfinite(ratio) and abs(ratio - 1.0) <= ANCHOR_TOL)

    print(f"\nrelic_coupling = {crossing:.4e}   anchor {ANCHOR_VALUE:.4e}   "
          f"ratio {ratio:.4f}   {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

"""Phase-3 follow-up, preflight (docs/phase3_followup.md).  Gates the
invalidation and the run.  ~10 min on >= 16 cores, cosmo_env active, from the
repo root:

    python GCE/scripts/phase3_followup_preflight.py          # all checks
    python GCE/scripts/phase3_followup_preflight.py -n 16    # cap the worker pool

Every check goes through the GENERATED follow-up drivers (omega_h2_joint, with
their PHASE3_KW), so it tests what the job will actually run.

  F1  forced Phase 3 keeps the validated anchor: alpha_relic at xi=2.5,
      eps=1e-9, mX=34.86, rv=0.950 must come back 8.386e-3 to 2 %.
  F2  the diagnosis reproduces: Omega/0.12 at the STORED (pre-follow-up)
      crossing of seven nodes must match the values measured on 2026-10-12 to
      1 %.  These are the corrections the follow-up applies (<sigma v> moves by
      exactly this factor); a mismatch means this solver or config is not the
      one the follow-up was planned with.
  F3  forced Phase 3 removes the mixed path at the worst xi=1 glitch node
      (mX=24.88, rv=0.760, eps=1e-9): at alpha=2.065e-3, where the campaign
      trigger stayed on the legacy path (Omega = 0.181), Omega must come back
      on the Phase-3 branch.

Exit code 0 iff F1, F2 and F3 all pass.  Prints a GATE line last.
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

ANCHOR = (20, 66, 8.3862e-3)  # (i, j, alpha_relic): mX=34.856, rv=0.950
ANCHOR_TOL = 0.02
F2_TOL = 0.01

# (driver, i, j, stored alpha_relic, Omega/0.12 there with Phase 3 forced),
# measured 2026-10-12 through solver 55550e14fe with the follow-up's PHASE3_KW.
F2_NODES = [
    ("d25_e09", 20, 26, 1.5090674965932695e-03, 1.12826),   # rv 0.159
    ("d25_e09", 20, 48, 3.0313732473365943e-03, 0.83275),   # rv 0.426
    ("d25_e10", 20, 26, 5.0040055616616480e-04, 1.40270),   # rv 0.159
    ("d25_e10", 20, 52, 1.4673228529848278e-03, 1.23277),   # rv 0.509
    ("d1_e09",  20, 26, 1.3707508141576406e-03, 0.91846),   # rv 0.159
    ("d1_e09",  20, 48, 2.0697999859483436e-03, 1.04070),   # rv 0.426
    ("d1_e10",  20, 66, 3.4156680113146584e-02, 1.12239),   # rv 0.950
]
F3_NODE = ("d1_e09", 12, 61, 3)          # alpha = ALPHAX[3] = 2.065e-3
F3_LEGACY = 0.18145                       # campaign trigger, legacy path
F3_FORCED = 0.308                         # forced trigger, Phase-3 path


def _load_driver(key, stem):
    """Import a generated driver under a dot-free module name (the xi_ini file
    names contain '.', which breaks pickling of their functions in mp.Pool)."""
    spec = importlib.util.spec_from_file_location(key, HERE / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[key] = mod
    spec.loader.exec_module(mod)
    mod.SOLVE_TIMEOUT = POOL_TIMEOUT
    mod._set_alphax(mod.EPS_JOINT[0])
    return mod


DRIVERS = {
    "d25_e09": _load_driver("d25_e09", "make_relic_grid_vP_xi_ini_2.5_fu_e09"),
    "d25_e10": _load_driver("d25_e10", "make_relic_grid_vP_xi_ini_2.5_fu_e10"),
    "d1_e09":  _load_driver("d1_e09",  "make_relic_grid_vP_fu_e09"),
    "d1_e10":  _load_driver("d1_e10",  "make_relic_grid_vP_fu_e10"),
}


def _omega(task):
    """Pool worker: (driver key, i, j, alpha) -> Omega h^2 through the driver."""
    key, i, j, alpha = task
    mod = DRIVERS[key]
    t0 = time.perf_counter()
    val = mod.omega_h2_joint(float(mod.mX_grid[i, j]), float(mod.mY_grid[i, j]),
                             alpha, mod.EPS_JOINT[0])
    return float(val), time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser(description="Phase-3 follow-up preflight")
    ap.add_argument("-n", "--ncores", type=int, default=min(32, mp.cpu_count()))
    args = ap.parse_args()

    for key, mod in DRIVERS.items():
        kw = mod.PHASE3_KW
        assert kw.get("enable_phase3") and kw.get("Sdot_switch_thresh") == 1e-3, \
            f"{key}: PHASE3_KW is not the follow-up's forced setting: {kw}"
        assert mod.ALPHAX.max() == 1.0, f"{key}: alpha_X axis is not (2e-5, 1.0)"

    d25 = DRIVERS["d25_e09"]
    i, j, _ = ANCHOR
    tasks = [("d25_e09", i, j, float(a)) for a in d25.ALPHAX]
    tasks += [(k, i, j, a) for k, i, j, a, _ in F2_NODES]
    k3, i3, j3, a3 = F3_NODE
    tasks.append((k3, i3, j3, float(DRIVERS[k3].ALPHAX[a3])))

    print(f"preflight: {len(tasks)} solves on {args.ncores} workers "
          f"(cap {POOL_TIMEOUT} s/solve)", flush=True)
    with mp.Pool(args.ncores) as pool:
        res = []
        for k, (om, secs) in enumerate(pool.imap(_omega, tasks, chunksize=1)):
            print(f"  [{k + 1:2d}/{len(tasks)}] {tasks[k]}: Omega={om:.5g} ({secs:.0f} s)",
                  flush=True)
            res.append(om)
    res = np.array(res)
    n1 = len(d25.ALPHAX)

    alpha = d25.relic_alphaX(res[:n1])
    ratio = alpha / ANCHOR[2]
    f1 = bool(np.isfinite(alpha) and abs(ratio - 1) <= ANCHOR_TOL)
    print(f"\nF1 anchor xi=2.5 eps=1e-9 rv=0.950: alpha_relic={alpha:.4e} "
          f"(anchor {ANCHOR[2]:.4e}, ratio {ratio:.4f}) {'PASS' if f1 else 'FAIL'}")

    f2 = True
    print("F2 Omega/0.12 at the stored crossing (= the <sigma v> correction):")
    for (key, i, j, a, want), om in zip(F2_NODES, res[n1:n1 + len(F2_NODES)]):
        got = om / 0.12
        ok = bool(np.isfinite(got) and abs(got / want - 1) <= F2_TOL)
        f2 &= ok
        mod = DRIVERS[key]
        print(f"   xi={getattr(mod, 'XI_INI', 1.0):<3} eps={mod.EPS_JOINT[0]:.0e} "
              f"rv={10 ** mod.log_rv[j]:.3f}: {got:.4f} (expected {want:.4f}) "
              f"{'PASS' if ok else 'FAIL'}")

    om3 = res[-1]
    f3 = bool(np.isfinite(om3) and abs(om3 / F3_FORCED - 1) <= 0.02)
    print(f"F3 glitch node xi=1 eps=1e-9 mX=24.88 rv=0.760 alpha={tasks[-1][3]:.4e}: "
          f"Omega={om3:.4g} (Phase-3 branch {F3_FORCED}, legacy {F3_LEGACY}) "
          f"{'PASS' if f3 else 'FAIL'}")

    go = f1 and f2 and f3
    print(f"\nGATE  invalidate + run: {'GO' if go else 'NO-GO'}")
    sys.exit(0 if go else 1)


if __name__ == "__main__":
    main()

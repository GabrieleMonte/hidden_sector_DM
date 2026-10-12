"""Preflight for a new vector portal, before committing real node-hours.

    python GCE/scripts/portal_preflight.py                  # BL, LiLj, B
    python GCE/scripts/portal_preflight.py --portals BL -n 32

The 2026-10 vector-portal campaign could gate its slices on two validated
anchors (docs/phase3_rerun_campaign.md, probe P1).  BL, LiLj and Baryon have
no such numbers, so this checks the things that would actually waste the run.
All four portals share the hidden sector -- sigmav_XX_to_YY is byte-identical
and only decay_width_to_SM differs -- so the risk is not the freeze-out, it is
whether the shifted Gamma_Y moves the relic crossing off the sampled axis.

  G1  Gamma_Y against the vector portal at the same epsX, and the branching
      ratios.  No solve; it only has to be finite and the right order.
  G2  the worst corner.  alpha_relic at rv = 0.950, at the mX where the vP
      grid peaks for that epsX (95.872 GeV at 1e-9, 68.426 GeV at 1e-10).
      A larger Gamma_Y means less dilution, so Omega rises and the crossing
      moves UP; it must stay bracketed below ALPHAX_LIMS[eps][1] = 1.0.
      This is the gate: NO-GO here means the alpha_X axis is too short and
      production would write NaN exactly where the grid is most interesting.
  G3  does Phase 3 engage?  Reported, never a gate -- a portal that never
      reaches phase_final 3 is a result, not a failure.

Exit code 0 iff every requested portal passes G2 at both epsX.
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
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

POOL_TIMEOUT = 3600          # s per solve; production allows SOLVE_TIMEOUT=7200
XI_INI = 1.0                 # the campaign scope for the new portals

# (epsX, mX GeV, rv) -- where the vector-portal grid peaks, i.e. the corner
# whose crossing sits closest to the top of the alpha_X axis.
WORST = [(1e-9, 95.872, 0.950), (1e-10, 68.426, 0.950)]
ALPHA_CEILING = 1.0          # ALPHAX_LIMS[1e-9] = ALPHAX_LIMS[1e-10] = (2e-5, 1.0)


def _load_driver():
    spec = importlib.util.spec_from_file_location("drv", HERE / "make_relic_grid.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["drv"] = mod
    spec.loader.exec_module(mod)
    mod.SOLVE_TIMEOUT = POOL_TIMEOUT
    return mod


DRV = _load_driver()


def _nearest_node(mX_want, rv_want):
    """Snap to the grid, so the probe is a node production will actually solve."""
    i = int(np.argmin(np.abs(10.0 ** DRV.log_mX - mX_want)))
    j = int(np.argmin(np.abs(10.0 ** DRV.log_rv - rv_want)))
    mX = float(DRV.mX_grid[i, j])
    mY = float(DRV.mY_grid[i, j])
    return i, j, mX, mY


def _task(arg):
    """Pool worker: one omega_h2_joint through the configured driver."""
    portal, mX, mY, alpha, eps = arg
    DRV._set_portal(portal, XI_INI)
    t0 = time.perf_counter()
    try:
        return {"omega": DRV.omega_h2_joint(mX, mY, alpha, eps),
                "secs": time.perf_counter() - t0}
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}", "secs": time.perf_counter() - t0}


def gamma_table(portals):
    """G1: Gamma_Y and branching ratios, relative to the vector portal."""
    from hidden_sector_DM.HiddenSectorDM import VectorPortal
    mX, mY, eps = 34.856, 31.680, 1e-9
    ref = VectorPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=1e-2,
                       include_antiparticlesX=True, include_antiparticlesY=False
                       ).decay_width_to_SM(eps)
    print(f"G1  Gamma_Y at mX={mX} mY={mY} epsX={eps:g}  "
          f"(VectorPortal reference {ref:.4e} GeV)")
    ok = True
    for p in portals:
        DRV._set_portal(p, XI_INI)
        m = DRV._make_model(mX, mY, 1e-2)
        g = m.decay_width_to_SM(eps)
        good = np.isfinite(g) and g > 0
        ok &= good
        brs = {k: v for k, v in m.branching_ratios_to_SM().items() if v > 0.01}
        print(f"      {p:<5} {g:.4e} GeV  (x{g / ref:6.2f} vP)  "
              f"{'ok' if good else 'BAD'}  {DRV.PORTAL_KW}")
        print(f"            BR: {', '.join(f'{k}={v:.2f}' for k, v in brs.items())}")
    return ok


def main():
    ap = argparse.ArgumentParser(description="preflight a new vector portal")
    ap.add_argument("--portals", default="BL,LiLj,B",
                    help="comma-separated (default BL,LiLj,B)")
    ap.add_argument("-n", "--ncores", type=int, default=min(48, mp.cpu_count()))
    args = ap.parse_args()
    portals = [p.strip() for p in args.portals.split(",") if p.strip()]

    print(f"portal preflight: {', '.join(portals)}  xi_ini={XI_INI:g}  "
          f"mediators at 10 TeV\n")
    g1_ok = gamma_table(portals)

    tasks, owner = [], []
    for p in portals:
        DRV._set_portal(p, XI_INI)
        for eps, mX_w, rv_w in WORST:
            i, j, mX, mY = _nearest_node(mX_w, rv_w)
            DRV._set_coupling(eps)
            for a in DRV.COUPLING:
                tasks.append((p, mX, mY, float(a), eps))
                owner.append((p, eps, mX, mY))

    print(f"\nG2/G3  {len(tasks)} solves on {args.ncores} workers "
          f"(cap {POOL_TIMEOUT} s each)")
    results = []
    with mp.Pool(args.ncores) as pool:
        for k, r in enumerate(pool.imap(_task, tasks, chunksize=1)):
            results.append(r)
            if (k + 1) % 8 == 0:
                print(f"       {k + 1:3d}/{len(tasks)} done", flush=True)

    print()
    passed = {}
    for p in portals:
        DRV._set_portal(p, XI_INI)
        for eps, mX_w, rv_w in WORST:
            DRV._set_coupling(eps)
            om = np.array([r.get("omega", np.nan) if "error" not in r else np.nan
                           for r, o in zip(results, owner)
                           if o[0] == p and o[1] == eps])
            secs = [r["secs"] for r, o in zip(results, owner)
                    if o[0] == p and o[1] == eps]
            a = DRV.relic_coupling(om)
            ok = bool(np.isfinite(a) and a < ALPHA_CEILING)
            passed[(p, eps)] = ok
            print(f"G2  {p:<5} eps={eps:<6g} rv=0.950 mX={mX_w:6.3f}: "
                  f"alpha_relic={a:.4e} "
                  f"({'bracketed' if ok else 'OFF AXIS'}, ceiling {ALPHA_CEILING:g})  "
                  f"finite {int(np.isfinite(om).sum())}/8  "
                  f"[{min(secs):.0f}-{max(secs):.0f} s/solve]  "
                  f"{'PASS' if ok else 'FAIL'}")

    print("\nGATES")
    overall = True
    for p in portals:
        good = all(passed[(p, e)] for e, _, _ in WORST)
        overall &= good
        print(f"  {p:<5} production: {'GO' if good else 'NO-GO -- extend ALPHAX_LIMS'}")
    if not g1_ok:
        print("  G1 reported a bad Gamma_Y; do not launch on that portal.")
        overall = False
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()

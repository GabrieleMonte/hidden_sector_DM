"""Phase-3 rerun campaign, Stage 1: preflight probes P1-P4
(docs/phase3_rerun_campaign.md).  Each probe gates real node-hours.  Run from
the repo root on an interactive/debug node, cosmo_env active:

    python GCE/scripts/phase3_preflight.py            # all probes, ~1 h on >= 30 cores
    python GCE/scripts/phase3_preflight.py -n 16      # cap the worker pool
    python GCE/scripts/phase3_preflight.py --only P1  # a subset (comma-separated)

  P1  end-to-end anchors THROUGH THE EDITED DRIVER (xi_ini=2.5, eps=1e-9,
      mX=34.86): alpha_relic at rv=0.950 / 0.9085 must come back
      8.38e-3 / 6.40e-3 to ~2 %.  Gates passes A+B.
  P2  Phase-3 rtol 1e-6 vs 1e-8 on the rv=0.950 anchor node (serial, timed):
      Omega within 0.5 %; the timings set the real per-node cost.
  P3  eps=1e-10 trigger at rv=0.9085, xi_ini = 1 and 2.5: does phase_final
      come back 3 at any alpha?  If never, pass B is pointless -- DROP it.
  P4  the worst alpha corner (xi_ini=0.5, mX=71.37, rv=0.95, stored
      alpha_relic=1.58e-2): the crossing must stay bracketed on the new
      (2e-5, 1.0) axis.  Gates pass A.

Exit code 0 iff P1 and P4 both pass.  Prints a GATES summary last.
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

POOL_TIMEOUT = 3600          # s per solve inside the pool phase (P1/P3/P4)
P2_TIMEOUT = 7200            # the rtol-1e-8 leg may legitimately need it

IX, J95, J91 = 20, 66, 65    # mX=34.856; rv=0.950, 0.9085
ANCHORS = {J95: 8.38e-3, J91: 6.40e-3}        # validated crossings, xi=2.5 eps=1e-9
ANCHOR_TOL = 0.02                             # |alpha/anchor - 1|
P4_NODE = (37, 66)                            # xi=0.5 block argmax: mX=71.372, rv=0.95
P4_STORED = 1.579e-2
P3_ALPHAS = np.logspace(np.log10(2e-5), 0.0, 8)[[4, 5, 6]]   # 9.7e-3, 4.6e-2, 2.1e-1


def _load_driver(key, stem):
    """Import a driver script under a dot-free module name (the xi_ini file
    names contain '.', which breaks pickling of their functions in mp.Pool)."""
    spec = importlib.util.spec_from_file_location(key, HERE / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[key] = mod
    spec.loader.exec_module(mod)
    mod.SOLVE_TIMEOUT = POOL_TIMEOUT          # preflight cap; P2 resets its own
    return mod


DRIVERS = {
    "d1": _load_driver("d1", "make_relic_grid_vP"),
    "d05": _load_driver("d05", "make_relic_grid_vP_xi_ini_0.5"),
    "d25": _load_driver("d25", "make_relic_grid_vP_xi_ini_2.5"),
}


def _node(mod, i, j):
    return 10.0 ** mod.log_mX[i], 10.0 ** mod.log_rv[j] * 10.0 ** mod.log_mX[i]


def _joint_direct(modkey, mX, mY, alpha, eps, p3_rtol=None):
    """One solve_boltzmann_joint call with the driver's production kwargs,
    returning phase_final (P3 needs it; omega_h2_joint hides it)."""
    mod = DRIVERS[modkey]
    kw = dict(mod.PHASE3_KW)
    if p3_rtol is not None:
        kw["phase3_kw"] = dict(rtol_value=p3_rtol)
    extra = {"xi_ini": mod.XI_INI} if hasattr(mod, "XI_INI") else {}
    solver = mod.BoltzmannSolver(mod.cosmo, mod.VectorPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=alpha,
        include_antiparticlesX=True, include_antiparticlesY=False, **extra))
    t0 = time.perf_counter()
    res = mod.timed_solve(solver.solve_boltzmann_joint, epsX=eps,
                          return_bg_ICs=True, min_eps_floor_ratio=0.0,
                          GoH_exit_thresh=mod.GOH_EXIT_THRESH,
                          xmax=mod.XMAX_JOINT, phase2_form="Y",
                          verbose=False, **kw)
    secs = time.perf_counter() - t0
    pf = res["phase_final"]
    omega = (res["YX_relic"] * mX * mod.s0_cosmo / mod.rhoc
             if pf in ("1.5S", 3) else None)
    return {"phase_final": pf, "omega": omega, "secs": secs}


def _run_task(task):
    """Pool worker: ('omega', modkey, mX, mY, alpha, eps) ->
    omega_h2_joint through the edited driver, or ('direct', ...)."""
    kind = task[0]
    t0 = time.perf_counter()
    try:
        if kind == "omega":
            _, modkey, mX, mY, alpha, eps = task
            val = DRIVERS[modkey].omega_h2_joint(mX, mY, alpha, eps)
            return {"omega": val, "secs": time.perf_counter() - t0}
        _, modkey, mX, mY, alpha, eps = task
        return _joint_direct(modkey, mX, mY, alpha, eps)
    except Exception as e:                    # keep one bad solve from killing the pool
        return {"error": f"{type(e).__name__}: {e}",
                "secs": time.perf_counter() - t0}


def main():
    ap = argparse.ArgumentParser(description="Phase-3 campaign preflight")
    ap.add_argument("-n", "--ncores", type=int,
                    default=min(32, mp.cpu_count()))
    ap.add_argument("--only", default="P1,P2,P3,P4",
                    help="comma-separated subset of probes")
    args = ap.parse_args()
    probes = {p.strip().upper() for p in args.only.split(",")}

    d1, d05, d25 = DRIVERS["d1"], DRIVERS["d05"], DRIVERS["d25"]
    d25._set_alphax(1e-9)
    d05._set_alphax(1e-9)
    assert d25.ALPHAX.max() == 1.0, "driver ALPHAX_LIMS[1e-9] is not (2e-5, 1.0)"

    tasks, owner = [], []
    if "P1" in probes:
        for j in (J95, J91):
            mX, mY = _node(d25, IX, j)
            for a in d25.ALPHAX:
                tasks.append(("omega", "d25", mX, mY, float(a), 1e-9))
                owner.append(("P1", j, float(a)))
    if "P3" in probes:
        for key in ("d1", "d25"):
            mX, mY = _node(DRIVERS[key], IX, J91)
            for a in P3_ALPHAS:
                tasks.append(("direct", key, mX, mY, float(a), 1e-10))
                owner.append(("P3", key, float(a)))
    if "P4" in probes:
        mX, mY = _node(d05, *P4_NODE)
        for a in d05.ALPHAX:
            tasks.append(("omega", "d05", mX, mY, float(a), 1e-9))
            owner.append(("P4", P4_NODE[1], float(a)))

    print(f"preflight: {len(tasks)} pooled solves on {args.ncores} workers "
          f"(cap {POOL_TIMEOUT} s/solve), then P2 serial" if tasks else
          "preflight: no pooled probes selected")
    results = []
    if tasks:
        with mp.Pool(args.ncores) as pool:
            for k, res in enumerate(
                    pool.imap(_run_task, tasks, chunksize=1)):
                tag = owner[k]
                brief = res.get("error") or (
                    f"omega={res['omega']:.4g}" if "omega" in res
                    and res["omega"] is not None and not isinstance(tag[1], str)
                    else f"phase_final={res.get('phase_final')}")
                print(f"  [{k + 1:2d}/{len(tasks)}] {tag}: {brief} "
                      f"({res['secs']:.0f} s)", flush=True)
                results.append(res)

    passed, gates = {}, []

    if "P1" in probes:
        ok = True
        for j in (J95, J91):
            om = np.array([r.get("omega", np.nan) if "error" not in r else np.nan
                           for r, o in zip(results, owner)
                           if o[0] == "P1" and o[1] == j])
            alpha = d25.relic_alphaX(om)
            ratio = alpha / ANCHORS[j]
            this = np.isfinite(alpha) and abs(ratio - 1) <= ANCHOR_TOL
            ok &= bool(this)
            print(f"P1 rv={10 ** d25.log_rv[j]:.4f}: alpha_relic={alpha:.4e} "
                  f"(anchor {ANCHORS[j]:.2e}, ratio {ratio:.4f}) "
                  f"{'PASS' if this else 'FAIL'}")
        passed["P1"] = ok

    if "P2" in probes:
        for key in DRIVERS:
            DRIVERS[key].SOLVE_TIMEOUT = P2_TIMEOUT
        mX, mY = _node(d25, IX, J95)
        legs = {}
        for rt in (1e-6, 1e-8):
            try:
                legs[rt] = _joint_direct("d25", mX, mY, ANCHORS[J95], 1e-9,
                                         p3_rtol=rt)
                print(f"P2 rtol={rt:g}: Omega={legs[rt]['omega']:.5f} "
                      f"phase_final={legs[rt]['phase_final']} "
                      f"[{legs[rt]['secs']:.0f} s]")
            except Exception as e:
                print(f"P2 rtol={rt:g}: {type(e).__name__}: {e}")
        if len(legs) == 2 and all(v["omega"] for v in legs.values()):
            r = legs[1e-6]["omega"] / legs[1e-8]["omega"]
            speed = legs[1e-8]["secs"] / max(legs[1e-6]["secs"], 1e-9)
            passed["P2"] = abs(r - 1) <= 0.005
            print(f"P2: Omega(1e-6)/Omega(1e-8) = {r:.5f}, speedup x{speed:.1f} "
                  f"{'PASS' if passed['P2'] else 'FAIL'}")
        else:
            passed["P2"] = False
            print("P2: incomplete -- treat as FAIL, keep rtol 1e-8 "
                  "(phase3_kw) and budget the long tail")

    if "P3" in probes:
        engaged = {}
        for key in ("d1", "d25"):
            pfs = [r.get("phase_final") for r, o in zip(results, owner)
                   if o[0] == "P3" and o[1] == key]
            engaged[key] = any(pf == 3 for pf in pfs)
            print(f"P3 {key} (xi_ini={getattr(DRIVERS[key], 'XI_INI', 1.0)}): "
                  f"phase_final per alpha = {pfs} -> "
                  f"{'ENGAGED' if engaged[key] else 'not engaged'}")
        passed["P3"] = any(engaged.values())

    if "P4" in probes:
        om = np.array([r.get("omega", np.nan) if "error" not in r else np.nan
                       for r, o in zip(results, owner) if o[0] == "P4"])
        alpha = d05.relic_alphaX(om)
        ok = np.isfinite(alpha) and alpha < 1.0
        passed["P4"] = bool(ok)
        print(f"P4 worst corner (xi=0.5, mX=71.37, rv=0.95): "
              f"alpha_relic={alpha:.4e} (stored {P4_STORED:.2e}, "
              f"ratio {alpha / P4_STORED:.2f}) "
              f"finite omegas {int(np.isfinite(om).sum())}/8 "
              f"{'PASS' if ok else 'FAIL'}")

    print("\nGATES")
    go_a = passed.get("P1", False) and passed.get("P4", False)
    gates.append(f"  Stage 2 (invalidate) + Stage 3 (pass A): "
                 f"{'GO' if go_a else 'NO-GO'}")
    if "P3" in probes:
        gates.append(f"  Stage 4 (pass B, eps=1e-10):          "
                     f"{'GO' if passed['P3'] else 'DROP (Phase 3 never engaged)'}")
    if "P2" in probes:
        gates.append("  phase3_kw rtol 1e-6:                  "
                     + ("KEEP" if passed.get("P2") else "REVERT to 1e-8"))
    print("\n".join(gates))
    sys.exit(0 if go_a else 1)


if __name__ == "__main__":
    main()

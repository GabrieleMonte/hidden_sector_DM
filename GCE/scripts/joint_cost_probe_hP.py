"""Why does every Higgs-portal joint solve near the eps=1e-10 crossing cost
more than SOLVE_TIMEOUT?

The 2026-10-20 rebuild finished eps=1e-9 in 1:14 on 16 nodes (306 timeouts,
nearly all above the crossing) but reached only 167 / 2083 nodes at eps=1e-10
in 2:00, logging 9517 timeouts -- 7275 of them on the five axis points
1.4e-2 .. 1.2e-1, which is exactly where the crossings sit.  At 600 s each
that is 1586 core-hours spent on solves that returned nothing.

This times the two halves of the driver separately (solve_boltzmann_joint,
then solve_background when the run does not exit in phase 1.5S or 3) so the
cost can be attributed, and repeats each point at a looser rtol to tell
"crawling, needs tolerance" from "genuinely long integration".

    cosmo_env
    sbatch GCE/scripts/joint_cost_probe_hP.sbatch

Writes GCE/output/diag/joint_cost_probe_hP.npz and prints a table.
"""

import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                    # noqa: E402
import make_relic_grid_hP as hp                                   # noqa: E402

from hidden_sector_DM.solver import BoltzmannSolver               # noqa: E402

# ---- configuration ---------------------------------------------------

PROBE_TIMEOUT = 1800     # s per half-solve; the production cap is 600
NODES = [(17.022, 0.695), (34.0, 0.40), (57.807, 0.869), (95.0, 0.20)]
LAM = np.array([8.377e-3, 1.425e-2, 2.424e-2, 4.125e-2, 7.017e-2, 1.194e-1])
EPS = (1e-9, 1e-10)
RTOL = (1e-6, 1e-4)      # 1e-6 is SOLVER_KW's; 1e-4 asks whether it crawls

DIAG_DIR = base.OUT_DIR.parent / "diag"
OUT = DIAG_DIR / "joint_cost_probe_hP.npz"


def one_task(arg):
    """(k, eps, lam, rtol) -> timings and Omega for one half-split solve."""
    k, mX, mY, eps, lam, rtol = arg
    base.SOLVE_TIMEOUT = PROBE_TIMEOUT

    kw = dict(base.SOLVER_KW)
    kw["rtol_value"] = rtol
    phase3_kw = dict(base.PHASE3_KW)
    phase3_kw["phase3_kw"] = dict(rtol_value=rtol)

    solver = BoltzmannSolver(base.cosmo, base._make_model(mX, mY, lam))
    rec = dict(k=k, eps=eps, lam=lam, rtol=rtol, t_joint=np.nan, t_bg=np.nan,
               omega=np.nan, phase_final="", converged=False, err="")

    t0 = time.perf_counter()
    try:
        joint = base.timed_solve(
            solver.solve_boltzmann_joint, epsX=eps, return_bg_ICs=True,
            min_eps_floor_ratio=0.0, GoH_exit_thresh=base.GOH_EXIT_THRESH,
            xmax=base.XMAX_JOINT, phase2_form="Y", verbose=False,
            **kw, **phase3_kw)
    except Exception as e:
        rec["t_joint"] = time.perf_counter() - t0
        rec["err"] = f"joint:{type(e).__name__}"
        return rec
    rec["t_joint"] = time.perf_counter() - t0
    rec["phase_final"] = str(joint["phase_final"])
    rec["converged"] = bool(joint["converged"])

    if joint["phase_final"] in ("1.5S", 3):
        rec["omega"] = joint["YX_relic"] * mX * base.s0_cosmo / base.rhoc
        return rec

    t1 = time.perf_counter()
    try:
        yield_X = base.timed_solve(
            solver.solve_background, epsX=eps, bg_ICs=joint["bg_ICs"],
            correct_Y_decay=True, traj_has_decay=True, stop_on_rhoY=True,
            verbose=False)["YX_final"]
    except Exception as e:
        rec["t_bg"] = time.perf_counter() - t1
        rec["err"] = f"bg:{type(e).__name__}"
        return rec
    rec["t_bg"] = time.perf_counter() - t1
    rec["omega"] = yield_X * mX * base.s0_cosmo / base.rhoc
    return rec


def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else base.NCORES
    hp._configure(1.0)

    tasks = [(k, mX, rv * mX, eps, float(lam), rtol)
             for k, (mX, rv) in enumerate(NODES)
             for eps in EPS for lam in LAM for rtol in RTOL]
    print(f"{len(tasks)} half-split solves, cap {PROBE_TIMEOUT} s, "
          f"{ncores} workers", flush=True)

    with mp.Pool(min(ncores, len(tasks))) as pool:
        recs = list(pool.imap_unordered(one_task, tasks))

    recs.sort(key=lambda r: (r["k"], r["eps"], r["lam"], r["rtol"]))
    print(f"\n{'mX':>7} {'rv':>5} {'eps':>7} {'lam':>9} {'rtol':>6} "
          f"{'t_joint':>8} {'t_bg':>8} {'phase':>6} {'cvg':>4} "
          f"{'Omega':>10}  err")
    for r in recs:
        mX, rv = NODES[r["k"]]
        print(f"{mX:7.3f} {rv:5.2f} {r['eps']:7.0e} {r['lam']:9.3e} "
              f"{r['rtol']:6.0e} {r['t_joint']:8.1f} {r['t_bg']:8.1f} "
              f"{r['phase_final']:>6} {str(r['converged']):>5} "
              f"{r['omega']:10.4g}  {r['err']}")

    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(OUT,
             mX=np.array([NODES[r["k"]][0] for r in recs]),
             rv=np.array([NODES[r["k"]][1] for r in recs]),
             eps=np.array([r["eps"] for r in recs]),
             lam=np.array([r["lam"] for r in recs]),
             rtol=np.array([r["rtol"] for r in recs]),
             t_joint=np.array([r["t_joint"] for r in recs]),
             t_bg=np.array([r["t_bg"] for r in recs]),
             omega=np.array([r["omega"] for r in recs]),
             converged=np.array([r["converged"] for r in recs]),
             phase_final=np.array([r["phase_final"] for r in recs]),
             err=np.array([r["err"] for r in recs]))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()

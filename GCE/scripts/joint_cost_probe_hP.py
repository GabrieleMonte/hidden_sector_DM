"""Why does every Higgs-portal joint solve near the eps=1e-10 crossing cost
more than SOLVE_TIMEOUT?

The 2026-10-20 rebuild finished eps=1e-9 in 1:14 on 16 nodes (306 timeouts,
nearly all above the crossing) but reached only 167 / 2083 nodes at eps=1e-10
in 2:00, logging 9517 timeouts -- 7275 of them on the five axis points
1.4e-2 .. 1.2e-1, which is exactly where the crossings sit.  The pre-fix
September runs timed out in the same band, so this is not the phase-entry
fix; the finer axis only stopped the coarse one from hiding it.

THE SUSPECT.  SOLVER_KW carries rtol_value=1e-6 with the note "1e-8 exceeds
any timeout at rv=0.95 (>45 min); 1e-6 solves those in ~1-6 min with
identical handoffs and Omega to 0.03%".  But omega_h2_secluded is the only
driver that passes SOLVER_KW: omega_h2_joint does not, so the joint
integration runs at solve_boltzmann_joint's own default of 1e-8.

Against that: timed_solve already retries a timed-out joint solve at
RETRY_RTOL=1e-6 (1e-6 < the 1e-8 in force), so each of those 9517 lines is a
600 s run at 1e-8 AND a 600 s run at 1e-6 that both ran out of time.  So 1e-6
alone is not obviously enough, and the question is whether it is slow or
unbounded.  This probe times the two halves of the driver separately at three
tolerances with the retry disabled, and records Omega at each: a tolerance is
only usable if it leaves the answer alone.

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

PROBE_TIMEOUT = 1800     # s per half-solve; production caps the pair at 600
NODES = [(17.022, 11.834), (22.866, 15.896), (57.807, 50.229), (95.000, 19.000)]
LAM = np.array([1.425e-2, 2.424e-2, 4.125e-2, 7.017e-2])
EPS = (1e-9, 1e-10)

# None = whatever solve_boltzmann_joint defaults to (1e-8), i.e. production.
# phase3_kw stays at its production 1e-6 throughout, so the only thing varying
# is the tolerance of the joint integration itself.
RTOL = (None, 1e-6, 1e-4)

DIAG_DIR = base.OUT_DIR.parent / "diag"
OUT = DIAG_DIR / "joint_cost_probe_hP.npz"


def one_task(arg):
    """(k, eps, lam, rtol) -> timings and Omega for one half-split solve.

    Mirrors omega_h2_joint exactly -- same keywords, same PHASE3_KW, no
    SOLVER_KW -- except for the rtol_value override and the disabled retry.
    """
    k, mX, mY, eps, lam, rtol = arg
    base.SOLVE_TIMEOUT = PROBE_TIMEOUT
    base.RETRY_RTOL = 0.0                    # never looser than what is in
                                             # force, so timed_solve re-raises
    kw = {} if rtol is None else dict(rtol_value=rtol)

    solver = BoltzmannSolver(base.cosmo, base._make_model(mX, mY, lam))
    rec = dict(k=k, eps=eps, lam=lam, rtol=(np.nan if rtol is None else rtol),
               t_joint=np.nan, t_bg=np.nan, omega=np.nan, phase_final="",
               converged=False, err="")

    t0 = time.perf_counter()
    try:
        joint = base.timed_solve(
            solver.solve_boltzmann_joint, epsX=eps, return_bg_ICs=True,
            min_eps_floor_ratio=0.0, GoH_exit_thresh=base.GOH_EXIT_THRESH,
            xmax=base.XMAX_JOINT, phase2_form="Y", verbose=False,
            **kw, **base.PHASE3_KW)
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

    tasks = [(k, mX, mY, eps, float(lam), rtol)
             for k, (mX, mY) in enumerate(NODES)
             for eps in EPS for lam in LAM for rtol in RTOL]
    print(f"{len(tasks)} half-split solves, cap {PROBE_TIMEOUT} s per half, "
          f"{ncores} workers", flush=True)

    with mp.Pool(min(ncores, len(tasks))) as pool:
        recs = list(pool.imap_unordered(one_task, tasks))

    order = {None: 0, 1e-6: 1, 1e-4: 2}
    recs.sort(key=lambda r: (r["k"], -r["eps"], r["lam"],
                             0 if np.isnan(r["rtol"]) else order[r["rtol"]]))

    print(f"\n{'mX':>7} {'mY':>7} {'eps':>7} {'lam':>9} {'rtol':>7} "
          f"{'t_joint':>8} {'t_bg':>8} {'phase':>6} {'cvg':>5} "
          f"{'Omega':>11}  err")
    for r in recs:
        mX, mY = NODES[r["k"]]
        rt = "1e-08*" if np.isnan(r["rtol"]) else f"{r['rtol']:.0e}"
        print(f"{mX:7.3f} {mY:7.3f} {r['eps']:7.0e} {r['lam']:9.3e} "
              f"{rt:>7} {r['t_joint']:8.1f} {r['t_bg']:8.1f} "
              f"{r['phase_final']:>6} {str(r['converged']):>5} "
              f"{r['omega']:11.5g}  {r['err']}")
    print("\n* = production: omega_h2_joint passes no rtol_value, so the joint "
          "solve runs at\n  solve_boltzmann_joint's own 1e-8 default.  "
          "phase3_kw is 1e-6 in every row.")

    # Does a looser tolerance move the answer?  That decides whether it is
    # usable at all, independently of what it costs.
    print(f"\n{'mX':>7} {'eps':>7} {'lam':>9} {'Omega(1e-8)':>12} "
          f"{'Omega(1e-6)':>12} {'Omega(1e-4)':>12} {'1e-6/1e-8':>10} "
          f"{'1e-4/1e-8':>10}")
    by = {(r["k"], r["eps"], r["lam"],
           None if np.isnan(r["rtol"]) else r["rtol"]): r["omega"] for r in recs}
    for k, (mX, _mY) in enumerate(NODES):
        for eps in EPS:
            for lam in LAM:
                o8 = by.get((k, eps, lam, None), np.nan)
                o6 = by.get((k, eps, lam, 1e-6), np.nan)
                o4 = by.get((k, eps, lam, 1e-4), np.nan)
                print(f"{mX:7.3f} {eps:7.0e} {lam:9.3e} {o8:12.5g} "
                      f"{o6:12.5g} {o4:12.5g} {o6/o8:10.4f} {o4/o8:10.4f}")

    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(OUT,
             mX=np.array([NODES[r["k"]][0] for r in recs]),
             mY=np.array([NODES[r["k"]][1] for r in recs]),
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

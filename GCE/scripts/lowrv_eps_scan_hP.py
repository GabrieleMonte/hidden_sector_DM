"""Why does only eps=1e-9 fail in the Higgs portal's rv < 0.2 corner?

The 2026-10-23 eps=1e-9 rebuild left 146 SOLVE_FAILED nodes, every one of them
at rv < 0.2, while eps=1e-8 and 1e-7 solve that same corner with none.  The
grid keeps only Omega, so it cannot say what the solver did differently.

This replays the FULL coupling axis -- no early stop, every point -- for a
handful of those nodes at all three eps, and records per point what the joint
solver actually did: Omega, phase_final, converged, wall time, and whether the
production gate would have kept it.  Controls are nodes in the same corner
that eps=1e-9 solved, so a difference cannot be put down to rv alone.

Also reports x_dec = mX / T_dec from Gamma_Y(eps) = H(T_dec), against
XMAX_JOINT: the joint integration wall is the obvious candidate for an
eps-dependent cliff, since x_dec ~ 1/eps.

    sbatch GCE/scripts/lowrv_eps_scan_hP.sbatch

Writes GCE/output/diag/lowrv_eps_scan_hP.npz and prints a table.
"""

import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid_hP as hp                                   # noqa: E402
import make_relic_grid as base                                    # noqa: E402

EPS_SCAN = (1e-7, 1e-8, 1e-9)
N_FAIL, N_CTRL = 4, 2                     # failing nodes, and solved controls
SCAN_TIMEOUT = 1200                       # per solve; the grid used 600
GRID = Path("/global/u1/g/gab97/secluded_DM/GCE/output/relic")
OUT = Path("/global/u1/g/gab97/secluded_DM/GCE/output/diag")


def x_dec(mX, mY, eps):
    """mX / T_dec, with Gamma_Y(eps) = H(T_dec).  NaN if unbracketed."""
    gamma = base._make_model(mX, mY, 0.1).decay_width_to_SM(eps)
    if not np.isfinite(gamma) or gamma <= 0.0:
        return np.nan

    def f(lnT):
        return np.log(base.cosmo.hubble(np.exp(lnT))) - np.log(gamma)

    lo, hi = np.log(1e-6), np.log(1e3)
    if f(lo) * f(hi) > 0:
        return np.nan
    return mX / np.exp(brentq(f, lo, hi, xtol=1e-10))


def one_point(task):
    """(node, eps, lam) -> what the joint solver did, with no early stop."""
    n, e, k, mX, mY, eps, lam = task
    solver = base.BoltzmannSolver(base.cosmo, base._make_model(mX, mY, lam))

    t0 = time.time()
    omega, phase, conv, kept, err = np.nan, -1.0, False, False, ""
    try:
        joint = base.timed_solve(solver.solve_boltzmann_joint, epsX=eps,
                                 return_bg_ICs=True, min_eps_floor_ratio=0.0,
                                 GoH_exit_thresh=base.GOH_EXIT_THRESH,
                                 xmax=base.XMAX_JOINT, phase2_form="Y",
                                 verbose=False, **base.PHASE3_KW)
        pf, conv = joint["phase_final"], bool(joint["converged"])
        phase = 1.5 if pf == "1.5S" else float(pf)   # 1.5S distinguished by SM_locked
        kept = conv or pf in base.ACCEPT_UNCONVERGED_PHASES

        if pf in ("1.5S", 3):
            omega = joint["YX_relic"] * mX * base.s0_cosmo / base.rhoc
        else:
            bg = base.timed_solve(solver.solve_background, epsX=eps,
                                  bg_ICs=joint["bg_ICs"], correct_Y_decay=True,
                                  traj_has_decay=True, stop_on_rhoY=True,
                                  verbose=False)
            omega = bg["YX_final"] * mX * base.s0_cosmo / base.rhoc
    except Exception as exc:
        err = f"{type(exc).__name__}: {exc}"[:80]

    return n, e, k, omega, phase, conv, kept, time.time() - t0, err


def pick_nodes():
    """Failing eps=1e-9 nodes spread over the corner, plus solved controls."""
    d = np.load(GRID / "relic_grid_hP_eps1e-09.npz")
    mX, mY, st = d["mX_grid"], d["mY_grid"], d["status"]
    ong = np.isfinite(mX)
    rv = np.where(ong, mY / np.where(ong, mX, 1.0), np.nan)

    fail = np.argwhere(ong & (st == 1))
    ctrl = np.argwhere(ong & (st == 0) & (rv < 0.2))

    def spread(cands, n):
        r = np.array([rv[i, j] for i, j in cands])
        want = np.linspace(r.min(), r.max(), n)
        return [tuple(cands[np.argmin(np.abs(r - w))]) for w in want]

    out = [(mX[i, j], mY[i, j], "fail") for i, j in spread(fail, N_FAIL)]
    out += [(mX[i, j], mY[i, j], "ctrl") for i, j in spread(ctrl, N_CTRL)]
    return out


def main():
    hp._configure(1)
    base.SOLVE_TIMEOUT = SCAN_TIMEOUT
    base._set_coupling(1e-9)                      # the hP lam axis; same for all eps
    lam_axis = np.array(base.COUPLING)

    nodes = pick_nodes()
    print(f"XMAX_JOINT = {base.XMAX_JOINT:g}   SOLVE_TIMEOUT = {SCAN_TIMEOUT} s")
    print(f"lam axis: {lam_axis.min():.4g} .. {lam_axis.max():.4g} "
          f"({len(lam_axis)} points), full sweep, no early stop\n")

    print(f"{'node':>4} {'kind':>5} {'mX':>8} {'mY':>8} {'rv':>7} | "
          + " ".join(f"x_dec({e:.0e})".rjust(14) for e in EPS_SCAN))
    for n, (mX, mY, kind) in enumerate(nodes):
        xs = " ".join(f"{x_dec(mX, mY, e):14.4g}" for e in EPS_SCAN)
        print(f"{n:4d} {kind:>5} {mX:8.3f} {mY:8.3f} {mY / mX:7.4f} | {xs}")

    tasks = [(n, e, k, mX, mY, eps, lam)
             for n, (mX, mY, _) in enumerate(nodes)
             for e, eps in enumerate(EPS_SCAN)
             for k, lam in enumerate(lam_axis)]
    print(f"\n{len(tasks)} solves on {mp.cpu_count()} cores", flush=True)

    shape = (len(nodes), len(EPS_SCAN), len(lam_axis))
    omega = np.full(shape, np.nan)
    phase = np.full(shape, np.nan)
    conv = np.zeros(shape, bool)
    kept = np.zeros(shape, bool)
    wall = np.full(shape, np.nan)
    errs = {}

    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else mp.cpu_count()
    with mp.Pool(n_workers) as pool:
        for done, r in enumerate(pool.imap_unordered(one_point, tasks), 1):
            n, e, k, om, ph, cv, kp, dt, err = r
            omega[n, e, k], phase[n, e, k] = om, ph
            conv[n, e, k], kept[n, e, k], wall[n, e, k] = cv, kp, dt
            if err:
                errs[(n, e, k)] = err
            if done % 20 == 0:
                print(f"  {done}/{len(tasks)}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / "lowrv_eps_scan_hP.npz",
             mX=np.array([m for m, _, _ in nodes]),
             mY=np.array([m for _, m, _ in nodes]),
             kind=np.array([k for _, _, k in nodes]),
             eps=np.array(EPS_SCAN), lam=lam_axis,
             omega=omega, phase=phase, converged=conv, kept=kept, wall=wall)

    # ---- per-node, per-eps: where the axis stops being cheap ----------
    for n, (mX, mY, kind) in enumerate(nodes):
        print(f"\n=== node {n} [{kind}] mX={mX:.3f} mY={mY:.3f} "
              f"rv={mY / mX:.4f} ===")
        print(f"{'lam':>10} | " + " | ".join(
            f"{'Om':>10} {'ph':>4} {'cv':>2} {'t[s]':>7}" for _ in EPS_SCAN))
        for k in range(len(lam_axis) - 1, -1, -1):
            cells = []
            for e in range(len(EPS_SCAN)):
                flag = "" if kept[n, e, k] else "x"
                cells.append(f"{omega[n, e, k]:10.3g} {phase[n, e, k]:4.1f} "
                             f"{str(conv[n, e, k])[0]:>1}{flag:1} "
                             f"{wall[n, e, k]:7.1f}")
            print(f"{lam_axis[k]:10.4g} | " + " | ".join(cells))

    if errs:
        print(f"\n{len(errs)} exceptions:")
        for (n, e, k), m in sorted(errs.items())[:40]:
            print(f"  node {n} eps={EPS_SCAN[e]:.0e} lam={lam_axis[k]:.4g}: {m}")

    # ---- the headline comparison --------------------------------------
    print("\nwall time over the whole axis, by node and eps (s):")
    print(f"{'node':>4} {'kind':>5} {'rv':>7} | "
          + " ".join(f"{e:>12.0e}" for e in EPS_SCAN))
    for n, (mX, mY, kind) in enumerate(nodes):
        tot = " ".join(f"{np.nansum(wall[n, e]):12.0f}" for e in range(len(EPS_SCAN)))
        print(f"{n:4d} {kind:>5} {mY / mX:7.4f} | {tot}")


if __name__ == "__main__":
    main()

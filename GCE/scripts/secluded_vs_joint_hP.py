"""Is the Higgs portal already secluded at eps = 1e-9 and 1e-10?

The campaign solves eps = 1e-9 and 1e-10 with solve_boltzmann_joint (SM and
hidden sector evolved together) and eps = 1e-11, 5e-12, 1e-12 with
solve_boltzmann_chempot_3phase + solve_background (hidden freeze-out first,
the SM background re-run per eps).  If the portal is already decoupled at
1e-9, the two routes must give the same Omega, and the joint route is then
just a much more expensive way to get the same number.

That matters for more than tidiness.  Gamma_phi ~ eps^2, so x_dec = mX/T_dec
scales as 1/eps: at eps=1e-9 the decay sits at x = 745-2594, inside
XMAX_JOINT = 5e3, and the joint solve exits on it in seconds.  At eps=1e-10
it sits at x = 7.4e3-2.6e4, BEYOND the wall, so the solve grinds to xmax at
rtol 1e-8 and 9517 of them exceeded the 600 s cap in the 2026-10-20 rebuild.
If the secluded route reproduces eps=1e-9, it also computes eps=1e-10 for
free.

Two comparisons per node, at the SAME lam axis:
  - crossing lam from the secluded route at eps = 1e-9 and 1e-10
  - crossing lam from the joint route at eps = 1e-9 (1e-10 is unaffordable)
plus a point check of both Omegas at the lam the delivered eps=1e-9 grid
already reports, where the joint answer is 0.12 by construction.

    cosmo_env
    sbatch GCE/scripts/secluded_vs_joint_hP.sbatch

Writes GCE/output/diag/secluded_vs_joint_hP.npz.
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

# ---- configuration ---------------------------------------------------

EPS_TEST = (1e-9, 1e-10)
PROBE_TIMEOUT = 1800     # s per solve; production caps at 600
N_NODES = 12             # spread over the delivered eps=1e-9 grid

GRID_1E9 = base.OUT_DIR / "relic_grid_hP_eps1e-09.npz"
DIAG_DIR = base.OUT_DIR.parent / "diag"
OUT = DIAG_DIR / "secluded_vs_joint_hP.npz"


def pick_nodes():
    """(mX, mY, lam_grid) at N_NODES solved points spread over the grid."""
    d = np.load(GRID_1E9, allow_pickle=True)
    ok = np.flatnonzero((d["status"] == 0).ravel())
    take = ok[np.linspace(0, ok.size - 1, N_NODES).astype(int)]
    mX = d["mX_grid"].ravel()[take]
    mY = d["mY_grid"].ravel()[take]
    lam = d["coupling_relic"].ravel()[take]
    return [(float(a), float(b), float(c)) for a, b, c in zip(mX, mY, lam)]


def one_task(arg):
    """(k, kind) -> crossing lam per eps, and Omega at the grid's own lam."""
    k, mX, mY, lam_grid, kind = arg
    base.SOLVE_TIMEOUT = PROBE_TIMEOUT
    base.EPS_SECLUDED = EPS_TEST          # re-aim the secluded driver
    t0 = time.perf_counter()

    if kind == "secluded":
        omega = base._sweep(lambda a: base.omega_h2_secluded(mX, mY, a),
                            len(EPS_TEST))
        cross = [base._crossing(omega[:, e])[0] for e in range(len(EPS_TEST))]
        point = base.omega_h2_secluded(mX, mY, lam_grid)
    else:
        eps = float(kind)
        omega = base._sweep(lambda a: base.omega_h2_joint(mX, mY, a, eps), 1)
        cross = [base._crossing(omega[:, 0])[0], np.nan]
        point = [base.omega_h2_joint(mX, mY, lam_grid, eps), np.nan]

    return k, kind, np.array(cross, float), np.array(point, float), \
        time.perf_counter() - t0


def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else base.NCORES
    hp._configure(1.0)

    nodes = pick_nodes()
    kinds = ["secluded", "1e-09"]
    tasks = [(k, mX, mY, lam, kind)
             for k, (mX, mY, lam) in enumerate(nodes) for kind in kinds]
    print(f"{len(nodes)} nodes x {len(kinds)} routes = {len(tasks)} sweeps, "
          f"{ncores} workers, cap {PROBE_TIMEOUT} s", flush=True)
    print(f"lam axis: {base.N_COUPLING} points "
          f"[{base.COUPLING[0]:.3g}, {base.COUPLING[-1]:.3g}]", flush=True)

    cross = {kd: np.full((len(nodes), len(EPS_TEST)), np.nan) for kd in kinds}
    point = {kd: np.full((len(nodes), len(EPS_TEST)), np.nan) for kd in kinds}
    secs = {kd: np.full(len(nodes), np.nan) for kd in kinds}

    with mp.Pool(min(ncores, len(tasks))) as pool:
        for k, kd, c, p, dt in pool.imap_unordered(one_task, tasks):
            cross[kd][k] = c
            point[kd][k] = p
            secs[kd][k] = dt

    # ---- Omega at the lam the delivered grid already reports -------
    print(f"\nOmega at the grid's own eps=1e-9 lam (joint should be ~0.12)")
    print(f"{'mX':>7} {'mY':>7} {'lam_grid':>10} {'Om_joint':>10} "
          f"{'Om_sec':>10} {'sec/joint':>10}")
    for k, (mX, mY, lam) in enumerate(nodes):
        oj, os_ = point["1e-09"][k, 0], point["secluded"][k, 0]
        print(f"{mX:7.3f} {mY:7.3f} {lam:10.4g} {oj:10.4g} {os_:10.4g} "
              f"{os_ / oj:10.4f}")

    # ---- the deliverable: crossing lam, both routes ----------------
    print(f"\ncrossing lam at eps=1e-9       "
          f"(t: secluded vs joint, seconds)")
    print(f"{'mX':>7} {'mY':>7} {'lam_grid':>10} {'lam_joint':>10} "
          f"{'lam_sec':>10} {'sec/joint':>10} {'t_sec':>8} {'t_joint':>8}")
    for k, (mX, mY, lam) in enumerate(nodes):
        lj, ls = cross["1e-09"][k, 0], cross["secluded"][k, 0]
        print(f"{mX:7.3f} {mY:7.3f} {lam:10.4g} {lj:10.4g} {ls:10.4g} "
              f"{ls / lj:10.4f} {secs['secluded'][k]:8.1f} "
              f"{secs['1e-09'][k]:8.1f}")

    r = cross["secluded"][:, 0] / cross["1e-09"][:, 0]
    good = np.isfinite(r)
    if good.any():
        print(f"\nlam_secluded / lam_joint at eps=1e-9: "
              f"min {r[good].min():.4f}  median {np.median(r[good]):.4f}  "
              f"max {r[good].max():.4f}  ({good.sum()}/{len(nodes)} nodes)")

    print(f"\ncrossing lam at eps=1e-10, secluded route only "
          f"(the joint route cannot finish)")
    print(f"{'mX':>7} {'mY':>7} {'lam_sec(1e-10)':>15}")
    for k, (mX, mY, _l) in enumerate(nodes):
        print(f"{mX:7.3f} {mY:7.3f} {cross['secluded'][k, 1]:15.4g}")

    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, eps=np.array(EPS_TEST),
             mX=np.array([n[0] for n in nodes]),
             mY=np.array([n[1] for n in nodes]),
             lam_grid=np.array([n[2] for n in nodes]),
             cross_secluded=cross["secluded"], cross_joint=cross["1e-09"],
             point_secluded=point["secluded"], point_joint=point["1e-09"],
             secs_secluded=secs["secluded"], secs_joint=secs["1e-09"],
             lam_axis=base.COUPLING)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()

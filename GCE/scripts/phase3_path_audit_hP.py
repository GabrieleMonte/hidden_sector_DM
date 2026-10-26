"""How much does the Higgs-portal grid move if every joint solve must end in Phase 3?

`make_relic_grid.py` sends any joint run that does not end in Phase 3 or SM-locked
(1.5S) through `solve_background`.  The vector-portal drivers instead re-solve it
with `Sdot_switch_thresh` lowered down PHASE3_RETRY_SDOT first (PHASE3_REQUIRE),
and keep the background value only when Phase 3 still never engages.  The two
routes do not agree, so the hP grids may carry a mix.

This replays whole coupling axes on the production path and on the retry ladder,
node by node, and reports the exit each lam sample took and the crossing each
route gives.  It writes nothing into the grids.

    python GCE/scripts/phase3_path_audit_hP.py --eps 1e-7 -n 32
"""

import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_relic_grid as base
import make_relic_grid_hP as hp

from hidden_sector_DM.HiddenSectorDM import BoltzmannSolver

# ---- configuration -------------------------------------------------------

RETRY_SDOT = (1e-4, 1e-5)       # the vP drivers' PHASE3_RETRY_SDOT, verbatim
N_NODES = 24                    # sampled nodes per eps, spread over the grid
SEED = 12345
OUT_DIR = ROOT / "GCE" / "output" / "relic"

PATH_P3, PATH_SMLOCK, PATH_LEGACY, PATH_FAILED = 3, 2, 1, 0
PATH_NAME = {PATH_P3: "phase3", PATH_SMLOCK: "1.5S", PATH_LEGACY: "legacy",
             PATH_FAILED: "failed"}


def _omega_from(joint, solver, mX, eps):
    """Omega h^2 from a finished joint solve, by whichever route its exit implies."""
    if joint["phase_final"] in ("1.5S", 3):
        return joint["YX_relic"] * mX * base.s0_cosmo / base.rhoc
    yield_X = base.timed_solve(solver.solve_background, epsX=eps,
                               bg_ICs=joint["bg_ICs"], correct_Y_decay=True,
                               traj_has_decay=True, stop_on_rhoY=True,
                               verbose=False)["YX_final"]
    return yield_X * mX * base.s0_cosmo / base.rhoc


def _solve(solver, eps, **kw):
    return base.timed_solve(solver.solve_boltzmann_joint, epsX=eps,
                            return_bg_ICs=True, min_eps_floor_ratio=0.0,
                            GoH_exit_thresh=base.GOH_EXIT_THRESH,
                            xmax=base.XMAX_JOINT, phase2_form="Y",
                            verbose=False, **kw)


def evaluate(mX, mY, lam, eps):
    """(Omega_production, Omega_phase3_required, path_production, path_required).

    The production value is exactly what `omega_h2_joint` would store.  The
    required value re-solves a legacy exit down RETRY_SDOT and falls back to the
    production value when Phase 3 still never engages, as PHASE3_REQUIRE does.
    """
    solver = BoltzmannSolver(base.cosmo, base._make_model(mX, mY, lam))
    try:
        first = _solve(solver, eps, **base.PHASE3_KW)
    except Exception:
        return np.nan, np.nan, PATH_FAILED, PATH_FAILED

    gate = (base.REQUIRE_CONVERGED and not first["converged"]
            and first["phase_final"] not in base.ACCEPT_UNCONVERGED_PHASES)
    if gate:
        return np.nan, np.nan, PATH_FAILED, PATH_FAILED

    pf = first["phase_final"]
    path = PATH_P3 if pf == 3 else PATH_SMLOCK if pf == "1.5S" else PATH_LEGACY
    try:
        om_prod = _omega_from(first, solver, mX, eps)
    except Exception:
        return np.nan, np.nan, PATH_FAILED, PATH_FAILED

    if path != PATH_LEGACY:                  # nothing for the ladder to improve
        return om_prod, om_prod, path, path

    for sdot in RETRY_SDOT:
        try:
            j = _solve(solver, eps, **{**base.PHASE3_KW, "Sdot_switch_thresh": sdot})
        except Exception:
            continue
        if j["phase_final"] == 3:
            om = j["YX_relic"] * mX * base.s0_cosmo / base.rhoc
            return om_prod, om, PATH_LEGACY, PATH_P3
    return om_prod, om_prod, PATH_LEGACY, PATH_LEGACY


def audit_node(task):
    """One (i, j): both crossings, and the exit every sampled lam took."""
    i, j, mX, mY, eps = task
    t0 = time.time()
    cols = base._sweep(lambda lam: np.array(evaluate(mX, mY, lam, eps)[:2]), 2)
    om_prod, om_req = cols[:, 0], cols[:, 1]
    both = np.isfinite(om_prod) & np.isfinite(om_req)
    moved = both & ~np.isclose(om_prod, om_req, rtol=1e-9, atol=0.0)
    lam_prod, st_prod = base._crossing(om_prod)
    lam_req, st_req = base._crossing(om_req)
    return dict(i=i, j=j, mX=mX, mY=mY, rv=mY / mX, eps=eps,
                lam_prod=lam_prod, lam_req=lam_req, st_prod=st_prod,
                st_req=st_req, n_moved=int(moved.sum()),
                n_solved=int(np.isfinite(om_prod).sum()), secs=time.time() - t0)


def pick_nodes(eps, n_nodes, rng):
    """`n_nodes` solved grid nodes, spread over rv, plus the reported bad node."""
    grid = np.load(OUT_DIR / f"relic_grid_hP_eps{eps:.0e}.npz")
    mX, mY, st = grid["mX_grid"], grid["mY_grid"], grid["status"]
    ok = np.argwhere((st == 0) & np.isfinite(mX))
    rv = np.array([mY[a, b] / mX[a, b] for a, b in ok])
    order = np.argsort(rv)
    take = order[np.linspace(0, len(order) - 1, n_nodes).astype(int)]
    sel = [tuple(ok[t]) for t in take]
    # the node the 2026-10-26 report names: m_DM = 81, r = 0.727
    d = np.abs(np.where(np.isfinite(mX), mX, np.inf) - 81.0) + np.abs(rv_grid(mX, mY) - 0.727)
    sel.insert(0, tuple(np.unravel_index(np.argmin(d), mX.shape)))
    out, seen = [], set()
    for i, j in sel:
        if (i, j) in seen:
            continue
        seen.add((i, j))
        out.append((int(i), int(j), float(mX[i, j]), float(mY[i, j]), eps))
    return out


def rv_grid(mX, mY):
    with np.errstate(invalid="ignore"):
        return np.where(np.isfinite(mX), mY / mX, np.inf)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=1e-7)
    ap.add_argument("-n", "--ncores", type=int, default=32)
    ap.add_argument("--nodes", type=int, default=N_NODES)
    args = ap.parse_args()

    hp._configure(1.0)
    base._set_coupling(args.eps)

    rng = np.random.default_rng(SEED)
    tasks = pick_nodes(args.eps, args.nodes, rng)
    print(f"eps={args.eps:g}  {len(tasks)} nodes  {base.N_COUPLING} lam points  "
          f"{args.ncores} cores", flush=True)

    rows = []
    with mp.Pool(args.ncores) as pool:
        for r in pool.imap_unordered(audit_node, tasks):
            rows.append(r)
            ratio = r["lam_req"] / r["lam_prod"] if np.isfinite(r["lam_prod"]) else np.nan
            print(f"  mX={r['mX']:7.3f} rv={r['rv']:.3f}  legacy-moved {r['n_moved']:2d}"
                  f"/{r['n_solved']:2d}  lam {r['lam_prod']:.5g} -> {r['lam_req']:.5g}"
                  f"  x{ratio:.4f}  [{r['secs']:.0f} s]", flush=True)

    fin = [r for r in rows if np.isfinite(r["lam_prod"]) and np.isfinite(r["lam_req"])]
    rat = np.array([r["lam_req"] / r["lam_prod"] for r in fin])
    moved = np.array([r["n_moved"] for r in rows])
    print(f"\nnodes {len(rows)}, both crossings finite {len(fin)}")
    print(f"lam samples whose Omega moved: {moved.sum()} over "
          f"{sum(r['n_solved'] for r in rows)} solved")
    print(f"nodes with >=1 moved sample: {int((moved > 0).sum())}")
    if len(fin):
        print(f"lam_required / lam_production: median {np.median(rat):.4f}  "
              f"min {rat.min():.4f}  max {rat.max():.4f}")
        print(f"  |change| > 1%: {int((np.abs(rat - 1) > 0.01).sum())} of {len(fin)}")
        print(f"  |change| > 5%: {int((np.abs(rat - 1) > 0.05).sum())} of {len(fin)}")

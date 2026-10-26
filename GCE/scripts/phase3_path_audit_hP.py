"""How much does the Higgs-portal relic grid move if every joint solve must end in Phase 3?

make_relic_grid.py used to send any joint run that did not end in Phase 3 or SM-locked
(1.5S) straight through solve_background.  The vector-portal drivers instead re-solve it
with Sdot_switch_thresh lowered down a retry ladder (PHASE3_REQUIRE) and keep the
background value only when Phase 3 still never engages.  The two routes do not agree, so
a grid built without the ladder carries a mix of them.

Stages:
  1. Pick N_NODES solved nodes per eps, spread over the mass ratio rv = mY/mX, with the
     node named in the 2026-10-26 report first.
  2. Replay each node's whole coupling axis twice, once per route, reusing the driver's
     own sweep and crossing code so the production value is exactly what it would store.
  3. Report the exit every lam sample took and the crossing each route gives.

Nothing is written into the grids.

Output: stdout only; redirect the job log to keep it.
"""

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

# ---- configuration --------------------------------------------------------

EPS_VALUES = (1e-7, 1e-8)     # the joint outputs; 1e-9 and 1e-10 are secluded, so immune
RETRY_SDOT = (1e-4, 1e-5)     # the vP drivers' PHASE3_RETRY_SDOT, verbatim
N_NODES = 24                  # sampled nodes per eps
NCORES = 32                   # worker processes
XI_INI = 1.0                  # dark-to-SM temperature ratio the grids were built at
REPORT_NODE = (81.0, 0.727)   # (mX [GeV], rv) the 2026-10-26 report names; audited first
OUT_DIR = ROOT/'GCE'/'output'/'relic'


def omega_pair(mX, mY, lam, eps):
    """(Omega_production, Omega_phase3_required) at one coupling [GeV masses, lam free].

    The production value is exactly what omega_h2_joint stored before the ladder existed.
    The required value re-solves a legacy exit down RETRY_SDOT and falls back to the
    production value when Phase 3 still never engages.  NaN where the solve fails or the
    driver's convergence gate rejects it.
    """
    solver = BoltzmannSolver(base.cosmo, base._make_model(mX, mY, lam))

    def solve(kw):
        return base.timed_solve(
            solver.solve_boltzmann_joint, epsX=eps, return_bg_ICs=True,
            min_eps_floor_ratio=0.0, GoH_exit_thresh=base.GOH_EXIT_THRESH,
            xmax=base.XMAX_JOINT, phase2_form='Y', verbose=False, **kw)

    try:
        first = solve(base.PHASE3_KW)
    except Exception:
        return np.nan, np.nan

    rejected = (base.REQUIRE_CONVERGED and not first['converged']
                and first['phase_final'] not in base.ACCEPT_UNCONVERGED_PHASES)
    if rejected:
        return np.nan, np.nan

    try:
        if first['phase_final'] in ('1.5S', 3):
            om_prod = first['YX_relic']*mX*base.s0_cosmo/base.rhoc
            return om_prod, om_prod                 # nothing for the ladder to improve
        yield_X = base.timed_solve(solver.solve_background, epsX=eps,
                                   bg_ICs=first['bg_ICs'], correct_Y_decay=True,
                                   traj_has_decay=True, stop_on_rhoY=True,
                                   verbose=False)['YX_final']
        om_prod = yield_X*mX*base.s0_cosmo/base.rhoc
    except Exception:
        return np.nan, np.nan

    for sdot in RETRY_SDOT:
        try:
            retry = solve({**base.PHASE3_KW, 'Sdot_switch_thresh': sdot})
        except Exception:
            continue                                # a failed retry = did not engage
        if retry['phase_final'] == 3:
            return om_prod, retry['YX_relic']*mX*base.s0_cosmo/base.rhoc
    return om_prod, om_prod


def audit_node(task):
    """One grid node -> both crossings and how many samples changed route."""
    i, j, mX, mY, eps = task
    t0 = time.time()

    omegas = base._sweep(lambda lam: np.array(omega_pair(mX, mY, lam, eps)), 2)
    om_prod, om_req = omegas[:, 0], omegas[:, 1]

    both = np.isfinite(om_prod) & np.isfinite(om_req)
    moved = both & ~np.isclose(om_prod, om_req, rtol=1e-9, atol=0.0)

    lam_prod, status_prod = base._crossing(om_prod)
    lam_req, status_req = base._crossing(om_req)

    return dict(i=i, j=j, mX=mX, mY=mY, rv=mY/mX, eps=eps,
                lam_prod=lam_prod, lam_req=lam_req,
                status_prod=status_prod, status_req=status_req,
                n_moved=int(moved.sum()), n_solved=int(np.isfinite(om_prod).sum()),
                secs=time.time() - t0)


def pick_nodes(eps, n_nodes):
    """n_nodes solved nodes spread over rv, with REPORT_NODE first."""
    grid = np.load(OUT_DIR/f'relic_grid_hP_eps{eps:.0e}.npz')
    mX, mY, status = grid['mX_grid'], grid['mY_grid'], grid['status']

    on_grid = np.isfinite(mX)
    with np.errstate(invalid='ignore'):
        rv_grid = np.where(on_grid, mY/mX, np.inf)

    solved = np.argwhere((status == 0) & on_grid)
    rvs = np.array([rv_grid[a, b] for a, b in solved])
    spread = np.argsort(rvs)[np.linspace(0, len(rvs) - 1, n_nodes).astype(int)]

    mX_v, rv_v = REPORT_NODE
    distance = np.abs(np.where(on_grid, mX, np.inf) - mX_v) + np.abs(rv_grid - rv_v)
    picked = [np.unravel_index(np.argmin(distance), mX.shape)]
    picked += [tuple(solved[s]) for s in spread]

    tasks, seen = [], set()
    for i, j in picked:
        if (i, j) in seen:
            continue
        seen.add((i, j))
        tasks.append((int(i), int(j), float(mX[i, j]), float(mY[i, j]), eps))
    return tasks


# ---- run ------------------------------------------------------------------

hp._configure(XI_INI)

for eps in EPS_VALUES:
    base._set_coupling(eps)
    tasks = pick_nodes(eps, N_NODES)
    print(f'\neps={eps:g}  {len(tasks)} nodes  {base.N_COUPLING} lam points  '
          f'{NCORES} cores', flush=True)

    rows = []
    with mp.Pool(NCORES) as pool:
        for row in pool.imap_unordered(audit_node, tasks):
            rows.append(row)
            ratio = row['lam_req']/row['lam_prod'] if np.isfinite(row['lam_prod']) else np.nan
            print(f"  mX={row['mX']:7.3f} rv={row['rv']:.3f}  moved {row['n_moved']:2d}"
                  f"/{row['n_solved']:2d}  lam {row['lam_prod']:.5g} -> {row['lam_req']:.5g}"
                  f"  x{ratio:.4f}  [{row['secs']:.0f} s]", flush=True)

    finite = [r for r in rows if np.isfinite(r['lam_prod']) and np.isfinite(r['lam_req'])]
    ratios = np.array([r['lam_req']/r['lam_prod'] for r in finite])
    moved = np.array([r['n_moved'] for r in rows])

    print(f'\nnodes {len(rows)}, both crossings finite {len(finite)}')
    print(f'lam samples whose Omega moved: {moved.sum()} over '
          f"{sum(r['n_solved'] for r in rows)} solved")
    print(f'nodes with at least one moved sample: {int((moved > 0).sum())}')
    if len(finite):
        print(f'lam_required/lam_production: median {np.median(ratios):.4f}  '
              f'min {ratios.min():.4f}  max {ratios.max():.4f}')
        print(f'  |change| > 1%: {int((np.abs(ratios - 1) > 0.01).sum())} of {len(finite)}')
        print(f'  |change| > 5%: {int((np.abs(ratios - 1) > 0.05).sum())} of {len(finite)}')

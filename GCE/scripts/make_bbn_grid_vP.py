"""
Minimum eps that satisfies the BBN hadronic bound (KKMT limit on m_Y*Y_Y vs
lifetime) at each (mX, mY, eps) already in GCE/output/relic/.  A benchmark eps
is BBN-excluded at a point where its own eps_bbn_min there exceeds that eps.

    cosmo_env
    cd hidden_sector_DM
    python GCE/scripts/make_bbn_grid_vP.py

alpha_X is already known per (mX, mY, eps) from make_relic_grid_vP.py, so this
is one freeze-out solve per grid node per eps (~640 solves total) feeding
find_epsilon_BBN_hadronic -- hours, resumable, writes GCE/output/relic/
bbn_grid_vP.npz after every mX column.
"""

import signal
import sys
from pathlib import Path
from tqdm import trange
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hidden_sector_DM.HiddenSectorDM import BoltzmannSolver, Cosmology, VectorPortal

RELIC_DIR = Path(__file__).resolve().parents[1] / "output" / "relic"
OUT_PATH  = RELIC_DIR / "bbn_grid_vP.npz"

EPS = (1e-9, 1e-10, 1e-11, 5e-12, 1e-12)
EPS_JOINT = (1e-9, 1e-10)                 # these use solve_boltzmann_joint
SOLVE_TIMEOUT = 240                       # s; a solve past this is abandoned
SOLVER_KW = dict(cannibal_switch_full=1, convergence_threshold=1e-2)

cosmo = Cosmology()

with np.load(RELIC_DIR / "relic_grid_vP_eps1e-09.npz") as d:
    log_mX, log_rv = d["log_mX"], d["log_rv"]
    mX_grid, mY_grid = d["mX_grid"], d["mY_grid"]
    alphaX_1e9 = d["alpha_relic"]
with np.load(RELIC_DIR / "relic_grid_vP_eps1e-10.npz") as d:
    alphaX_1e10 = d["alpha_relic"]
with np.load(RELIC_DIR / "relic_grid_vP_eps_secluded.npz") as d:
    alphaX_secluded = d["alpha_relic"]    # (3, N_MX, N_RV): eps = 1e-11, 5e-12, 1e-12

ALPHAX = {1e-9: alphaX_1e9, 1e-10: alphaX_1e10, 1e-11: alphaX_secluded[0],
          5e-12: alphaX_secluded[1], 1e-12: alphaX_secluded[2]}
N_MX, N_RV = mX_grid.shape


def _raise_timeout(signum, frame):
    raise TimeoutError


signal.signal(signal.SIGALRM, _raise_timeout)


def yy_final(mX, mY, alphaX, eps):
    """Freeze-out solve at the already-known relic alpha_X; (solver, YY_final)."""
    solver = BoltzmannSolver(cosmo, VectorPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=alphaX,
        include_antiparticlesX=True, include_antiparticlesY=False))
    signal.alarm(SOLVE_TIMEOUT)
    try:
        if eps in EPS_JOINT:
            sol = solver.solve_boltzmann_joint(epsX=eps, return_bg_ICs=True,
                                               min_eps_floor_ratio=0.0,
                                               phase2_form="Y", verbose=False)
        else:
            sol = solver.solve_boltzmann_chempot_3phase(return_bg_ICs=True, verbose=False,
                                                        **SOLVER_KW)
    finally:
        signal.alarm(0)
    return solver, sol["YY_final"]


eps_bbn_min = np.full((len(EPS), N_MX, N_RV), np.nan)
start = 0
if OUT_PATH.exists():
    with np.load(OUT_PATH) as saved:
        if np.array_equal(saved["log_mX"], log_mX) and np.array_equal(saved["log_rv"], log_rv):
            eps_bbn_min = saved["eps_bbn_min"].copy()
            start = int(saved["done"])

for i in trange(start, N_MX):
    print(f"-- column {i + 1}/{N_MX}, mX={mX_grid[i, 0]:.1f} GeV --", flush=True)
    for j in range(N_RV):
        mX, mY = mX_grid[i, j], mY_grid[i, j]
        for e, eps in enumerate(EPS):
            alphaX = ALPHAX[eps][i, j]
            if not np.isfinite(alphaX):
                continue
            try:
                solver, YY = yy_final(mX, mY, alphaX, eps)
                eps_bbn_min[e, i, j] = solver.find_epsilon_BBN_hadronic(
                    YY, log10eps_range=(-20, -8.5), root_rtol=1e-3)
            except Exception as exc:
                print(f"  mY={mY:.1f} eps={eps:.0e}  FAILED: {exc}", flush=True)
    np.savez(OUT_PATH, log_mX=log_mX, log_rv=log_rv, mX_grid=mX_grid, mY_grid=mY_grid,
             eps=np.array(EPS), eps_bbn_min=eps_bbn_min, done=i + 1)

n_valid = int(np.isfinite(mX_grid).sum())
print(f"eps_bbn_min: {np.isfinite(eps_bbn_min).sum()}/{n_valid * len(EPS)} solved -> {OUT_PATH}")

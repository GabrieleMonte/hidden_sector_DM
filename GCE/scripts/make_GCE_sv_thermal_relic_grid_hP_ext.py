"""Thermal (SM-locked) relic coupling and <sigma v> for the Higgs portal over
the mX = 100-130 GeV extension grid.

The extension half of make_GCE_sv_thermal_relic_grids.py, which stops at the
main grid's mX = 100 GeV edge.  The Higgs portal only: the vector portals' GCE
contours close inside the main grid, so nothing above 100 GeV is needed for
them.  Nothing in GCE_sv_thermal_relic_grid_{hP,vP,bP,blP}.npz is read or
rewritten; this writes its own file.

    cosmo_env
    python GCE/scripts/make_GCE_sv_thermal_relic_grid_hP_ext.py

Output: GCE/output/sv_relic_fits/GCE_sv_thermal_relic_grid_hP.npz, APPENDED to
rather than written beside, so one npz is the portal's whole grid.  The main driver's
16384 points keep their order and the 578 extension points follow them:
    mX, mY                   [GeV]    the points, the extension in the order
                                      make_channel_spectra.extension_grid gives them
    coupling_relic                    lam giving Omega_c h^2 = OCH2_TARGET
    sigmav_relic             [cm^3/s] the s-wave <sigma v> there
    mX_coarse_ext, mY_coarse_ext      this driver's own coarse mesh, kept beside the
    coupling_coarse_ext               main driver's under a _ext suffix

The target points are imported from make_channel_spectra rather than rebuilt,
so the thermal relic lands on exactly the (mX, mY) the spectra cache holds.

The relic root-find below mirrors make_GCE_sv_thermal_relic_grids.py's: same
solver, same target, same tolerance, same scaled-guess / log-log secant /
brentq ladder.  The two must be changed together.

COARSE MESH
-----------
As in the main driver the coupling is solved on a coarse (log10 mX, log10 rv)
mesh and log10 of it is interpolated onto the target points.  The mesh is
padded by one step on either side in mX and one step below in rv, so that
every target point falls inside the interpolation hull rather than on its
boundary: the extension's masses are rounded to 2 decimals, and an unpadded
hull loses edge points to that round-off.

rv is capped at RV_MAX because the solve fails at mY = mX.  The extension's top
row is rv = 1 exactly, where the s-wave closes and <sigma v> = 0 whatever the
coupling, so nothing there needs one.

Above rv ~ 0.93 the rows are placed at a uniform step in log10 beta, beta =
sqrt(1 - rv^2), rather than in log10 rv, for the reason given on coarse_rv_rows in
make_GCE_sv_thermal_relic_grids.py: the s-wave closes as beta -> 0, and a ladder
uniform in log10 rv leaves one interval spanning several fine rows.  The two
drivers' meshes must be changed together, as their root-finds must.
"""

import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import griddata
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from hidden_sector_DM.HiddenSectorDM import (                  # noqa: E402
    Cosmology, BoltzmannSolver, HiggsPortal, s0_cosmo, rhoc,
)
from GCE.spectrum import GEVM2_TO_CM3S                         # noqa: E402

from make_channel_spectra import N_MAIN, extension_grid        # noqa: E402

# ---- configuration -------------------------------------------------------

OCH2_TARGET = 0.12
RV_MAX      = 0.995          # the solve fails at mY = mX
RELIC_RTOL  = 3e-3
NCORES      = len(os.sched_getaffinity(0))   # as in the main driver: one node per task

N_RV_BETA   = 4              # top rv rows placed by beta, not by log10 rv
N_MX_COARSE = 8              # coarse mesh before padding; the extension spans only
N_RV_COARSE = 8              # 0.1139 dex in mX and 0.2218 dex in rv, so 8 x 8 keeps
                             # the main driver's coarse density (0.0183 / 0.0394 dex)

COUPLING_RANGE = (2e-5, 1.0)   # lam, the axis make_relic_grid_hP.py samples
COUPLING_GUESS = 0.1
OMEGA_POWER    = 4             # sigmav ~ lam_s^2 lam_p^2 = lam^4, so Omega ~ lam^-4

GX_MAJORANA = 2              # Majorana fermion DM: 2 spin states, no antiparticle
GY_REAL_SCALAR = 1           # real scalar mediator: one degree of freedom

OUT_DIR = ROOT / 'GCE' / 'output' / 'sv_relic_fits'
OUT_NAME = 'GCE_sv_thermal_relic_grid_hP.npz'   # appended to, not written beside

cosmo = Cosmology()

# ---- the relic coupling at one point -------------------------------------

def _model(mX, mY, lam):
    """The Higgs-portal model at one point.  mX, mY in GeV; lam dimensionless."""
    return HiggsPortal(mX=mX, mY=mY, gX=GX_MAJORANA, gY=GY_REAL_SCALAR, lam=lam,
                       include_antiparticlesX=False, include_antiparticlesY=False)


def _och2(mX, mY, lam):
    """Omega_c h^2 from the SM-thermalised (hidden sector locked to the SM) solve."""
    sol = BoltzmannSolver(cosmo, _model(mX, mY, lam)).solve_boltzmann_thermSM(
        verbose=False, include_hidden_in_gstar=False)
    return sol['YX_relic']*mX*s0_cosmo/rhoc


def _coupling_relic(mX, mY):
    """lam for Omega h^2 = OCH2_TARGET (nan if outside COUPLING_RANGE);
    scaled guess + log-log secant, brentq fallback."""
    a_lo, a_hi = COUPLING_RANGE
    a0 = COUPLING_GUESS
    o0 = _och2(mX, mY, a0)
    a1 = float(np.clip(a0*(o0/OCH2_TARGET)**(1.0/OMEGA_POWER), a_lo, a_hi))
    o1 = _och2(mX, mY, a1)
    if abs(o1/OCH2_TARGET - 1.0) < RELIC_RTOL:
        return a1

    la = (np.log(a1) + (np.log(OCH2_TARGET) - np.log(o1))
          *(np.log(a1) - np.log(a0))/(np.log(o1) - np.log(o0)))
    a2 = float(np.clip(np.exp(la), a_lo, a_hi))
    if abs(_och2(mX, mY, a2)/OCH2_TARGET - 1.0) < RELIC_RTOL:
        return a2

    f = lambda l: np.log(_och2(mX, mY, np.exp(l))/OCH2_TARGET)
    for lo, hi in ((a2/3.0, a2*3.0), (a_lo, a_hi)):
        try:
            return float(np.exp(brentq(f, np.log(max(lo, a_lo)),
                                       np.log(min(hi, a_hi)), rtol=1e-3)))
        except ValueError:
            continue
    return np.nan


# ---- the extension grid --------------------------------------------------

MX_GRID, MY_GRID = extension_grid()          # (2, 578), rounded to 2 decimals
n_pts = MX_GRID.size
RV_GRID = MY_GRID/MX_GRID

lmx = np.log10(MX_GRID)
lrv = np.log10(RV_GRID)

# One padding step on either side in mX and one below in rv; see COARSE MESH.
dmx = (lmx.max() - lmx.min())/(N_MX_COARSE - 1)
drv = (lrv.max() - lrv.min())/(N_RV_COARSE - 1)
def coarse_rv_rows(lrv_lo, rv_max=RV_MAX, n_rv=N_RV_COARSE + 1, n_beta=N_RV_BETA):
    """log10 of the coarse mY/mX rows, dense where the s-wave closes.

    make_GCE_sv_thermal_relic_grids.coarse_rv_rows, on this driver's own ladder:
    the uniform log10 rv rows up to the second-to-last, then `n_beta` rows at a
    uniform step in log10 beta.  See that one for why.
    """
    ladder = np.linspace(lrv_lo, np.log10(rv_max), n_rv)
    beta = np.sqrt(1.0 - (10.0 ** ladder[-2:]) ** 2)
    log_beta = np.linspace(np.log10(beta[0]), np.log10(beta[1]), n_beta + 1)[1:]
    return np.concatenate([ladder[:-1], 0.5 * np.log10(1.0 - (10.0 ** log_beta) ** 2)])


lmx_c = np.linspace(lmx.min() - dmx, lmx.max() + dmx, N_MX_COARSE + 2)
lrv_c = coarse_rv_rows(lrv.min() - drv)

LMX_C, LRV_C = np.meshgrid(lmx_c, lrv_c, indexing='ij')

print(f'{n_pts} extension points, mX {MX_GRID.min():.2f}-{MX_GRID.max():.2f} GeV, '
      f'mY {MY_GRID.min():.2f}-{MY_GRID.max():.2f} GeV, '
      f'rv {RV_GRID.min():.4f}-{RV_GRID.max():.4f}')
print(f'coarse mesh {LMX_C.shape[0]} x {LMX_C.shape[1]} = {LMX_C.size} nodes, '
      f'mX {10.0**lmx_c[0]:.2f}-{10.0**lmx_c[-1]:.2f} GeV, '
      f'rv {10.0**lrv_c[0]:.4f}-{10.0**lrv_c[-1]:.4f}')

# ---- solve the coarse mesh, interpolate, evaluate <sigma v> --------------

def _solve_node(task):
    """One coarse node, for the pool: (i, j, mX, mY) -> (i, j, coupling, message).

    The failure is carried back rather than raised, so one bad node still cannot kill
    the sweep now that it runs in a worker.
    """
    i, j, mX, mY = task
    try:
        return i, j, _coupling_relic(mX, mY), ''
    except Exception as exc:
        return i, j, np.nan, f'  FAILED mX={mX:.2f} mY={mY:.2f}: {exc}'


tasks = []
for i in range(LMX_C.shape[0]):
    for j in range(LMX_C.shape[1]):
        mX = 10.0**LMX_C[i, j]
        mY = mX*10.0**LRV_C[i, j]
        if mY > RV_MAX*mX:
            continue
        tasks.append((i, j, mX, mY))

# See the note on relic_grid in make_GCE_sv_thermal_relic_grids.py: the nodes are
# independent, so the sweep is a map and the numbers do not depend on NCORES.
coupling_c = np.full(LMX_C.shape, np.nan)
n_failed = 0
with mp.get_context('fork').Pool(NCORES) as pool:
    for i, j, coupling, message in tqdm(pool.imap_unordered(_solve_node, tasks,
                                                            chunksize=1),
                                        total=len(tasks), desc='hP ext'):
        coupling_c[i, j] = coupling
        if message:
            n_failed += 1
            print(message, flush=True)

ok_c = np.isfinite(coupling_c)
print(f'coarse: {ok_c.sum()}/{coupling_c.size} solved, {n_failed} raised')

coupling = 10.0**griddata(
    np.column_stack([LMX_C[ok_c], LRV_C[ok_c]]), np.log10(coupling_c[ok_c]),
    np.column_stack([lmx, lrv]), method='linear')

sigmav = np.full(n_pts, np.nan)
for k in np.flatnonzero(np.isfinite(coupling)):
    m = _model(float(MX_GRID[k]), float(MY_GRID[k]), float(coupling[k]))
    sigmav[k] = m.sigmav_XX_to_YY_swave()*GEVM2_TO_CM3S

# the s-wave closes at mY = mX: nothing reaches the GCE there, whatever the coupling
sigmav[MY_GRID >= MX_GRID] = 0.0

print(f'targets: {np.isfinite(sigmav).sum()}/{n_pts} with a <sigma v>')

# ---- fold into the portal's grid -----------------------------------------
# One npz per portal, so the notebook reads one file and never has to know where the
# seam is.  Idempotent: the main driver's own N_MAIN points are kept and anything
# after them is replaced, so re-running this appends the extension once, not twice.
# The coarse mesh is the two drivers' own diagnostic and they solve different ones,
# so this one is carried under a _ext suffix rather than overwriting the main's.

out_path = OUT_DIR / OUT_NAME
if not out_path.exists():
    raise FileNotFoundError(f'{out_path} not found; run '
                            'make_GCE_sv_thermal_relic_grids.py first')

with np.load(out_path) as d:
    if d['mX'].size < N_MAIN:
        raise ValueError(f'{OUT_NAME} holds {d["mX"].size} points, fewer than the '
                         f'{N_MAIN} of the main grid; rerun the main driver')
    saved = {k: d[k] for k in d.files}

appended = {'mX': MX_GRID, 'mY': MY_GRID,
            'coupling_relic': coupling, 'sigmav_relic': sigmav}
for key, column in appended.items():
    saved[key] = np.concatenate([saved[key][:N_MAIN], column])
saved.update(mX_coarse_ext=10.0**LMX_C, mY_coarse_ext=10.0**(LMX_C + LRV_C),
             coupling_coarse_ext=coupling_c)

np.savez(out_path, **saved)
print(f'-> {out_path}: {N_MAIN} + {n_pts} = {saved["mX"].size} points')

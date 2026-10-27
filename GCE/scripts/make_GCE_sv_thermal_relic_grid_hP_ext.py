"""Thermal (SM-locked) relic coupling and <sigma v> for the Higgs portal over
the mX = 100-130 GeV extension grid.

The extension half of make_GCE_sv_thermal_relic_grids.py, which stops at the
main grid's mX = 100 GeV edge.  The Higgs portal only: the vector portals' GCE
contours close inside the main grid, so nothing above 100 GeV is needed for
them.  Nothing in GCE_sv_thermal_relic_grid_{hP,vP,bP,blP}.npz is read or
rewritten; this writes its own file.

    cosmo_env
    python GCE/scripts/make_GCE_sv_thermal_relic_grid_hP_ext.py

Output: GCE/output/sv_relic_fits/GCE_sv_thermal_relic_grid_hP_ext.npz
    mX, mY               [GeV]    the 578 extension points, in the order
                                  make_channel_spectra.extension_grid gives them
    coupling_relic                lam giving Omega_c h^2 = OCH2_TARGET, interpolated
    sigmav_relic         [cm^3/s] the s-wave <sigma v> there
    mX_coarse, mY_coarse [GeV]    the solved coarse mesh
    coupling_coarse               lam solved on it (nan where not solved)

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
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import griddata
from tqdm import trange

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from hidden_sector_DM.HiddenSectorDM import (                  # noqa: E402
    Cosmology, BoltzmannSolver, HiggsPortal, s0_cosmo, rhoc,
)
from GCE.spectrum import GEVM2_TO_CM3S                         # noqa: E402

from make_channel_spectra import extension_grid                # noqa: E402

# ---- configuration -------------------------------------------------------

OCH2_TARGET = 0.12
RV_MAX      = 0.995          # the solve fails at mY = mX
RELIC_RTOL  = 3e-3

N_MX_COARSE = 8              # coarse mesh before padding; the extension spans only
N_RV_COARSE = 8              # 0.1139 dex in mX and 0.2218 dex in rv, so 8 x 8 keeps
                             # the main driver's coarse density (0.0183 / 0.0394 dex)

COUPLING_RANGE = (2e-5, 1.0)   # lam, the axis make_relic_grid_hP.py samples
COUPLING_GUESS = 0.1
OMEGA_POWER    = 4             # sigmav ~ lam_s^2 lam_p^2 = lam^4, so Omega ~ lam^-4

GX_MAJORANA = 2              # Majorana fermion DM: 2 spin states, no antiparticle
GY_REAL_SCALAR = 1           # real scalar mediator: one degree of freedom

OUT_DIR = ROOT / 'GCE' / 'output' / 'sv_relic_fits'
OUT_NAME = 'GCE_sv_thermal_relic_grid_hP_ext.npz'

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
lmx_c = np.linspace(lmx.min() - dmx, lmx.max() + dmx, N_MX_COARSE + 2)
lrv_c = np.linspace(lrv.min() - drv, np.log10(RV_MAX), N_RV_COARSE + 1)

LMX_C, LRV_C = np.meshgrid(lmx_c, lrv_c, indexing='ij')

print(f'{n_pts} extension points, mX {MX_GRID.min():.2f}-{MX_GRID.max():.2f} GeV, '
      f'mY {MY_GRID.min():.2f}-{MY_GRID.max():.2f} GeV, '
      f'rv {RV_GRID.min():.4f}-{RV_GRID.max():.4f}')
print(f'coarse mesh {LMX_C.shape[0]} x {LMX_C.shape[1]} = {LMX_C.size} nodes, '
      f'mX {10.0**lmx_c[0]:.2f}-{10.0**lmx_c[-1]:.2f} GeV, '
      f'rv {10.0**lrv_c[0]:.4f}-{10.0**lrv_c[-1]:.4f}')

# ---- solve the coarse mesh, interpolate, evaluate <sigma v> --------------

coupling_c = np.full(LMX_C.shape, np.nan)
n_failed = 0

for i in trange(LMX_C.shape[0], desc='hP ext'):
    for j in range(LMX_C.shape[1]):
        mX = 10.0**LMX_C[i, j]
        mY = mX*10.0**LRV_C[i, j]
        if mY > RV_MAX*mX:
            continue
        try:
            coupling_c[i, j] = _coupling_relic(mX, mY)
        except Exception as exc:              # one bad node must not kill the sweep
            n_failed += 1
            print(f'  FAILED mX={mX:.2f} mY={mY:.2f}: {exc}', flush=True)

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

OUT_DIR.mkdir(parents=True, exist_ok=True)
np.savez(OUT_DIR / OUT_NAME, mX=MX_GRID, mY=MY_GRID,
         coupling_relic=coupling, sigmav_relic=sigmav,
         mX_coarse=10.0**LMX_C, mY_coarse=10.0**(LMX_C + LRV_C),
         coupling_coarse=coupling_c)
print(f'-> {OUT_DIR / OUT_NAME}')

"""Free-<sigma v> GCE fit of every portal over the (mX, mY) grid.

For each grid point the cascade template is built from the channel-spectra cache and
<sigma v> is minimised out against the GCE measurement by a generalised least squares
fit, so <sigma v> is read off the gamma rays alone and no relic constraint enters.

    cosmo_env
    python GCE/scripts/make_GCE_sv_grids.py

Set PORTALS to pick which portals run.  Above the main grid's mX = 100 GeV edge the
Higgs portal continues over `make_channel_spectra.EXT_GRID`, written to its own file.

Output: GCE/output/sv_fits/
    GCE_sv_grid_{portal}_{REGION}.npz          the 16384 main-grid points
    GCE_sv_grid_{portal}_ext_{REGION}.npz      the 578 extension points, hP only
        mX, mY       [GeV]     the points, in the grid's own order
        chi2                   chi^2 at the best-fit <sigma v>
        sigmav       [cm^3/s]  the GLS <sigma v>, clamped at >= 0
        sigmav_err   [cm^3/s]  its 1-sigma error
        ts                     the 14-bin template at <sigma v> = 1 cm^3/s
"""

import sys
from pathlib import Path
from tqdm import trange
import numpy as np
ROOT = Path.cwd().resolve()
while ROOT.name and not (ROOT / 'hidden_sector_DM').is_dir():
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
from hidden_sector_DM.HiddenSectorDM import VectorPortal, BaryonPortal, HiggsPortal,BLPortal  

from GCE.likelihood import cascade_template, gls_fit, measurement
from GCE.pythia_runner import DEFAULT_SEED
from GCE.halo import HALO_2112, RHO_LOCAL_PRIOR
from GCE.spectrum import GEVM2_TO_CM3S
import time
from scipy.optimize import brentq
from scipy.interpolate import griddata
from hidden_sector_DM.HiddenSectorDM import Cosmology, BoltzmannSolver, s0_cosmo, rhoc

sys.path.insert(0, str(ROOT / 'GCE' / 'scripts'))    # sibling scripts
from make_channel_spectra import EXT_GRID               # the hP mX = 100-130 block

# ---- configuration -------------------------------------------------------
MX_LIM = (15.0, 100.0)   # GeV
MY_LIM = (5.0, 100.0)    # GeV, capped from above by mX
N_POINTS = 16384
N_COLS = 124
DECIMALS = 2

def triangle_grid(n_points: int = N_POINTS, n_cols: int = N_COLS,
                  mx_lim=MX_LIM, my_lim=MY_LIM, decimals: int = DECIMALS):
    """The grid as a (2, `n_points`) array of (mX, mY) in GeV."""
    u = np.linspace(*np.log10(mx_lim), n_cols)          # log10 mX per column
    v_lo = np.log10(my_lim[0])
    v_hi = np.minimum(u, np.log10(my_lim[1]))           # the mY <= mX edge
    height = v_hi - v_lo
    if not (height > 0).all():
        raise ValueError("some column has no room in mY; check the limits")

    # Points per column, proportional to column height.  Largest-remainder
    # apportionment, so the total is exactly `n_points`.
    raw = height / height.sum() * n_points
    n = np.floor(raw).astype(int)
    short = n_points - n.sum()
    if short:
        n[np.argsort(n - raw)[:short]] += 1
    if n.min() < 2:
        raise ValueError(f"column of {n.min()} point(s): lower N_COLS")

    mX = np.repeat(10.0 ** u, n)
    mY = np.concatenate([10.0 ** np.linspace(v_lo, hi, k)
                         for hi, k in zip(v_hi, n)])
    return np.round(np.vstack([mX, mY]), decimals)

STRIDE = 1

def thin(grid, stride: int):
    """Every `stride`-th point in both directions, trapezoid edges kept."""
    if stride == 1:
        return grid
    mX = grid[0]
    starts = np.flatnonzero(np.r_[True, np.diff(mX) != 0.0])   # column starts
    ends = np.r_[starts[1:], mX.size]
    cols = np.unique(np.r_[np.arange(0, starts.size, stride), starts.size - 1])
    keep = np.concatenate([
        np.unique(np.r_[np.arange(starts[c], ends[c], stride), ends[c] - 1])
        for c in cols])
    return grid[:, keep]


GRID = thin(triangle_grid(), STRIDE)    
MX_GRID, MY_GRID = GRID     
print(f"{GRID.shape[1]} of {N_POINTS} points, "
      f"{np.unique(MX_GRID).size} of {N_COLS} mX columns  (STRIDE={STRIDE})")

REGION = "south"        # 'both' (40x40), 'south' or 'north'
N_EVENTS = 400_000      # what the channel_spectra cache holds
SEED = DEFAULT_SEED
ALPHAX = 1e-2           # the dataclass wants a default value
N_GCE_BINS = 14         # bins in the measurement, so columns in a template

PORTAL_NAMES = ('hP', 'vP', 'bP', 'blP')

# Which portals to fit, and which of them continue past the main grid's mX = 100 GeV
# edge.  Set PORTALS to ('hP',) to refresh the Higgs portal alone: the other three are
# then left as they are rather than rewritten with the same numbers.  Only the Higgs
# portal has an extension -- make_channel_spectra generates that block's channels for
# it alone -- so EXT_PORTALS cannot name any other.
PORTALS = ('hP', 'vP', 'bP', 'blP')
EXT_PORTALS = ('hP',)

if not set(PORTALS) <= set(PORTAL_NAMES):
    raise ValueError(f'PORTALS {PORTALS} must come from {PORTAL_NAMES}')
if not set(EXT_PORTALS) <= set(PORTALS) & {'hP'}:
    raise ValueError(f'EXT_PORTALS {EXT_PORTALS} must be fitted portals, and only the '
                     'Higgs portal has cached spectra above mX = 100 GeV')


def portal_model(name, mX, mY):
    """The portal `name` at (`mX`, `mY`) in GeV, at the overlay coupling ALPHAX.

    The coupling never reaches the template -- `cascade_template` divides the
    normalisation out -- so it is only here because the dataclasses want a value.
    """
    if name == 'hP':
        return HiggsPortal(mX=mX, mY=mY, gX=2, gY=1, lam=ALPHAX)
    elif name == 'vP':
        return VectorPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=ALPHAX)
    elif name == 'bP':
        return BaryonPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=ALPHAX)
    elif name == 'blP':
        return BLPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=ALPHAX)
    raise ValueError(f'unknown portal {name!r}; choose one of {PORTAL_NAMES}')


def sigmav_chi2_GCE(model, phi_GCE, Cinv_GCE):
    t = cascade_template(model, n_events=N_EVENTS, seed=SEED)
    sigmav, sigmav_err, chi2= gls_fit(t, phi_GCE, Cinv_GCE)
    return t, sigmav, sigmav_err, chi2

# ---- the fit -------------------------------------------------------------

phi, Cinv = measurement(region=REGION)      # once; independent of the model

def fit_grid(name, grid, desc):
    """GLS fit of portal `name` over `grid`, a (2, n) array of (mX, mY) in GeV.

    Returns the save dictionary: `chi2`, `sigmav` [cm^3/s], `sigmav_err` [cm^3/s] and
    the templates `ts`, one row per grid point in the grid's own order.
    """
    n_pts = grid.shape[1]
    chi2 = np.full(n_pts, np.nan)
    sv = np.full(n_pts, np.nan)
    sv_err = np.full(n_pts, np.nan)
    ts = np.full((n_pts, N_GCE_BINS), np.nan)

    for k in trange(n_pts, desc=desc):
        mX, mY = (float(v) for v in grid[:, k])
        model = portal_model(name, mX, mY)
        ts[k], sv[k], sv_err[k], chi2[k] = sigmav_chi2_GCE(model, phi, Cinv)

    return {'chi2': chi2, 'sigmav': sv, 'sigmav_err': sv_err, 'ts': ts}


#save them in GCE/output/sv_fits/
OUT_DIR = ROOT / 'GCE' / 'output' / 'sv_fits'
OUT_DIR.mkdir(parents=True, exist_ok=True)

for name in PORTALS:
    dic = fit_grid(name, GRID, f'{name} main')
    np.savez(OUT_DIR / f'GCE_sv_grid_{name}_{REGION}.npz', mX=MX_GRID, mY=MY_GRID, **dic)

    if name in EXT_PORTALS:
        dic = fit_grid(name, EXT_GRID, f'{name} ext')
        np.savez(OUT_DIR / f'GCE_sv_grid_{name}_ext_{REGION}.npz',
                 mX=EXT_GRID[0], mY=EXT_GRID[1], **dic)

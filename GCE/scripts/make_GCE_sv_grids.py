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

def sigmav_chi2_GCE(model, phi_GCE, Cinv_GCE):
    t = cascade_template(model, n_events=N_EVENTS, seed=SEED)
    sigmav, sigmav_err, chi2= gls_fit(t, phi_GCE, Cinv_GCE)
    return t, sigmav, sigmav_err, chi2

# ---- the fit -------------------------------------------------------------

phi, Cinv = measurement(region=REGION)      # once; independent of the model

n_pts = GRID.shape[1]

chi2_hP, chi2_vP, chi2_bP, chi2_blP = np.full(n_pts, np.nan), np.full(n_pts, np.nan), np.full(n_pts, np.nan), np.full(n_pts, np.nan)
ts_hP, ts_vP, ts_bP, ts_blP = np.full((n_pts, 14), np.nan), np.full((n_pts, 14), np.nan), np.full((n_pts, 14), np.nan), np.full((n_pts, 14), np.nan)
sv_hP, sv_vP, sv_bP, sv_blP = np.full(n_pts, np.nan), np.full(n_pts, np.nan), np.full(n_pts, np.nan), np.full(n_pts, np.nan)
sv_err_hP, sv_err_vP, sv_err_bP, sv_err_blP = np.full(n_pts, np.nan), np.full(n_pts, np.nan), np.full(n_pts, np.nan), np.full(n_pts, np.nan)

for k in trange(n_pts):
    mX, mY = (float(v) for v in GRID[:, k])
    mod_hP = HiggsPortal(mX=mX, mY=mY, gX=2, gY=1, lam=ALPHAX)
    mod_vP = VectorPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=ALPHAX)
    mod_bP = BaryonPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=ALPHAX)
    mod_blP = BLPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=ALPHAX)
    ts_hP[k], sv_hP[k], sv_err_hP[k],  chi2_hP[k]  = sigmav_chi2_GCE(mod_hP, phi, Cinv)
    ts_vP[k], sv_vP[k], sv_err_vP[k],  chi2_vP[k]  =sigmav_chi2_GCE(mod_vP, phi, Cinv)
    ts_bP[k], sv_bP[k], sv_err_bP[k],  chi2_bP[k] = sigmav_chi2_GCE(mod_bP, phi, Cinv)
    ts_blP[k], sv_blP[k], sv_err_blP[k],  chi2_blP[k] =sigmav_chi2_GCE(mod_blP, phi, Cinv)

dic_hP={'chi2':chi2_hP, 'sigmav': sv_hP, 'sigmav_err': sv_err_hP, 'ts':ts_hP}
dic_vP={'chi2':chi2_vP, 'sigmav': sv_vP, 'sigmav_err': sv_err_vP, 'ts':ts_vP}
dic_bP={'chi2':chi2_bP, 'sigmav': sv_bP, 'sigmav_err': sv_err_bP, 'ts':ts_bP}
dic_blP={'chi2':chi2_blP, 'sigmav': sv_blP, 'sigmav_err': sv_err_blP, 'ts':ts_blP}

#save them in GCE/output/sv_fits/
OUT_DIR = ROOT / 'GCE' / 'output' / 'sv_fits'
OUT_DIR.mkdir(parents=True, exist_ok=True)
for name, dic in [('hP', dic_hP), ('vP', dic_vP), ('bP', dic_bP), ('blP', dic_blP)]:
    np.savez(OUT_DIR / f'GCE_sv_grid_{name}_{REGION}.npz', mX=MX_GRID, mY=MY_GRID, **dic)

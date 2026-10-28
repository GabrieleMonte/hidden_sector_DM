import sys
from pathlib import Path
from tqdm import trange
import numpy as np
ROOT = Path.cwd().resolve()
while ROOT.name and not (ROOT / 'hidden_sector_DM').is_dir():
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
from hidden_sector_DM.HiddenSectorDM import VectorPortal, HiggsPortal

from GCE.spectrum import GEVM2_TO_CM3S
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

# ---- coupling that gives Omega_c h^2 = 0.12 in the SM-thermalised regime ----
OCH2_TARGET  = 0.12
N_MX_COARSE  = 46            # coarse (log10 mX, log10 mY/mX) grid, interpolated onto GRID
N_RV_COARSE  = 34            # rows the uniform log10 rv ladder carries, up to RV_MAX
N_RV_BETA    = 4             # rows replacing its top one, uniform in log10 beta
RV_MAX       = 0.995         # the solve fails at mY = mX; 0.995 covers every GRID point below it
RELIC_RTOL   = 3e-3

# Per relic computation: the coupling searched, its starting guess, and n in
# Omega ~ coupling^-n, which sets the first scaled step.
COUPLING_RANGE = {'hP': (2e-5, 1.0),    # lam, the axis make_relic_grid_hP.py samples
                  'vP': (1e-5, 1.0)}    # perturbative ceiling alpha_X = gX^2 / 4pi < 1
COUPLING_GUESS = {'hP': 0.1, 'vP': 1e-2}
OMEGA_POWER    = {'hP': 4,              # sigmav ~ lam_s^2 lam_p^2 = lam^4
                  'vP': 2}              # sigmav ~ alpha_X^2

# Which portals to solve.  'vP' also writes the bP and blP files, for the reason
# given at the save below.  Set to ('hP',) to refresh the Higgs-portal grid alone:
# that halves the run, and the three vector-portal files are then left untouched
# rather than rewritten with the same numbers.
PORTALS = ('hP', 'vP')

if not set(PORTALS) <= set(COUPLING_RANGE):
    raise ValueError(f'PORTALS {PORTALS} must come from {tuple(COUPLING_RANGE)}; '
                     'bP and blP are not solved, they are written from the vP solve')

cosmo = Cosmology()

def _model(portal, mX, mY, coupling):
    if portal == 'hP':
        return HiggsPortal(mX=mX, mY=mY, gX=2, gY=1, lam=coupling,
                           include_antiparticlesX=False, include_antiparticlesY=False)
    return VectorPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=coupling,
                        include_antiparticlesX=True, include_antiparticlesY=False)

def _och2(portal, mX, mY, coupling):
    sol = BoltzmannSolver(cosmo, _model(portal, mX, mY, coupling)).solve_boltzmann_thermSM(
        verbose=False, include_hidden_in_gstar=False)
    return sol['YX_relic'] * mX * s0_cosmo / rhoc

def _coupling_relic(portal, mX, mY):
    """Coupling for Omega h^2 = OCH2_TARGET (nan if outside COUPLING_RANGE);
    scaled guess + log-log secant, brentq fallback."""
    a_lo, a_hi = COUPLING_RANGE[portal]
    a0 = COUPLING_GUESS[portal]
    o0 = _och2(portal, mX, mY, a0)
    a1 = float(np.clip(a0 * (o0 / OCH2_TARGET) ** (1.0 / OMEGA_POWER[portal]), a_lo, a_hi))
    o1 = _och2(portal, mX, mY, a1)
    if abs(o1 / OCH2_TARGET - 1.0) < RELIC_RTOL:
        return a1
    la = (np.log(a1) + (np.log(OCH2_TARGET) - np.log(o1))
          * (np.log(a1) - np.log(a0)) / (np.log(o1) - np.log(o0)))
    a2 = float(np.clip(np.exp(la), a_lo, a_hi))
    if abs(_och2(portal, mX, mY, a2) / OCH2_TARGET - 1.0) < RELIC_RTOL:
        return a2
    f = lambda l: np.log(_och2(portal, mX, mY, np.exp(l)) / OCH2_TARGET)
    for lo, hi in ((a2 / 3.0, a2 * 3.0), (a_lo, a_hi)):
        try:
            return float(np.exp(brentq(f, np.log(max(lo, a_lo)),
                                       np.log(min(hi, a_hi)), rtol=1e-3)))
        except ValueError:
            continue
    return np.nan

# ---- the relic -----------------------------------------------------------

def coarse_rv_rows(rv_min, rv_max=RV_MAX, n_rv=N_RV_COARSE, n_beta=N_RV_BETA):
    """log10 of the coarse mY/mX rows, placed densely where the s-wave closes.

    Everything below rv ~ 0.91 is the old uniform ladder in log10 rv, unchanged.  Its
    single top row is replaced by `n_beta` rows at a uniform step in log10 beta, with
    beta = sqrt(1 - rv^2).

    The s-wave shuts off as beta -> 0, so a ladder uniform in log10 rv leaves ONE
    interval from 0.909 to `rv_max` covering six fine-grid rows, and the <sigma v>
    interpolated across it is wrong by up to 65% (vP) or 9% (hP) against an exact
    solve.  The same number of rows placed by beta brings both under 4%, at every mX
    tried.  Interpolating in log beta while leaving the rows alone was also tried and
    is WORSE for hP -- its lam barely moves with rv and turns over near `rv_max` -- so
    it is the row placement that matters, not the interpolation variable.  Measured in
    the scratchpad script rv_mesh_compare.py.
    """
    ladder = np.linspace(np.log10(rv_min), np.log10(rv_max), n_rv)
    beta = np.sqrt(1.0 - (10.0 ** ladder[-2:]) ** 2)
    log_beta = np.linspace(np.log10(beta[0]), np.log10(beta[1]), n_beta + 1)[1:]
    return np.concatenate([ladder[:-1], 0.5 * np.log10(1.0 - (10.0 ** log_beta) ** 2)])


lmx_c = np.linspace(np.log10(MX_GRID.min()), np.log10(MX_GRID.max()), N_MX_COARSE)
lrv_c = coarse_rv_rows((MY_GRID / MX_GRID).min())
LMX_C, LRV_C = np.meshgrid(lmx_c, lrv_c, indexing='ij')
RV_STEP = 10.0 ** (lrv_c[1] - lrv_c[0])      # one coarse row in mY/mX, at the bottom

n_pts = GRID.shape[1]

def relic_grid(portal):
    """Solve the coupling on the coarse grid, interpolate log10 of it onto GRID,
    and evaluate the s-wave <sigma v> there [cm^3/s]."""
    coupling_c = np.full(LMX_C.shape, np.nan)
    for i in trange(N_MX_COARSE, desc=portal):
        for j in range(lrv_c.size):
            mX = 10.0 ** LMX_C[i, j]
            mY = mX * 10.0 ** LRV_C[i, j]
            # one row below the mY floor, so the interpolation covers mY = MY_LIM[0]
            if MY_LIM[0] / RV_STEP <= mY <= RV_MAX * mX:
                coupling_c[i, j] = _coupling_relic(portal, mX, mY)

    ok_c = np.isfinite(coupling_c)
    coupling = 10.0 ** griddata(
        np.column_stack([LMX_C[ok_c], LRV_C[ok_c]]), np.log10(coupling_c[ok_c]),
        np.column_stack([np.log10(MX_GRID), np.log10(MY_GRID / MX_GRID)]), method='linear')

    sv = np.full(n_pts, np.nan)
    for k in np.flatnonzero(np.isfinite(coupling)):
        m = _model(portal, float(MX_GRID[k]), float(MY_GRID[k]), float(coupling[k]))
        sv[k] = m.sigmav_XX_to_YY_swave() * GEVM2_TO_CM3S
    # the s-wave closes at mY = mX: nothing reaches the GCE there, whatever the coupling
    sv[MY_GRID >= MX_GRID] = 0.0
    return {'coupling_relic': coupling, 'sigmav_relic': sv,
            'mX_coarse': 10.0 ** LMX_C, 'mY_coarse': 10.0 ** (LMX_C + LRV_C),
            'coupling_coarse': coupling_c}

dics = {portal: relic_grid(portal) for portal in PORTALS}
# The vector portals share the hidden sector, and solve_boltzmann_thermSM sees only
# mX, mY, gX, gY and sigmav_XX_to_YY -- never the portal -- so their thermal relic
# is identical (checked: same Omega h^2 to 4 digits).  Solved once, saved three times.
if 'vP' in dics:
    dics['bP'] = dics['vP']
    dics['blP'] = dics['vP']

#save them in GCE/output/sv_relic_fits/
OUT_DIR = ROOT / 'GCE' / 'output' / 'sv_relic_fits'
OUT_DIR.mkdir(parents=True, exist_ok=True)
for name, dic in dics.items():
    np.savez(OUT_DIR / f'GCE_sv_thermal_relic_grid_{name}.npz', mX=MX_GRID, mY=MY_GRID, **dic)

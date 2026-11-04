"""Free-<sigma v> GCE fit of every portal over the (mX, mY) grid.

For each grid point the cascade template is built from the channel-spectra cache and
<sigma v> is minimised out against the GCE measurement by a generalised least squares
fit, so <sigma v> is read off the gamma rays alone and no relic constraint enters.
Loads only: every spectrum must already be in GCE/channel_spectra/
(make_channel_spectra.py), and a missing one raises naming the file.

    cosmo_env
    python GCE/scripts/make_gce_fits.py

The (mX, mY) points are make_channel_spectra's, so the fit lands on exactly the points
the cache holds.  Set PORTALS to pick which portals run.  Above the main grid's
mX = 100 GeV edge the Higgs portal continues over `make_channel_spectra.EXT_GRID`,
APPENDED to the main grid rather than written to its own file, so one npz is one
portal's whole grid and the reader never has to know where the seam is.  hP therefore
has 16384 + 578 = 16962 points where the other three have 16384.

Output: GCE/output/fits/gce_fit_{portal}_{REGION}.npz
        mX, mY       [GeV]     the points, in the grid's own order
        chi2                   chi^2 at the best-fit <sigma v>
        sigmav       [cm^3/s]  the GLS <sigma v>, clamped at >= 0
        sigmav_err   [cm^3/s]  its 1-sigma error
        ts                     the 14-bin template at <sigma v> = 1 cm^3/s
"""

import sys
from pathlib import Path

import numpy as np
from tqdm import trange

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))   # sibling scripts

from hidden_sector_DM.HiddenSectorDM import (                  # noqa: E402
    BaryonPortal, BLPortal, HiggsPortal, VectorPortal,
)
from GCE.source.likelihood import cascade_template, gls_fit, measurement   # noqa: E402
from GCE.source.pythia_runner import DEFAULT_SEED              # noqa: E402

from make_channel_spectra import EXT_GRID, MAIN_GRID           # noqa: E402

# ---- configuration -------------------------------------------------------

REGION = 'south'        # 'both' (40x40), 'south' or 'north'
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

OUT_DIR = ROOT / 'GCE' / 'output' / 'fits'


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


def portal_grid(name):
    """The (2, n) array of (mX, mY) in GeV that portal `name` is fitted on.

    An extended portal carries the mX = 100-130 GeV block appended to the main grid,
    so its arrays are 578 points longer than the other three.
    """
    if name in EXT_PORTALS:
        return np.hstack([MAIN_GRID, EXT_GRID])
    else:
        return MAIN_GRID


def sigmav_chi2_GCE(model, phi_GCE, Cinv_GCE):
    t = cascade_template(model, n_events=N_EVENTS, seed=SEED)
    sigmav, sigmav_err, chi2 = gls_fit(t, phi_GCE, Cinv_GCE)
    return t, sigmav, sigmav_err, chi2


def fit_grid(name, grid, phi, Cinv):
    """GLS fit of portal `name` over `grid`, a (2, n) array of (mX, mY) in GeV.

    Returns the save dictionary: `chi2`, `sigmav` [cm^3/s], `sigmav_err` [cm^3/s] and
    the templates `ts`, one row per grid point in the grid's own order.
    """
    n_pts = grid.shape[1]
    chi2 = np.full(n_pts, np.nan)
    sv = np.full(n_pts, np.nan)
    sv_err = np.full(n_pts, np.nan)
    ts = np.full((n_pts, N_GCE_BINS), np.nan)

    for k in trange(n_pts, desc=name):
        mX, mY = (float(v) for v in grid[:, k])
        model = portal_model(name, mX, mY)
        ts[k], sv[k], sv_err[k], chi2[k] = sigmav_chi2_GCE(model, phi, Cinv)

    return {'chi2': chi2, 'sigmav': sv, 'sigmav_err': sv_err, 'ts': ts}


# ---- the fit -------------------------------------------------------------

phi, Cinv = measurement(region=REGION)      # once; independent of the model
OUT_DIR.mkdir(parents=True, exist_ok=True)

for name in PORTALS:
    grid = portal_grid(name)
    dic = fit_grid(name, grid, phi, Cinv)
    np.savez(OUT_DIR / f'gce_fit_{name}_{REGION}.npz', mX=grid[0], mY=grid[1], **dic)

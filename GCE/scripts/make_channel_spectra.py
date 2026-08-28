"""
Build the library of per-channel photon spectra over the (mX, mY) grid,
one process per core.

    cosmo_env
    cd hidden_sector_DM
    python GCE/scripts/make_channel_spectra.py       # every core it can see
    python GCE/scripts/make_channel_spectra.py -n 128

No portal appears below.  A cached spectrum is counts PER DECAY of one channel;
branching ratios only weight them later, in `cascade_spectrum`.  So every
channel runs at every point with BR = 1, and one cache serves any model --
including channels a given portal never opens.

One task is one (mX, mY), writing one .npz per open channel into `OUT_DIR`.
Existing files are loaded rather than regenerated, so a job killed at the wall
clock is simply relaunched.
"""

import argparse
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))   # sibling scripts

from GCE.pythia_runner import (                             
    CHANNEL_TO_PDG, MASSES, PythiaRunner, channels_for_species, run_spectrum,
)
from GCE.spectrum import default_ebins                      

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


GRID = triangle_grid()          # (2, 16384) masses in GeV
N_EVENTS = 200_000
SEED = 12345
SPECIES = "gamma"               # picks the channels; the .npz holds all targets
N_BINS = 180                    # `default_ebins` resolution, per point
OUT_DIR = Path(__file__).resolve().parents[1] / "channel_spectra"

# Every channel Pythia can shower, with BR = 1 so none is dropped for being
# small.  Neutrino channels still go, since they make no photons at all.
CHANNELS = [ch for ch, _br in channels_for_species({c: 1.0 for c in CHANNEL_TO_PDG}, SPECIES)]
_RUNNER = None                  # one Pythia per worker process, never pickled

def run_point(i: int) -> str:
    """Run every open channel of grid point `i`.  Returns '' on success."""
    global _RUNNER
    mX, mY = (float(v) for v in GRID[:, i])
    ebins = default_ebins(mX, n_bins=N_BINS)
    ch = ""
    try:
        if _RUNNER is None:                 # 128 of these start at once
            _RUNNER = PythiaRunner(seed=SEED)
        for ch in CHANNELS:
            if sum(MASSES[abs(p)] for p in CHANNEL_TO_PDG[ch]) >= mY:
                continue                    # Y is too light to decay this way
            run_spectrum(_RUNNER, mX, mY, ch, N_EVENTS, ebins, OUT_DIR,
                         verbose=False)
    except Exception as exc:                # one bad point must not kill the pool
        return f"FAILED mX={mX:g} mY={mY:g} [{ch}]: {exc}"
    return ""


ap = argparse.ArgumentParser()
ap.add_argument("-n", "--ncores", type=int,
                    default=len(os.sched_getaffinity(0)))
args = ap.parse_args()
n = GRID.shape[1]
OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"{n} points, <={len(CHANNELS)} channels each, "
      f"{args.ncores} cores -> {OUT_DIR}", flush=True)
t0 = time.time()

with mp.Pool(args.ncores) as pool:
    for k, msg in enumerate(pool.imap_unordered(run_point, range(n), 1), 1):
        if msg:
            print(msg, flush=True)
        if k % 200 == 0 or k == n:
            dt = time.time() - t0
            print(f"  {k}/{n}  {dt/60:.0f} min elapsed, "
                  f"{dt/k*(n-k)/60:.0f} min left", flush=True)

print(f"done in {(time.time()-t0)/60:.1f} min", flush=True)

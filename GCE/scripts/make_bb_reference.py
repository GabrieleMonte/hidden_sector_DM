"""
Build the direct-annihilation photon spectra: DM DM -> b bbar at sqrt(s) = 2 m_DM,
one Pythia run per DM mass.  The 1-D counterpart of `make_channel_spectra.py`,
which does the two-mass cascade grid.

    cosmo_env
    cd hidden_sector_DM
    python GCE/scripts/make_bb_reference.py       # every core it can see
    python GCE/scripts/make_bb_reference.py -n 16

A direct annihilation is a mediator AT REST decaying to the channel, so it is
cached at (m_X, m_Y) = (m_DM, 2 m_DM): `run_channel` clamps the space-like
momentum to zero and the daughters go back-to-back with E = m_DM each.  This is
exactly what `spectrum.direct_spectrum` / `likelihood.direct_template` load, so
one cache serves both.  One .npz per (mass, channel); existing files are skipped,
so a job killed at the wall clock is simply relaunched.

No portal and no <sigma v>: a cached spectrum is counts PER ANNIHILATION, and
kappa / <sigma v> / halo all enter later as cheap reweightings.  The favoured
region itself is built from these in a notebook, not here.
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

from GCE.pythia_runner import CHANNEL_TO_PDG, MASSES, PythiaRunner, run_spectrum
from GCE.spectrum import default_ebins

# ---- configuration -------------------------------------------------------

M_LIM = (15.0, 80.0)          # DM mass range [GeV]
N_MASS = 50                    # log-spaced points across M_LIM
DECIMALS = 2                   # rounding; the consumer must ask for the same masses

MASS_GRID = np.round(np.logspace(*np.log10(M_LIM), N_MASS), DECIMALS)

N_EVENTS = 200_000
SEED = 12345
N_BINS = 180                   # `default_ebins` resolution, per mass
CHANNELS = ("bb",)             # the arXiv:2112.09706 Fig. 18 reference; extend as needed
OUT_DIR = Path(__file__).resolve().parents[1] / "channel_spectra"

_RUNNER = None                 # one Pythia per worker process, never pickled


def run_point(i: int) -> str:
    """Run every channel at DM mass `MASS_GRID[i]`.  Returns '' on success."""
    global _RUNNER
    m_DM = float(MASS_GRID[i])
    m_X, m_Y = m_DM, 2.0 * m_DM          # mediator at rest == direct annihilation
    ebins = default_ebins(m_X, n_bins=N_BINS)
    ch = ""
    try:
        if _RUNNER is None:
            _RUNNER = PythiaRunner(seed=SEED)
        for ch in CHANNELS:
            if sum(MASSES[abs(p)] for p in CHANNEL_TO_PDG[ch]) >= m_Y:
                continue                 # sqrt(s) too low for this channel
            run_spectrum(_RUNNER, m_X, m_Y, ch, N_EVENTS, ebins, OUT_DIR,
                         verbose=False)
    except Exception as exc:             # one bad mass must not kill the pool
        return f"FAILED m_DM={m_DM:g} [{ch}]: {exc}"
    return ""


ap = argparse.ArgumentParser()
ap.add_argument("-n", "--ncores", type=int,
                default=len(os.sched_getaffinity(0)))
args = ap.parse_args()
n = len(MASS_GRID)
OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"{n} DM masses over {M_LIM[0]:g}-{M_LIM[1]:g} GeV, {len(CHANNELS)} "
      f"channel(s) each, {args.ncores} cores -> {OUT_DIR}", flush=True)
t0 = time.time()

with mp.Pool(args.ncores) as pool:
    for k, msg in enumerate(pool.imap_unordered(run_point, range(n), 1), 1):
        if msg:
            print(msg, flush=True)
        if k % 20 == 0 or k == n:
            dt = time.time() - t0
            print(f"  {k}/{n}  {dt/60:.0f} min elapsed, "
                  f"{dt/k*(n-k)/60:.0f} min left", flush=True)

print(f"done in {(time.time()-t0)/60:.1f} min", flush=True)

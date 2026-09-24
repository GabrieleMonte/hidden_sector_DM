"""
Upgrade the channel-spectra cache to the decay-in-flight convention.

`pythia_runner.PythiaRunner` now decays the long-lived species (mu+-, pi+-,
K+-, K_L: `LONG_LIVED_DECAYED`) that Pythia keeps detector-stable by default,
because over kpc of propagation they all decay and their photons arrive with
the prompt spectrum.  A fresh `make_channel_spectra.py` run produces the new
convention directly.  This script instead UPGRADES the already-computed cache
(200k events x ~160k files, weeks of CPU) to the same convention at ~1% of
that cost, exploiting that the extra photons factorise:

    the correction per Y decay depends only on (channel, sqrt(s) = mY) in the
    mediator rest frame; the mX dependence is the isotropic-decay boost.

Stages (run in this order; each is resumable / refuses to clobber):

  --move      GCE/channel_spectra -> GCE/legacy/channel_spectra (one rename),
              then recreates an empty GCE/channel_spectra.
  --library   For every (channel, sqrt(s)) needed by the legacy inventory, run
              N_LIB events at rest with the NEW runner and store only the
              CORRECTION photons on a fine energy grid -> LIB_FILE.
              A photon belongs to the correction iff its ancestry contains a
              decay product (|status| 91-99) whose mother is long-lived: those
              photons did not exist under the old convention.  Photons merely
              RADIATED off a long-lived particle (FSR at production, e.g. off
              a primary muon) exist in the legacy spectra too and are not
              double-counted.
  --migrate   new gamma = legacy gamma + boosted correction; pbar copied
              unchanged (mu/pi/K decays make no antibaryons).  Same filename,
              provenance keys added ('decayed_long_lived', 'decay_correction').
  --validate  At NV_NODES random legacy nodes, generate the full spectrum from
              scratch with the new runner (fresh seed) and compare with the
              migrated file: totals and per-bin pulls against the combined MC
              noise of the two.

Usage (from the repo root, cosmo_env active):
    python GCE/scripts/upgrade_spectra_decay_convention.py --move
    python GCE/scripts/upgrade_spectra_decay_convention.py --library -n 32
    python GCE/scripts/upgrade_spectra_decay_convention.py --migrate -n 32
    python GCE/scripts/upgrade_spectra_decay_convention.py --validate -n 8
"""

import argparse
import multiprocessing as mp
import os
import re
import sys
import time
from bisect import bisect_right
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from GCE.pythia_runner import (
    CHANNEL_TO_PDG, LONG_LIVED_DECAYED, MASSES, PythiaRunner, SPECTRA_DIR,
    _pythia_version,
)

# ---- configuration -------------------------------------------------------

LEGACY_DIR = SPECTRA_DIR.parent / "legacy" / "channel_spectra"
OUT_DIR    = SPECTRA_DIR                       # migrated files land here
LIB_FILE   = SPECTRA_DIR.parent / "legacy" / "decay_correction_library.npz"

N_LIB      = 400_000    # events per (channel, sqrt(s)) library point: the
                        # correction then carries at least the statistics of a
                        # 200k-event node on the component it adds, and each
                        # point is shared by many nodes (shared noise smooths,
                        # it does not scatter, the chi^2 surface)
N_SQRTS    = 30         # library sqrt(s) points per channel (log-spaced over
                        # that channel's mY range in the legacy inventory)
LIB_SEED   = 20260924   # never a cache seed: library events are not cache events

# Fine rest-frame grid the correction is stored on.  Extends below every cache
# ebins floor (default_ebins reaches ~0.1 GeV) and above sqrt(s)/2 for the
# largest sqrt(s) in the cache -- the direct b bbar reference points sit at
# mY = 2 m_DM, up to 200 GeV, so rest-frame photons reach 100 GeV.
FINE_EDGES = np.geomspace(1e-3, 110.0, 451)

# Channels whose correction is identically zero (no long-lived decay products).
NO_CORRECTION = {"ee", "gamgam", "nunu_e", "nunu_mu", "nunu_tau"}

NV_NODES   = 12         # --validate sample size
NV_EVENTS  = 200_000    # fresh-run statistics for validation
NV_SEED    = 777        # fresh seed: an independent MC draw, not a replay

_LONGLIVED_ABS = {abs(p) for p in LONG_LIVED_DECAYED}
_FNAME = re.compile(r"^(?P<ch>.+)_mX(?P<mX>[\d.eE+-]+)_mY(?P<mY>[\d.eE+-]+)"
                    r"_N(?P<N>\d+)_s(?P<s>\d+)\.npz$")


# ---- correction extraction (rest frame, decay-vertex tagged) ---------------

def _is_correction_photon(ev, i) -> bool:
    """Ancestry contains a decay product (|status| 91-99) of a long-lived
    particle.  Walks mother chains once, cycle-safe."""
    seen = set()
    stack = [i]
    while stack:
        j = stack.pop()
        if j <= 0 or j in seen:
            continue
        seen.add(j)
        prt = ev[j]
        if 91 <= abs(prt.status()) <= 99:
            for k in prt.motherList():
                if k > 0 and abs(ev[k].id()) in _LONGLIVED_ABS:
                    return True
        stack.extend(prt.motherList())
    return False


def correction_at_rest(runner, channel: str, sqrt_s: float,
                       n_events: int) -> np.ndarray:
    """Correction-photon counts per decay on FINE_EDGES, Y at rest."""
    counts = np.zeros(len(FINE_EDGES) - 1)
    edges = FINE_EDGES.tolist()
    lo, hi = FINE_EDGES[0], FINE_EDGES[-1]
    n_ok = 0
    for _ in range(int(n_events)):
        ev = runner.p.event
        ev.reset()
        runner._inject_decay(ev, CHANNEL_TO_PDG[channel], sqrt_s, 0.0)
        if not runner.p.next():
            continue
        n_ok += 1
        for i in range(ev.size()):
            prt = ev[i]
            if prt.id() != 22 or not prt.isFinal():
                continue
            e = prt.e()
            if lo <= e < hi and _is_correction_photon(ev, i):
                counts[bisect_right(edges, e) - 1] += 1.0
    return counts / float(n_ok if n_ok else n_events)


# ---- boost: fine rest-frame counts -> counts on a file's ebins -------------

def boost_counts(fine_counts: np.ndarray, m_X: float, m_Y: float,
                 ebins: np.ndarray) -> np.ndarray:
    """Isotropic-decay boost of per-decay photon counts onto `ebins`.

    Uses the same momentum clamp as `PythiaRunner.run_channel`:
    p = sqrt(max(mX^2 - mY^2, 0)), so mY >= mX means "mediator at rest"
    (the direct-annihilation cache points).  A photon of rest energy E' from
    an isotropically moving source is uniform in lab energy over
    [gamma (1-beta) E', gamma (1+beta) E'] -- each fine bin's count is spread
    over that box and integrated against the target bins exactly.
    """
    p = np.sqrt(max(m_X * m_X - m_Y * m_Y, 0.0))
    E_Y = np.sqrt(p * p + m_Y * m_Y)
    gamma, beta = E_Y / m_Y, p / E_Y
    centers = np.sqrt(FINE_EDGES[:-1] * FINE_EDGES[1:])
    out = np.zeros(len(ebins) - 1)
    if beta < 1e-12:                      # at rest: plain rebin by centre
        idx = np.searchsorted(ebins, centers, side="right") - 1
        ok = (idx >= 0) & (idx < len(out))
        np.add.at(out, idx[ok], fine_counts[ok])
        return out
    a = gamma * (1.0 - beta) * centers    # box edges per fine bin
    b = gamma * (1.0 + beta) * centers
    w = b - a
    lo = np.maximum(a[:, None], ebins[None, :-1])
    hi = np.minimum(b[:, None], ebins[None, 1:])
    frac = np.clip(hi - lo, 0.0, None) / w[:, None]
    return frac.T @ fine_counts


# ---- the library ------------------------------------------------------------

def scan_legacy():
    """{channel: sorted unique mY list} and the full file inventory."""
    files = sorted(LEGACY_DIR.glob("*.npz"))
    if not files:
        sys.exit(f"no legacy cache at {LEGACY_DIR} -- run --move first")
    inv, by_ch = [], {}
    for f in files:
        m = _FNAME.match(f.name)
        if not m:
            continue
        ch, mY = m["ch"], float(m["mY"])
        inv.append((f, ch, float(m["mX"]), mY))
        by_ch.setdefault(ch, set()).add(mY)
    return {ch: np.array(sorted(v)) for ch, v in by_ch.items()}, inv


def library_tasks():
    by_ch, _ = scan_legacy()
    tasks = []
    for ch, mYs in sorted(by_ch.items()):
        if ch in NO_CORRECTION:
            continue
        # exact span of the inventory, so --migrate never extrapolates
        grid = np.geomspace(mYs.min(), mYs.max(), N_SQRTS)
        grid[0], grid[-1] = mYs.min(), mYs.max()
        thr = sum(MASSES[abs(p)] for p in CHANNEL_TO_PDG[ch])
        tasks += [(ch, float(s)) for s in grid if s > thr * 1.001]
    return tasks


_RUNNER = None


def _lib_worker(task):
    global _RUNNER
    ch, s = task
    if _RUNNER is None:
        _RUNNER = PythiaRunner(seed=LIB_SEED)
    t0 = time.time()
    c = correction_at_rest(_RUNNER, ch, s, N_LIB)
    return ch, s, c, time.time() - t0


def build_library(ncores: int) -> None:
    tasks = library_tasks()
    done = {}
    if LIB_FILE.exists():                 # resume: keep finished points
        with np.load(LIB_FILE, allow_pickle=False) as z:
            if np.array_equal(z["fine_edges"], FINE_EDGES) \
                    and int(z["n_lib"]) == N_LIB:
                for ch, s, c in zip(z["channel"], z["sqrt_s"], z["counts"]):
                    done[(str(ch), float(s))] = c
    todo = [t for t in tasks if t not in done]
    print(f"library: {len(tasks)} (channel, sqrt_s) points, "
          f"{len(done)} cached, {len(todo)} to run at {N_LIB} events each")
    t0 = time.time()

    def save():
        keys = sorted(done)
        np.savez_compressed(
            LIB_FILE, fine_edges=FINE_EDGES, n_lib=N_LIB, seed=LIB_SEED,
            pythia_version=_pythia_version(),
            decayed_long_lived=np.array(LONG_LIVED_DECAYED, dtype=int),
            channel=np.array([k[0] for k in keys]),
            sqrt_s=np.array([k[1] for k in keys]),
            counts=np.array([done[k] for k in keys]))

    with mp.Pool(ncores) as pool:
        for k, (ch, s, c, dt) in enumerate(
                pool.imap_unordered(_lib_worker, todo, 1), 1):
            done[(ch, s)] = c
            print(f"  [{k}/{len(todo)}] {ch:8s} sqrt_s={s:7.2f}  "
                  f"{c.sum():.3f} corr-gamma/decay  ({dt:.0f} s)", flush=True)
            if k % 8 == 0:
                save()
    save()
    print(f"library -> {LIB_FILE}  ({(time.time()-t0)/60:.1f} min)")


class Library:
    """Correction counts interpolated in sqrt(s), per channel."""

    def __init__(self):
        with np.load(LIB_FILE, allow_pickle=False) as z:
            if not np.array_equal(z["fine_edges"], FINE_EDGES):
                sys.exit(f"{LIB_FILE} was built on a different fine grid")
            self.by_ch = {}
            for ch, s, c in zip(z["channel"], z["sqrt_s"], z["counts"]):
                self.by_ch.setdefault(str(ch), []).append((float(s), c))
        for ch in self.by_ch:
            self.by_ch[ch].sort(key=lambda t: t[0])

    def at(self, channel: str, mY: float) -> np.ndarray:
        """Rest-frame correction counts at sqrt(s)=mY (linear in log sqrt(s));
        zero for channels with no correction or below their library floor."""
        pts = self.by_ch.get(channel)
        if pts is None:
            if channel in NO_CORRECTION:
                return np.zeros(len(FINE_EDGES) - 1)
            raise KeyError(f"channel '{channel}' missing from the library")
        ss = [p[0] for p in pts]
        if mY <= ss[0]:
            return pts[0][1] if np.isclose(mY, ss[0], rtol=1e-6) \
                else np.zeros(len(FINE_EDGES) - 1)   # below threshold point
        if mY >= ss[-1]:
            return pts[-1][1]
        k = bisect_right(ss, mY)
        w = (np.log(mY) - np.log(ss[k - 1])) / (np.log(ss[k]) - np.log(ss[k - 1]))
        return (1.0 - w) * pts[k - 1][1] + w * pts[k][1]


# ---- migrate ----------------------------------------------------------------

_LIB = None


def _migrate_worker(args):
    global _LIB
    path_str, ch, mX, mY = args
    if _LIB is None:
        _LIB = Library()
    src = Path(path_str)
    dst = OUT_DIR / src.name
    if dst.exists():
        return ""
    with np.load(src, allow_pickle=False) as z:
        if "decayed_long_lived" in z.files:
            return f"SKIP {src.name}: already new-convention"
        d = {k: z[k] for k in z.files}
    delta = boost_counts(_LIB.at(ch, mY), mX, mY, d["ebins"])
    d["gamma"] = d["gamma"] + delta            # pbar and everything else copied
    d["decayed_long_lived"] = np.array(LONG_LIVED_DECAYED, dtype=int)
    d["decay_correction"] = np.array(
        f"legacy+{LIB_FILE.name} (N_lib={N_LIB}, boosted)")
    tmp = dst.with_suffix(".tmp.npz")
    np.savez(tmp, **d)
    tmp.rename(dst)
    return ""


def migrate(ncores: int) -> None:
    if not LIB_FILE.exists():
        sys.exit(f"missing {LIB_FILE} -- run --library first")
    _, inv = scan_legacy()
    todo = [(str(f), ch, mX, mY) for f, ch, mX, mY in inv
            if not (OUT_DIR / f.name).exists()]
    print(f"migrate: {len(inv)} legacy files, {len(inv) - len(todo)} done, "
          f"{len(todo)} to write -> {OUT_DIR}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with mp.Pool(ncores) as pool:
        for k, msg in enumerate(
                pool.imap_unordered(_migrate_worker, todo, 64), 1):
            if msg:
                print(f"  {msg}", flush=True)
            if k % 5000 == 0 or k == len(todo):
                dt = time.time() - t0
                print(f"  {k}/{len(todo)}  {dt/60:.1f} min elapsed", flush=True)
    print(f"migrate done in {(time.time()-t0)/60:.1f} min")


# ---- move / validate --------------------------------------------------------

def move() -> None:
    if LEGACY_DIR.exists() and any(LEGACY_DIR.iterdir()):
        sys.exit(f"{LEGACY_DIR} already exists and is not empty -- refusing")
    if not SPECTRA_DIR.exists():
        sys.exit(f"nothing to move: {SPECTRA_DIR} does not exist")
    LEGACY_DIR.parent.mkdir(parents=True, exist_ok=True)
    SPECTRA_DIR.rename(LEGACY_DIR)
    SPECTRA_DIR.mkdir()
    print(f"moved {SPECTRA_DIR} -> {LEGACY_DIR}\n"
          f"recreated empty {SPECTRA_DIR}")


def _validate_worker(args):
    path_str, ch, mX, mY = args
    runner = PythiaRunner(seed=NV_SEED)
    mig = Path(path_str)
    with np.load(mig, allow_pickle=False) as z:
        ebins, g_mig = z["ebins"], z["gamma"]
        n_ev = int(z["n_events"])
    fresh = runner.run_channel(mX, mY, ch, NV_EVENTS, ebins)["gamma"]
    # combined MC sigma per bin (counts/decay from n_ev resp. NV_EVENTS events)
    var = g_mig / n_ev + fresh / NV_EVENTS
    ok = var > 0
    pulls = (g_mig[ok] - fresh[ok]) / np.sqrt(var[ok])
    tot_m, tot_f = g_mig.sum(), fresh.sum()
    return (f"{ch:8s} mX={mX:7.2f} mY={mY:7.2f}: total {tot_m:8.3f} vs "
            f"{tot_f:8.3f} ({tot_m/tot_f - 1:+.2%}), "
            f"pull rms {pulls.std():.2f}, max|pull| {np.abs(pulls).max():.1f}")


def validate(ncores: int) -> None:
    _, inv = scan_legacy()
    rng = np.random.default_rng(0)
    with_corr = [t for t in inv if t[1] not in NO_CORRECTION]
    sample = [with_corr[i] for i in
              rng.choice(len(with_corr), size=NV_NODES, replace=False)]
    tasks = []
    for f, ch, mX, mY in sample:
        dst = OUT_DIR / f.name
        if not dst.exists():
            print(f"  not migrated yet, skipping {f.name}")
            continue
        tasks.append((str(dst), ch, mX, mY))
    print(f"validate: {len(tasks)} nodes, fresh {NV_EVENTS} events each "
          f"(seed {NV_SEED}) vs migrated files")
    with mp.Pool(min(ncores, max(len(tasks), 1))) as pool:
        for line in pool.imap_unordered(_validate_worker, tasks):
            print(f"  {line}", flush=True)
    print("expect: totals within ~1% and pull rms ~1 if the migration is "
          "statistically indistinguishable from a fresh run")


# ---- main -------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--move", action="store_true")
    ap.add_argument("--library", action="store_true")
    ap.add_argument("--migrate", action="store_true")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("-n", "--ncores", type=int,
                    default=len(os.sched_getaffinity(0)))
    args = ap.parse_args()
    if not (args.move or args.library or args.migrate or args.validate):
        ap.error("pick a stage: --move / --library / --migrate / --validate")
    if args.move:
        move()
    if args.library:
        build_library(args.ncores)
    if args.migrate:
        migrate(args.ncores)
    if args.validate:
        validate(args.ncores)


if __name__ == "__main__":
    main()

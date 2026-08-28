"""Pythia 8 spectra for one Y -> channel decay, cached on disk.

`load_spectrum` and `run_spectrum` are the entry points; everything above them
exists to serve those two.  Hadron-Level Standalone recipe (Pythia manual):
process level off, daughters injected by hand, shower forced with
forceTimeShower before pythia.next() hadronises and decays.

NOTE: never import this from a directory containing a `python/` subdirectory. 
The pythia8mc wheel execs every .py in a *relative* `python/` at import.
"""

from __future__ import annotations
import math
import warnings
from bisect import bisect_right
from pathlib import Path

import numpy as np
import pythia8mc as pythia8


# ---- 1. particle bookkeeping ---------------------------------------------
"""
 Particles are named by their the Particle Data Group (PDG) code
 integer per species, which Pythia and every other collider tool speaks.
 1-6 are the quarks d, u, s, c, b, t; 11-16 the leptons e, nu_e, mu, nu_mu,
 tau, nu_tau; 21-25 the bosons g, gamma, Z, W, h; baryons get four digits
 (2212 proton, 2112 neutron).  A negative code is the antiparticle.

 Final-state species we count, and the PDG codes each collects.  To add one
 (positrons, -11) put it here and in `_NEGLIGIBLE`: it then needs a fresh
 cache, since existing files hold only what was asked for when they were made.
Note:'pbar' counts anti-p and anti-n together, the cosmic-ray convention of 
PPPC4DMID (n-bar decays to p-bar on propagation timescales).
"""
TARGETS = {"gamma": (22,), "pbar": (-2212, -2112)}

# Masses [GeV] of every daughter CHANNEL_TO_PDG can inject, by PDG code.
MASSES = {
    1: 0.005, 2: 0.002, 3: 0.095, 4: 1.27, 5: 4.18, 6: 172.5,  # d u s c b t
    11: 0.000511, 12: 0.0, 13: 0.10566, 14: 0.0,               # e nu_e mu nu_mu
    15: 1.77686, 16: 0.0,                                      # tau nu_tau
    21: 0.0, 22: 0.0, 23: 91.1876, 24: 80.379, 25: 125.0,      # g gamma Z W h
}

# The two-body final state each annihilation channel injects into Pythia.
CHANNEL_TO_PDG: dict[str, tuple[int, int]] = {
    "ee":      (11, -11),
    "mumu":    (13, -13),
    "tautau":  (15, -15),
    "uu":      (2,  -2),
    "dd":      (1,  -1),
    "ss":      (3,  -3),
    "cc":      (4,  -4),
    "bb":      (5,  -5),
    "tt":      (6,  -6),
    "nunu_e":  (12, -12),
    "nunu_mu": (14, -14),
    "nunu_tau":(16, -16),
    "gg":      (21,  21),
    "gamgam":  (22,  22),
    "WW":      (24, -24),
    "ZZ":      (23,  23),
    "hh":      (25,  25),
}

# Channels that produce none of a given species, so are not worth showering.
_NEGLIGIBLE = {
    "gamma":    {"nunu_e", "nunu_mu", "nunu_tau"},
    "pbar":     {"nunu_e", "nunu_mu", "nunu_tau", "ee", "mumu", "gamgam"},
}

# Channels already warned about.  A grid scan calls `channels_for_species` once
# per point, so without this an unmapped channel would warn thousands of times.
_warned_channels: set[str] = set()


def channels_for_species(brs: dict, species: str) -> list[tuple[str, float]]:
    """(channel, BR) pairs worth showering: BR > 0, mapped, and not negligible."""
    skip = _NEGLIGIBLE.get(species, set())
    out: list[tuple[str, float]] = []
    for ch, br in brs.items():
        if br <= 0.0 or ch in skip:
            continue
        if ch not in CHANNEL_TO_PDG:
            if ch not in _warned_channels:
                _warned_channels.add(ch)
                warnings.warn(
                    f"channel '{ch}' (BR={br:.3g}) has no PDG mapping; dropped.",
                    RuntimeWarning, stacklevel=2)
            continue
        out.append((ch, float(br)))
    return out


# ---- 2. the Pythia engine ------------------------------------------------

def _sample_decay(rng, m_Y: float, m1: float, m2: float, p_Y_lab: float):
    """Sample Y -> d1 d2 isotropically in Y rest frame, boost to lab (Y along +z)."""
    if m1 + m2 >= m_Y:
        raise ValueError(f"channel mass {m1+m2:.3g} >= m_Y={m_Y:.3g}")
    E1r = (m_Y * m_Y + m1 * m1 - m2 * m2) / (2.0 * m_Y)
    E2r = m_Y - E1r
    p_r = math.sqrt(max(E1r * E1r - m1 * m1, 0.0))
    cos_t = rng.uniform(-1.0, 1.0)
    sin_t = math.sqrt(max(1.0 - cos_t * cos_t, 0.0))
    phi = rng.uniform(0.0, 2.0 * math.pi)
    px = p_r * sin_t * math.cos(phi)
    py = p_r * sin_t * math.sin(phi)
    pz_r = p_r * cos_t

    E_Y_lab = math.sqrt(p_Y_lab * p_Y_lab + m_Y * m_Y)
    gamma = E_Y_lab / m_Y
    bg = p_Y_lab / m_Y
    pz1 = gamma * pz_r + bg * E1r
    E1 = gamma * E1r + bg * pz_r
    pz2 = gamma * (-pz_r) + bg * E2r
    E2 = gamma * E2r + bg * (-pz_r)
    return (px, py, pz1, E1, -px, -py, pz2, E2)


def _colour_lines(id1: int, id2: int):
    """Colour-flow tags (col1, acol1, col2, acol2) for the injected pair.

    String fragmentation needs to know which coloured partons are connected.
    Pythia expresses that by giving each colour line an arbitrary index, carried
    as `col` by the parton with the colour and as `acol` by the one with the
    matching anticolour; 0 means neither.  A q qbar pair is joined by one line,
    gg by two, and anything colourless gets no tags.
    """
    if 1 <= abs(id1) <= 6 and id1 == -id2:          # q qbar
        return (101, 0, 0, 101) if id1 > 0 else (0, 101, 101, 0)
    if id1 == 21 and id2 == 21:                     # gg
        return (101, 102, 102, 101)
    return (0, 0, 0, 0)


class PythiaRunner:
    """Hadron-level-standalone Pythia 8 driver."""

    def __init__(self, seed: int = 12345, quiet: bool = True):
        self.seed = int(seed)   # kept: the cache filename encodes it
        self.p = pythia8.Pythia("", False)
        self.p.readString("ProcessLevel:all = off")
        self.p.readString("Random:setSeed = on")
        self.p.readString(f"Random:seed = {int(seed)}")
        for flag in ("TimeShower:QEDshowerByQ", "TimeShower:QEDshowerByL",
                     "TimeShower:QEDshowerByOther"):
            self.p.readString(f"{flag} = on")
        if quiet:
            self.p.readString("Print:quiet = on")
            self.p.readString("Next:numberCount = 0")
            for s in ("Init:showProcesses", "Init:showMultipartonInteractions",
                      "Init:showChangedSettings", "Init:showChangedParticleData"):
                self.p.readString(f"{s} = off")
        if not self.p.init():
            raise RuntimeError("Pythia init() failed")
        # Decay angles come from here, not the global numpy RNG, so that `seed`
        # fixes the whole run rather than just Pythia's half of it.
        self._rng = np.random.default_rng(int(seed))

    def _inject_decay(self, ev, pdg_pair, m_Y: float, p_Y: float):
        id1, id2 = int(pdg_pair[0]), int(pdg_pair[1])
        m1, m2 = MASSES.get(abs(id1), 0.0), MASSES.get(abs(id2), 0.0)
        px1, py1, pz1, E1, px2, py2, pz2, E2 = _sample_decay(
            self._rng, m_Y, m1, m2, p_Y)
        c1, ac1, c2, ac2 = _colour_lines(id1, id2)
        i1 = ev.append(id1, 23, 0, 0, 0, 0, c1, ac1, px1, py1, pz1, E1, m1)
        i2 = ev.append(id2, 23, 0, 0, 0, 0, c2, ac2, px2, py2, pz2, E2, m2)
        pTmax = 0.5 * m_Y
        # Without these scale() calls forceTimeShower is a silent no-op.
        ev[i1].scale(pTmax)
        ev[i2].scale(pTmax)
        self.p.forceTimeShower(i1, i2, pTmax)

    def run_channel(self, m_X: float, m_Y: float, channel: str, n_events: int,
                    ebins, want=tuple(TARGETS)) -> dict[str, np.ndarray]:
        """Run n_events Y -> `channel` decays. Returns {species: counts_per_decay}."""
        if channel not in CHANNEL_TO_PDG:
            raise KeyError(f"channel '{channel}' has no PDG mapping")
        ebins = np.asarray(ebins, dtype=float)
        nb = len(ebins) - 1
        if nb <= 0:
            raise ValueError("ebins must have at least 2 entries")
        for k in want:
            if k not in TARGETS:
                raise ValueError(f"unknown species '{k}'")

        # The clamp at 0 is what makes m_Y = 2 m_X mean "at rest", which
        # is how `spectrum.direct_spectrum` gets a direct annihilation.
        p_Y = math.sqrt(max(m_X * m_X - m_Y * m_Y, 0.0))
        e_lo, e_hi = float(ebins[0]), float(ebins[-1])

        # PDG code -> the species it counts towards, flattened from TARGETS.
        pid_species = {pid: k for k in want for pid in TARGETS[k]}

        edges = ebins.tolist()                  # bisect needs a plain list
        counts = {k: [0.0] * nb for k in want}

        n_ok = 0
        for _ in range(int(n_events)):
            ev = self.p.event
            ev.reset()
            self._inject_decay(ev, CHANNEL_TO_PDG[channel], m_Y, p_Y)
            if not self.p.next():
                continue
            n_ok += 1
            for i in range(ev.size()):
                prt = ev[i]
                species = pid_species.get(prt.id())
                if species is None or not prt.isFinal():
                    continue
                e = prt.e()
                if e_lo <= e < e_hi:
                    counts[species][bisect_right(edges, e) - 1] += 1.0

        denom = float(n_ok if n_ok > 0 else n_events)
        return {k: np.array(v) / denom for k, v in counts.items()}


# ---- 3. cached spectra on disk ------------------------------------------
"""
 One file is the per-decay histogram for one channel at one (m_X, m_Y). It
 depends on exactly (m_X, m_Y, channel, ebins, n_events, seed, Pythia
 version) and nothing else -- not the portal, the couplings, <sigma v>, the
 halo or the ROI, which all enter later as cheap reweightings.

 The filename carries everything but `ebins`, so a directory listing tells
 you what you have. `ebins` is written into the .npz and checked on load: a
 spectrum binned on a different grid raises instead of quietly standing in.

 `load_spectrum` never runs Pythia -- it takes no runner to run one with.
 `run_spectrum` does, and you have to hand it the runner. Nothing starts a
 multi-minute job implicitly.
"""
SPECTRA_DIR = Path(__file__).resolve().parent / "channel_spectra"
DEFAULT_SEED = 12345


def _pythia_version() -> str:
    """Installed pythia8mc version, without importing the module."""
    try:
        from importlib.metadata import version
        return version("pythia8mc")
    except Exception:
        return "unknown"


def spectrum_path(m_X, m_Y, channel, n_events, seed=DEFAULT_SEED,
                  directory=None) -> Path:
    """Cache path, e.g. 'spectra/bb_mX40_mY80_N100000_s12345.npz'."""
    d = Path(directory) if directory is not None else SPECTRA_DIR
    return d / (f"{channel}_mX{m_X:g}_mY{m_Y:g}"
                f"_N{int(n_events)}_s{int(seed)}.npz")


def _binned_differently(d, ebins) -> bool:
    """Does an open .npz disagree with `ebins`? False if it stores none.
    `ebins` is the one part of the cache key the filename cannot carry, so it
    is the one thing to verify on load.
    """
    if "ebins" not in d.files:
        return False
    eb = d["ebins"]
    return eb.shape != ebins.shape or not np.allclose(eb, ebins, rtol=1e-10)


def load_spectrum(m_X, m_Y, channel, n_events, ebins, seed=DEFAULT_SEED,
                  directory=None) -> dict:
    """One cached spectrum, {species: counts_per_decay}.  Never runs Pythia.

    Raises FileNotFoundError naming the file if it has not been generated, and
    ValueError if it exists but was binned on a different grid.
    """
    ebins = np.asarray(ebins, dtype=float)
    path = spectrum_path(m_X, m_Y, channel, n_events, seed, directory)
    if not path.exists():
        raise FileNotFoundError(
            f"no cached spectrum for channel '{channel}' at m_X={m_X:g}, "
            f"m_Y={m_Y:g} with n_events={int(n_events)}, seed={int(seed)}\n"
            f"  expected: {path}\n"
            f"  generate it with run_spectrum(PythiaRunner(seed={int(seed)}), "
            f"{m_X:g}, {m_Y:g}, '{channel}', {int(n_events)}, ebins), or pass "
            f"runner= to whatever asked for it")
    with np.load(path) as d:
        if _binned_differently(d, ebins):
            eb = d["ebins"]
            raise ValueError(
                f"{path.name} was binned differently: file has {len(eb) - 1} "
                f"bins over [{eb[0]:g}, {eb[-1]:g}] GeV, you asked for "
                f"{len(ebins) - 1} over [{ebins[0]:g}, {ebins[-1]:g}] GeV. "
                f"Pass on_mismatch='regenerate' (with a runner) to rebuild it, "
                f"or ask with the file's ebins.")
        return {k: d[k].copy() for k in TARGETS if k in d.files}


def run_spectrum(runner, m_X, m_Y, channel, n_events, ebins, directory=None,
                 overwrite=False, on_mismatch="raise", verbose=True) -> dict:
    """Run Pythia for one channel and cache it, or load it if already there.

    Safe to call in a loop: an existing file is loaded rather than regenerated
    unless `overwrite`.  The seed comes off `runner`, so it cannot disagree with
    the filename.

    A cached file can match the filename and still be binned differently, since
    the name cannot carry `ebins`.  `on_mismatch='raise'` (default) refuses,
    naming both grids; 'regenerate' reruns and overwrites.  'regenerate' fires
    only on an actual mismatch, so it is safe to leave on across a scan -- but
    it does replace a file another analysis may want, which is unavoidable while
    two grids share one name.
    """
    if on_mismatch not in ("raise", "regenerate"):
        raise ValueError(f"on_mismatch must be 'raise' or 'regenerate', "
                         f"got {on_mismatch!r}")
    ebins = np.asarray(ebins, dtype=float)
    path = spectrum_path(m_X, m_Y, channel, n_events, runner.seed, directory)

    if path.exists() and not overwrite:
        with np.load(path) as d:
            stale = _binned_differently(d, ebins)
        if not stale or on_mismatch == "raise":
            if verbose and not stale:
                print(f"  {channel:8s} cached, skipping ({path.name})", flush=True)
            # A stale file under 'raise' falls through to load_spectrum, which
            # reports both grids.
            return load_spectrum(m_X, m_Y, channel, n_events, ebins,
                                 runner.seed, directory)
        if verbose:
            print(f"  {channel:8s} rebinned, regenerating ({path.name})", flush=True)

    hist = runner.run_channel(m_X, m_Y, channel, n_events, ebins)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path,
             ebins=ebins,
             m_X=float(m_X), m_Y=float(m_Y),
             channel=str(channel),
             n_events=int(n_events), seed=int(runner.seed),
             pythia_version=_pythia_version(),
             **hist)
    if verbose:
        print(f"  {channel:8s} wrote {path.name}", flush=True)
    return hist

"""
Per-annihilation spectra, and the Galactic-Centre intensity built on them.

Two annihilation topologies, sharing everything downstream of the spectrum:

  cascade   X Xbar -> Y Y, each Y decaying to the SM through the portal.  The
            branching ratios pick the channels and the cached per-channel
            histograms are summed BR-weighted, times two mediators.
  direct    DM DM -> f fbar at sqrt(s) = 2 m_DM.  The reference arm, and what
            the 2112.09706 contours in `data/contours_2112.09706/` were fit with.

    ebins     = default_ebins(m_DM)                 # 0.02 GeV -> m_DM
    E, dNdE   = direct_spectrum(m_DM, ebins, channel='bb')
    E, dPhidE = direct_flux(m_DM, ebins, sigmav)    # + halo -> intensity

    E, dNdE   = cascade_spectrum(model, ebins)
    E, dPhidE = cascade_flux(model, ebins, sigmav=1.0)

Everything here LOADS ONLY.  Pass `runner=PythiaRunner(seed=...)` to generate
what is missing: Pythia never starts unless a call site hands over the object
that runs it, so nothing kicks off a multi-minute job by accident.
"""

from __future__ import annotations

import numpy as np

from .halo import HALO_2112, j_factor
# Re-exported: the channel bookkeeping lives with the rest of the particle
# bookkeeping in `pythia_runner`, but this stayed its import site.
from .pythia_runner import (  # noqa: F401
    CHANNEL_TO_PDG, DEFAULT_SEED, MASSES, channels_for_species, load_spectrum,
    run_spectrum,
)

# ---- unit conversions ----------------------------------------------------
#
# Everything is assembled in natural units and converted once, on the last line of the flux functions.

# hbar*c [GeV cm] -- the only dimensionful input needed to leave natural units.
HBARC = 1.9732698e-14
# <sigma v>: GeV^-2 -> cm^3/s, i.e. (hbar c)^2 * c.
GEVM2_TO_CM3S = HBARC ** 2 * 2.99792458e10
# Jbar: GeV^2 cm^-5 -> GeV^7, i.e. (hbar c)^5.
JBAR_TO_NATURAL = HBARC ** 5
# Fermi's lowest GCE bin in 2112.09706 sits at ~0.28 GeV; go an order below.
EMIN_DEF = 0.02  # GeV


def default_ebins(m_DM, e_min: float = EMIN_DEF, n_bins: int = 180):
    """Log bins spanning `e_min` to `m_DM`, the kinematic endpoint per annihilation."""
    return np.logspace(np.log10(e_min), np.log10(m_DM), n_bins + 1)


# ---- 1. per-annihilation dN/dE -------------------------------------------

def _dnde(ebins, counts):
    """Counts per bin -> (geometric bin centres, dN/dE)."""
    return np.sqrt(ebins[:-1] * ebins[1:]), counts / np.diff(ebins)


def cascade_spectrum(model, ebins, n_events: int = 200_000,
                     species: str = "gamma", runner=None, seed=DEFAULT_SEED,
                     directory=None, on_mismatch: str = "raise",
                     verbose: bool = False):
    """Per-annihilation dN/dE for `species` ('gamma', 'pbar', 'positron').

    One cached spectrum per open channel, BR-weighted and summed.  With
    `runner=None` an ungenerated channel raises FileNotFoundError naming the
    file; pass a runner to generate it instead, using the runner's seed.

    `ebins` must match the grid the cache was written on, or it raises an error 
    rather than silently rebinning.  `on_mismatch='regenerate'` (needs a runner)
    rebuilds only the channels that actually disagree.

    Returns (E_centers, dNdE), including the factor of 2 for the two mediators.
    """
    ebins = np.asarray(ebins, dtype=float)
    acc = np.zeros(len(ebins) - 1)
    for ch, br in channels_for_species(model.branching_ratios_to_SM(), species):
        # Mediator too light to decay this way -- HDECAY leaves a ~1e-10 off-shell
        # WW*/ZZ* BR at all masses, and the cache (rightly) never holds these.
        if (ch in CHANNEL_TO_PDG
                and sum(MASSES[abs(p)] for p in CHANNEL_TO_PDG[ch]) >= model.mY):
            continue
        if runner is None:
            spec = load_spectrum(model.mX, model.mY, ch, n_events, ebins,
                                  seed, directory)
        else:
            spec = run_spectrum(runner, model.mX, model.mY, ch, n_events,
                                 ebins, directory, on_mismatch=on_mismatch,
                                 verbose=verbose)
        if species not in spec:
            raise KeyError(f"species {species!r} not in the {ch} cache; "
                           f"have {list(spec)}")
        acc += br * spec[species]
        if verbose:
            print(f"  [{species}] {ch:8s} BR={br:.3f}  "
                  f"n_per_decay={float(spec[species].sum()):.2f}", flush=True)
    return _dnde(ebins, 2.0 * acc)


def direct_spectrum(m_DM, ebins, channel: str = "bb", n_events: int = 200_000,
                    species: str = "gamma", runner=None, seed=DEFAULT_SEED,
                    directory=None, on_mismatch: str = "raise",
                    verbose: bool = False):
    """Per-annihilation dN/dE [GeV^-1] for DM DM -> `channel` at rest.

    Cached at (m_X, m_Y) = (m_DM, 2 m_DM).  That is not a physical mediator: a
    particle of mass 2 m_DM at sqrt(s) = 2 m_DM has E = m_DM, less than its own
    mass, and `run_channel` clamps the imaginary momentum to zero.  What comes
    out is a mediator at rest, whose daughters go back-to-back with E = m_DM
    each -- exactly a direct annihilation.  Do not "fix" that clamp without
    fixing this.

    One decay is one annihilation here, so unlike `cascade_spectrum` there is no
    factor of 2.  Returns (E_centers, dNdE).
    """
    ebins = np.asarray(ebins, dtype=float)
    m_Y = 2.0 * m_DM
    if runner is None:
        spec = load_spectrum(m_DM, m_Y, channel, n_events, ebins, seed,
                              directory)
    else:
        spec = run_spectrum(runner, m_DM, m_Y, channel, n_events, ebins,
                             directory, on_mismatch=on_mismatch, verbose=verbose)
    if species not in spec:
        raise KeyError(f"species {species!r} not cached; have {list(spec)}")
    if verbose:
        print(f"  [{species}] {channel:8s} direct  "
              f"n_per_annihilation={float(spec[species].sum()):.2f}", flush=True)
    return _dnde(ebins, spec[species])


# ---- 2. the Galactic-Centre intensity ------------------------------------

def _intensity(E, dNdE, m_DM, sigmav, kappa, halo, roi, units, per_sr):
    """The flux assembly shared by both arms.  `sigmav` in cm^3/s."""
    if units not in ("cm", "natural"):
        raise ValueError(f"units must be 'cm' or 'natural', got {units!r}")
    kw = dict(HALO_2112 if halo is None else halo)
    if roi is not None:                 # the one piece call sites routinely vary
        kw["roi"] = roi
    jbar, domega = j_factor(return_solid_angle=True, **kw)
    if not per_sr:
        jbar = jbar * domega
    # <sigma v> [GeV^-2], m_DM [GeV], dN/dE [GeV^-1], Jbar [GeV^2 cm^-5 sr^-1].
    base = (sigmav / GEVM2_TO_CM3S) / (kappa * np.pi * m_DM ** 2) * jbar * dNdE
    return E, base * (GEVM2_TO_CM3S if units == "cm" else JBAR_TO_NATURAL)


def cascade_flux(model, ebins, sigmav=None, n_events: int = 200_000,
                 species: str = "gamma", majorana=None, units: str = "cm",
                 per_sr: bool = True, halo=None, roi=None, runner=None,
                 seed=DEFAULT_SEED, directory=None, on_mismatch: str = "raise",
                 verbose: bool = False):
    """Galactic-Centre intensity from a secluded-DM model.

        dPhi/dE = <sigma v> / (kappa pi mX^2) * Jbar * dN/dE

    kappa is 8 for self-conjugate DM and 16 for a Dirac pair, taken from
    `model.include_antiparticlesX` unless `majorana` overrides it.

    `sigmav` is in cm^3/s; None falls back to the model's own s-wave value
    (halo velocities are v ~ 1e-3, so the p-wave piece is irrelevant).  Every
    fit in this package passes `sigmav=1.0` and fits the normalisation instead.

    `halo` defaults to `halo.HALO_2112` and `roi` overrides just its window.
    `units` is 'cm' for GeV^-1 cm^-2 s^-1 [sr^-1] or 'natural' for GeV^2.
    `per_sr=True` uses the ROI-averaged Jbar and returns an intensity, as in
    2112.09706; False multiplies by the solid angle for the ROI-integrated J of
    arXiv:2509.08043 Eq. A.2, whose flux carries no sr^-1.

    Returns (E_centers, dPhi_dE).  Only 'gamma' travels in straight lines -- for
    'pbar'/'positron' this is the source term, not an observed flux.
    """
    if majorana is None:
        majorana = not model.include_antiparticlesX
    if sigmav is None:
        sigmav = model.sigmav_XX_to_YY_swave() * GEVM2_TO_CM3S
    E, dNdE = cascade_spectrum(model, ebins, n_events=n_events, species=species,
                               runner=runner, seed=seed, directory=directory,
                               on_mismatch=on_mismatch, verbose=verbose)
    return _intensity(E, dNdE, model.mX, sigmav, 8.0 if majorana else 16.0,
                      halo, roi, units, per_sr)


def direct_flux(m_DM, ebins, sigmav, channel: str = "bb",
                n_events: int = 200_000, species: str = "gamma",
                majorana: bool = True, units: str = "cm", per_sr: bool = True,
                halo=None, roi=None, runner=None, seed=DEFAULT_SEED,
                directory=None, on_mismatch: str = "raise",
                verbose: bool = False):
    """Galactic-Centre intensity for direct DM DM -> `channel`.

    As `cascade_flux`, but with no model to read kappa off: `majorana` is an
    explicit flag, defaulting to True because that is the convention of
    2112.09706.  `sigmav` is required, in cm^3/s
    """
    E, dNdE = direct_spectrum(m_DM, ebins, channel=channel, n_events=n_events,
                              species=species, runner=runner, seed=seed,
                              directory=directory, on_mismatch=on_mismatch,
                              verbose=verbose)
    return _intensity(E, dNdE, m_DM, sigmav, 8.0 if majorana else 16.0,
                      halo, roi, units, per_sr)

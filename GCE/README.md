# GCE — the Galactic-Centre-excess favoured region

An add-on analysis layered on `hidden_sector_DM`.

**The question.** Which dark-matter masses and annihilation cross sections
reproduce the GeV gamma-ray excess towards the Galactic Centre — for direct
annihilation `DM DM -> b bbar`, and for a secluded cascade `X Xbar -> Y Y` with
each `Y` decaying to the SM through a portal.

The prediction is a Pythia-generated spectrum, compared with the measured excess
of [arXiv:2112.09706](https://arxiv.org/abs/2112.09706) over its 14 energy bins.
`<sigma v>` enters the prediction *linearly*, so it profiles out in closed form
and a favoured region is linear algebra on a 14x14 covariance — no sampler, no
minimiser anywhere in the package.

## Layout

```
GCE/
├── pythia_runner.py   Pythia 8 -> one cached per-decay histogram per channel
├── spectrum.py        BR-weighted dN/dE, and the GC intensity built on it
├── halo.py            gNFW profile, ROI-averaged J-factor, halo systematic on <sigma v>
├── gce_data.py        the measured excess and its covariance (2112.09706)
├── likelihood.py      14-bin templates, and the closed-form <sigma v> fit
├── __init__.py        path constants and the curated re-exports
│
├── scripts/           the Pythia-cache builders (run from the repository root)
├── data/              third-party inputs, as downloaded
├── channel_spectra/   the Pythia cache: one .npz per (channel, mX, mY); regeneratable
├── notebooks/         the favoured-region fits and figures, one per portal
└── figs/              .png outputs
```

Modules are meant to be read top to bottom; each uses only the ones above it:
`pythia_runner -> halo -> spectrum -> gce_data -> likelihood`.

## Environment

Python 3.9+ with NumPy, SciPy, tqdm, and Pythia 8 supplied by the `pythia8mc`
wheel (`pip install pythia8mc`; version 8.317). It imports as `pythia8mc` and is
aliased internally to `pythia8`. Put the repository root on `PYTHONPATH` so
`import GCE` and `import hidden_sector_DM` resolve, and **never run from a
directory that contains a `python/` subdirectory** — the `pythia8mc` wheel execs
every `.py` in a relative `python/` at import time.

## The pipeline

**`pythia_runner.py`** — `PythiaRunner.run_channel(mX, mY, channel, n_events,
ebins)` injects one `Y -> channel` two-body decay by hand, boosts `Y` by
`p = sqrt(mX^2 - mY^2)` (clamped to zero, so `mY = 2 mX` means "at rest"), forces
the time shower, hadronises, and histograms the final-state photons and
antiprotons **per decay**. `run_spectrum` / `load_spectrum` cache the result as
one `.npz` per `(channel, mX, mY, n_events, seed)` under `channel_spectra/`.
`load_spectrum` never starts Pythia; `run_spectrum` does, and only if you hand
it a runner.

**`spectrum.py`** — per-annihilation `dN/dE`, then the GC intensity.
`cascade_spectrum(model, ebins)` sums the cached per-channel histograms
BR-weighted (`model.branching_ratios_to_SM()`), times two mediators.
`direct_spectrum(m_DM, ebins, channel='bb')` is the reference arm,
`DM DM -> channel` at rest, cached at `(m_DM, 2 m_DM)`. `cascade_flux` /
`direct_flux` multiply by `<sigma v> / (kappa pi mX^2) * Jbar`.
`default_ebins(m_DM)` is 180 log bins from 0.02 GeV to `m_DM`.

**`halo.py`** — the gNFW profile of 2112.09706 Eq. 6, its ROI-averaged J-factor
`j_factor(...)`, and `sigmav_halo_range(...)`, which turns the `rho_local`
uncertainty (`RHO_LOCAL_PRIOR = 0.44 +/- 0.13`) into a multiplicative range on
the inferred `<sigma v>`. `HALO_2112` is the fixed profile every template is
built at; `ROIS` holds the 40x40 deg window and its north / south halves, with
`|b| <= 2 deg` masked.

**`gce_data.py`** — `load_sed(model, region)` and `load_covariance(...)` return
the excess `E^2 dPhi/dE` on the 14 Table III bins and the Eq. 18 covariance
(statistical diagonal plus truncated interstellar-emission-model systematics).
`region` is `'both'`, `'north'` or `'south'`.

**`likelihood.py`** — the fit.
```python
phi, Cinv = measurement(region="south")     # once; reuse across points
t         = cascade_template(model)         # or direct_template(m_DM)
sigmav_hat, sigmav_err, chi2_min = gls_fit(t, phi, Cinv)
```
`cascade_template` / `direct_template` give the 14-bin `E^2 dPhi/dE` at
`<sigma v> = 1`, built at `HALO_2112`. `gls_fit` profiles `<sigma v>` out
analytically and clamps it at `>= 0`; for an interior fit
`chi2(s) = chi2_min + ((s - sigmav_hat) / sigmav_err)^2`.

**`notebooks/`** — `GCE_fit_<portal>.ipynb` (`vP`, `bP`, `blP`, `hP`) scans
`cascade_template` x `gls_fit` over the `(mX, mY)` triangle grid, then folds the
`rho_local` prior into the `(mX, <sigma v>)` contour.
`GCE_fit_validation_bb_annihilation.ipynb` is the closure test against the
`b bbar` contour of 2112.09706 Fig. 18; `spectrum_validation.ipynb` reproduces
the Fig. 4 spectra of [arXiv:1912.08821](https://arxiv.org/abs/1912.08821).

## Building the cache

Both scripts write into `channel_spectra/`, run one process per core, and are
resumable — an existing `.npz` is loaded, not regenerated:

```bash
# from the repository root
python GCE/scripts/make_channel_spectra.py     # (mX, mY) triangle grid, every channel at BR = 1
python GCE/scripts/make_bb_reference.py        # direct b bbar at 50 DM masses
```

Edit the module-level configuration blocks (`MX_LIM`, `N_POINTS`, `N_EVENTS`,
`CHANNELS`, ...) rather than adding CLI flags. The favoured region itself is
built in the notebooks, not here.

## Using it as a library

```python
import numpy as np
from hidden_sector_DM.HiddenSectorDM import VectorPortal
from GCE.pythia_runner import PythiaRunner
from GCE.spectrum import cascade_spectrum, default_ebins

model = VectorPortal(mX=76, mY=60, gX=2, gY=3, alphaX=1e-2)
ebins = default_ebins(model.mX)

# Loads only. A channel that is not cached raises, naming the file it wanted.
E, dNdE = cascade_spectrum(model, ebins, n_events=200_000)

# Pass a runner to generate whatever is missing, then carry on as before.
E, dNdE = cascade_spectrum(model, ebins, n_events=200_000,
                           runner=PythiaRunner(seed=12345))
```

The filename carries everything but the bin edges; `ebins` lives inside the
`.npz` and is checked on load, so a spectrum binned on a different grid raises
rather than quietly standing in. With a runner in hand,
`on_mismatch='regenerate'` rebuilds the stale entries instead.

## Conventions that bite

- **`'pbar'` counts both anti-protons (-2212) and anti-neutrons (-2112)** — the
  cosmic-ray convention of 1912.08821 and PPPC4DMID.
- **`kappa`** is 8 for self-conjugate DM (the 2112.09706 convention, used by the
  direct arm) and 16 for a Dirac pair (the portal models), read off
  `model.include_antiparticlesX`. A portal region therefore sits a factor two
  above a Majorana comparison at equal flux — say so wherever the two overlap.
- **`alphaX` does not enter the fit.** Branching ratios depend on `mY` alone and
  templates are built at `<sigma v> = 1`, so no coupling reaches the prediction.
- **The halo is fixed inside a contour.** Templates are built at `HALO_2112`;
  `rho_local` uncertainty is applied afterwards, to `<sigma v>` only. Present it
  as a contour plus a shaded band, not a merged region.
- **`default_ebins` must reach below the lowest GCE bin** (0.275 GeV) once
  rebinned, or `cascade_template` / `direct_template` raise on the empty bin.
  `e_min = 0.02 GeV` has the margin.
- **Sub-GeV hadronic `HiggsPortal` channels** (`pi0pi0`, `K0K0`, ...) have no
  PDG mapping in `CHANNEL_TO_PDG` and are dropped with a one-time `RuntimeWarning`.

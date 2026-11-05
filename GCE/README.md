# GCE — the Galactic-Centre-excess favoured regions

An analysis layered on the `hidden_sector_DM` toolkit. It asks which dark-matter masses and
annihilation cross sections reproduce the GeV gamma-ray excess towards the Galactic Centre, for the
secluded cascade $X\bar X \to YY$ with each mediator $Y$ decaying to the Standard Model through
one of four portals — hypercharge (`vP`), baryonic (`bP`), $B-L$ (`blP`) and Higgs (`hP`) — and
where those regions meet the observed relic abundance $\Omega_c h^2 = 0.12$.

The prediction is a Pythia 8 spectrum, compared with the excess measured by
[arXiv:2112.09706](https://arxiv.org/abs/2112.09706) (south ROI) over its 14 energy bins against
the full Eq. 18 covariance. $\langle\sigma v\rangle$ enters the prediction linearly, so it is
profiled out in closed form and a favoured region is linear algebra on a 14×14 matrix: there is
no sampler and no minimiser anywhere in the package.

## Layout

```
GCE/
├── source/            the library
│   ├── pythia_runner.py   Pythia 8 -> one cached per-decay histogram per channel
│   ├── spectrum.py        BR-weighted dN/dE, and the GC intensity built on it
│   ├── halo.py            gNFW profile, ROI-averaged J-factor, the rho_sun prior
│   ├── gce_data.py        the measured excess and its covariance (2112.09706)
│   └── likelihood.py      14-bin templates and the closed-form <sigma v> fit
├── scripts/           the production drivers, run from the repository root
│   ├── make_channel_spectra.py   Pythia cache over the (mX, mY) grid
│   ├── make_bb_reference.py      direct DM DM -> b bbar spectra, for the closure test
│   ├── make_gce_fits.py          free-<sigma v> fit of every portal over the grid
│   └── make_relic_grid.py        relic coupling and <sigma v>, --portal vP|bP|blP|hP
├── notebooks/
│   ├── paper_figs.ipynb                           every portal figure of the paper
│   └── GCE_fit_validation_bb_annihilation.ipynb   closure test against 2112.09706 Fig. 18
├── data/              third-party inputs (below)
├── channel_spectra/   the Pythia cache: one .npz per (channel, mX, mY)
├── output/
│   ├── fits/          gce_fit_<portal>_south.npz
│   └── relic/         relic_grid_<portal>[_xi_ini_<x>]_{thermal,eps<eps>,eps_secluded}.npz
└── paper_figs/        the figures, as PNG and PDF
```

The modules are meant to be read top to bottom; each uses only the ones above it:
`pythia_runner -> halo -> spectrum -> gce_data -> likelihood`.

## Environment

Python ≥ 3.9 with NumPy, SciPy, Matplotlib (with LaTeX) and tqdm, plus Pythia 8 from the
`pythia8mc` wheel (`pip install pythia8mc`, version 8.317), which imports as `pythia8mc`. The relic
grids also need the `hidden_sector_DM` package and, for the Higgs portal, its compiled HDECAY
(`bash hidden_sector_DM/phi_decay/setup_phi_decay.sh`).

Run scripts **from the repository root**; they put it on `sys.path` themselves. **Never run from a
directory that contains a `python/` subdirectory**: the `pythia8mc` wheel executes every `.py` file
in a relative `python/` at import time.

## Reproducing the paper

Every product the notebooks read ships with the repository, so the figures can be remade
straight away. To rebuild a product from scratch, run its stage; each stage only needs the ones
above it.

| Stage | Command | Writes | Cost |
|---|---|---|---|
| 1. Pythia cache | `python GCE/scripts/make_channel_spectra.py` | `channel_spectra/` (162,373 files) | ~5000 core-hours |
| 2. Direct $b\bar b$ | `python GCE/scripts/make_bb_reference.py` | `channel_spectra/bb_mX<m>_mY<2m>_*` | 50 Pythia runs |
| 3. Fits | `python GCE/scripts/make_gce_fits.py` | `output/fits/` | loads only |
| 4. Relic grids | `python GCE/scripts/make_relic_grid.py --portal <p> [--xi-ini <x>]` | `output/relic/` | cluster scale |
| 5. Figures | `notebooks/paper_figs.ipynb`, `notebooks/GCE_fit_validation_bb_annihilation.ipynb` | `paper_figs/` | minutes |

- **Stage 1** is resumable: an existing `.npz` is loaded, never regenerated, so a job killed at the
  wall clock is just relaunched. Across cluster ranks, rank `w` of `W` takes every `W`-th grid
  column: `--worker $RANK --nworkers $NRANKS -n 128`. Set `CHANNEL_SPECTRA_DIR` to write somewhere
  faster than the repository, and `--check` lists what is still missing.
- **Stage 4** runs one portal at a time and, without `--only`, all of its outputs: the thermal relic
  (hidden sector locked to the SM), the joint-solver eps, and the secluded eps. The vector portal is
  also run at `--xi-ini 0.5` and `2.5`. Each output can be its own job (`--only thermal`,
  `--only secluded`, `--only 1e-9`, ...) and split across jobs with `--nshards/--shard/--merge`;
  `--dry-run` reports what is left. The joint and secluded files double as checkpoints.
- Scripts are configured by the constants at the top of each file (grids, event counts, seeds,
  portal settings). Edit those rather than adding flags; the flags say how to run, not what to
  compute.

### Provenance of the shipped relic grids

Every thermal and every secluded grid, and the Higgs-portal joint grids, are
`make_relic_grid.py`'s output. The joint grids of the three Z' portals were solved by earlier
revisions of the relic code (October 2026) and are kept as they are; their mesh and coupling axis
are the script's, so they load and resume, but a rerun will not reproduce them bit for bit. The
docstring of `make_relic_grid.py` records exactly what differs.

## The pipeline

**`pythia_runner.py`** — `PythiaRunner.run_channel(mX, mY, channel, n_events, ebins)` injects one
two-body decay `Y -> channel` by hand, boosts `Y` by $p = \sqrt{m_X^2 - m_Y^2}$ (clamped to zero,
so `mY = 2 mX` means "at rest", i.e. a direct annihilation), forces the time shower, hadronises,
and histograms the final-state photons and antinucleons **per decay**. `run_spectrum` /
`load_spectrum` cache the result as `<channel>_mX<mX>_mY<mY>_N<n_events>_s<seed>.npz`; the energy
bins are stored inside and checked on load. `load_spectrum` never starts Pythia; `run_spectrum`
does, and only with a runner handed to it.

**`spectrum.py`** — `cascade_spectrum(model, ebins)` sums the cached per-channel histograms
weighted by `model.branching_ratios_to_SM()`, times two mediators; `direct_spectrum(m_DM, ebins)`
is the same chain with the mediator at rest. `cascade_flux` / `direct_flux` multiply by
$\langle\sigma v\rangle/(\kappa\pi m^2)\,\bar J$. `default_ebins(m)` is 180 log bins from
0.02 GeV to `m`.

**`halo.py`** — the gNFW profile of 2112.09706 Eq. 6, its ROI-averaged J-factor `j_factor`, and
`RHO_LOCAL_PRIOR = (0.44, 0.13)` GeV/cm³. Every template is built at the fixed halo `HALO_2112`
($\gamma = 1.2$, $r_c = 20$ kpc, $\rho_\odot = 0.44$ GeV/cm³, $r_\odot = 8.5$ kpc). The one
exception is the closure test, which compares with 2112.09706 at their own
`RHO_LOCAL_2112 = 0.4` GeV/cm³.

**`gce_data.py`** — `load_sed(model, region)` and `load_covariance(...)` give the excess
$E^2\,d\Phi/dE$ on the 14 Table III bins and the Eq. 18 covariance: the statistical diagonal plus
the three leading principal components of $\Sigma_{ij,\rm mod}$.

**`likelihood.py`** — the fit:

```python
phi, Cinv = measurement(region='south')     # once; reuse across points
t = cascade_template(model)                  # or direct_template(m_DM)
sigmav_hat, sigmav_err, chi2_min = gls_fit(t, phi, Cinv)
```

Templates are the 14-bin $E^2\,d\Phi/dE$ at $\langle\sigma v\rangle = 1$ cm³/s. `gls_fit` profiles
$\langle\sigma v\rangle$ out analytically and clamps it at ≥ 0; at an interior fit
$\chi^2(s) = \chi^2_{\min} + ((s - \hat s)/\sigma_s)^2$.

**`notebooks/paper_figs.ipynb`** reads `output/fits/` and `output/relic/`, builds the favoured
regions in $(m_X, m_Y)$ and $(m_X, \langle\sigma v\rangle)$ with $\rho_\odot$ fixed and profiled
over its prior, overlays the relic condition at each eps, and draws the best-fit spectra.

## Using it as a library

```python
import sys
sys.path.insert(0, '.')                       # the repository root

from hidden_sector_DM.HiddenSectorDM import VectorPortal
from GCE.source.pythia_runner import PythiaRunner
from GCE.source.spectrum import cascade_spectrum, default_ebins

model = VectorPortal(mX=76, mY=60, gX=2, gY=3, alphaX=1e-2)
ebins = default_ebins(model.mX)

# Loads only: a channel that is not cached raises, naming the file it wanted.
E, dNdE = cascade_spectrum(model, ebins, n_events=400_000)

# Pass a runner to generate whatever is missing, then carry on as before.
E, dNdE = cascade_spectrum(model, ebins, n_events=400_000, runner=PythiaRunner(seed=12345))
```

## Data

- `data/zenodo_2112.09706/` — the files of [Zenodo 6423495](https://zenodo.org/records/6423495)
  (Cholis et al., arXiv:2112.09706) this code reads: the best-fit GCE spectrum in the full window
  and its north and south halves, the systematic covariance `cov_mat_21Dec02.npy`, and the
  statistical errors. The other emission models and the diffuse-emission maps are on Zenodo;
  copy a model's spectrum into `Figures_12_and_14_GCE_Spectra/` to use it with `load_sed`.
- `data/contours_2112.09706/` — the $b\bar b$ region of 2112.09706 Fig. 18 (south ROI), digitised,
  for the closure test.

## Conventions that bite

- **`'pbar'` counts antiprotons and antineutrons** (PDG −2212 and −2112): the antineutron decays to
  an antiproton on propagation timescales, the convention of 1912.08821 and PPPC4DMID.
- **Long-lived species decay in flight.** μ±, π±, K± and K_L are decayed by Pythia
  (`LONG_LIVED_DECAYED`), since over kpc every one of them does; a cached file records this in its
  `decayed_long_lived` key.
- **`forceTimeShower` is a silent no-op** unless every injected daughter first gets
  `ev[i].scale(pTmax)`.
- **Heavy-quark channels open at the open-flavour threshold**, $2m_B$ for $b\bar b$, not at
  $2m_b$. `HADRON_MASSES` here and `_OPEN_FLAVOUR` in `hidden_sector_DM/model.py` must stay in step.
- **$\kappa$** is 8 for self-conjugate DM (the 2112.09706 convention, used by the direct arm) and 16
  for a Dirac pair, read from `model.include_antiparticlesX`. A Dirac portal region therefore sits a
  factor two above a Majorana comparison at equal flux.
- **The hemisphere covariance is a prescription.** Zenodo ships $\Sigma_{ij,\rm mod}$ for the full
  window only. `sigma_mod_scaling='paper'` (the default) scales it by the paper's 0.55² (north) or
  0.45² (south); `'flux'` rescales it bin by bin to the hemisphere's spectrum and runs 1.5–4×
  wider in the south.
- **`'BestFit'` is a different emission model in each region**: XLIX in the full window, XV in the
  north, XLVIII in the south.
- **The coupling never reaches the fit.** Branching ratios depend on `mY` alone and templates are
  built at $\langle\sigma v\rangle = 1$, so $\alpha_X$ (or $\lambda$) is an overlay on a finished
  region.
- **`n_events` is part of the cache key**, so a cache mixed across event counts cannot be read by
  one call. The grid is homogeneous at 400k events, seed 12345.
- **Sub-GeV hadronic Higgs-portal channels** (`pi0pi0`, `K0K0`, ...) have no PDG mapping and are
  dropped with a one-time `RuntimeWarning`.

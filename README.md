# Hidden Sector Dark Matter

Toolkit for computing the relic abundance of hidden-sector dark matter, where the DM particle $X$ annihilates into unstable mediators $Y$ that subsequently decay to the Standard Model through a small portal coupling $\varepsilon$. The defining hierarchy is $g_D \gg \varepsilon^2$, such that the hidden-sector freeze-out is controlled by the dark coupling $g_D$, while $\varepsilon$ only sets the mediator lifetime and hence the late-time entropy dilution.

The code solves the coupled Boltzmann equations for the DM yield $Y_X$, the mediator yield $Y_Y$, and the dark-to-SM temperature ratio $\xi = T_h / T_{\rm SM}$, using a three-phase strategy that exploits the natural hierarchy of freeze-out timescales:
1. **Phase 1 (equilibrium)** — both $X$ and $Y$ track their equilibrium abundances; only $\ln\xi$ evolves.
2. **Phase 1.5 (mediator-in-equilibrium)** — $X$ departs from equilibrium while cannibal reactions ($3\to2$) keep $Y$ pinned at its equilibrium abundance ($\bar\mu_Y = 0$); solver evolves $(\ln\xi,\bar\mu_X)$.
3. **Phase 2 (full)** — cannibal processes decouple and the full three-variable system $(\ln\xi,\bar\mu_X,\bar\mu_Y)$ is integrated to freeze-out.

Once the Boltzmann system has converged, a separate background solver evolves the mediator decay $Y\to\text{SM}$ in $e$-folds, tracking entropy injection and the resulting dilution of the DM yield. This **secluded** workflow ($\varepsilon \ll \varepsilon_{\rm therm}$, with $Y$ treated as stable through freeze-out) is the primary product of the code and what almost all bounds in the documentation are derived from.

For larger portal couplings, where $\Gamma_Y$ becomes comparable to Hubble during freeze-out and the secluded approximation breaks down, a complementary **joint** solver (`solve_boltzmann_joint`) is provided. It evolves the dark sector and the portal-induced $Y\leftrightarrow\text{SM}$ exchanges together throughout freeze-out, splitting Phase 1.5 into a dark-locked branch (1.5-D, $\xi$ free) and an SM-locked branch (1.5-S, $\xi=1$ pinned, entered automatically when $\Gamma_Y\,n_Y^{\rm eq}(T_{\rm SM})/(H\,n_Y) \geq 1$). The joint solver hands off to the same background evolution and covers the full range up to and including the thermalisation floor, smoothly recovering the secluded result in the small-$\varepsilon$ limit.

A full technical description — Boltzmann derivation, cross sections, numerics, and validity — is available in [`docs/hidden_sector_DM_toolkit.pdf`](docs/hidden_sector_DM_toolkit.pdf).

## Models

Five portal scenarios are implemented:

| Model | DM | Mediator | Portal coupling | Reference |
|-------|-----|----------|----------------|-----------|
| `VectorPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with hypercharge | [arXiv:1602.08490](https://arxiv.org/abs/1602.08490) |
| `BLPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with $U(1)_{B-L}$ | [arXiv:1912.08821](https://arxiv.org/abs/1912.08821) |
| `LiLjPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with $U(1)_{L_i - L_j}$ | [arXiv:1912.08821](https://arxiv.org/abs/1912.08821) |
| `BaryonPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with $U(1)_B$ | [arXiv:1912.08821](https://arxiv.org/abs/1912.08821) |
| `HiggsPortal` | Majorana fermion | Dark scalar (spin-0) | Scalar mixing with SM Higgs | [arXiv:1609.02555](https://arxiv.org/abs/1609.02555) |

All models inherit from a common `HiddenSectorModel` interface providing:
- $2\to2$ annihilation cross section ($X\bar X \to YY$), split into $s$- and $p$-wave coefficients
- $3\to2$ cannibal cross sections ($YYY \to YY$, $YYX \to YX$, $YXX \to XX$)
- Mediator decay width $\Gamma(Y \to \text{SM})$
- Spin-independent DM–nucleon cross section $\sigma_{\rm SI}$

The DM symmetry factor (`include_antiparticlesX = True/False`) is honoured throughout: Dirac models pick up a $1/2$ for $X\bar X \to YY$ annihilations, Majorana models do not (HiggsPortal is Majorana by default; the vector / $B$-$L$ / $B$ / $L_iL_j$ portals are Dirac).

For the `BLPortal`, `BaryonPortal`, and `LiLjPortal`, the gauge couplings (`g_BL`, `g_B`, `g_LiLj`) default to $1$, so `eps` represents the effective portal coupling product $\varepsilon \cdot g$.

## Usage

### Secluded workflow (primary)

```python
from hidden_sector_DM.HiddenSectorDM import Cosmology, VectorPortal, BoltzmannSolver

# Set up cosmology and model
cosmo = Cosmology(gstar_choice="standard", gstarpath="./notebooks/")
model = VectorPortal(mX=50.0, mY=30.0, alphaX=1e-2)

# Solve Boltzmann equations (three-phase QSSA in chemical-potential variables)
solver = BoltzmannSolver(cosmo, model)
sol = solver.solve_boltzmann_chempot_3phase(return_bg_ICs=True, verbose=True)
print(f"YX_frozen = {sol['YX_final']:.4e}")

# Late-time Y -> SM decay and entropy injection
bg = solver.solve_background(epsX=1e-9, bg_ICs=sol['bg_ICs'],
                             correct_Y_decay=True)
```

### Joint workflow (large portal coupling, optional)

For $\varepsilon$ near or above the thermalisation floor, use the joint solver and its mass-finder counterpart `find_mX_DMRelic_joint` (which fixes $\varepsilon$ and root-finds $m_X$ on the relic abundance):

```python
from hidden_sector_DM.HiddenSectorDM import find_mX_DMRelic_joint

factory = lambda mX, mY: VectorPortal(mX=mX, mY=mY, alphaX=1e-2)
res = find_mX_DMRelic_joint(cosmo, factory, epsX=1e-7, r=2.5, verbose=True)
print(f"mX = {res['mX']:.4e} GeV   GoH_handoff = {res['GoH_handoff']:.2e}")
```

The reported `GoH_handoff` ($\Gamma_Y/H$ at the joint $\to$ background handoff) is the natural diagnostic for when the joint workflow is no longer needed: once it drops well below $1$, $Y$ decay is sub-Hubble during freeze-out and the secluded workflow above is the appropriate (and more efficient) tool. For the fully thermalised limit ($\varepsilon \gg \varepsilon_{\rm therm}$, $\xi=1$ throughout), the standalone solver `solve_boltzmann_thermSM` is also exposed; it assumes $Y$ in detailed balance with the SM from the outset and integrates a single ODE for the comoving DM yield $Y_X$.

### Mapping the portal coupling $\varepsilon$

Once the hidden-sector evolution is solved, the portal coupling $\varepsilon$ can be mapped against three independent constraints (relic abundance, BBN hadronic injection, direct detection) plus a heuristic thermalisation-floor estimate:

```python
# 1. Relic abundance: eps that gives Omega_DM h^2 = 0.12
eps_relic = solver.find_epsilon_DMRelic(sol['bg_ICs'])

# 2. BBN hadronic injection: minimum eps satisfying Kawasaki et al. (2017)
eps_bbn = solver.find_epsilon_BBN_hadronic(sol['YY_final'])

# 3. Direct detection: eps that saturates a SI cross-section limit
eps_LZ  = solver.find_epsilon_DD(constraint='LZ')          # LZ 90% CL
eps_nu  = solver.find_epsilon_DD(constraint='nufloor_Xe')  # Xe neutrino floor

# 4. Thermalization floor: minimum eps to keep hidden sector decoupled consistently
eps_therm = solver.find_epsilon_thermal_floor()            # default x_FO = 20
```

**Relic abundance.** `find_epsilon_DMRelic` runs the Boltzmann solver once, then root-finds on `log10(eps)` (Brent's method over $[-20, -6]$) for the value that reproduces $\Omega_c h^2 = 0.12$ via the background evolution. A unique solution exists only when the frozen yield over-produces; otherwise no amount of entropy dilution can recover the observed abundance. The complementary entry point `find_mX_DMRelic_joint` (joint workflow) fixes $\varepsilon$ and root-finds $m_X$ instead.

**BBN.** Uses 2D log-log interpolation over the Kawasaki, Kohri, Moroi & Takaesu (2017) upper bounds on $m_Y \cdot Y_Y$ as a function of lifetime, combining `uu` and `bb` injection channels weighted by the model's branching ratios.

**Direct detection.** Root-finds the $\varepsilon$ for which `model.sigma_SI(eps)` equals the experimental limit at the given $m_X$. Limit curves (LZ, XENONnT, Xe neutrino floor) are loaded as log-log interpolators from tabulated data. Returns `np.inf` for models with no tree-level SI cross section (e.g. `LiLjPortal`), and `np.nan` if $m_X$ is outside the data range.

**Thermalisation floor.** *Heuristic estimate*, not a solver result. Returns the analytic value $\varepsilon_{\rm therm} = \sqrt{H(T_{\rm FO})/\Gamma(\varepsilon{=}1)}$ from the criterion $\Gamma(Y\to\text{SM}) \geq H(T_{\rm FO})$ at $T_{\rm FO} = m_X / x_{\rm FO}$. Use as a seed/anchor; for an accurate boundary in the dynamical regime, run the joint solver and inspect `GoH_handoff`.

## Structure

```
hidden_sector_DM/         # Python package
  HiddenSectorDM.py       # Public API (re-exports)
  model.py                # HiddenSectorModel + VectorPortal, BLPortal, LiLjPortal, BaryonPortal, HiggsPortal
  solver.py               # BoltzmannSolver: secluded (3phase / 2phase legacy), joint, thermSM,
                          #   solve_background, find_epsilon_*, find_mX_DMRelic_joint
  cosmology.py            # Cosmology (g_star tables, thermodynamics, Bessel moments)
  constants.py            # SM parameters, fermion tables
  bbn.py                  # BBN hadronic-injection constraints (Kawasaki et al. 2017)
  direct_detection.py     # Direct-detection limit interpolators (LZ, XENONnT, Xe neutrino floor)
  phi_decay/              # Mediator decay widths (2HDECAY + scalar_portal)
  data/
    bbn_constraints_decay/  # BBN bound tables (uu, bb channels, reference masses)
    direct_detection_data/  # LZ / XENONnT SI 90% CL, neutrino floors (Xe SI, Ar SI, Xe SD)
notebooks/                # Scripts to reproduce the figures of our pre-print
  letter_vp_figures.ipynb # VectorPortal figures
  bP_figure.ipynb         # BaryonPortal figures
  blP_figure.ipynb        # BLPortal figures
  hP_figure.ipynb         # HiggsPortal figures
scripts/
  run_relic.py              # Single-point relic-abundance runner
  get_vPDM_bounds_all.py    # VectorPortal parameter scans (joint + secluded stages)
  get_bLPDM_bounds_all.py   # BLPortal parameter scans
  get_BPDM_bounds_all.py    # BaryonPortal parameter scans
  get_hPDM_bounds_all.py    # HiggsPortal parameter scans
  branching_ratios_pie.py   # Mediator branching-ratio plots
mathematica/                # FeynCalc notebooks deriving the 2->2 and 3->2 cross sections
docs/
  hidden_sector_DM_toolkit.pdf  # Full technical documentation
```

## Dependencies

- Python $\geq$ 3.9
- NumPy, SciPy, tqdm
- [2HDECAY](https://github.com/marcel-krause/2HDECAY) (bundled in `hidden_sector_DM/phi_decay/` for HiggsPortal above 5 GeV)
- [scalar_portal](https://github.com/JLTastet/scalar_portal) (for HiggsPortal below 5 GeV)

## Citation

If you use this code, please cite:

```bibtex
@article{Hooper:2026iga,
    author = "Hooper, Dan and Krnjaic, Gordan and Montefalcone, Gabriele",
    title = "{WIMP-like Dark Matter Without Thermalization At Freeze-Out}",
    eprint = "2605.xxxxx",
    archivePrefix = "arXiv",
    primaryClass = "hep.ph",
    reportNumber = "UTWI-16-2026",
    month = "5",
    year = "2026"
}
```

## AI Disclosure

The preparation of this public repository was aided by Claude (Anthropic). Specifically, AI assistance was used to organize the directory structure for public release, add descriptive comments and markdown sections to the notebooks, and produce an initial draft of this README. All AI-generated content was reviewed and verified by the authors.

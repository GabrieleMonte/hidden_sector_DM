# Hidden Sector Dark Matter

Toolkit for computing the relic abundance of hidden-sector dark matter, where the DM particle $X$ annihilates into unstable mediators $Y$ that subsequently decay to the Standard Model through a small portal coupling $\varepsilon$. The defining hierarchy is $\alpha_X \gg \varepsilon^2$: hidden-sector freeze-out is controlled by the dark coupling, while $\varepsilon$ only sets the mediator lifetime and hence the late-time entropy dilution.

The code solves the coupled Boltzmann equations for the DM yield $Y_X$, the mediator yield $Y_Y$, and the dark-to-SM temperature ratio $\xi = T_h / T_{\rm SM}$, using a three-phase strategy that exploits the natural hierarchy of freeze-out timescales:
1. **Phase 1 (equilibrium)** — both $X$ and $Y$ track their equilibrium abundances; only $\ln\xi$ evolves.
2. **Phase 1.5 (QSSA)** — $X$ departs from equilibrium while cannibal reactions ($3\to2$) keep $Y$ in chemical equilibrium ($\bar\mu_Y = 0$); solver evolves $(\ln\xi,\bar\mu_X)$.
3. **Phase 2 (full)** — cannibal processes decouple and the full three-variable system $(\ln\xi,\bar\mu_X,\bar\mu_Y)$ is integrated to freeze-out.

Once the Boltzmann system has converged, a separate background solver evolves the mediator decay $Y\to\text{SM}$ in $e$-folds, tracking entropy injection and the resulting dilution of the DM yield.

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

For the `BLPortal`, `BaryonPortal`, and `LiLjPortal`, the gauge couplings (`g_BL`, `g_B`, `g_LiLj`) default to $1$, so `eps` represents the effective portal coupling product $\varepsilon \cdot g$.

## Usage

```python
from hidden_sector_DM.HiddenSectorDM import Cosmology, VectorPortal, BoltzmannSolver, find_mXstar

# Set up cosmology and model
cosmo = Cosmology(gstar_choice="standard", gstarpath="./notebooks/")
model = VectorPortal(mX=50.0, mY=30.0, alphaX=1e-2)

# Solve Boltzmann equations (three-phase QSSA in chemical-potential variables)
solver = BoltzmannSolver(cosmo, model)
sol = solver.solve_boltzmann_chempot_3phase(return_bg_ICs=True, verbose=True)
print(f"YX_frozen = {sol['YX_final']:.4e}")

# Find critical mass mX* (minimum DM mass below which no epsilon reproduces Omega_DM h^2)
factory = lambda mX, mY: VectorPortal(mX=mX, mY=mY, alphaX=1e-2)
res = find_mXstar(cosmo, factory, r=1.5, verbose=True)
print(f"mXstar = {res['mXstar']:.4e} GeV")
```

### Mapping the portal coupling $\varepsilon$

Once the hidden-sector evolution is solved, the portal coupling $\varepsilon$ can be determined from four independent constraints using the background solver:

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

**Relic abundance.** `find_epsilon_DMRelic` runs the Boltzmann solver once, then root-finds on `log10(eps)` (Brent's method over $[-20, -6]$) for the value that reproduces $\Omega_c h^2 = 0.12$ via the background evolution. A unique solution exists only when the frozen yield over-produces ($m_X \geq m_X^\star$); otherwise no amount of entropy dilution can recover the observed abundance. `find_mXstar` locates this critical mass.

**BBN.** Uses 2D log-log interpolation over the Kawasaki, Kohri, Moroi & Takaesu (2017) upper bounds on $m_Y \cdot Y_Y$ as a function of lifetime, combining `uu` and `bb` injection channels weighted by the model's branching ratios.

**Direct detection.** Root-finds the $\varepsilon$ for which `model.sigma_SI(eps)` equals the experimental limit at the given $m_X$. Limit curves (LZ and Xe neutrino floor) are loaded as log-log interpolators from tabulated data. Returns `np.inf` for models with no tree-level SI cross section (e.g. `LiLjPortal`), and `np.nan` if $m_X$ is outside the data range.

**Thermalization floor.** Estimates the minimum $\varepsilon$ such that $Y$ decay/inverse-decay keeps the dark and SM sectors in thermal contact at freeze-out, requiring $\Gamma(Y\to\text{SM}) \geq H(T_{\rm FO})$. Since $\Gamma \sim \varepsilon^2$, this gives $\varepsilon = \sqrt{H_{\rm FO}/\Gamma(\varepsilon{=}1)}$ analytically. Approximate lower bound; may underestimate the true thermalization coupling for very light mediators ($m_Y \lesssim 0.1\,m_X$).

## Structure

```
hidden_sector_DM/         # Python package
  HiddenSectorDM.py       # Public API (re-exports)
  model.py                # HiddenSectorModel, VectorPortal, BLPortal, LiLjPortal, BaryonPortal, HiggsPortal
  solver.py               # BoltzmannSolver (3phase, 2phase, hybrid, Yeq), solve_background, find_epsilon_*, find_mXstar
  cosmology.py            # Cosmology (g_star tables, thermodynamics, Bessel moments)
  constants.py            # SM parameters, fermion tables
  bbn.py                  # BBN hadronic-injection constraints (Kawasaki et al. 2017)
  direct_detection.py     # Direct-detection limit interpolators (LZ, neutrino floors)
  phi_decay/              # Mediator decay widths (2HDECAY + scalar_portal)
  data/
    bbn_constraints_decay/  # BBN bound tables (uu, bb channels, reference masses)
    direct_detection_data/  # LZ SI 90% CL, neutrino floors (Xe SI, Ar SI, Xe SD)
notebooks/
  vP_example.ipynb        # VectorPortal walkthrough
  bP_example.ipynb        # BLPortal / BaryonPortal walkthrough
  hP_example.ipynb        # HiggsPortal walkthrough
  rates_test.ipynb        # Cross-check of hidden-sector rates
python/
  run_relic.py              # Single-point relic-abundance runner
  get_vPDM_bounds_all.py    # VectorPortal parameter scans
  get_bLPDM_bounds_all.py   # BLPortal parameter scans
  get_BPDM_bounds_all.py    # BaryonPortal parameter scans
  get_hPDM_bounds_all.py    # HiggsPortal parameter scans
  compare_bbn.py            # BBN-constraint comparison utility
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

## AI Disclosure

The preparation of this public repository was aided by Claude (Anthropic). Specifically, AI assistance was used to organize the directory structure for public release, add descriptive comments and markdown sections to the notebooks, and produce an initial draft of this README. All AI-generated content was reviewed and verified by the authors.

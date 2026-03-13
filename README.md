# Secluded Dark Matter

Toolkit for computing the relic abundance of secluded dark matter, where dark matter annihilates into unstable mediators that subsequently decay to the Standard Model through a portal interaction.

The code solves the coupled Boltzmann equations for the dark matter yield and dark sector chemical potential using a three-phase QSSA (quasi-static steady-state approximation) method:
1. **Phase 1 (equilibrium)** -- both X and Y track their equilibrium abundances
2. **Phase 1.5 (QSSA)** -- X freezes out while Y is maintained in chemical equilibrium by cannibal processes
3. **Phase 2 (full)** -- cannibal processes decouple, full coupled ODE system is solved

## Models

Five portal scenarios are implemented:

| Model | DM | Mediator | Portal coupling | Reference |
|-------|-----|----------|----------------|-----------|
| `VectorPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with hypercharge | [arXiv:1602.08490](https://arxiv.org/abs/1602.08490) |
| `BLPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with U(1)_{B-L} | [arXiv:1912.08821](https://arxiv.org/abs/1912.08821) |
| `LiLjPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with U(1)_{L_i - L_j} | [arXiv:1912.08821](https://arxiv.org/abs/1912.08821) |
| `BaryonPortal` | Dirac fermion | Dark photon (spin-1) | Kinetic mixing with U(1)_B | [arXiv:1912.08821](https://arxiv.org/abs/1912.08821) |
| `HiggsPortal` | Majorana fermion | Dark scalar (spin-0) | Scalar mixing with SM Higgs | [arXiv:1609.02555](https://arxiv.org/abs/1609.02555) |

All models implement a common `DarkSectorModel` interface providing:
- 2->2 annihilation cross sections (XX -> YY)
- 3->2 cannibal cross sections (YYY -> YY, YYX -> YX, YXX -> XX)
- Mediator decay width to SM final states
- Spin-independent DM-nucleon cross section (sigma_SI)

## Usage

```python
from source.SecludedDM import Cosmology, VectorPortal, BoltzmannSolver, find_mXstar

# Set up cosmology and model
cosmo = Cosmology(gstar_choice="standard", gstarpath="./notebooks/")
model = VectorPortal(mX=50.0, mY=30.0, alphaX=1e-2)

# Solve Boltzmann equations
solver = BoltzmannSolver(cosmo, model)
sol = solver.solve_boltzmann_chempot_QSSA(return_bg_ICs=True, verbose=True)
print(f"YX_relic = {sol['YX_relic']:.4e}")

# Find critical mass (minimum epsilon)
factory = lambda mX, mY: VectorPortal(mX=mX, mY=mY, alphaX=1e-2)
res = find_mXstar(cosmo, factory, r=1.5, verbose=True)
print(f"mXstar = {res['mXstar']:.4e} GeV")
```

### Finding portal coupling constraints

Once the secluded-sector evolution is solved, the portal coupling `eps` can be determined from three independent constraints:

```python
# 1. Relic abundance: eps that gives Omega_DM h^2 = 0.12
eps_relic = solver.find_epsilon_DMRelic(sol['bg_ICs'])

# 2. BBN hadronic injection: minimum eps satisfying Kawasaki et al. (2017)
eps_bbn = solver.find_epsilon_BBN_hadronic(sol['YY_final'])

# 3. Direct detection: eps that saturates a SI cross-section limit
eps_LZ = solver.find_epsilon_DD(constraint='LZ')          # LZ 90% CL
eps_nu = solver.find_epsilon_DD(constraint='nufloor_Xe')   # Xe neutrino floor
```

The BBN constraint uses 2D log-log interpolation over the Kawasaki, Kohri, Moroi & Takaesu (2017) upper bounds on `m_Y * Y_Y` as a function of lifetime, combining `uu` and `bb` injection channels weighted by the model's branching ratios.

The direct detection constraint root-finds the `eps` for which `model.sigma_SI(eps)` equals the experimental limit at the given `mX`. Limit curves (LZ and Xe neutrino floor) are loaded as log-log interpolators from tabulated data. Returns `np.inf` for models with no tree-level SI cross section (e.g. `LiLjPortal`), and `np.nan` if `mX` is outside the data range.

For the `BLPortal`, `BaryonPortal`, and `LiLjPortal`, the gauge couplings (`g_BL`, `g_B`, `g_LiLj`) default to 1, so `eps` directly represents the effective portal coupling product `eps * g`.

## Structure

```
source/
  SecludedDM.py       # Public API (re-exports)
  model.py            # DarkSectorModel, VectorPortal, BLPortal, LiLjPortal, BaryonPortal, HiggsPortal
  solver.py           # BoltzmannSolver, find_mXstar, find_epsilon_*
  cosmology.py        # Cosmology (g_star tables, thermodynamics)
  constants.py        # SM parameters, fermion tables
  bbn.py              # BBN hadronic-injection constraints (Kawasaki et al. 2017)
  direct_detection.py # Direct-detection limit interpolators (LZ, neutrino floor)
  phi_decay/          # Mediator decay widths (HDECAY + scalar_portal)
  data/
    bbn_constraints_decay/   # BBN bound tables (uu, bb channels, 6 reference masses)
    direct_detection_data/   # LZ SI 90% CL, neutrino floors (Xe SI, Ar SI, Xe SD)
notebooks/
  example1.ipynb      # VectorPortal examples
  example2.ipynb      # Mass scans and parameter space
  example3.ipynb      # HiggsPortal examples
python/
  get_vPDM_bounds_all.py   # VectorPortal parameter scans
  get_bLPDM_bounds_all.py  # BLPortal parameter scans
  get_BPDM_bounds_all.py   # BaryonPortal parameter scans
  get_hPDM_bounds_all.py   # HiggsPortal parameter scans
output/
  vP/   # VectorPortal results
  bLP/  # BLPortal results
  BLP/  # BaryonPortal results
  hP/   # HiggsPortal results
```

## Dependencies

- Python >= 3.9
- NumPy, SciPy, tqdm
- HDECAY (bundled in `source/phi_decay/` for HiggsPortal above 5 GeV)
- scalar_portal (for HiggsPortal below 5 GeV)

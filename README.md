# Secluded Dark Matter

Toolkit for computing the relic abundance of secluded dark matter, where dark matter annihilates into unstable mediators that subsequently decay to the Standard Model through a portal interaction.

The code solves the coupled Boltzmann equations for the dark matter yield and dark sector chemical potential using a three-phase QSSA (quasi-static steady-state approximation) method:
1. **Phase 1 (equilibrium)** -- both X and Y track their equilibrium abundances
2. **Phase 1.5 (QSSA)** -- X freezes out while Y is maintained in chemical equilibrium by cannibal processes
3. **Phase 2 (full)** -- cannibal processes decouple, full coupled ODE system is solved

## Models

Three portal scenarios are implemented:

| Model | Mediator | Portal coupling | Reference |
|-------|----------|----------------|-----------|
| `VectorPortal` | Dark photon (spin-1) | Kinetic mixing with hypercharge | -- |
| `BLPortal` | Dark photon (spin-1) | Kinetic mixing with U(1)_{B-L} | [arXiv:1912.08821](https://arxiv.org/abs/1912.08821) |
| `HiggsPortal` | Dark scalar (spin-0) | Scalar mixing with SM Higgs | -- |

All models implement a common `DarkSectorModel` interface providing:
- 2->2 annihilation cross sections (XX -> YY)
- 3->2 cannibal cross sections (YYY -> YY, YYX -> YX, YXX -> XX)
- Mediator decay width to SM final states

## Usage

```python
from source.SecludedDM import Cosmology, VectorPortal, BoltzmannSolver, find_mXstar

# Set up cosmology and model
cosmo = Cosmology(gstar_choice="standard", gstarpath="./notebooks/")
model = VectorPortal(mX=50.0, mY=30.0, alphaX=1e-2,
                     Delta1=1e-2, Delta2=1e-2, Delta3=1e-2)

# Solve Boltzmann equations
solver = BoltzmannSolver(cosmo, model)
sol = solver.solve_boltzmann_chempot_QSSA(return_bg_ICs=True, verbose=True)
print(f"YX_relic = {sol['YX_relic']:.4e}")

# Find critical mass (minimum epsilon)
factory = lambda mX, mY: VectorPortal(mX=mX, mY=mY, alphaX=1e-2,
                                       Delta1=1e-2, Delta2=1e-2, Delta3=1e-2)
res = find_mXstar(cosmo, factory, r=1.5, verbose=True)
print(f"mXstar = {res['mXstar']:.4e} GeV")
```

## Structure

```
source/
  SecludedDM.py    # Public API (re-exports)
  model.py         # DarkSectorModel, VectorPortal, BLPortal, HiggsPortal
  solver.py        # BoltzmannSolver, find_mXstar
  cosmology.py     # Cosmology (g_star tables, thermodynamics)
  constants.py     # SM parameters, fermion tables
  phi_decay/       # Mediator decay widths (HDECAY + scalar_portal)
notebooks/
  example1.ipynb   # VectorPortal examples
  example2.ipynb   # Mass scans and parameter space
  example3.ipynb   # HiggsPortal examples
```

## Dependencies

- Python >= 3.9
- NumPy, SciPy, tqdm
- HDECAY (bundled in `source/phi_decay/` for HiggsPortal above 5 GeV)
- scalar_portal (for HiggsPortal below 5 GeV)

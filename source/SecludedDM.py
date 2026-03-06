"""
SecludedDM.py — backward-compatible convenience wrapper.

Re-exports the public API so that existing scripts can do::

    from source.SecludedDM import Cosmology, VectorPortal, BoltzmannSolver, ...

without knowing the internal module layout.
"""

# --- Physical constants ---
from .constants import (
    Me, Mmu, Mtau, Mu, Mc, Mt, Md, Ms, Mb,
    Mpip, MKp, MZ, MW,
    Mpl,
    sW2, sW, cW, tW, s2W,
    alphaS, gW, gY, WZ, alphaEM, eC,
    s0_cosmo, rhoc, Och2, mY_relic,
    SM_MYQ,
)

# --- Core classes ---
from .cosmology import Cosmology
from .model import DarkSectorModel, VectorPortal, BLPortal, HiggsPortal
from .solver import BoltzmannSolver, find_mXstar

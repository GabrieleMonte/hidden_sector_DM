"""
DM_Relic.py — backward-compatibility shim.

The monolithic module has been refactored into:
    constants.py, cosmology.py, model.py, solver.py

This file re-exports the public API so that old imports still work::

    from source.DM_Relic import Cosmology, VectorPortal
"""

from .SecludedDM import *          # noqa: F401,F403

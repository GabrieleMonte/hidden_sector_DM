"""
direct_detection.py
===================
Direct-detection limit interpolators for spin-independent DM-nucleon
cross sections.

Provides lazy-loaded, log-log interpolators for:
  - LZ 90% CL upper limit  (SI, Xe target)
  - XENONnT 2024 combined   (SI, Xe target)
  - Xenon SI neutrino floor

Data files live in ``source/direct_detection_data/``.
"""

import os
import numpy as np
from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
#  Data directory
# ---------------------------------------------------------------------------
_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'direct_detection_data')

# ---------------------------------------------------------------------------
#  Lazy cache:  filename -> (interp1d, mX_min, mX_max)
# ---------------------------------------------------------------------------
_cache = {}


def _load_limit(filename):
    """Load a 2-column (mX [GeV], sigma [cm^2]) file and build a
    log10-log10 interpolator.  Cached after first call."""
    if filename in _cache:
        return _cache[filename]

    fpath = os.path.join(_DATA_DIR, filename)
    data = np.loadtxt(fpath)
    mX = data[:, 0]
    sigma = data[:, 1]

    idx = np.argsort(mX)
    mX, sigma = mX[idx], sigma[idx]

    interp = interp1d(
        np.log10(mX), np.log10(sigma),
        kind='linear', bounds_error=True,
    )
    mX_min, mX_max = mX[0], mX[-1]
    _cache[filename] = (interp, mX_min, mX_max)
    return interp, mX_min, mX_max


# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------

def sigma_SI_LZ(mX):
    """LZ 90% CL upper limit on SI DM-nucleon cross section [cm^2].

    Parameters
    ----------
    mX : float
        DM mass [GeV].  Valid range: 9–10000 GeV.

    Returns
    -------
    sigma : float
        Cross-section limit [cm^2].

    Raises
    ------
    ValueError
        If mX is outside the data range.
    """
    interp, mX_min, mX_max = _load_limit('LZ_SI_90CL.txt')
    if mX < mX_min or mX > mX_max:
        raise ValueError(
            f"mX = {mX} GeV is outside the LZ data range "
            f"[{mX_min}, {mX_max}] GeV")
    return 10.0 ** float(interp(np.log10(mX)))


def sigma_SI_XENONnT(mX):
    """XENONnT 2024 combined 90% CL upper limit on SI DM-nucleon cross section [cm^2].

    Parameters
    ----------
    mX : float
        DM mass [GeV].

    Returns
    -------
    sigma : float
        Cross-section limit [cm^2].

    Raises
    ------
    ValueError
        If mX is outside the data range.
    """
    interp, mX_min, mX_max = _load_limit('XENONnT_2024_Combined.txt')
    if mX < mX_min or mX > mX_max:
        raise ValueError(
            f"mX = {mX} GeV is outside the XENONnT data range "
            f"[{mX_min}, {mX_max}] GeV")
    return 10.0 ** float(interp(np.log10(mX)))


def sigma_SI_nufloor_Xe(mX):
    """Xenon SI neutrino floor [cm^2].

    Parameters
    ----------
    mX : float
        DM mass [GeV].  Valid range: 0.1–10000 GeV.

    Returns
    -------
    sigma : float
        Neutrino-floor cross section [cm^2].

    Raises
    ------
    ValueError
        If mX is outside the data range.
    """
    interp, mX_min, mX_max = _load_limit(
        os.path.join('NuFloors', 'Xe_SI.txt'))
    if mX < mX_min or mX > mX_max:
        raise ValueError(
            f"mX = {mX} GeV is outside the Xe neutrino floor data range "
            f"[{mX_min}, {mX_max}] GeV")
    return 10.0 ** float(interp(np.log10(mX)))

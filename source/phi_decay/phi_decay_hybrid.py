"""
phi_decay_hybrid.py
===================
Decay width of a scalar mediator phi mixing with the SM Higgs,
valid from m_phi = 0.1 GeV to 1 TeV.

Physics
-------
All SM-like partial widths scale as sin^2(eps) relative to the SM Higgs.
The module returns normalised widths Gamma/sin^2(eps) in GeV; multiply by
sin^2(eps) to get the physical width.

The only exception is phi->hh, whose coupling is model-dependent
(see width_to_hh).  By default we use the singlet-extension result:
    lambda_{phi_hh} ~ sin(eps) * (2 m_phi^2 + m_h^2) / (2v)

Method
------
Three mass regimes are stitched together:

    m_phi <  2   GeV : Winkler dispersive analysis         [scalar_portal]
    2 <= m_phi < 5 GeV : perturbative QCD Nf=4             [scalar_portal]
    m_phi >=  5   GeV : NNLO QCD + NLO EW + Prophecy4f     [HDECAY]

On top of the SM-like width from any of these methods, the
model-dependent phi->hh channel is added for m_phi > 2 m_h.

References
----------
[1] Djouadi, Kalinowski, Muhlleitner, Spira, CPC 238 (2019) 214
    (HDECAY)
[2] Winkler, 1809.01876 (dispersive scalar widths below 2 GeV)
[3] Boiarska et al., 1904.10447 (scalar_portal package)
[4] Berlin, Hooper, Krnjaic, 1609.02555 (singlet-extension coupling)
"""

import os
import sys
import warnings
import numpy as np

# ---------------------------------------------------------------------------
#  Physical constants
# ---------------------------------------------------------------------------
_V_EW  = 246.0        # Higgs VEV [GeV]
_M_H   = 125.10       # SM Higgs mass [GeV]
_HBAR  = 6.582119514e-25  # hbar [GeV s]

# Regime boundaries
_M_DISPERSIVE_MAX = 2.0   # switch dispersive -> pQCD Nf=4
_M_SCALAR_PORTAL_MAX = 5.0  # switch scalar_portal -> HDECAY

# ---------------------------------------------------------------------------
#  (1)  Low-mass regime:  scalar_portal  (m_phi < 5 GeV)
# ---------------------------------------------------------------------------
_sp_light = None       # scalar_portal Model for dispersive (< 2 GeV)
_sp_heavy = None       # scalar_portal Model for pQCD Nf=4  (2-5 GeV)
_sp_ready = None       # None = not tried, True/False = outcome


def _init_scalar_portal():
    """Lazy initialisation of scalar_portal models."""
    global _sp_light, _sp_heavy, _sp_ready
    if _sp_ready is not None:
        return _sp_ready
    try:
        # Look for scalar_portal adjacent to phi_decay package or in PYTHONPATH
        here = os.path.dirname(os.path.abspath(__file__))
        for candidate in [os.path.join(here, 'scalar_portal'),
                          os.path.join(here, '..', 'scalar_portal'),
                          os.path.join(here, '..', '..', 'scalar_portal')]:
            parent = os.path.dirname(os.path.abspath(candidate))
            if os.path.isdir(candidate) and parent not in sys.path:
                sys.path.insert(0, parent)

        from scalar_portal import Model

        _sp_light = Model()
        _sp_light.decay.enable('LightScalar')

        _sp_heavy = Model()
        _sp_heavy.decay.enable('HeavyScalar')

        _sp_ready = True
    except ImportError as exc:
        warnings.warn(f"scalar_portal unavailable ({exc}); "
                      f"m_phi < {_M_SCALAR_PORTAL_MAX} GeV widths will be zero.")
        _sp_ready = False
    return _sp_ready


def _scalar_portal_width(mS):
    """Normalised total width Gamma/sin^2(eps) [GeV] via scalar_portal (m < 5 GeV)."""
    if not _init_scalar_portal():
        return 0.0
    model = _sp_light if mS < _M_DISPERSIVE_MAX else _sp_heavy
    return model.compute_branching_ratios(mS, theta=1.0).total_width


# ---------------------------------------------------------------------------
#  (2)  High-mass regime:  HDECAY  (m_phi >= 5 GeV)
# ---------------------------------------------------------------------------
from .hdecay_interface import (hdecay_branching_ratios as _hdecay_br,
                               hdecay_total_width      as _hdecay_width)

# Internal cache: mass (rounded to 0.01 GeV) -> full result dict
_hdecay_cache = {}


def _hdecay_result(mS):
    """Return cached HDECAY result dict for a single mass."""
    key = round(mS, 2)
    if key not in _hdecay_cache:
        _hdecay_cache[key] = _hdecay_br(key)
    return _hdecay_cache[key]


def _hdecay_total_width(mS):
    """SM total width Gamma_SM(m_H = m_phi) [GeV] from HDECAY."""
    return _hdecay_result(mS)['total_width']


# ---------------------------------------------------------------------------
#  (3)  Model-dependent channel:  phi -> hh
# ---------------------------------------------------------------------------

def width_to_hh(mS, coupling='singlet'):
    """
    Normalised partial width Gamma(phi->hh) / sin^2(eps)  [GeV].

    Parameters
    ----------
    mS : float
        Scalar mass [GeV].
    coupling : str
        'singlet'  - default singlet-extension coupling.
        'zero'     - turn off hh channel entirely.
    """
    if coupling == 'zero' or mS <= 2.0 * _M_H:
        return 0.0

    beta = np.sqrt(1.0 - (2.0 * _M_H / mS) ** 2)
    lam2 = (2.0 * mS**2 + _M_H**2) ** 2 / (4.0 * _V_EW**2)  # |lambda_{phi_hh}/sin(eps)|^2
    return lam2 / (32.0 * np.pi * mS) * beta


# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------

def phi_total_width_normalised(mS, hh_coupling='singlet'):
    """
    Normalised total decay width  Gamma_phi / sin^2(eps)   [GeV].

    Multiply by sin^2(eps) to obtain the physical width.

    Parameters
    ----------
    mS : float
        Scalar mediator mass [GeV].  Valid: 0.01 - 1000.
    hh_coupling : str
        'singlet' (default) or 'zero'.
    """
    mS = float(mS)
    if mS < 0.01:
        return 0.0

    # SM-like width (all channels except hh)
    if mS < _M_SCALAR_PORTAL_MAX:
        w_sm = _scalar_portal_width(mS)
    else:
        w_sm = _hdecay_total_width(mS)

    # Model-dependent hh on top
    w_hh = width_to_hh(mS, coupling=hh_coupling)

    return w_sm + w_hh


# American-spelling alias
phi_total_width_normalized = phi_total_width_normalised


def phi_total_width(mS, sin_epsilon, hh_coupling='singlet'):
    """Physical total decay width  Gamma_phi  [GeV]."""
    return sin_epsilon**2 * phi_total_width_normalised(mS, hh_coupling)


def phi_lifetime(mS, sin_epsilon, hh_coupling='singlet'):
    """Proper lifetime tau_phi  [seconds]."""
    gamma = phi_total_width(mS, sin_epsilon, hh_coupling)
    return _HBAR / gamma if gamma > 0 else np.inf


def phi_partial_widths(mS, hh_coupling='singlet'):
    """
    Dictionary of normalised partial widths  Gamma_i / sin^2(eps)  [GeV].

    For m_phi >= 5 GeV the SM channels come from HDECAY; for lower
    masses the breakdown is not available (only the total from
    scalar_portal is returned under the key 'total_sp').
    """
    mS = float(mS)
    w_hh = width_to_hh(mS, coupling=hh_coupling)

    if mS < _M_SCALAR_PORTAL_MAX:
        w_sp = _scalar_portal_width(mS)
        return {'total_sp': w_sp, 'hh': w_hh, 'total': w_sp + w_hh}

    r = _hdecay_result(mS)
    w_tot_sm = r['total_width']   # HDECAY total SM width [GeV]

    # Convert BRs to partial widths
    out = {}
    for ch in ('bb', 'tautau', 'mumu', 'ss', 'cc', 'tt',
               'gg', 'gamgam', 'Zgam', 'WW', 'ZZ'):
        out[ch] = r[ch] * w_tot_sm

    out['hh'] = w_hh
    out['total'] = w_tot_sm + w_hh
    return out


def phi_branching_ratios(mS, hh_coupling='singlet'):
    """
    Dictionary of branching ratios  (including phi->hh).
    """
    pw = phi_partial_widths(mS, hh_coupling)
    total = pw['total']
    if total <= 0:
        return {k: 0.0 for k in pw}
    return {k: v / total for k, v in pw.items()}

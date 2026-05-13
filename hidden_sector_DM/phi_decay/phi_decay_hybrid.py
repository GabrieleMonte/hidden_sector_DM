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
from ..constants import vH as _V_EW, Mh as _M_H
_HBAR  = 6.582119514e-25  # hbar [GeV s]

# Regime boundaries
_M_DISPERSIVE_MAX = 2.0   # switch dispersive -> pQCD Nf=4
_M_SCALAR_PORTAL_MAX = 5.0  # scalar_portal returns NaN above this

# scalar_portal <-> HDECAY blend window (smooths the ~2x step at mY=5).
_M_BLEND_LO = 4.0
_M_BLEND_HI = 5.0


def _smoothstep(t):
    if t <= 0.0: return 0.0
    if t >= 1.0: return 1.0
    return t * t * (3.0 - 2.0 * t)

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


# Channel name mapping: scalar_portal verbose -> standardised
_SP_CHANNEL_MAP = {
    "S -> e+ e-":     "ee",
    "S -> mu+ mu-":   "mumu",
    "S -> tau+ tau-":  "tautau",
    "S -> s sbar":    "ss",
    "S -> c cbar":    "cc",
    "S -> b bbar":    "bb",
    "S -> t tbar":    "tt",
    "S -> g g":       "gg",
    "S -> pi0 pi0":   "pi0pi0",
    "S -> pi+ pi-":   "pipi",
    "S -> K0 Kbar0":  "K0K0",
    "S -> K+ K-":     "KK",
    "S -> mesons...":  "mesons",
}


def _scalar_portal_partial_widths(mS):
    """Per-channel normalised widths [GeV] from scalar_portal (m < 5 GeV)."""
    if not _init_scalar_portal():
        return {}
    model = _sp_light if mS < _M_DISPERSIVE_MAX else _sp_heavy
    res = model.compute_branching_ratios(mS, theta=1.0)
    out = {}
    for sp_name, w in res.decay.widths.items():
        key = _SP_CHANNEL_MAP.get(sp_name, sp_name)
        if float(w) > 0:
            out[key] = float(w)
    return out


# ---------------------------------------------------------------------------
#  (2)  High-mass regime:  HDECAY  (5 GeV <= m_phi <= 500 GeV)
# ---------------------------------------------------------------------------
from .hdecay_interface import (hdecay_branching_ratios as _hdecay_br,
                               hdecay_total_width      as _hdecay_width)

_M_HDECAY_MAX = 500.0  # Above ~500 GeV, HDECAY's un-resummed EW Sudakov
                       # logs inflate WW/ZZ; tree-level is more reliable.

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
#  (2b)  Very high mass: tree-level with exact kinematics (m_phi > 500 GeV)
#
#  HDECAY's NLO EW corrections contain un-resummed Sudakov logarithms
#  ~log^n(m_H/M_W) that become uncontrolled above ~1 TeV (producing
#  negative partial widths by 1.5 TeV).  For m_phi > 500 GeV the tree-level
#  formulas with full mass dependence are the appropriate baseline.
# ---------------------------------------------------------------------------
_GF = 1.1663787e-5   # Fermi constant [GeV^-2]
_MW = 80.385
_MZ = 91.1876
_MT = 173.2
_MB = 4.18
_MC = 1.28
_MTAU = 1.777
_MMU = 0.10566
_ALPHA_S_MZ = 0.118


def _alpha_s_1loop(mu, nf=6):
    """One-loop running alpha_s."""
    b0 = (33 - 2*nf) / (12*np.pi)
    return _ALPHA_S_MZ / (1 + _ALPHA_S_MZ * b0 * np.log(mu / 91.1876))


def _tree_level_width(mS):
    """
    SM-like Gamma/sin^2(eps) [GeV] from tree-level formulas with exact
    kinematics and LO QCD corrections.  Used for m_phi > 500 GeV.

    Channels: ff-bar (all flavours), WW, ZZ, gg.
    Does NOT include hh (added separately).
    """
    w = 0.0
    pi = np.pi

    # ff-bar channels: Gamma = n_c G_F m_f^2 m_phi / (4 sqrt(2) pi) * beta^3 * (1 + Delta_QCD)
    for nc, mf, is_q in [(3, _MT, True), (3, _MB, True), (3, _MC, True),
                          (1, _MTAU, False), (1, _MMU, False)]:
        if mS > 2*mf:
            beta = np.sqrt(1 - (2*mf/mS)**2)
            gf = nc * _GF * mf**2 * mS / (4*np.sqrt(2)*pi) * beta**3
            if is_q:
                aS = _alpha_s_1loop(max(mS, 2*mf))
                gf *= (1 + 5.67 * aS/pi)
            w += gf

    # WW: Gamma = G_F m^3/(8 sqrt(2) pi) * sqrt(1-x)(1-x+3x^2/4)
    xW = (2*_MW/mS)**2
    if xW < 1:
        w += _GF*mS**3/(8*np.sqrt(2)*pi) * np.sqrt(1-xW) * (1-xW+0.75*xW**2)

    # ZZ: Gamma = G_F m^3/(16 sqrt(2) pi) * sqrt(1-x)(1-x+3x^2/4)
    xZ = (2*_MZ/mS)**2
    if xZ < 1:
        w += _GF*mS**3/(16*np.sqrt(2)*pi) * np.sqrt(1-xZ) * (1-xZ+0.75*xZ**2)

    # gg (top-loop, high-mass limit): Gamma = G_F alpha_s^2 m^3/(36 sqrt(2) pi^3)
    aS = _alpha_s_1loop(mS)
    w += _GF * aS**2 * mS**3 / (36*np.sqrt(2)*pi**3)

    return w


def _tree_level_partial_widths(mS):
    """Partial width dict from tree-level for m_phi > 500 GeV."""
    out = {}
    pi = np.pi

    def _gff(nc, mf, is_q):
        if mS <= 2*mf: return 0.0
        beta = np.sqrt(1 - (2*mf/mS)**2)
        g = nc * _GF * mf**2 * mS / (4*np.sqrt(2)*pi) * beta**3
        if is_q:
            aS = _alpha_s_1loop(max(mS, 2*mf))
            g *= (1 + 5.67 * aS/pi)
        return g

    out['tt']  = _gff(3, _MT, True)
    out['bb']  = _gff(3, _MB, True)
    out['cc']  = _gff(3, _MC, True)
    out['ss']  = 0.0  # negligible at TeV scale
    out['tautau'] = _gff(1, _MTAU, False)
    out['mumu']   = _gff(1, _MMU, False)

    xW = (2*_MW/mS)**2
    out['WW'] = _GF*mS**3/(8*np.sqrt(2)*pi)*np.sqrt(1-xW)*(1-xW+0.75*xW**2) if xW<1 else 0.0
    xZ = (2*_MZ/mS)**2
    out['ZZ'] = _GF*mS**3/(16*np.sqrt(2)*pi)*np.sqrt(1-xZ)*(1-xZ+0.75*xZ**2) if xZ<1 else 0.0

    aS = _alpha_s_1loop(mS)
    out['gg'] = _GF * aS**2 * mS**3 / (36*np.sqrt(2)*pi**3)
    out['gamgam'] = 0.0
    out['Zgam']   = 0.0

    return out


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
        Scalar mediator mass [GeV].  Valid: 0.01 - 10^5.
    hh_coupling : str
        'singlet' (default) or 'zero'.
    """
    mS = float(mS)
    if mS < 0.01:
        return 0.0

    # SM-like width (all channels except hh)
    if mS < _M_BLEND_LO:
        w_sm = _scalar_portal_width(mS)
    elif mS <= _M_BLEND_HI:
        w_sp = _scalar_portal_width(mS)
        w_hd = _hdecay_total_width(mS)
        w = _smoothstep((mS - _M_BLEND_LO) / (_M_BLEND_HI - _M_BLEND_LO))
        w_sm = (1.0 - w) * w_sp + w * w_hd
    elif mS <= _M_HDECAY_MAX:
        w_sm = _hdecay_total_width(mS)
    else:
        w_sm = _tree_level_width(mS)

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
        out = _scalar_portal_partial_widths(mS)
        if not out:
            # scalar_portal unavailable — fall back to total only
            w_sp = _scalar_portal_width(mS)
            return {'total_sp': w_sp, 'hh': w_hh, 'total': w_sp + w_hh}
        out['hh'] = w_hh
        out['total'] = sum(out.values())
        return out

    if mS > _M_HDECAY_MAX:
        out = _tree_level_partial_widths(mS)
        out['hh'] = w_hh
        out['total'] = sum(out.values())
        return out

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

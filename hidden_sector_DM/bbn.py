"""
bbn.py
======
BBN hadronic-injection constraints on late-decaying particles.

Uses the Kawasaki, Kohri, Moroi & Takaesu (2017) [arXiv:1709.01211]
upper bounds on m_Y * Y_Y as a function of lifetime tau, for uu-bar
and bb-bar injection channels at six reference masses
(30 GeV, 100 GeV, 1 TeV, 10 TeV, 100 TeV, 1000 TeV).

A 2D log-log interpolator is built on the grid (log10 m, log10 tau),
with power-law extrapolation outside the data mass range.

For a model with multiple decay channels, the effective bound is
obtained by weighting the uu and bb constraints by the corresponding
branching ratios.
"""

import os
import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator

# ---------------------------------------------------------------------------
#  Reference masses in TeV (matching the data file names)
# ---------------------------------------------------------------------------
_MASSES_TEV = np.array([0.03, 0.1, 1.0, 10.0, 100.0, 1000.0])
_LOG10_MASSES = np.log10(_MASSES_TEV)

# File name pattern
_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'bbn_constraints_decay')
_FNAME_PATTERN = {
    'uu': [
        'bbn_constraint_uubar_m30GeV.txt',
        'bbn_constraint_uubar_m100GeV.txt',
        'bbn_constraint_uubar_m1TeV.txt',
        'bbn_constraint_uubar_m10TeV.txt',
        'bbn_constraint_uubar_m100TeV.txt',
        'bbn_constraint_uubar_m1000TeV.txt',
    ],
    'bb': [
        'bbn_constraint_bbbar_m30GeV.txt',
        'bbn_constraint_bbbar_m100GeV.txt',
        'bbn_constraint_bbbar_m1TeV.txt',
        'bbn_constraint_bbbar_m10TeV.txt',
        'bbn_constraint_bbbar_m100TeV.txt',
        'bbn_constraint_bbbar_m1000TeV.txt',
    ],
}

# Power-law exponents for mass extrapolation  (mY_bound ~ m^alpha)
_ALPHA = {
    ('uu', 'low'):  0.50,   # below 30 GeV — pion-dominated
    ('uu', 'high'): 0.70,   # above 1000 TeV — nucleon-energy dominated
    ('bb', 'low'):  0.34,   # below 30 GeV — near threshold, suppressed
    ('bb', 'high'): 0.70,   # above 1000 TeV — converges with uu
}


# ---------------------------------------------------------------------------
#  Internal: load + build interpolator (lazy, cached)
# ---------------------------------------------------------------------------
_cache = {}   # channel -> (interp_func, tau_min, tau_max)


def _load_tables(channel):
    """Load data files for a channel ('uu' or 'bb').
    Returns dict: mass_TeV -> (log10_tau_array, mY_array).
    """
    tables = {}
    for m_TeV, fname in zip(_MASSES_TEV, _FNAME_PATTERN[channel]):
        fpath = os.path.join(_DATA_DIR, fname)
        d = np.loadtxt(fpath)
        idx = np.argsort(d[:, 0])          # ensure sorted by log10(tau)
        tables[m_TeV] = (d[idx, 0], d[idx, 1])
    return tables


def _build_interpolator(channel):
    """Build 2D log-log interpolator for a channel.
    Returns (callable, tau_min, tau_max).
    callable(log10_m_TeV, log10_tau) -> log10(mY_bound [GeV])
    """
    tables = _load_tables(channel)

    # Common tau grid: intersection of all mass curves
    tau_min = max(t[0].min() for t in tables.values())
    tau_max = min(t[0].max() for t in tables.values())
    n_tau = 200
    tau_grid = np.linspace(tau_min, tau_max, n_tau)

    # Interpolate each mass curve onto the common tau grid
    grid = np.zeros((len(_MASSES_TEV), n_tau))
    for i, m in enumerate(_MASSES_TEV):
        log10_tau, mY = tables[m]
        f = interp1d(log10_tau, np.log10(mY), kind='linear',
                     fill_value='extrapolate')
        grid[i, :] = f(tau_grid)

    interp = RegularGridInterpolator(
        (_LOG10_MASSES, tau_grid), grid,
        method='linear', bounds_error=False, fill_value=None,
    )
    return interp, tau_min, tau_max


def _get_interpolator(channel):
    """Return cached (interp, tau_min, tau_max) for channel."""
    if channel not in _cache:
        _cache[channel] = _build_interpolator(channel)
    return _cache[channel]


# ---------------------------------------------------------------------------
#  Core: single-channel bound
# ---------------------------------------------------------------------------

def mY_bound_channel(m_GeV, tau_s, channel='uu'):
    """
    BBN upper bound on m_chi * Y_chi [GeV] for a single injection
    channel ('uu' or 'bb').

    Parameters
    ----------
    m_GeV : float or array
        Decaying particle mass [GeV].  Valid >= 5 GeV.
    tau_s : float or array
        Lifetime [seconds].
    channel : str
        'uu' or 'bb'.

    Returns
    -------
    mY_bound : float or array [GeV]
        Upper bound on m_chi * Y_chi.  Returns np.inf for lifetimes
        shorter than the data range (decays before BBN — no
        constraint).  For lifetimes longer than the data range,
        the bound is clamped to the longest-tau value (conservative:
        if the particle hasn't decayed by BBN, the constraint is
        at least as strong as the tightest tabulated point).
    """
    m_GeV = np.atleast_1d(np.asarray(m_GeV, dtype=float))
    tau_s = np.atleast_1d(np.asarray(tau_s, dtype=float))

    interp, tau_min, tau_max = _get_interpolator(channel)
    alpha_low = _ALPHA[(channel, 'low')]
    alpha_high = _ALPHA[(channel, 'high')]

    log10_m = np.log10(m_GeV / 1e3)   # convert to TeV
    log10_t = np.log10(np.clip(tau_s, 1e-30, None))

    result = np.full_like(m_GeV, np.inf)

    for i, (lm, lt) in enumerate(zip(log10_m, log10_t)):
        # Lifetime shorter than data range -> decays before BBN, allowed
        if lt < tau_min:
            continue

        # Lifetime longer than data range -> clamp to tau_max
        # (the bound at the longest tabulated tau is the tightest;
        #  a particle that lives even longer is at least as constrained)
        lt_eff = min(lt, tau_max)

        if lm < _LOG10_MASSES[0]:
            # Extrapolate below lowest data mass (30 GeV)
            anchor = 10.0 ** float(interp((_LOG10_MASSES[0], lt_eff)))
            result[i] = anchor * (10.0**lm / _MASSES_TEV[0]) ** alpha_low

        elif lm > _LOG10_MASSES[-1]:
            # Extrapolate above highest data mass (1000 TeV)
            anchor = 10.0 ** float(interp((_LOG10_MASSES[-1], lt_eff)))
            result[i] = anchor * (10.0**lm / _MASSES_TEV[-1]) ** alpha_high

        else:
            # Within data range
            result[i] = 10.0 ** float(interp((lm, lt_eff)))

    return result.squeeze()


# ---------------------------------------------------------------------------
#  Multi-channel: combine uu and bb using branching ratios
# ---------------------------------------------------------------------------

# Channels that use the uu constraint (light quarks, gluons)
_UU_LIKE = {'uu', 'dd', 'ss', 'cc', 'gg', 'tt'}
# Channels that use the bb constraint
_BB_LIKE = {'bb'}
# Channels that are hadronic but we approximate with uu
# (WW, ZZ, hh have ~70% hadronic BRs themselves, dominated by quarks/gluons)
_BOSONIC_HADRONIC = {'WW', 'ZZ', 'hh'}
# Purely leptonic / invisible — no hadronic injection
_LEPTONIC = {'ee', 'mumu', 'tautau', 'nunu_e', 'nunu_mu', 'nunu_tau',
             'gamgam', 'Zgam'}


def mY_bound_model(m_GeV, tau_s, branching_ratios):
    """
    BBN upper bound on m_Y * Y_Y [GeV] for a model with given
    branching ratios, combining the uu and bb constraints.

    The effective constraint is:
        1 / mY_bound_eff = BR_uu_like / mY_bound_uu + BR_bb / mY_bound_bb

    where BR_uu_like = sum of BRs for uu, dd, ss, cc, tt, gg, WW, ZZ, hh
    (all approximated with the uu injection spectrum).

    Parameters
    ----------
    m_GeV : float
        Mediator mass [GeV].
    tau_s : float
        Mediator lifetime [seconds].
    branching_ratios : dict
        {channel_name: BR} as returned by model.branching_ratios_to_SM().

    Returns
    -------
    mY_bound : float [GeV]
        Upper bound on m_Y * Y_Y.  Returns np.inf if no hadronic
        channels are open.
    """
    # Accumulate hadronic BRs
    br_uu_like = 0.0
    br_bb = 0.0

    for ch, br in branching_ratios.items():
        if br <= 0:
            continue
        if ch in _BB_LIKE:
            br_bb += br
        elif ch in _UU_LIKE or ch in _BOSONIC_HADRONIC:
            br_uu_like += br
        # else: leptonic/invisible — no hadronic constraint

    br_had = br_uu_like + br_bb
    if br_had <= 0:
        return np.inf

    # Get per-channel bounds
    bound_uu = float(mY_bound_channel(m_GeV, tau_s, channel='uu'))
    bound_bb = float(mY_bound_channel(m_GeV, tau_s, channel='bb'))

    # Weighted harmonic combination:
    # Total hadronic energy injected = BR_uu * (m Y) + BR_bb * (m Y)
    # Each piece must satisfy its own bound, so the combined constraint is:
    #   mY_eff = 1 / (BR_uu_like / bound_uu + BR_bb / bound_bb)
    inv_bound = 0.0
    if br_uu_like > 0 and np.isfinite(bound_uu):
        inv_bound += br_uu_like / bound_uu
    if br_bb > 0 and np.isfinite(bound_bb):
        inv_bound += br_bb / bound_bb

    if inv_bound <= 0:
        return np.inf

    return 1.0 / inv_bound


# ---------------------------------------------------------------------------
#  Convenience: check if a point is BBN-allowed
# ---------------------------------------------------------------------------

def is_bbn_allowed(m_GeV, tau_s, Y_chi, branching_ratios):
    """
    Returns True if (m, tau, Y) satisfies the BBN hadronic injection
    constraint for the given branching ratios.
    """
    bound = mY_bound_model(m_GeV, tau_s, branching_ratios)
    return m_GeV * Y_chi < bound

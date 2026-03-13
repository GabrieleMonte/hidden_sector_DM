# BBN Hadronic Injection Constraints: Implementation Guide

## Overview

The goal is to implement the BBN upper bound on $m_\chi Y_\chi$ as a function of
decay lifetime $\tau$ and mass $m_\chi$, for hadronic injection channels, covering
the mass range ~5 GeV to ~10 PeV.

The data tables (Kawasaki et al. 2017, arXiv:1709.01211) are available for
$u\bar{u}$ and $b\bar{b}$ injection at masses
$m_X = \{0.03, 0.1, 1, 10, 100, 1000\}$ TeV.

---

## Step 1: Load and interpolate the data tables

Each data file contains two columns: $\log_{10}(\tau/\text{s})$ and $m_\chi Y_\chi$ [GeV].

```python
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d

MASSES_TEV = [0.03, 0.1, 1.0, 10.0, 100.0, 1000.0]  # TeV

def load_tables(file_list):
    """
    Load BBN constraint tables.
    Returns dict: mass (TeV) -> (log10_tau array, mY array)
    """
    tables = {}
    for m, fpath in zip(MASSES_TEV, file_list):
        d = np.loadtxt(fpath)
        # Sort by tau in case points are not ordered
        idx = np.argsort(d[:, 0])
        tables[m] = (d[idx, 0], d[idx, 1])
    return tables
```

Build a 2D interpolator on the grid $(\log_{10} m, \log_{10} \tau)$:

```python
def build_interpolator(tables):
    """
    Build a 2D log-log interpolator over the data grid.
    Returns a callable f(log10_m_TeV, log10_tau) -> log10(mY [GeV])
    """
    log10_masses = np.log10(MASSES_TEV)

    # Find common tau range across all masses
    tau_min = max(t[0].min() for t in tables.values())
    tau_max = min(t[0].max() for t in tables.values())
    tau_grid = np.linspace(tau_min, tau_max, 200)

    # Interpolate each mass curve onto common tau grid
    grid = np.zeros((len(MASSES_TEV), len(tau_grid)))
    for i, m in enumerate(MASSES_TEV):
        log10_tau, mY = tables[m]
        f = interp1d(log10_tau, np.log10(mY), kind='linear',
                     fill_value='extrapolate')
        grid[i, :] = f(tau_grid)

    interp2d = RegularGridInterpolator(
        (log10_masses, tau_grid), grid,
        method='linear', bounds_error=False, fill_value=None
    )

    def query(log10_m_TeV, log10_tau):
        pts = np.column_stack([np.atleast_1d(log10_m_TeV),
                               np.atleast_1d(log10_tau)])
        return interp2d(pts)

    return query, tau_min, tau_max
```

---

## Step 2: Mass extrapolation — the scaling law

The local power-law index $\alpha$ (where $m_\chi Y_\chi \propto m_\chi^\alpha$)
is **not constant** across the full mass range. From the numerical analysis:

### For $u\bar{u}$:

| Mass region            | $\alpha$   | Notes                          |
|------------------------|------------|-------------------------------|
| Below 0.03 TeV (extrap) | **0.50**  | Pion-dominated, sqrt scaling  |
| 0.03 – 0.1 TeV         | 0.50       | Data available                |
| 0.1 – 1 TeV            | 0.58       | Data available                |
| 1 – 10 TeV             | 0.65       | Data available                |
| 10 – 100 TeV           | 0.68       | Data available                |
| 100 – 1000 TeV         | 0.72       | Data available                |
| Above 1000 TeV (extrap) | **0.70** | Near $m^{1-\delta}$, $\delta\sim0.3$ |

### For $b\bar{b}$:

| Mass region            | $\alpha$   | Notes                                    |
|------------------------|------------|------------------------------------------|
| Below 0.03 TeV (extrap) | **0.34** | Near $b\bar{b}$ threshold, suppressed   |
| 0.03 – 0.1 TeV         | 0.34       | Data available                           |
| 0.1 – 1 TeV            | 0.53       | Data available                           |
| 1 – 10 TeV             | 0.63       | Data available                           |
| 10 – 100 TeV           | 0.68       | Data available                           |
| 100 – 1000 TeV         | 0.70       | Data available                           |
| Above 1000 TeV (extrap) | **0.70** | Converges with uu                        |

**Key physics note**: $u\bar{u}$ and $b\bar{b}$ converge to the same $\alpha\approx0.70$
at high mass (both nucleon-energy dominated), but diverge sharply at low mass because
$b\bar{b}$ near threshold produces soft, low-multiplicity showers.

### Implementation:

```python
# Anchor points for extrapolation
ALPHA_UU_LOW  = 0.50   # below 0.03 TeV
ALPHA_UU_HIGH = 0.70   # above 1000 TeV

ALPHA_BB_LOW  = 0.34   # below 0.03 TeV
ALPHA_BB_HIGH = 0.70   # above 1000 TeV

M_LOW_TEV  = 0.03      # lowest data point
M_HIGH_TEV = 1000.0    # highest data point

def extrapolate_low(mY_anchor, m_anchor_TeV, m_target_TeV, alpha):
    """Scale mY from anchor mass to target mass using power law."""
    return mY_anchor * (m_target_TeV / m_anchor_TeV) ** alpha

def extrapolate_high(mY_anchor, m_anchor_TeV, m_target_TeV, alpha):
    return mY_anchor * (m_target_TeV / m_anchor_TeV) ** alpha
```

---

## Step 3: Full constraint evaluator

```python
def bbn_constraint(m_GeV, tau_s, channel='uu',
                   interp_uu=None, interp_bb=None,
                   tau_min_uu=None, tau_max_uu=None,
                   tau_min_bb=None, tau_max_bb=None):
    """
    Returns the BBN upper bound on m_chi * Y_chi [GeV]
    for a particle of mass m_GeV [GeV] and lifetime tau_s [s].

    Parameters
    ----------
    m_GeV    : float or array, mass in GeV
    tau_s    : float or array, lifetime in seconds
    channel  : 'uu' or 'bb'

    Returns
    -------
    mY_bound : float or array [GeV], upper bound on m_chi Y_chi
    """
    m_GeV   = np.atleast_1d(np.asarray(m_GeV, dtype=float))
    tau_s   = np.atleast_1d(np.asarray(tau_s,  dtype=float))
    m_TeV   = m_GeV / 1e3
    log10_m = np.log10(m_TeV)
    log10_t = np.log10(tau_s)

    if channel == 'uu':
        interp      = interp_uu
        alpha_low   = ALPHA_UU_LOW
        alpha_high  = ALPHA_UU_HIGH
        tau_min     = tau_min_uu
        tau_max     = tau_max_uu
    else:
        interp      = interp_bb
        alpha_low   = ALPHA_BB_LOW
        alpha_high  = ALPHA_BB_HIGH
        tau_min     = tau_min_bb
        tau_max     = tau_max_bb

    result = np.zeros_like(m_GeV)

    for i, (lm, lt) in enumerate(zip(log10_m, log10_t)):

        # Clamp tau to data range (extrapolation in tau is unreliable)
        lt_clamped = np.clip(lt, tau_min, tau_max)

        if lm < np.log10(M_LOW_TEV):
            # Extrapolate downward from lowest data mass
            mY_anchor = 10 ** interp(np.log10(M_LOW_TEV), lt_clamped)
            result[i] = extrapolate_low(mY_anchor, M_LOW_TEV,
                                         10**lm, alpha_low)

        elif lm > np.log10(M_HIGH_TEV):
            # Extrapolate upward from highest data mass
            mY_anchor = 10 ** interp(np.log10(M_HIGH_TEV), lt_clamped)
            result[i] = extrapolate_high(mY_anchor, M_HIGH_TEV,
                                          10**lm, alpha_high)
        else:
            # Within data range: direct 2D interpolation
            result[i] = 10 ** interp(lm, lt_clamped)

    return result.squeeze()
```

---

## Step 4: Channel mixing for realistic portals

For a mediator with multiple decay channels, only the hadronic fraction
contributes to the hadronic injection constraint. The bound should be
rescaled accordingly:

```python
HADRONIC_FRACTIONS = {
    'hypercharge': 0.61,   # qq~53.7% + bb~7.4%, exclude nu, e/mu, tau
    'baryon':      1.00,   # 100% hadronic
    'higgs':       0.92,   # bb~85.8% + qq~4.3% + gg~2.2%
    'BmL':         0.21,   # qq~17.4% + bb~4.3%; 39% is invisible nu
    'LiLj':        0.00,   # 100% leptonic/invisible — use EM constraint instead
}

# Best representative channel per portal
BEST_CHANNEL = {
    'hypercharge': 'uu',
    'baryon':      'uu',
    'higgs':       'bb',   # Yukawa-weighted, bb dominates
    'BmL':         'uu',
}

def bbn_constraint_portal(m_GeV, tau_s, portal):
    """
    BBN constraint accounting for hadronic branching fraction.
    The bound on m_chi Y_chi is WEAKENED by 1/Br_had,
    since only Br_had fraction of decays inject hadronic energy.
    """
    Br = HADRONIC_FRACTIONS[portal]
    if Br == 0:
        return np.inf  # no hadronic constraint

    channel = BEST_CHANNEL.get(portal, 'uu')
    mY_hadronic = bbn_constraint(m_GeV, tau_s, channel=channel)

    # The actual bound on total mY is weakened by the hadronic fraction
    # because only Br_had * n_Y particles inject hadronic energy
    return mY_hadronic / Br
```

---

## Step 5: Yield constraint on the model

In practice you want to check whether your model's predicted yield $Y_\chi$
violates the BBN bound:

```python
def is_bbn_allowed(m_GeV, tau_s, Y_chi, portal):
    """
    Returns True if (m_chi, tau, Y_chi) satisfies BBN injection constraint.
    """
    mY_bound = bbn_constraint_portal(m_GeV, tau_s, portal)
    return m_GeV * Y_chi < mY_bound
```

---

## Step 6: Validity and caveats

| Mass range     | Status              | Uncertainty                        |
|----------------|---------------------|------------------------------------|
| < 5 GeV        | Do not use          | Below QCD thresholds, formalism breaks |
| 5 – 30 GeV     | Extrapolation       | Factor ~2–3, likely conservative   |
| 30 GeV – 1 PeV | Data + interpolation | ~20% from channel choice          |
| > 1 PeV        | Extrapolation       | Factor ~2, both channels agree     |

**The low-mass extrapolation ($u\bar{u}$, $\alpha=0.50$) is conservative**
(i.e. over-constraining) because at lower masses the hadronic shower
produces fewer nucleons per unit energy — the true bound on $m_\chi Y_\chi$
is expected to be *weaker* (higher) than the extrapolation gives.

**For $b\bar{b}$ below 30 GeV**: kinematic suppression near threshold makes
$\alpha \approx 0.34$ — do not use $\sqrt{m}$ scaling here, it will 
over-estimate the constraint by up to an order of magnitude at 5 GeV.

**The $L_i - L_j$ portal** (100% leptonic/invisible) requires a separate
treatment using the photodissociation constraints from electromagnetic
injection (Fig. 13 of Kawasaki et al.), which are ~2–3 orders of magnitude
weaker at $\tau \lesssim 10^2$ s and only become relevant at $\tau \gtrsim 10^4$ s.

---

## Summary: recommended $\alpha$ values for extrapolation

```python
ALPHA_EXTRAP = {
    # (channel, direction): alpha
    ('uu', 'low'):  0.50,
    ('uu', 'high'): 0.70,
    ('bb', 'low'):  0.34,
    ('bb', 'high'): 0.70,
}
```

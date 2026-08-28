"""
chi^2 of one model point against the measured GCE spectrum of arXiv:2112.09706.

A point is a mass -- m_DM for direct annihilation, or the (mX, mY) carried by a
portal model -- plus a cross section <sigma v>.  Its chi^2 over the 14 measured
bins is

    chi2 = (phi - <sigma v> t)^T C^-1 (phi - <sigma v> t)

with phi and C from `gce_data`, and t the point's spectrum at <sigma v> = 1:

    phi, Cinv  = measurement(region="south")    # once; reuse across points
    t          = cascade_template(model)        # or direct_template(m_DM)
    sigmav_hat, sigmav_err, chi2_min = gls_fit(t, phi, Cinv)

<sigma v> enters linearly, so `gls_fit` minimises over it in closed form, and
the three numbers it returns give chi^2 at any other value without refitting.

Templates are built at the fixed halo `halo.HALO_2112`, whose gamma has to be
the profile the measured spectrum was extracted with.  Halo uncertainty is
applied afterwards, to <sigma v> only, by `halo.sigmav_halo_range`.
"""

from __future__ import annotations

import numpy as np

from . import gce_data as gd
from .halo import HALO_2112
from .pythia_runner import DEFAULT_SEED
from .spectrum import cascade_flux, default_ebins, direct_flux


# ---- the measured side ---------------------------------------------------

def measurement(region: str = "both", model: str = "BestFit",
                n_pc: int | None = 3):
    """The measured spectrum and the inverse of its covariance, `(phi, Cinv)`.

    `phi` is E^2 dPhi/dE on the 14 bins [GeV cm^-2 s^-1 sr^-1], `Cinv` its inverse.

    Call once and reuse: the inverse is the only real cost here, and it does not
    depend on the model point.  Arguments go straight to `gce_data`.
    """
    phi = gd.load_sed(model, region)
    return phi, np.linalg.inv(gd.load_covariance(n_pc=n_pc, model=model,
                                                 region=region))


# ---- the predicted side --------------------------------------------------

def rebin_to_gce(E, dPhidE):
    """Bin-average E^2 dPhi/dE onto the 14 bins of `gce_data`.

    A bin running past the top of the energy grid is integrated over the part
    that exists and divided by the *full* bin width.  That is the correct
    average, not an approximation: no annihilation puts a photon above E = m.
    Without it the top bin (23.7-51.9 GeV) would be unusable below m = 51.9 GeV,
    and that is the bin which discriminates against heavy dark matter.

    A bin below the *start* of the grid has no flux to average and returns nan.
    """
    E = np.asarray(E, dtype=float)
    y = np.asarray(dPhidE, dtype=float) * E ** 2
    # Log-log interpolation, as a falling spectrum wants.  Empty bins go to
    # -inf and come back as zero flux, which is what they are.
    with np.errstate(divide="ignore"):
        lny, lnE = np.log(y), np.log(E)
    out = np.full(gd.N_EBINS, np.nan)
    for i, (lo, hi) in enumerate(zip(gd.EBIN_LO, gd.EBIN_HI)):
        if lo < E[0]:
            continue                      # grid does not reach down this far
        if lo >= E[-1]:
            out[i] = 0.0                  # entirely above the kinematic endpoint
            continue
        g = np.linspace(np.log(lo), np.log(min(hi, E[-1])), 33)
        out[i] = np.trapz(np.exp(np.interp(g, lnE, lny)), g) / np.log(hi / lo)
    return out


def _check_bin_coverage(t, where: str):
    """Raise unless the spectrum reached every GCE bin.
    `rebin_to_gce` returns nan for any bin starting below the bottom of the
    energy grid, since there is no flux there to average. 
    """
    if not np.isfinite(t).all():
        raise ValueError(
            f"template at {where} does not reach GCE bins "
            f"{np.flatnonzero(~np.isfinite(t)).tolist()}; lower e_min")
    return t


def direct_template(m_DM, channel: str = "bb", n_events: int = 200_000,
                    majorana: bool = True, seed=DEFAULT_SEED, n_bins: int = 180,
                    e_min=None, halo=HALO_2112, runner=None,
                    on_mismatch: str = "raise"):
    """14-bin E^2 dPhi/dE for DM DM -> `channel`, at <sigma v> = 1 cm^3/s.

    `majorana=True` sets kappa = 8, the convention of 2112.09706, so this is
    directly comparable with their contours.  `n_events` defaults to the 200k
    the direct-arm cache holds.
    """
    kw = {} if e_min is None else dict(e_min=e_min)
    E, dPhidE = direct_flux(m_DM, default_ebins(m_DM, n_bins=n_bins, **kw),
                            sigmav=1.0, channel=channel, n_events=n_events,
                            majorana=majorana, halo=halo, runner=runner,
                            seed=seed, on_mismatch=on_mismatch)
    return _check_bin_coverage(rebin_to_gce(E, dPhidE), f"m_DM = {m_DM:g} GeV")


def cascade_template(model, n_events: int = 100_000, seed=DEFAULT_SEED,
                     n_bins: int = 180, e_min=None, halo=HALO_2112,
                     runner=None, on_mismatch: str = "raise"):
    """14-bin E^2 dPhi/dE for the cascade of `model`, at <sigma v> = 1 cm^3/s.

    Depends on (mX, mY) alone -- branching ratios from mY, boost from mX/mY,
    normalisation divided out -- so no coupling reaches it and alphaX is an
    overlay on a finished region, never a fit parameter.  kappa comes from
    `model.include_antiparticlesX`: 16 for a Dirac X against 8 for the direct
    reference, so say which is which when overlaying the two.  `n_events`
    defaults to the 100k the cascade cache holds.
    """
    kw = {} if e_min is None else dict(e_min=e_min)
    E, dPhidE = cascade_flux(model, default_ebins(model.mX, n_bins=n_bins, **kw),
                             sigmav=1.0, n_events=n_events, halo=halo,
                             runner=runner, seed=seed, on_mismatch=on_mismatch)
    return _check_bin_coverage(rebin_to_gce(E, dPhidE),
                               f"mX = {model.mX:g}, mY = {model.mY:g} GeV")


# ---- the comparison ------------------------------------------------------

def gls_fit(t, phi, Cinv):
    """<sigma v> minimised out: `(sigmav_hat, sigmav_err, chi2_min)`,
    with `sigmav_hat` clamped at >= 0.

    Generalised least squares: least squares against the full covariance rather
    than its diagonal, which matters because C is strongly correlated here.

        sigmav_hat = (t C^-1 phi) / (t C^-1 t),   sigmav_err = (t C^-1 t)^-1/2

    chi^2 is exactly parabolic in <sigma v>, so for an interior fit these three
    numbers give it everywhere without refitting:

        chi2(s) = chi2_min + ((s - sigmav_hat) / sigmav_err)**2

    The unconstrained solution goes negative where the template anti-correlates
    with the residual under C^-1.  A negative annihilation cross section is not
    a hypothesis, so it is clamped: `sigmav_hat = 0` and `chi2_min` becomes the
    chi^2 with no signal at all.  Such a point is a one-sided limit, and the
    parabola above does NOT hold there -- for chi^2 at some other s on a clamped
    point, form the residual directly: `r = phi - s * t; float(r @ Cinv @ r)`.
    """
    tC = np.asarray(t) @ Cinv
    sigmav_hat = (tC @ phi) / (tC @ t)
    sigmav_err = float(1.0 / np.sqrt(tC @ t))
    r = phi - sigmav_hat * t
    chi2_min = float(r @ Cinv @ r)
    if sigmav_hat < 0.0:
        chi2_min += (sigmav_hat / sigmav_err) ** 2      # the chi^2 at sv = 0
        sigmav_hat = 0.0
    return float(sigmav_hat), sigmav_err, chi2_min

"""
Generalised-NFW halo profile, ROI-averaged J-factors, and the halo systematic.

The profile follows arXiv:2112.09706 Eq. 6,

    rho(r) = rho_0 / [ (r/r_c)^gamma * (1 + r/r_c)^(3-gamma) ]

with rho_0 fixed by rho(r_sun) = rho_local.  The J-factor is the ROI-*averaged*
line-of-sight integral of rho^2,

    Jbar = (1/dOmega) Int_ROI dOmega Int ds rho(r(s,psi))^2  [GeV^2 cm^-5 sr^-1]

which is what multiplies <sigma v> dN/dE to give an intensity.  `HALO_2112` is
the fixed profile everything here is built at, and `sigmav_halo_range` turns
uncertainty on it into a range on the inferred <sigma v>.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.special import erfinv
from scipy.interpolate import interp1d

KPC_CM = 3.0856775814913673e21  # cm per kpc

# ---- gNFW defaults of 2112.09706 -----------------------------------------

GAMMA_DEF = 1.2
RC_DEF = 20.0        # kpc
RHO_LOCAL_DEF = 0.4  # GeV/cm^3
R_SUN_DEF = 8.5      # kpc; 2112.09706 keeps 8.5 rather than 8.3, as GALPROP assumes it
R_HALO_DEF = 100.0   # kpc, outer truncation of the line-of-sight integral

# The reference halo, as a kwargs dict for `j_factor` and the flux functions.
# Every template and contour is built at exactly this profile: Fig. 18 of
# 2112.09706 was derived without varying the halo, so comparing against it is
# only meaningful at the same fixed choice.  Uncertainty is applied afterwards,
# by `sigmav_halo_range`, and only to <sigma v>.
HALO_2112 = dict(roi="both", gamma=GAMMA_DEF, r_c=RC_DEF,
                 rho_local=RHO_LOCAL_DEF, r_sun=R_SUN_DEF)

# The 40 deg x 40 deg window of 2112.09706 and its two halves; |b| <= 2 deg is
# their galactic-disk mask, excluded from every fit.  The entries differ only in
# `hemisphere`: same window, same mask.  The keys are the region names
# `gce_data`/`likelihood` use, so one string names the ROI throughout.
ROIS: dict[str, dict] = {
    "both":  dict(l_max=20.0, b_max=20.0, b_mask=2.0, hemisphere="both"),
    "south": dict(l_max=20.0, b_max=20.0, b_mask=2.0, hemisphere="south"),
    "north": dict(l_max=20.0, b_max=20.0, b_mask=2.0, hemisphere="north"),
}


def roi_window(roi):
    """Look up an ROI name in `ROIS`, raising on anything unrecognised."""
    try:
        return dict(ROIS[roi])
    except (KeyError, TypeError):
        raise ValueError(
            f"unknown ROI {roi!r}; choose one of {sorted(ROIS)}") from None


# ---- the profile and its line-of-sight integral --------------------------

def rho_gnfw(r, gamma=GAMMA_DEF, r_c=RC_DEF, rho_local=RHO_LOCAL_DEF,
             r_sun=R_SUN_DEF):
    """gNFW density [GeV/cm^3] at galactocentric radius `r` [kpc]."""
    r = np.asarray(r, dtype=float)

    def shape(x):
        return 1.0 / ((x / r_c) ** gamma * (1.0 + x / r_c) ** (3.0 - gamma))

    return rho_local * shape(r) / shape(r_sun)


def _los_integral(psi, gamma, r_c, rho_local, r_sun, r_halo):
    """Int rho^2 ds along a line of sight at angle `psi` from the GC [GeV^2/cm^5]."""
    cpsi = np.cos(psi)
    # distance at which the line of sight exits the halo sphere of radius r_halo
    s_max = r_sun * cpsi + np.sqrt(r_halo ** 2 - (r_sun * np.sin(psi)) ** 2)

    def integrand(s):
        r = np.sqrt(r_sun ** 2 + s ** 2 - 2.0 * r_sun * s * cpsi)
        return rho_gnfw(max(r, 1e-8), gamma, r_c, rho_local, r_sun) ** 2

    # split at the point of closest approach, where the integrand is peaked
    s_peak = r_sun * cpsi
    pts = [0.0] + ([s_peak] if 0.0 < s_peak < s_max else []) + [s_max]
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        tot += quad(integrand, a, b, limit=200, epsabs=0.0, epsrel=1e-8)[0]
    return tot * KPC_CM  # ds was in kpc


@lru_cache(maxsize=None)
def j_factor(roi="both", gamma=GAMMA_DEF, r_c=RC_DEF,
             rho_local=RHO_LOCAL_DEF, r_sun=R_SUN_DEF, r_halo=R_HALO_DEF,
             n_grid=201, return_solid_angle=False):
    """ROI-averaged J-factor [GeV^2 cm^-5 sr^-1].

    Jbar is the same for all three ROIs: the profile is spherical and the window
    symmetric about b = 0, so a hemisphere only halves the solid angle, which
    cancels in an average.  `return_solid_angle` also returns dOmega [sr], where
    that halving does show up.

    Memoised: the line-of-sight quad is the dominant cost of a template, and a
    scan over (mX, mY) asks for the same fixed halo at every point -- ~8x on a
    16k-point grid.  Every argument is a scalar or a string and the result is a
    float or a tuple of them, so there is nothing here to be unhashable or to
    mutate under a caller.  Pass floats, not 0-d arrays.
    """
    win = roi_window(roi)
    l_max, b_max, b_mask = win["l_max"], win["b_max"], win["b_mask"]

    # J depends on the line of sight only through psi, so tabulate it once.
    psi_min = np.deg2rad(b_mask) if b_mask > 0 else 1e-4
    psi_max = np.deg2rad(np.hypot(l_max, b_max)) * 1.001
    psi_tab = np.logspace(np.log10(psi_min * 0.99), np.log10(psi_max), 160)
    j_tab = [_los_integral(p, gamma, r_c, rho_local, r_sun, r_halo)
             for p in psi_tab]
    log_j = interp1d(np.log(psi_tab), np.log(j_tab), kind="cubic")

    # Average over the l > 0, b > 0 quadrant; symmetry supplies the rest.
    lon = np.linspace(0.0, np.deg2rad(l_max), n_grid)
    lat = np.linspace(np.deg2rad(b_mask), np.deg2rad(b_max), n_grid)
    L, B = np.meshgrid(lon, lat, indexing="ij")
    psi = np.arccos(np.clip(np.cos(B) * np.cos(L), -1.0, 1.0))
    J = np.exp(log_j(np.log(np.clip(psi, psi_tab[0], psi_tab[-1]))))

    w = np.cos(B)                                   # dOmega = cos(b) db dl
    dO = np.trapz(np.trapz(w, lat, axis=1), lon)
    jbar = np.trapz(np.trapz(J * w, lat, axis=1), lon) / dO

    if not return_solid_angle:
        return jbar
    # x2 for l -> -l, x2 again for b -> -b unless one hemisphere is used.  The
    # factor cancels in `jbar`, so it is applied only here.
    return jbar, dO * (2.0 if win["hemisphere"] != "both" else 4.0)


# ---- halo uncertainty on the inferred <sigma v> --------------------------

# Local dark-matter density [GeV/cm^3], central value and 1 sigma.  Note this
# is NOT HALO_2112's 0.4: the templates are built at the reference, so a prior
# centred elsewhere shifts the allowed <sigma v> as well as widening it.
RHO_LOCAL_PRIOR = (0.44, 0.13)


def sigmav_halo_range(cl=95, rho_local=RHO_LOCAL_PRIOR, gamma=None,
                      ref=HALO_2112):
    """Multiplicative range the halo uncertainty puts on the inferred <sigma v>.

    A fit at the reference halo gives <sigma v>; the true value then lies in
    `<sigma v> * [R_lo, R_hi]` at confidence `cl`, in percent.

    The halo reaches the flux only as Jbar = rho_local^2 j(gamma), one overall
    factor, so <sigma v> ~ 1/Jbar: it slides the region along <sigma v> and
    cannot move a mass contour at all.

    `rho_local` and `gamma` are `(mean, sigma)` pairs, or None to hold that one
    at the reference.  rho_local is log-normal, which keeps Jbar positive and
    makes the propagation exact in log space; both push Jbar monotonically, so
    the range comes from their endpoints.  gamma is None by default because the
    measured spectrum was extracted with the spatial template already fixed.
    """
    # One parameter (the Jbar factor), so cl -> n sigma is the normal quantile.
    n = np.sqrt(2.0) * erfinv(cl / 100.0)
    lo = hi = 1.0
    if rho_local is not None:
        mu, sig = rho_local
        centre = (mu / ref["rho_local"]) ** 2
        lo *= centre * np.exp(-2 * n * sig / mu)
        hi *= centre * np.exp(+2 * n * sig / mu)
    if gamma is not None:
        mu, sig = gamma
        def j(g):                       # `ref`, with gamma overridden
            return j_factor(**{**ref, "gamma": g})
        j_ref = j(ref["gamma"])
        lo *= j(mu - n * sig) / j_ref
        hi *= j(mu + n * sig) / j_ref
    return 1.0 / hi, 1.0 / lo          # <sigma v> ~ 1 / Jbar

"""Measured GCE spectrum and covariance, from arXiv:2112.09706 (Zenodo 6423495).

The excess is tabulated as phi_i = (E^2 dPhi/dE)_i, in GeV cm^-2 s^-1 sr^-1, on
the 14 energy bins of their Table III.  It is always measured in the same
40 deg x 40 deg window around the Galactic Centre, with the |b| <= 2 deg disk
masked out.  `region` chooses how much of that window to use: 'both' takes all
of it, 'north' and 'south' take its b > 0 and b < 0 halves.  The halves share
the window, the mask and the bins, so they are not separate regions of interest;
they are tabulated only for the five best-fitting emission models, not all 80.

Two errors come with each bin, and Eq. 18 of the paper combines them into the
covariance the fit uses:

    C_ij = sigma_i^2 delta_ij + Sigma_ij,mod^trunc

sigma_i is statistical and uncorrelated between bins.  Sigma_ij,mod covers the
choice of interstellar-emission model, and is both larger and strongly
correlated: its mean off-diagonal correlation is 0.82, and its leading principal
component, which carries 86% of its variance, has the same sign in all 14 bins.
A shared mode of that size is what a fit with one free normalisation is least
able to tell apart from a signal, so the full C is inverted rather than just its
diagonal.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np

_DATA = Path(__file__).resolve().parent / "data" / "zenodo_2112.09706"
_SPECTRA = _DATA / "Figures_12_and_14_GCE_Spectra"
_COV = _DATA / "Covariance_Matrix_Information"

# Table III bin edges in GeV
EBIN_LO = np.array([274.698, 357.014, 463.995, 603.034, 783.737, 1018.59,
                    1323.82, 1720.51, 2236.07, 2906.12, 3776.96, 4908.75,
                    10776.0, 23656.1]) / 1e3
EBIN_HI = np.array([357.014, 463.995, 603.034, 783.737, 1018.59, 1323.82,
                    1720.51, 2236.07, 2906.12, 3776.96, 4908.75, 10776.0,
                    23656.1, 51931.2]) / 1e3
# Bin centres as tabulated in GeV
EBIN_CTR = np.array([313.162910358, 407.004243322, 528.965751061, 687.473829541,
                     893.47989989, 1161.21704886, 1509.18340158, 1961.42016848,
                     2549.17266733, 3313.04908164, 4305.82610508, 7273.00221156,
                     15966.1266302, 35049.7899155]) / 1e3
N_EBINS = len(EBIN_CTR)

# Tolerance for recognising a file as living on these bins.
_BIN_RTOL = 1e-5

_RANKED = {"BestFit": "BestFitModel", "2ndBestFit": "2ndBestFitModel",
           "3rdBestFit": "3rdBestFitModel", "4thBestFit": "4thBestFitModel",
           "5thBestFit": "5thBestFitModel"}
_REGION_TAG = {"both": "", "north": "_North", "south": "_South"}


def _spectrum_path(model: str, region: str) -> Path:
    try:
        tag = _REGION_TAG[region]
    except KeyError:
        raise ValueError(f"unknown region {region!r}; "
                         f"choose one of {sorted(_REGION_TAG)}") from None
    token = _RANKED.get(model, f"Model{model}")
    path = _SPECTRA / f"GCE_{token}{tag}_flux_Inner40x40_masked_disk.dat"
    if not path.exists():
        extra = ("; north/south files exist only for the five ranked models"
                 if region != "both" and model not in _RANKED else "")
        raise FileNotFoundError(f"no GCE spectrum for model={model!r}, "
                                f"region={region!r} at {path}{extra}")
    return path


def _check_bins(E: np.ndarray, what: str) -> None:
    """Refuse a file that is not on the Table III bins."""
    if not np.allclose(E, EBIN_CTR, rtol=_BIN_RTOL):
        raise ValueError(f"{what} is not on the Table III bins")


def load_sed(model: str = "BestFit", region: str = "both") -> np.ndarray:
    """Measured E^2 dPhi/dE on the 14 bins [GeV cm^-2 s^-1 sr^-1].

    `model` is a rank ('BestFit' ... '5thBestFit') or a Roman numeral ('XLIX').
    The file's own error bars are dropped here; use `load_covariance`.
    """
    E, phi, _lo, _hi = np.loadtxt(_spectrum_path(model, region)).T
    _check_bins(E, f"model {model!r}")
    return phi


def _sigma_stat(model: str, region: str) -> np.ndarray:
    """Per-bin statistical errors sigma_i on E^2 dPhi/dE.
    """
    if (model, region) == ("BestFit", "both"):
        tab = np.loadtxt(_COV / "GCE_Statistical_errors.dat")
        _check_bins(tab[:, 0], "statistical-error file")
        return tab[:, 1]
    E, _phi, lo, hi = np.loadtxt(_spectrum_path(model, region)).T
    _check_bins(E, f"model {model!r}")
    return 0.5 * (hi - lo)


def load_covariance(n_pc: int | None = 3, include_stat: bool = True,
                    model: str = "BestFit", region: str = "both") -> np.ndarray:
    """Fit covariance C_ij of Eq. 18 [(GeV cm^-2 s^-1 sr^-1)^2].

    `n_pc` says how many principal components of the systematics matrix to
    keep.  Sigma_ij,mod is replaced by its `n_pc` largest ones and the rest are
    discarded; the paper keeps 3, which retains 98% of its variance.  Pass None
    to keep Sigma_ij,mod whole.

    `include_stat=True` (the default) then adds the statistical variances along
    the diagonal, which gives the full C of Eq. 18.  Set it to False to get the
    systematics on their own.

    ASSUMPTION for a hemisphere.  Zenodo ships Sigma_ij,mod for the full window
    only, and it cannot be rebuilt from the spectra (the halves exist for five
    models, not 80).  A hemisphere therefore reuses the full-window correlation
    matrix and fractional variances, rescaled to its own spectrum:

        Sigma^reg_ij = Sigma_ij * (phi^reg_i phi^reg_j) / (phi_i phi_j)

    i.e. the systematic is taken to be multiplicative in the flux -- a
    leading-order guess, not a measurement.  The statistical part stays exact.
    """
    Sigma = np.load(_COV / "cov_mat_21Dec02.npy")
    if Sigma.shape != (N_EBINS, N_EBINS):
        raise ValueError(f"expected a {N_EBINS}x{N_EBINS} covariance, "
                         f"got {Sigma.shape}")

    if (model, region) != ("BestFit", "both"):
        scale = load_sed(model, region) / load_sed("BestFit", "both")
        Sigma = Sigma * np.outer(scale, scale)

    if n_pc is not None:
        if not 0 < n_pc <= N_EBINS:
            raise ValueError(f"n_pc must be in 1..{N_EBINS}, got {n_pc}")
        # Sigma is symmetric positive-semidefinite, so the SVD of the Zenodo
        # script is an eigendecomposition; eigh gets there without the round trip.
        w, V = np.linalg.eigh(Sigma)
        keep = np.argsort(w)[::-1][:n_pc]
        Sigma = (V[:, keep] * w[keep]) @ V[:, keep].T

    if not include_stat:
        return Sigma
    return Sigma + np.diag(_sigma_stat(model, region) ** 2)

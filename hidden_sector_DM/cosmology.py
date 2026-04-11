"""
Cosmology — SM thermal background.

Loads g*(T) and g*_S(T) tables, provides Hubble rate, radiation
energy density, entropy density, and full species thermodynamics
(integral / Maxwell-Boltzmann hybrid).
"""

import os
import numpy as np
import warnings
from dataclasses import dataclass
from typing import Optional
from pathlib import Path
from scipy import special as scisp
from scipy.integrate import quad
from scipy.special import kv as Knu, zeta
from scipy.integrate._quadpack_py import IntegrationWarning

from .constants import Mpl


@dataclass
class Cosmology:
    """
    Standard-Model thermal background.

    Loads g*(T) and g*_S(T) tables, provides Hubble rate, radiation
    energy density, entropy density, and full species thermodynamics
    (integral / Maxwell-Boltzmann hybrid).

    Parameters
    ----------
    gstar_choice : str
        Which g* table to load (e.g. "standard", "HP_A", ...).
    gstarpath : str
        Base path prepended to the filename from the table map.
    gstar_dir : str
        Directory that contains the .tab files.
    """
    gstar_choice: str = "standard"
    gstarpath: str = "./"
    gstar_dir: Optional[str] = None
    def __post_init__(self):
        # Resolve repo root relative to THIS source file: main_directory/source/this_file.py
        source_dir = Path(__file__).resolve().parent          # .../source
        repo_root  = source_dir.parent                        # .../main_directory

        if self.gstar_dir is None:
            self.gstar_dir = str(repo_root / "gstar")
        self._load_gstar(self.gstar_choice, self.gstarpath)

    # ----------------------------------------------------------
    #  g* table loading and interpolation
    # ----------------------------------------------------------
    def _load_gstar(self, choice: str, gstarpath: str):
        files = {
            "standard": "std.tab",
            "HP_A":     "HP_A.tab",
            "HP_B":     "HP_B.tab",
            "HP_B2":    "HP_B2.tab",
            "HP_B3":    "HP_B3.tab",
            "HP_C":     "HP_C.tab",
        }

        # Build the final path robustly.
        filepath = Path(gstarpath) / Path(self.gstar_dir) / files[choice]
        filepath = filepath.resolve()

        if not filepath.exists():
            raise FileNotFoundError(
                f"g* table not found: {filepath}\n"
                f"cwd={Path.cwd()}\n"
                f"gstar_dir={self.gstar_dir}, gstarpath={gstarpath}, choice={choice}"
            )

        T, gS, g = np.loadtxt(
            filepath,
            comments="#", usecols=(0, 1, 2), unpack=True,
        )

        o = np.argsort(T)
        T, gS, g = T[o], gS[o], g[o]

        dg  = np.gradient(g,  T, edge_order=2)
        dgS = np.gradient(gS, T, edge_order=2)
        eps = 1e-300

        self.Tvec             = T
        self.gvec             = g
        self.gSvec            = gS
        self._dlngdlnT_vec   = (T / np.maximum(g,  eps)) * dg
        self._dlngSdlnT_vec  = (T / np.maximum(gS, eps)) * dgS

    def _interp(self, T, y):
        T = np.asarray(T, float)
        return np.interp(T, self.Tvec, y, left=y[0], right=y[-1])

    # --- public interpolators ---
    def gstar(self, T):
        return self._interp(T, self.gvec)

    def gstarS(self, T):
        return self._interp(T, self.gSvec)

    def dlngdlnT(self, T):
        return self._interp(T, self._dlngdlnT_vec)

    def dlngSdlnT(self, T):
        return self._interp(T, self._dlngSdlnT_vec)

    # ----------------------------------------------------------
    #  Background thermodynamics
    # ----------------------------------------------------------
    def rho_rad(self, T):
        """SM radiation energy density  rho_rad = (pi^2/30) g*(T) T^4."""
        T = np.asarray(T, float)
        return (np.pi**2 / 30.0) * self.gstar(T) * T**4

    def s_entropy(self, T):
        """SM entropy density  s = (2 pi^2/45) g*_S(T) T^3."""
        T = np.asarray(T, float)
        return (2.0 * np.pi**2 / 45.0) * self.gstarS(T) * T**3

    def hubble(self, T):
        """Hubble rate from SM radiation only."""
        T = np.asarray(T, float)
        return np.sqrt((np.pi**2 / 90.0) * self.gstar(T)) * T**2 / Mpl

    def g_tilde(self, T):
        """g_tilde = 1 + (1/3) d ln g*_S / d ln T."""
        return 1.0 + self.dlngSdlnT(T) / 3.0

    # ----------------------------------------------------------
    #  Bessel helpers
    # ----------------------------------------------------------
    @staticmethod
    def K1_over_K2(z: float) -> float:
        """K_1(z)/K_2(z) with asymptotic expansion for large z."""
        if z < 50.0:
            return float(Knu(1, z) / Knu(2, z))
        return 1.0 - 1.5 / z + 1.875 / z**2 - 6.5625 / z**3

    # ----------------------------------------------------------
    #  B-factors  (energy-weighted MB moments)
    # ----------------------------------------------------------
    @staticmethod
    def Bfac1(z, m, z_switch=50.0):
        """B_1 = m (K_1/K_2 + 3/z)."""
        if z < z_switch:
            R = scisp.kv(1, z) / scisp.kv(2, z)
        else:
            R = 1.0 - 3.0 / (2 * z) + 15.0 / (8 * z**2)
        return m * (R + 3.0 / z)

    @staticmethod
    def Bfac2(z, m, z_switch=50.0):
        """B_2 = m (K_1/K_2 + 4/z)."""
        if z < z_switch:
            R = scisp.kv(1, z) / scisp.kv(2, z)
        else:
            R = 1.0 - 3.0 / (2 * z) + 15.0 / (8 * z**2)
        return m * (R + 4.0 / z)

    @staticmethod
    def dBfac1dT(z, m, z_switch=50.0):
        """dB_1/dT_dark (derivative w.r.t. dark temperature)."""
        if z < z_switch:
            R = scisp.kv(1, z) / scisp.kv(2, z)
            G = R**2 - 1.0 + 3.0 / z * R
        else:
            G = 3.0 / (2 * z**2) - 15.0 / (4 * z**3)
        return 3.0 - z**2 * G

    @staticmethod
    def dBfac2dT(z, m, z_switch=50.0):
        """dB_2/dT_dark."""
        if z < z_switch:
            R = scisp.kv(1, z) / scisp.kv(2, z)
            G = R**2 - 1.0 + 3.0 / z * R
        else:
            G = 3.0 / (2 * z**2) - 15.0 / (4 * z**3)
        return 4.0 - z**2 * G

    # ----------------------------------------------------------
    #  Equilibrium number density (Maxwell-Boltzmann)
    # ----------------------------------------------------------
    @staticmethod
    def neq_MB(z, g, m, include_antiparticles=False):
        """n_eq in the MB approximation:  g m^3/(2 pi^2) K_2(z)/z."""
        if z > 500:
            return 0.0
        neq = g * m**3 / (2.0 * np.pi**2) * scisp.kv(2, z) / z
        if include_antiparticles:
            neq *= 2
        return neq

    @staticmethod
    def dneqdT_MB(z, g, m, include_antiparticles=False):
        """dn_eq/dT in the MB approximation."""
        if z > 500:
            return 0.0
        val = g * m**2 / (2.0 * np.pi**2) * (z * scisp.kv(1, z) + 3 * scisp.kv(2, z))
        if include_antiparticles:
            val *= 2
        return val

    @staticmethod
    def ln_Yeq(z, g, m, include_antiparticles, s):
        """
        ln(n_eq / s), safe for any z.
        Uses the asymptotic Bessel expansion for z >= 300.
        """
        if z < 300:
            neq = Cosmology.neq_MB(z, g, m, include_antiparticles)
            if neq > 0 and s > 0:
                return np.log(neq / s)
        ln_neq = (np.log(g * m**3 / (2.0 * np.pi**2))
                  + 0.5 * np.log(np.pi / (2.0 * z)) - z
                  + np.log(1.0 + 15.0 / (8 * z) + 105.0 / (128 * z**2))
                  - np.log(z))
        if include_antiparticles:
            ln_neq += np.log(2)
        return ln_neq - np.log(s) if s > 0 else -1e10

    @staticmethod
    def lambda_i(z):
        """d ln n_i^eq / d ln T_dark = z K_1(z)/K_2(z) + 3."""
        if z < 50:
            return z * scisp.kv(1, z) / scisp.kv(2, z) + 3
        return z * (1.0 - 3.0 / (2 * z) + 15.0 / (8 * z**2)) + 3

    # ----------------------------------------------------------
    #  Full species thermodynamics (integral / MB hybrid)
    # ----------------------------------------------------------
    @staticmethod
    def _eta(statistics: str) -> float:
        s = statistics.lower()
        if s not in ("fermion", "boson"):
            raise ValueError("statistics must be 'fermion' or 'boson'.")
        return +1.0 if s == "fermion" else -1.0

    def species_thermo(
        self,
        m: float,
        T: float,
        g: float = 1.0,
        statistics: str = "fermion",
        include_antiparticles: bool = False,
        z_rel: float = 1e-3,
        z_switch: float = 1.0,
        what: str = "all",
        quad_epsrel: float = 1e-9,
        warn_relerr: float = 1e-6,
    ):
        """
        Full species thermodynamics (rho, P, n) with three regimes:
            z < z_rel    -> ultra-relativistic analytic
            z > z_switch -> Maxwell-Boltzmann
            otherwise    -> numerical integral
        """
        if T <= 0:
            raise ValueError("T must be positive.")
        if m < 0:
            raise ValueError("m must be non-negative.")
        what = what.lower()
        if what not in ("n", "rho", "p", "all"):
            raise ValueError("what must be 'n', 'rho', 'P', or 'all'.")

        eta  = self._eta(statistics)
        z    = m / T
        mult = 2.0 if include_antiparticles else 1.0

        # --- ultra-relativistic analytic ---
        if z <= z_rel:
            if statistics.lower() == "boson":
                rho = g * (np.pi**2 / 30.0) * T**4
                n   = g * (zeta(3) / np.pi**2) * T**3
            else:
                rho = g * (7.0 / 8.0) * (np.pi**2 / 30.0) * T**4
                n   = g * (3.0 / 4.0) * (zeta(3) / np.pi**2) * T**3
            P = rho / 3.0
            rho *= mult; P *= mult; n *= mult
            if what == "rho": return rho
            if what == "p":   return P
            if what == "n":   return n
            return {"rho": rho, "P": P, "n": n, "z": z, "regime": "UR"}

        # --- Maxwell-Boltzmann ---
        if z >= z_switch:
            K2 = Knu(2, z)
            n  = g * T**3 / (2.0 * np.pi**2) * z**2 * K2 * mult
            if what == "n":
                return n
            P = T * n
            R = self.K1_over_K2(z)
            B1 = z * R + 3.0
            rho = B1 * n
            if what == "p":   return P
            if what == "rho": return rho
            return {"rho": rho, "P": P, "n": n, "z": z, "regime": "Maxwell-Boltzmann"}

        # --- Integral region ---
        def denom_fn(x):
            if x > 50.0:
                return np.exp(-x)
            return 1.0 / (np.exp(x) + eta)

        def I_rho(u):
            x = np.sqrt(u * u + z * z)
            return u * u * x * denom_fn(x)

        def I_P(u):
            x = np.sqrt(u * u + z * z)
            return (u**4 / x) * denom_fn(x)

        def I_n(u):
            x = np.sqrt(u * u + z * z)
            return u * u * denom_fn(x)

        need_rho = what in ("rho", "all")
        need_P   = what in ("p",   "all")
        need_n   = what in ("n",   "all")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", IntegrationWarning)
            if need_rho:
                Irho, erho = quad(I_rho, 0.0, np.inf, epsabs=0.0,
                                  epsrel=quad_epsrel, limit=300)
                rrho = erho / max(abs(Irho), 1e-300)
            if need_P:
                IP, eP = quad(I_P, 0.0, np.inf, epsabs=0.0,
                              epsrel=quad_epsrel, limit=300)
                rP = eP / max(abs(IP), 1e-300)
            if need_n:
                In, en = quad(I_n, 0.0, np.inf, epsabs=0.0,
                              epsrel=quad_epsrel, limit=300)
                rn = en / max(abs(In), 1e-300)

            rels = []
            if need_rho: rels.append(rrho)
            if need_P:   rels.append(rP)
            if need_n:   rels.append(rn)
            if (rels
                and any(isinstance(ww.message, IntegrationWarning) for ww in w)
                and max(rels) > warn_relerr):
                warnings.warn(
                    f"quad IntegrationWarning: rel errors {rels} at z=m/T={z:.3g}",
                    IntegrationWarning,
                )

        prefac = mult * g * T**4 / (2.0 * np.pi**2)
        if need_rho:
            rho = prefac * Irho
        if need_P:
            P = prefac / 3.0 * IP
        if need_n:
            n = mult * g * T**3 / (2.0 * np.pi**2) * In

        if what == "rho": return rho
        if what == "p":   return P
        if what == "n":   return n
        return {"rho": rho, "P": P, "n": n, "z": z, "regime": "Integral"}

    # ----------------------------------------------------------
    #  Numerical derivative helper
    # ----------------------------------------------------------
    @staticmethod
    def dndT_nearest_stable(Ts, ns):
        """
        Stable numerical dn/dT on a log-spaced positive grid.
        Uses nearest-neighbour finite differences in log-log space.
        """
        Ts = np.asarray(Ts)
        ns = np.asarray(ns)
        tiny = np.finfo(float).tiny
        ns_pos = np.maximum(ns, tiny)

        x = np.log(Ts)
        y = np.log(ns_pos)

        dlogn = np.empty_like(ns_pos)
        dlogn[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
        dlogn[0]    = (y[1]  - y[0])   / (x[1]  - x[0])
        dlogn[-1]   = (y[-1] - y[-2])  / (x[-1] - x[-2])

        return (ns_pos / Ts) * dlogn

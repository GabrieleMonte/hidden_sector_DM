"""
Dark-sector model interface and concrete implementations.

DarkSectorModel : abstract base class defining the rate interface
VectorPortal    : secluded dark sector with kinetic mixing (hypercharge)
BLPortal        : secluded dark sector with kinetic mixing (B-L)
LiLjPortal      : secluded dark sector with kinetic mixing (Li-Lj)
BaryonPortal    : secluded dark sector with kinetic mixing (baryon number)
HiggsPortal     : secluded dark sector with scalar mixing
"""

import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

from .constants import (SM_MYQ, BL_FERMIONS, LILJ_FERMIONS, BARYON_FERMIONS,
                        MZ, gW, sW, cW, gY as gY_SM, eC,
                        Mh, vH, mN, HBARC2)

# ---------------------------------------------------------------------------
#  Direct detection helpers
# ---------------------------------------------------------------------------
_Q_u, _T3_u = 2.0 / 3.0, 0.5    # up-quark charge and weak isospin
_Q_d, _T3_d = -1.0 / 3.0, -0.5  # down-quark charge and weak isospin
_A_Xe, _Z_Xe = 131, 54           # Xenon-131 (LZ / XENON default)


# =====================================================================
#  Abstract base
# =====================================================================

class DarkSectorModel(ABC):
    """
    Interface that any secluded dark-sector model must implement.

    The BoltzmannSolver reads masses, DOF, and rate methods from this
    interface — all model-specific physics lives here.
    """

    # --- Required attributes (set by subclass __init__) ---
    mX: float                      # DM mass  [GeV]
    mY: float                      # mediator mass  [GeV]
    gX: int                        # DM internal DOF
    gY: int                        # mediator internal DOF
    include_antiparticlesX: bool
    include_antiparticlesY: bool
    xi_infl: float                 # T_dark / T_SM at reheating

    # --- 2 -> 2 annihilation ---

    @abstractmethod
    def sigmav_XX_to_YY(self, TD: float) -> float:
        """Full thermally-averaged <sigma v> for XX -> YY at dark temperature TD."""

    @abstractmethod
    def sigmav_XX_to_YY_swave(self) -> float:
        """s-wave-only <sigma v> for XX -> YY (temperature-independent)."""

    # --- 3 -> 2 cannibal processes ---

    @abstractmethod
    def sigmav2_YYY_to_YY(self, TD: float) -> float:
        """<sigma^2 v^2> for YYY -> YY at dark temperature TD."""

    @abstractmethod
    def sigmav2_YXX_to_XX(self) -> float:
        """<sigma^2 v^2> for YXX -> XX (temperature-independent)."""

    @abstractmethod
    def sigmav2_YYX_to_YX(self) -> float:
        """<sigma^2 v^2> for YYX -> YX (temperature-independent)."""

    # --- Decay ---

    @abstractmethod
    def decay_width_to_SM(self, epsX: float) -> float:
        """Total decay width Gamma(Y -> SM) for portal coupling epsX."""

    @abstractmethod
    def branching_ratios_to_SM(self) -> dict:
        """Branching ratios for mediator decay to SM final states.

        Returns dict {channel: BR} where channel names follow the
        convention: 'ee', 'mumu', 'tautau', 'uu', 'dd', 'ss', 'cc',
        'bb', 'tt', 'nunu_e', 'nunu_mu', 'nunu_tau', 'gg', 'gamgam',
        'Zgam', 'WW', 'ZZ', 'hh', 'pipi', 'KK', etc.

        Branching ratios are independent of the portal coupling epsX.
        """

    # --- Direct detection ---

    @abstractmethod
    def sigma_SI(self, epsX: float, A_N: int = _A_Xe, Z_N: int = _Z_Xe) -> float:
        """Spin-independent DM-nucleon cross section [cm^2].

        Parameters
        ----------
        epsX : float
            Portal coupling (kinetic mixing parameter).
        A_N : int
            Atomic mass number of target (default: Xe-131).
        Z_N : int
            Atomic number of target (default: Xe, Z=54).
        """


# Fermion key -> pair channel name
_PAIR_NAME = {
    "e": "ee", "mu": "mumu", "tau": "tautau",
    "u": "uu", "d": "dd", "s": "ss", "c": "cc", "b": "bb", "t": "tt",
    "nu_e": "nunu_e", "nu_mu": "nunu_mu", "nu_tau": "nunu_tau",
}


# =====================================================================
#  VectorPortal
# =====================================================================

@dataclass
class VectorPortal(DarkSectorModel):
    """
    Secluded dark sector coupled to the SM via kinetic mixing.

    Provides Z'-fermion couplings, decay widths, and all cross
    sections needed by the BoltzmannSolver.

    Parameters
    ----------
    mX, mY : float
        Dark-matter and mediator masses (GeV).
    gX, gY : int
        Internal degrees of freedom.
    alphaX : float
        Dark fine-structure constant.
    include_antiparticlesX, include_antiparticlesY : bool
        Whether to double-count for Xbar or Ybar.
    Delta1, Delta2 : float
        Phenomenological coefficients for YYY->YY and YXX->XX.
    Delta3 : float or None
        Coefficient for YYX->YX.  If None (default), the exact
        tree-level result is computed from mX/mY.
    xi_infl : float
        Initial dark-to-SM temperature ratio.
    """
    mX: float = 100.0
    mY: float = 10.0
    gX: int   = 2
    gY: int   = 3
    alphaX: float = 1e-2
    include_antiparticlesX: bool = True
    include_antiparticlesY: bool = False
    Delta1: float = 0.5
    Delta2: float = 0.5
    Delta3: Optional[float] = None
    xi_infl: float = 1.0

    def __post_init__(self):
        if self.Delta3 is None:
            self.Delta3 = self._Delta3_exact(self.mX / self.mY)

    # ----------------------------------------------------------
    #  Exact Delta3 coefficient
    # ----------------------------------------------------------
    @staticmethod
    def _Delta3_exact(r: float) -> float:
        """
        Exact tree-level Delta3 for Z'Z'X -> Z'X.

        Parameters
        ----------
        r : float
            Mass ratio mX / mY.
        """
        return (
            (4.0 * np.pi**2
             * r**5
             * np.sqrt(4.0 / 3.0 * r * (r + 2) + 1)
             * (2 * r * (r * (r * (r * (2 * r * (r * (8 * r * (10 * r + 71)
               + 1683) + 2798) + 6107) + 4722) + 2335) + 578) + 195))
            / ((r + 2)**3 * (2 * r + 1)**4
               * (2 * (r + 2) * r**2 + r - 1)**2)
        )

    # ----------------------------------------------------------
    #  DarkSectorModel interface
    # ----------------------------------------------------------

    def sigmav_XX_to_YY(self, TD: float) -> float:
        return self._sigmav_XX_to_YY_full(self.mX, self.mY, self.alphaX, TD)

    def sigmav_XX_to_YY_swave(self) -> float:
        return self._sigmav_XX_to_YY_swave(self.mX, self.alphaX)

    def sigmav2_YYY_to_YY(self, TD: float) -> float:
        return self._sigmav2_YYY_to_YY(self.mX, self.alphaX, TD, self.Delta1)

    def sigmav2_YXX_to_XX(self) -> float:
        return self._sigmav2_YXX_to_XX(self.mX, self.alphaX, self.Delta2)

    def sigmav2_YYX_to_YX(self) -> float:
        return self._sigmav2_YYX_to_YX(self.mX, self.alphaX, self.Delta3)

    def decay_width_to_SM(self, epsX: float) -> float:
        return self.decay_width_to_all(self.mY, epsX)

    def branching_ratios_to_SM(self) -> dict:
        partials = {}
        total = 0.0
        for f, (mf, *_rest) in SM_MYQ.items():
            if self.mY > 2.0 * mf:
                w = self.decay_width_to_ff(f, self.mY, 1.0)
                partials[_PAIR_NAME[f]] = w
                total += w
        if total <= 0:
            return {ch: 0.0 for ch in partials}
        return {ch: w / total for ch, w in partials.items()}

    def sigma_SI(self, epsX: float, A_N: int = _A_Xe, Z_N: int = _Z_Xe) -> float:
        """
        Spin-independent DM-nucleon cross section via Z-Z' mixing.

        Full formula from Eq. 3.13 of arXiv:2509.08043, including both
        Z and Z' exchange.  Significant isospin violation (sigma_n/sigma_p ~ 9).

        Parameters
        ----------
        epsX : float
            Kinetic mixing parameter.
        A_N : int
            Atomic mass number (default: Xe-131).
        Z_N : int
            Atomic number (default: Xe, Z=54).

        Returns
        -------
        sigma : float
            SI cross section [cm^2].
        """
        gD = np.sqrt(4.0 * np.pi * self.alphaX)
        mZp = self.mY

        # Z-Z' mixing angle
        tan2xi = -2.0 * epsX * sW * MZ**2 / (MZ**2 - mZp**2)
        xi = 0.5 * np.arctan(tan2xi)
        s_xi, c_xi = np.sin(xi), np.cos(xi)

        # DM couplings  (Eqs. 3.16-3.17)
        V_chi_Z = gD * s_xi
        V_chi_Zp = gD * c_xi

        # Quark couplings  (Eqs. 3.14-3.15)
        eswcw = eC / (sW * cW)
        V_u_Z = (eC * epsX * s_xi * _Q_u
                 + eswcw * (c_xi + epsX * sW * s_xi) * (_T3_u - sW**2 * _Q_u))
        V_d_Z = (eC * epsX * s_xi * _Q_d
                 + eswcw * (c_xi + epsX * sW * s_xi) * (_T3_d - sW**2 * _Q_d))
        V_u_Zp = (eC * epsX * c_xi * _Q_u
                  + eswcw * (s_xi - epsX * sW * c_xi) * (_T3_u - sW**2 * _Q_u))
        V_d_Zp = (eC * epsX * c_xi * _Q_d
                  + eswcw * (s_xi - epsX * sW * c_xi) * (_T3_d - sW**2 * _Q_d))

        mu = self.mX * mN / (self.mX + mN)

        bracket = 0.0
        for m_i, Vchi, Vu, Vd in [(MZ,  V_chi_Z,  V_u_Z,  V_d_Z),
                                   (mZp, V_chi_Zp, V_u_Zp, V_d_Zp)]:
            bracket += Vchi * (Z_N * (2*Vu + Vd)
                               + (A_N - Z_N) * (Vu + 2*Vd)) / (A_N * m_i**2)

        return mu**2 / np.pi * bracket**2 * HBARC2

    # ----------------------------------------------------------
    #  Z'-fermion couplings
    # ----------------------------------------------------------
    @staticmethod
    def gfv_gfa(f: str, mZp: float, eps: float):
        """
        Return (g_fV, g_fA) for SM charged fermion *f*.

        Uses g_{R,L} = eps (mZp^2 gY Y_{R,L} - mZ^2 g sW cW Q) / (mZ^2 - mZp^2)
        and  gV/A = (gR +/- gL) / 2.
        """
        _, Q, YL, YR, _ = SM_MYQ[f]
        den    = MZ**2 - mZp**2
        common = MZ**2 * gW * sW * cW * Q
        gL = eps * (mZp**2 * gY_SM * YL - common) / den
        gR = eps * (mZp**2 * gY_SM * YR - common) / den
        return 0.5 * (gR + gL), 0.5 * (gR - gL)

    # ----------------------------------------------------------
    #  Decay widths
    # ----------------------------------------------------------
    @staticmethod
    def decay_width_to_ff(f: str, mA: float, epsX: float) -> float:
        """Partial width Gamma(A' -> f fbar)."""
        gVf, gAf = VectorPortal.gfv_gfa(f, mA, epsX)
        NCf = SM_MYQ[f][-1]
        mf  = SM_MYQ[f][0]
        r2  = 4.0 * mf**2 / mA**2
        return (NCf * mA / (12.0 * np.pi)
                * np.sqrt(1.0 - r2)
                * (gVf**2 * (1.0 + 0.5 * r2) + gAf**2 * (1.0 - r2)))

    @staticmethod
    def decay_width_to_all(mA: float, epsX: float) -> float:
        """Total visible width Gamma(A' -> SM), summing over open channels."""
        total = 0.0
        for f in SM_MYQ:
            if mA > 2.0 * SM_MYQ[f][0]:
                total += VectorPortal.decay_width_to_ff(f, mA, epsX)
        return total

    # ----------------------------------------------------------
    #  Cross-section formulas (static, used by the interface)
    # ----------------------------------------------------------
    @staticmethod
    def _sigmav_XX_to_YY_swave(mX: float, alphaX: float) -> float:
        """<sigma v> for X Xbar -> Y Y  (s-wave, mY -> 0 limit)."""
        return np.pi * alphaX**2 / mX**2

    @staticmethod
    def _sigmav_XX_to_YY_full(mX: float, mY: float, alphaX: float, TX: float) -> float:
        """<sigma v> for X Xbar -> Y Y  including finite mY and p-wave."""
        rv = mY / mX
        xX = mX / TX
        av = (2.0 * np.pi * alphaX**2 / mX**2
              * (1.0 - rv**2)**1.5 / (2.0 - rv**2)**2)
        bv = (np.pi * alphaX**2 / (12.0 * mX**2)
              * np.sqrt(1.0 - rv**2)
              * (24.0 + 28.0 * rv**2 - 36.0 * rv**4 + 17.0 * rv**6)
              / (2.0 - rv**2)**4)
        return 2.0 * (av + 6.0 * bv / xX)

    @staticmethod
    def _sigmav2_YYY_to_YY(mX: float, alphaX: float, Th: float,
                            Delta1: float = 0.5) -> float:
        """3 -> 2 cannibal:  Y Y Y -> Y Y."""
        return Delta1 * alphaX**5 * Th**7 / mX**12

    @staticmethod
    def _sigmav2_YXX_to_XX(mX: float, alphaX: float,
                            Delta2: float = 0.5) -> float:
        """3 -> 2:  Y X X -> X X."""
        return Delta2 * alphaX**3 / mX**5

    @staticmethod
    def _sigmav2_YYX_to_YX(mX: float, alphaX: float,
                            Delta3: float = 0.5) -> float:
        """3 -> 2:  Y Y X -> Y X."""
        return Delta3 * alphaX**3 / mX**5


# =====================================================================
#  BLPortal
# =====================================================================

@dataclass
class BLPortal(DarkSectorModel):
    """
    Secluded dark sector coupled to the SM via B-L kinetic mixing.

    The hidden-sector Z' mixes kinetically with the U(1)_{B-L} gauge
    boson Z_{B-L}.  After restoring canonical kinetic terms the Z'
    acquires purely vectorial couplings to SM fermions proportional
    to their (B-L) charge (Eq. 8 of arXiv:1912.08821).

    Dark-sector cross sections (XX->YY, 3->2 cannibal) are identical
    to the VectorPortal — only the decay Z'->SM changes.

    Parameters
    ----------
    mX, mY : float
        Dark-matter and mediator masses (GeV).
    alphaX : float
        Dark fine-structure constant.
    g_BL : float
        U(1)_{B-L} gauge coupling.
    m_ZBL : float
        Mass of the Z_{B-L} gauge boson (GeV).
    Delta1, Delta2 : float
        Phenomenological coefficients for YYY->YY and YXX->XX.
    Delta3 : float or None
        Coefficient for YYX->YX.  If None, exact tree-level result.
    """
    mX: float = 100.0
    mY: float = 10.0
    gX: int   = 2
    gY: int   = 3
    alphaX: float = 1e-2
    g_BL: float = 1.0
    m_ZBL: float = 1000.0
    include_antiparticlesX: bool = True
    include_antiparticlesY: bool = False
    Delta1: float = 0.5
    Delta2: float = 0.5
    Delta3: Optional[float] = None
    xi_infl: float = 1.0

    def __post_init__(self):
        if self.Delta3 is None:
            self.Delta3 = VectorPortal._Delta3_exact(self.mX / self.mY)

    # ----------------------------------------------------------
    #  B-L fermion couplings  (Eq. 8 of 1912.08821)
    # ----------------------------------------------------------
    @staticmethod
    def gfv_gfa_BL(f: str, mZp: float, eps: float,
                   g_BL: float, m_ZBL: float):
        """
        Return (g_fV, g_fA) for fermion *f* in the B-L portal.

        g_fV = eps * g_BL * (B-L)_f * |(m_ZBL^2 + m_Zp^2) / (m_ZBL^2 - m_Zp^2)|
        g_fA = 0
        """
        _, BL_charge, _ = BL_FERMIONS[f]
        mixing = abs((m_ZBL**2 + mZp**2) / (m_ZBL**2 - mZp**2))
        gfV = eps * g_BL * BL_charge * mixing
        return gfV, 0.0

    # ----------------------------------------------------------
    #  Decay widths
    # ----------------------------------------------------------
    @staticmethod
    def decay_width_to_ff_BL(f: str, mZp: float, epsX: float,
                             g_BL: float, m_ZBL: float) -> float:
        """Partial width Gamma(Z' -> f fbar) via B-L portal."""
        mf, _, Nc = BL_FERMIONS[f]
        if mZp <= 2.0 * mf:
            return 0.0
        gfV, _ = BLPortal.gfv_gfa_BL(f, mZp, epsX, g_BL, m_ZBL)
        r2 = 4.0 * mf**2 / mZp**2
        return (Nc * mZp / (12.0 * np.pi)
                * np.sqrt(1.0 - r2)
                * gfV**2 * (1.0 + 0.5 * r2))

    # ----------------------------------------------------------
    #  DarkSectorModel interface
    # ----------------------------------------------------------

    def sigmav_XX_to_YY(self, TD: float) -> float:
        return VectorPortal._sigmav_XX_to_YY_full(
            self.mX, self.mY, self.alphaX, TD)

    def sigmav_XX_to_YY_swave(self) -> float:
        return VectorPortal._sigmav_XX_to_YY_swave(self.mX, self.alphaX)

    def sigmav2_YYY_to_YY(self, TD: float) -> float:
        return VectorPortal._sigmav2_YYY_to_YY(
            self.mX, self.alphaX, TD, self.Delta1)

    def sigmav2_YXX_to_XX(self) -> float:
        return VectorPortal._sigmav2_YXX_to_XX(
            self.mX, self.alphaX, self.Delta2)

    def sigmav2_YYX_to_YX(self) -> float:
        return VectorPortal._sigmav2_YYX_to_YX(
            self.mX, self.alphaX, self.Delta3)

    def decay_width_to_SM(self, epsX: float) -> float:
        total = 0.0
        for f in BL_FERMIONS:
            total += self.decay_width_to_ff_BL(
                f, self.mY, epsX, self.g_BL, self.m_ZBL)
        return total

    def branching_ratios_to_SM(self) -> dict:
        partials = {}
        total = 0.0
        for f, (mf, _, Nc) in BL_FERMIONS.items():
            if self.mY > 2.0 * mf:
                w = self.decay_width_to_ff_BL(
                    f, self.mY, 1.0, self.g_BL, self.m_ZBL)
                partials[_PAIR_NAME[f]] = w
                total += w
        if total <= 0:
            return {ch: 0.0 for ch in partials}
        return {ch: w / total for ch, w in partials.items()}

    def sigma_SI(self, epsX: float, A_N: int = _A_Xe, Z_N: int = _Z_Xe) -> float:
        """
        SI DM-nucleon cross section for the B-L portal.

        Purely vectorial, isospin-universal: all quarks have (B-L) = 1/3,
        so proton and neutron couplings are equal and the result is
        independent of (A_N, Z_N).
        """
        gD = np.sqrt(4.0 * np.pi * self.alphaX)
        mixing = abs((self.m_ZBL**2 + self.mY**2)
                     / (self.m_ZBL**2 - self.mY**2))
        gqV = epsX * self.g_BL * (1.0 / 3.0) * mixing
        mu = self.mX * mN / (self.mX + mN)
        amplitude = gD * 3.0 * gqV / self.mY**2
        return mu**2 / np.pi * amplitude**2 * HBARC2


# =====================================================================
#  LiLjPortal
# =====================================================================

_VALID_FLAVORS = ("e-mu", "mu-tau", "e-tau")

@dataclass
class LiLjPortal(DarkSectorModel):
    """
    Secluded dark sector coupled to the SM via L_i - L_j kinetic mixing.

    The hidden-sector Z' mixes kinetically with the gauge boson of a
    U(1)_{L_i - L_j} symmetry.  The Z' acquires purely vectorial couplings
    to leptons and neutrinos of the two involved families (Eq. 10 of
    arXiv:1912.08821).  Quarks are uncharged — the Z' decays entirely
    to leptons at tree level.

    Parameters
    ----------
    mX, mY : float
        Dark-matter and mediator masses (GeV).
    alphaX : float
        Dark fine-structure constant.
    flavor : str
        Which lepton-flavor combination: 'e-mu', 'mu-tau', or 'e-tau'.
    g_LiLj : float
        U(1)_{L_i - L_j} gauge coupling.
    m_ZLiLj : float
        Mass of the Z_{L_i - L_j} gauge boson (GeV).
    Delta1, Delta2 : float
        Phenomenological coefficients for YYY->YY and YXX->XX.
    Delta3 : float or None
        Coefficient for YYX->YX.  If None, exact tree-level result.
    """
    mX: float = 100.0
    mY: float = 10.0
    gX: int   = 2
    gY: int   = 3
    alphaX: float = 1e-2
    flavor: str = "mu-tau"
    g_LiLj: float = 1.0
    m_ZLiLj: float = 1000.0
    include_antiparticlesX: bool = True
    include_antiparticlesY: bool = False
    Delta1: float = 0.5
    Delta2: float = 0.5
    Delta3: Optional[float] = None
    xi_infl: float = 1.0

    def __post_init__(self):
        if self.flavor not in _VALID_FLAVORS:
            raise ValueError(
                f"flavor must be one of {_VALID_FLAVORS}, got '{self.flavor}'")
        if self.Delta3 is None:
            self.Delta3 = VectorPortal._Delta3_exact(self.mX / self.mY)
        self._fermion_table = LILJ_FERMIONS[self.flavor]

    # ----------------------------------------------------------
    #  Li-Lj fermion couplings  (Eq. 10 of 1912.08821)
    # ----------------------------------------------------------
    @staticmethod
    def gfv_gfa_LiLj(f: str, mZp: float, eps: float,
                     g_LiLj: float, m_ZLiLj: float,
                     fermion_table: dict):
        """
        Return (g_fV, g_fA) for fermion *f* in the Li-Lj portal.

        g_fV = eps * g_{Li-Lj} * (Li-Lj)_f * |(m²_{ZLiLj} + m²_{Z'}) / (m²_{ZLiLj} - m²_{Z'})|
        g_fA = 0
        """
        _, charge, _ = fermion_table[f]
        mixing = abs((m_ZLiLj**2 + mZp**2) / (m_ZLiLj**2 - mZp**2))
        gfV = eps * g_LiLj * charge * mixing
        return gfV, 0.0

    # ----------------------------------------------------------
    #  DarkSectorModel interface
    # ----------------------------------------------------------

    def sigmav_XX_to_YY(self, TD: float) -> float:
        return VectorPortal._sigmav_XX_to_YY_full(
            self.mX, self.mY, self.alphaX, TD)

    def sigmav_XX_to_YY_swave(self) -> float:
        return VectorPortal._sigmav_XX_to_YY_swave(self.mX, self.alphaX)

    def sigmav2_YYY_to_YY(self, TD: float) -> float:
        return VectorPortal._sigmav2_YYY_to_YY(
            self.mX, self.alphaX, TD, self.Delta1)

    def sigmav2_YXX_to_XX(self) -> float:
        return VectorPortal._sigmav2_YXX_to_XX(
            self.mX, self.alphaX, self.Delta2)

    def sigmav2_YYX_to_YX(self) -> float:
        return VectorPortal._sigmav2_YYX_to_YX(
            self.mX, self.alphaX, self.Delta3)

    def decay_width_to_SM(self, epsX: float) -> float:
        total = 0.0
        for f, (mf, charge, Nc) in self._fermion_table.items():
            if self.mY <= 2.0 * mf:
                continue
            gfV, _ = self.gfv_gfa_LiLj(
                f, self.mY, epsX, self.g_LiLj, self.m_ZLiLj,
                self._fermion_table)
            r2 = 4.0 * mf**2 / self.mY**2
            total += (Nc * self.mY / (12.0 * np.pi)
                      * np.sqrt(1.0 - r2)
                      * gfV**2 * (1.0 + 0.5 * r2))
        return total

    def branching_ratios_to_SM(self) -> dict:
        partials = {}
        total = 0.0
        for f, (mf, charge, Nc) in self._fermion_table.items():
            if self.mY <= 2.0 * mf:
                continue
            gfV, _ = self.gfv_gfa_LiLj(
                f, self.mY, 1.0, self.g_LiLj, self.m_ZLiLj,
                self._fermion_table)
            r2 = 4.0 * mf**2 / self.mY**2
            w = (Nc * self.mY / (12.0 * np.pi)
                 * np.sqrt(1.0 - r2)
                 * gfV**2 * (1.0 + 0.5 * r2))
            partials[_PAIR_NAME[f]] = w
            total += w
        if total <= 0:
            return {ch: 0.0 for ch in partials}
        return {ch: w / total for ch, w in partials.items()}

    def sigma_SI(self, epsX: float, A_N: int = _A_Xe, Z_N: int = _Z_Xe) -> float:
        """
        SI DM-nucleon cross section for L_i - L_j portal.

        Quarks are uncharged under L_i - L_j, so the tree-level
        DM-nucleon cross section is zero.
        """
        return 0.0


# =====================================================================
#  BaryonPortal
# =====================================================================

@dataclass
class BaryonPortal(DarkSectorModel):
    """
    Secluded dark sector coupled to the SM via baryon-number kinetic mixing.

    The hidden-sector Z' mixes kinetically with the gauge boson of a
    U(1)_B symmetry.  The Z' acquires purely vectorial couplings to
    quarks proportional to their baryon number (Eq. 9 of arXiv:1912.08821).
    Leptons are uncharged — the Z' decays entirely to quarks at tree level.

    Parameters
    ----------
    mX, mY : float
        Dark-matter and mediator masses (GeV).
    alphaX : float
        Dark fine-structure constant.
    g_B : float
        U(1)_B gauge coupling.
    m_ZB : float
        Mass of the Z_B gauge boson (GeV).
    Delta1, Delta2 : float
        Phenomenological coefficients for YYY->YY and YXX->XX.
    Delta3 : float or None
        Coefficient for YYX->YX.  If None, exact tree-level result.
    """
    mX: float = 100.0
    mY: float = 10.0
    gX: int   = 2
    gY: int   = 3
    alphaX: float = 1e-2
    g_B: float = 1.0
    m_ZB: float = 1000.0
    include_antiparticlesX: bool = True
    include_antiparticlesY: bool = False
    Delta1: float = 0.5
    Delta2: float = 0.5
    Delta3: Optional[float] = None
    xi_infl: float = 1.0

    def __post_init__(self):
        if self.Delta3 is None:
            self.Delta3 = VectorPortal._Delta3_exact(self.mX / self.mY)

    # ----------------------------------------------------------
    #  Baryon-number couplings  (Eq. 9 of 1912.08821)
    # ----------------------------------------------------------
    @staticmethod
    def gfv_gfa_B(f: str, mZp: float, eps: float,
                  g_B: float, m_ZB: float):
        """
        Return (g_fV, g_fA) for quark *f* in the baryon portal.

        g_fV = eps * g_B * B_f * |(m_ZB^2 + m_Zp^2) / (m_ZB^2 - m_Zp^2)|
        g_fA = 0
        """
        _, B_charge, _ = BARYON_FERMIONS[f]
        mixing = abs((m_ZB**2 + mZp**2) / (m_ZB**2 - mZp**2))
        gfV = eps * g_B * B_charge * mixing
        return gfV, 0.0

    # ----------------------------------------------------------
    #  DarkSectorModel interface
    # ----------------------------------------------------------

    def sigmav_XX_to_YY(self, TD: float) -> float:
        return VectorPortal._sigmav_XX_to_YY_full(
            self.mX, self.mY, self.alphaX, TD)

    def sigmav_XX_to_YY_swave(self) -> float:
        return VectorPortal._sigmav_XX_to_YY_swave(self.mX, self.alphaX)

    def sigmav2_YYY_to_YY(self, TD: float) -> float:
        return VectorPortal._sigmav2_YYY_to_YY(
            self.mX, self.alphaX, TD, self.Delta1)

    def sigmav2_YXX_to_XX(self) -> float:
        return VectorPortal._sigmav2_YXX_to_XX(
            self.mX, self.alphaX, self.Delta2)

    def sigmav2_YYX_to_YX(self) -> float:
        return VectorPortal._sigmav2_YYX_to_YX(
            self.mX, self.alphaX, self.Delta3)

    def decay_width_to_SM(self, epsX: float) -> float:
        total = 0.0
        for f, (mf, B_charge, Nc) in BARYON_FERMIONS.items():
            if self.mY <= 2.0 * mf:
                continue
            gfV, _ = self.gfv_gfa_B(
                f, self.mY, epsX, self.g_B, self.m_ZB)
            r2 = 4.0 * mf**2 / self.mY**2
            total += (Nc * self.mY / (12.0 * np.pi)
                      * np.sqrt(1.0 - r2)
                      * gfV**2 * (1.0 + 0.5 * r2))
        return total

    def branching_ratios_to_SM(self) -> dict:
        partials = {}
        total = 0.0
        for f, (mf, B_charge, Nc) in BARYON_FERMIONS.items():
            if self.mY <= 2.0 * mf:
                continue
            gfV, _ = self.gfv_gfa_B(
                f, self.mY, 1.0, self.g_B, self.m_ZB)
            r2 = 4.0 * mf**2 / self.mY**2
            w = (Nc * self.mY / (12.0 * np.pi)
                 * np.sqrt(1.0 - r2)
                 * gfV**2 * (1.0 + 0.5 * r2))
            partials[_PAIR_NAME[f]] = w
            total += w
        if total <= 0:
            return {ch: 0.0 for ch in partials}
        return {ch: w / total for ch, w in partials.items()}

    def sigma_SI(self, epsX: float, A_N: int = _A_Xe, Z_N: int = _Z_Xe) -> float:
        """
        SI DM-nucleon cross section for the Baryon portal.

        Purely vectorial, isospin-universal: all quarks have B = 1/3,
        so proton and neutron couplings are equal and the result is
        independent of (A_N, Z_N).
        """
        gD = np.sqrt(4.0 * np.pi * self.alphaX)
        mixing = abs((self.m_ZB**2 + self.mY**2)
                     / (self.m_ZB**2 - self.mY**2))
        gqV = epsX * self.g_B * (1.0 / 3.0) * mixing
        mu = self.mX * mN / (self.mX + mN)
        amplitude = gD * 3.0 * gqV / self.mY**2
        return mu**2 / np.pi * amplitude**2 * HBARC2


# =====================================================================
#  HiggsPortal
# =====================================================================

@dataclass
class HiggsPortal(DarkSectorModel):
    """
    Secluded dark sector coupled to the SM via Higgs-portal mixing.

    The mediator phi is a real scalar that mixes with the SM Higgs
    through the mixing angle eps (small-angle limit: sin(eps) ~ eps).

    Parameters
    ----------
    mX, mY : float
        Dark-matter and mediator masses (GeV).
    gX, gY : int
        Internal degrees of freedom.
    lam : float or None
        Shorthand coupling: sets lam_s = lam_p = lam.
    lam_s, lam_p : float or None
        Scalar and pseudo-scalar couplings (specify both, or use lam).
    include_antiparticlesX, include_antiparticlesY : bool
        Whether to double-count for Xbar or Ybar.
    Delta1, Delta3 : float
        Multiplicative rescaling of YYY->YY and YYX->YX cannibal rates.
        Set to 0 to turn off a channel; default 1.0 (no rescaling).
    xi_infl : float
        Initial dark-to-SM temperature ratio.
    hh_coupling : str
        Coupling scheme for phi->hh: 'singlet' or 'zero'.
    """
    mX: float = 100.0
    mY: float = 10.0
    gX: int   = 2       # Dirac fermion DM
    gY: int   = 1       # real scalar mediator
    lam: Optional[float] = 1e-2
    lam_s: Optional[float] = None
    lam_p: Optional[float] = None
    include_antiparticlesX: bool = True
    include_antiparticlesY: bool = False
    Delta1: float = 1.0
    Delta3: float = 1.0
    xi_infl: float = 1.0
    hh_coupling: str = 'singlet'

    def __post_init__(self):
        # --- Resolve couplings ---
        if self.lam is not None:
            if self.lam_s is not None or self.lam_p is not None:
                raise ValueError(
                    "Specify either 'lam' (sets lam_s = lam_p = lam) "
                    "or both 'lam_s' and 'lam_p', not both.")
            self.lam_s = self.lam
            self.lam_p = self.lam
        else:
            if self.lam_s is None or self.lam_p is None:
                raise ValueError(
                    "Must specify either 'lam' or both 'lam_s' and 'lam_p'.")

        # --- Ensure HDECAY is ready ---
        from .phi_decay import ensure_hdecay_ready
        ensure_hdecay_ready()

    # ----------------------------------------------------------
    #  DarkSectorModel interface
    # ----------------------------------------------------------

    def sigmav_XX_to_YY(self, TD: float) -> float:
        return self._sigmav_XX_to_YY_full(
            self.mX, self.mY, self.lam_s, self.lam_p, TD)

    def sigmav_XX_to_YY_swave(self) -> float:
        return self._sigmav_XX_to_YY_swave(
            self.mX, self.mY, self.lam_s, self.lam_p)

    def sigmav2_YYY_to_YY(self, TD: float) -> float:
        return self.Delta1 * self._sigmav2_YYY_to_YY(
            self.mX, self.mY, self.lam_s, self.lam_p)

    def sigmav2_YXX_to_XX(self) -> float:
        return 0.0   # subdominant, no simplified form available

    def sigmav2_YYX_to_YX(self) -> float:
        return self.Delta3 * self._sigmav2_YYX_to_YX(
            self.mX, self.mY, self.lam_s, self.lam_p)

    def decay_width_to_SM(self, epsX: float) -> float:
        from .phi_decay import phi_total_width_normalised
        return epsX**2 * phi_total_width_normalised(
            self.mY, hh_coupling=self.hh_coupling)

    # ----------------------------------------------------------
    #  Cross-section formulas
    # ----------------------------------------------------------

    @staticmethod
    def _sigmav_XX_to_YY_full(mX, mY, lam_s, lam_p, TX):
        """<sigma v> for X Xbar -> phi phi  (s-wave + p-wave)."""
        r = mY / mX
        xX = mX / TX
        a = (2.0 * np.sqrt(1.0 - r**2) * lam_p**2 * lam_s**2
             ) / (mX**2 * np.pi * (r**2 - 2.0)**2)
        b = (-2.0 * (r**2 - 1.0)**3 * lam_p**4
             + 3.0 * (r**6 - 8.0*r**4 + 20.0*r**2 - 12.0) * lam_s**2 * lam_p**2
             + 2.0 * (-2.0*r**6 + 10.0*r**4 - 17.0*r**2 + 9.0) * lam_s**4
             ) / (12.0 * mX**2 * np.pi * np.sqrt(1.0 - r**2) * (r**2 - 2.0)**4)
        if xX < 1.0:
            return 2.0 * a
        return 2.0 * a + 6.0 * b / xX

    @staticmethod
    def _sigmav_XX_to_YY_swave(mX, mY, lam_s, lam_p):
        """s-wave-only <sigma v> for X Xbar -> phi phi."""
        r = mY / mX
        return (4.0 * np.sqrt(1.0 - r**2) * lam_p**2 * lam_s**2
                ) / (mX**2 * np.pi * (r**2 - 2.0)**2)

    @staticmethod
    def _sigmav2_YYY_to_YY(mX, mY, lam_s, lam_p):
        """<sigma^2 v^2> for phi phi phi -> phi phi."""
        return ((np.sqrt(5.0) / (3.0 * np.pi**5))
                * lam_s**2
                * (3.0*lam_s**4 + 10.0*lam_s**2*lam_p**2
                   + 15.0*lam_p**4)**2
                / (mX**2 * mY**3))

    @staticmethod
    def _sigmav2_YYX_to_YX(mX, mY, lam_s, lam_p):
        """<sigma^2 v^2> for phi phi X -> phi X."""
        r = mX / mY
        Fac = (9.0 * r**5 * np.sqrt(12.0*r*(r + 2.0) + 9.0) * (
            3.0*lam_p**6 * (1.0 - 2.0*r*(r + 1.0))**2
            + lam_p**4 * lam_s**2 * (2.0*r + 1.0) * (
                2.0*r*(2.0*r*(r*(r*(2.0*r*(16.0*r*(r + 5.0) + 149.0)
                + 271.0) + 142.0) + 40.0) - 5.0) + 9.0)
            + lam_s**4 * (4.0*r*(r*(4.0*r*(2.0*r*(2.0*r + 9.0) + 25.0)
                + 45.0) + 4.0) + 9.0) * (lam_p*(2.0*r + 1.0))**2
            + lam_s**6 * (2.0*r + 1.0)**5 * (2.0*r + 3.0)
            )) / (64.0 * np.pi * (r + 2.0)**3 * (2.0*r + 1.0)**4
                  * (2.0*r**2*(r + 2.0) + r - 1.0)**2)
        return Fac / mX**5

    def branching_ratios_to_SM(self) -> dict:
        from .phi_decay import phi_branching_ratios
        br = phi_branching_ratios(self.mY, hh_coupling=self.hh_coupling)
        # Remove internal 'total' key — return only channel BRs
        br.pop('total', None)
        br.pop('total_sp', None)
        return br

    def sigma_SI(self, epsX: float, A_N: int = _A_Xe, Z_N: int = _Z_Xe,
                 f_N: float = 0.30) -> float:
        """
        SI DM-nucleon cross section for the Higgs portal.

        Majorana fermion chi scattering via t-channel exchange of two
        mixed scalars (h, S) in the small-mixing limit (sin(theta) = epsX).
        From arXiv:1609.02555, Sec. VI.B.

        Parameters
        ----------
        epsX : float
            Higgs-singlet mixing angle (sin(theta)).
        A_N, Z_N : int
            Target nucleus (unused — isospin-universal via Higgs coupling).
        f_N : float
            Higgs-nucleon form factor (default 0.30).
        """
        mu = self.mX * mN / (self.mX + mN)
        propagator = 1.0 / self.mY**2 - 1.0 / Mh**2
        a_N = self.lam_s * epsX * f_N * mN / vH * propagator

        return 4.0 * mu**2 / np.pi * a_N**2 * HBARC2

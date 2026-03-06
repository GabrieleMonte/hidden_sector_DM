"""
Dark-sector model interface and concrete implementations.

DarkSectorModel : abstract base class defining the rate interface
VectorPortal    : secluded dark sector with kinetic mixing (hypercharge)
BLPortal        : secluded dark sector with kinetic mixing (B-L)
HiggsPortal     : secluded dark sector with scalar mixing
"""

import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

from .constants import SM_MYQ, BL_FERMIONS, MZ, gW, sW, cW, gY as gY_SM


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
                * (gVf**2 * (1.0 - r2) + gAf**2 * (1.0 + 2.0 * mf**2 / mA**2)))

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
    g_BL: float = 0.1
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

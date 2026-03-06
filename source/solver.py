"""
BoltzmannSolver — model-independent Boltzmann equation machinery.

All dark-sector rates are read from a DarkSectorModel instance.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Dict, Any
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy import special as scisp
from tqdm import tqdm

from .constants import Mpl, rhoc, s0_cosmo, Och2, mY_relic
from .cosmology import Cosmology
from .model import DarkSectorModel


class BoltzmannSolver:
    """
    Model-independent Boltzmann solver for secluded dark sectors.

    Parameters
    ----------
    cosmo : Cosmology
        Shared SM-background instance.
    model : DarkSectorModel
        Dark-sector model providing masses, DOF, and rate methods.
    """

    def __init__(self, cosmo: Cosmology, model: DarkSectorModel):
        self.cosmo = cosmo
        self.model = model

    # =============================================================
    #  Boltzmann solver — chemical-potential formulation (2-phase)
    # =============================================================
    def solve_boltzmann_chempot(
        self,
        xmin: float = 1e-3,
        n_points: int = 500,
        rtol_value: float = 1e-8,
        atol_value: float = 1e-15,
        Gamma_switch_threshold: float = 1e7,
        convergence_threshold: float = 1e-2,
        convergence_mode: str = "dlnxi_NR",
        Tinf: float = 1e14,
        return_bg_ICs: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Two-phase Boltzmann solver using chemical-potential variables
        (ln xi, mubar_X, mubar_Y) with Y_i = Y_i^eq exp(mubar_i).

        Phase 1 (equilibrium): mubar = 0, evolve only ln xi
        Phase 2 (full):        evolve all three via decoupled A,B,C,D,E,F,G
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        include_antiparticlesX = model.include_antiparticlesX
        include_antiparticlesY = model.include_antiparticlesY

        # x-grid
        xfList = np.concatenate([
            np.array([1e-2, 1e-1, 0.99]),
            np.arange(1, 100, 1),
            np.array([2e2, 5e2, 1e3, 2e3, 5e3, 1e4]),
        ])

        # Temperature-independent cross sections
        sv_YXX_XX = model.sigmav2_YXX_to_XX()
        sv_YYX_YX = model.sigmav2_YYX_to_YX()

        r = mX / mY
        state = {"x_FO": None, "x_switch": None, "x_converged": None}

        # --- local helpers ---

        def _compute_thermo(x, lnxi, muX, muY):
            xi = np.exp(np.clip(lnxi, -20, 20))
            T  = mX / x
            TD = xi * T
            s  = cosmo.s_entropy(T)
            H  = cosmo.hubble(T)
            dlngSdlnT_val = cosmo.dlngSdlnT(T)
            g_tilde = 1.0 + dlngSdlnT_val / 3.0

            zX = mX / TD;  zY = mY / TD

            BX1 = cosmo.Bfac1(zX, mX);  BX2 = cosmo.Bfac2(zX, mX);  dBX1 = cosmo.dBfac1dT(zX, mX)
            BY1 = cosmo.Bfac1(zY, mY);  BY2 = cosmo.Bfac2(zY, mY);  dBY1 = cosmo.dBfac1dT(zY, mY)

            lamX = cosmo.lambda_i(zX);  lamY = cosmo.lambda_i(zY)

            lnYeqX = cosmo.ln_Yeq(zX, gX, mX, include_antiparticlesX, s)
            lnYeqY = cosmo.ln_Yeq(zY, gY, mY, include_antiparticlesY, s)
            lnYX = lnYeqX + muX;  lnYY = lnYeqY + muY

            YX = np.exp(np.clip(lnYX, -700, 700))
            YY = np.exp(np.clip(lnYY, -700, 700))
            nX = YX * s;  nY = YY * s

            rho_h = max(BX1 * nX + BY1 * nY, 0.0)
            Htot  = np.sqrt(H**2 + rho_h / (3.0 * Mpl**2))

            return dict(
                xi=xi, T=T, TD=TD, s=s, H=H, Htot=Htot,
                dlngSdlnT=dlngSdlnT_val, g_tilde=g_tilde, zX=zX, zY=zY,
                BX1=BX1, BX2=BX2, dBX1=dBX1,
                BY1=BY1, BY2=BY2, dBY1=dBY1,
                lnYeqX=lnYeqX, lnYeqY=lnYeqY, lnRXY=lnYeqX - lnYeqY,
                YX=YX, YY=YY, lnYX=lnYX, lnYY=lnYY, nX=nX, nY=nY,
                lamX=lamX, lamY=lamY,
                sv_YYY_YY=model.sigmav2_YYY_to_YY(TD),
                sv_XXYY=model.sigmav_XX_to_YY(TD),
            )

        # --- Phase 1: equilibrium ---

        def dlnxi_dx_equilibrium(x, lnxi):
            th = _compute_thermo(x, lnxi, 0.0, 0.0)
            neqX, neqY = th['nX'], th['nY']

            Num = neqY * th['BY2'] + neqX * th['BX2']
            dneqX = cosmo.dneqdT_MB(th['zX'], gX, mX, include_antiparticlesX)
            dneqY = cosmo.dneqdT_MB(th['zY'], gY, mY, include_antiparticlesY)
            Den = (th['BX1'] * dneqX + th['dBX1'] * neqX
                   + th['BY1'] * dneqY + th['dBY1'] * neqY)
            if np.abs(Den) < 1e-300:
                Den = np.copysign(1e-300, Den)

            dlnxi = (1.0 / x) * (1.0 - (3.0 + th['dlngSdlnT']) * Num / (th['TD'] * Den))
            Gamma_ann = (th['s'] / th['Htot']
                         * (th['sv_XXYY'] / 2.0) * th['YX'])
            return dlnxi, Gamma_ann

        def BEQs_eq(t, y):
            dlnxi, _ = dlnxi_dx_equilibrium(t, y[0])
            return [np.clip(dlnxi, -10.0 / t, 10.0 / t)]

        # --- Phase 2: full system ---

        def compute_derivatives_full(x, lnxi, muX, muY):
            th = _compute_thermo(x, lnxi, muX, muY)
            g_tilde = th['g_tilde']
            YX, YY  = th['YX'], th['YY']
            nX, nY  = th['nX'], th['nY']

            Pfac = th['s'] * g_tilde / (th['Htot'] * x) * (th['sv_XXYY'] / 2.0)
            RXY2_YY2 = np.exp(np.clip(2.0 * th['lnRXY'] + 2.0 * th['lnYY'], -700, 700))

            coll_X = Pfac * (YX - RXY2_YY2 / YX) if YX > 1e-300 else 0.0
            A = -coll_X + (th['lamX'] - 3.0 * g_tilde) / x
            B = -th['lamX']

            coll_Y_2to2 = Pfac * (YX**2 - RXY2_YY2) / YY if YY > 1e-300 else 0.0

            S_3to2 = (YY**2 * th['sv_YYY_YY']
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            Gamma_can = th['s']**2 * g_tilde / (th['Htot'] * x) * S_3to2

            if muY > 50:
                one_m_expnmuY = 1.0
            elif muY < -50:
                one_m_expnmuY = -np.exp(-muY)
            else:
                one_m_expnmuY = 1.0 - np.exp(-muY)

            C = coll_Y_2to2 - Gamma_can * one_m_expnmuY + (th['lamY'] - 3.0 * g_tilde) / x
            D = -th['lamY']

            D_cal  = nX * th['dBX1'] + nY * th['dBY1']
            E_cal  = th['BX1'] * nX * th['lamX'] + th['BY1'] * nY * th['lamY']
            Dtilde = th['TD'] * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            rho_plus_P = nX * th['BX2'] + nY * th['BY2']
            E = (1.0 / x) * (1.0 - 3.0 * g_tilde * rho_plus_P / Dtilde)
            F = -th['BX1'] * nX / Dtilde
            G = -th['BY1'] * nY / Dtilde

            denom = 1.0 - F * B - G * D
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)

            dlnxi_dx = np.clip((E + F * A + G * C) / denom, -10.0 / x, 10.0 / x)
            dmuX_dx  = A + B * dlnxi_dx
            dmuY_dx  = C + D * dlnxi_dx

            if muX > 0.05 and state["x_FO"] is None:
                state["x_FO"] = x

            return dlnxi_dx, dmuX_dx, dmuY_dx

        def BEQs_full(t, y):
            return compute_derivatives_full(t, *y)

        # --- Initial conditions ---

        x0    = xmin
        lnxi0 = (1.0 / 3.0) * np.log(
            np.round(cosmo.gstarS(mX / xmin) / cosmo.gstarS(Tinf), 4)) + np.log(model.xi_infl)
        lnYX_check_prev = cosmo.ln_Yeq(
            mX / (mX / x0), gX, mX, include_antiparticlesX, cosmo.s_entropy(mX / x0))
        x_check_prev = x0

        x_all, lnxi_all = [x0], [lnxi0]
        muX_all, muY_all = [0.0], [0.0]
        in_equilibrium = True
        y0_eq   = [0.0]
        y0_full = None

        if verbose:
            print("=" * 60)
            print("Boltzmann Solver -- Chemical Potential Formulation")
            print("=" * 60)
            print(f"  mX = {mX:.2e} GeV, mY = {mY:.2e} GeV, r = {r:.2f}")
            print(f"  Phase 1->2 switch: Gamma_ann < {Gamma_switch_threshold:.0e}")
            print(f"  Convergence mode: {convergence_mode}")
            print(f"  Convergence threshold: {convergence_threshold}")
            print("-" * 60)

        # --- Main loop ---

        converged = False
        iterator = tqdm(range(len(xfList)), disable=not verbose, desc="Evolving")

        for j in iterator:
            xf = xfList[j]
            if x0 >= xf:
                continue

            if x0 < 1:
                xs = np.logspace(np.log10(x0 * 1.001), np.log10(xf * 0.999), n_points)
            else:
                xs = np.linspace(x0 * 1.001, xf * 0.999, n_points)

            if in_equilibrium:
                sol = solve_ivp(BEQs_eq, (x0, xf), y0_eq, t_eval=xs,
                                rtol=rtol_value, atol=1e-12, method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (eq): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break

                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend([0.0] * len(sol.t))
                muY_all.extend([0.0] * len(sol.t))
                x0 = sol.t[-1];  y0_eq = [sol.y[0, -1]]

                _, Gamma_now = dlnxi_dx_equilibrium(x0, y0_eq[0])
                if Gamma_now < Gamma_switch_threshold:
                    in_equilibrium = False
                    state["x_switch"] = x0
                    y0_full = [y0_eq[0], 0.0, 0.0]
                    if verbose:
                        print(f"\n  -> Full system at x = {x0:.2f}"
                              f" (Gamma_ann = {Gamma_now:.1e})")

            else:
                sol = solve_ivp(BEQs_full, (x0, xf), y0_full, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (full): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break

                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend(sol.y[1].tolist())
                muY_all.extend(sol.y[2].tolist())
                x0 = sol.t[-1];  y0_full = sol.y[:, -1].tolist()

                # convergence check
                if state["x_FO"] is not None and x0 > 2.0 * state["x_FO"]:

                    if convergence_mode == 'dlnYX':
                        xi_now  = np.exp(np.clip(sol.y[0, -1], -20, 20))
                        muX_now = sol.y[1, -1]
                        T_now   = mX / x0
                        zX_now  = mX / (xi_now * T_now)
                        lnYX_now = muX_now + cosmo.ln_Yeq(
                            zX_now, gX, mX, include_antiparticlesX,
                            cosmo.s_entropy(T_now))

                        if (np.isfinite(lnYX_check_prev) and np.isfinite(lnYX_now)
                                and x0 > x_check_prev * 1.5):
                            dlnYX_dlnx = ((lnYX_now - lnYX_check_prev)
                                          / np.log(x0 / x_check_prev))
                            if verbose:
                                iterator.set_postfix({
                                    'x': f'{x0:.0f}', 'muX': f'{muX_now:.1f}',
                                    '|dlnYX/dlnx|': f'{np.abs(dlnYX_dlnx):.1e}'})
                            if np.abs(dlnYX_dlnx) < convergence_threshold:
                                converged = True
                                state["x_converged"] = x0
                                if verbose:
                                    print(f"\n  Converged (dlnYX) at x = {x0:.1f}")
                                break
                            lnYX_check_prev = lnYX_now
                            x_check_prev = x0

                    elif convergence_mode == 'dlnxi_NR':
                        lnxi_now = sol.y[0, -1]
                        muX_now  = sol.y[1, -1]
                        muY_now  = sol.y[2, -1]
                        dlnxi_num, _, _ = compute_derivatives_full(
                            x0, lnxi_now, muX_now, muY_now)

                        g_tilde_now = _compute_thermo(
                            x0, lnxi_now, muX_now, muY_now)['g_tilde']
                        dlnxi_NR = (1.0 - 2.0 * g_tilde_now) / x0
                        rel_dev = np.abs(dlnxi_num - dlnxi_NR) * x0

                        if verbose:
                            iterator.set_postfix({
                                'x': f'{x0:.0f}', 'muX': f'{muX_now:.1f}',
                                'dlnxi_dev': f'{rel_dev:.1e}'})
                        if rel_dev < convergence_threshold:
                            converged = True
                            state["x_converged"] = x0
                            if verbose:
                                print(f"\n  Converged (dlnxi_NR) at x = {x0:.1f}")
                            break

        # --- Post-process ---

        x_arr    = np.array(x_all)
        lnxi_arr = np.array(lnxi_all)
        muX_arr  = np.array(muX_all)
        muY_arr  = np.array(muY_all)

        xi_arr = np.exp(np.clip(lnxi_arr, -20, 20))
        T_arr  = mX / x_arr
        TD_arr = xi_arr * T_arr
        zX_arr = mX / TD_arr;  zY_arr = mY / TD_arr

        lnYX_arr   = np.zeros_like(x_arr)
        lnYY_arr   = np.zeros_like(x_arr)
        lnYXeq_arr = np.zeros_like(x_arr)
        lnYYeq_arr = np.zeros_like(x_arr)
        for i in range(len(x_arr)):
            s_i = cosmo.s_entropy(T_arr[i])
            lnYXeq_arr[i] = cosmo.ln_Yeq(zX_arr[i], gX, mX, include_antiparticlesX, s_i)
            lnYYeq_arr[i] = cosmo.ln_Yeq(zY_arr[i], gY, mY, include_antiparticlesY, s_i)
            lnYX_arr[i] = muX_arr[i] + lnYXeq_arr[i]
            lnYY_arr[i] = muY_arr[i] + lnYYeq_arr[i]

        YX_arr   = np.exp(np.clip(lnYX_arr,   -700, 700))
        YY_arr   = np.exp(np.clip(lnYY_arr,   -700, 700))
        YXeq_arr = np.exp(np.clip(lnYXeq_arr, -700, 700))
        YYeq_arr = np.exp(np.clip(lnYYeq_arr, -700, 700))

        if verbose:
            print("-" * 60)
            print("RESULTS:")
            if state['x_switch']:
                print(f"  x_switch (eq -> full): {state['x_switch']:.2f}")
            if state['x_FO']:
                print(f"  x_FO (muX > 0.05): {state['x_FO']:.1f}")
            print(f"  Converged: {converged}"
                  + (f" at x = {state['x_converged']:.1f}" if converged else ""))
            print(f"  YX_relic = {YX_arr[-1]:.6e}")
            print("=" * 60)

        result = {
            'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
            'YXeq': YXeq_arr, 'YYeq': YYeq_arr,
            'mubar_X': muX_arr, 'mubar_Y': muY_arr,
            'YX_relic': YX_arr[-1], 'YY_final': YY_arr[-1], 'xi_final': xi_arr[-1],
            'x_FO': state['x_FO'], 'x_switch': state['x_switch'],
            'x_converged': state['x_converged'], 'converged': converged,
        }
        if return_bg_ICs:
            result['bg_ICs'] = {
                'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
                'T': T_arr, 'TD': TD_arr,
            }
        return result

    # =============================================================
    #  Boltzmann solver — QSSA three-phase
    # =============================================================
    def solve_boltzmann_chempot_QSSA(
        self,
        xmin: float = 1e-3,
        xmax: float = None,
        n_points: int = 500,
        rtol_value: float = 1e-8,
        atol_value: float = 1e-15,
        Gamma_switch_QSSA: float = 5e9,
        Gamma_switch_full: float = 1e4,
        cannibal_switch_full: float = 1,
        convergence_threshold: float = 1e-2,
        convergence_mode: str = "dlnxi_NR",
        Tinf: float = 1e14,
        return_bg_ICs: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Three-phase QSSA Boltzmann solver using chemical-potential variables.

        Phase 1   : equilibrium (muX=muY=0, evolve ln xi only)
        Phase 1.5 : QSSA / Y-in-eq (muY=0, evolve ln xi & muX)
        Phase 2   : full system (evolve ln xi, muX, muY)
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        include_antiparticlesX = model.include_antiparticlesX
        include_antiparticlesY = model.include_antiparticlesY

        r = mX / mY
        xmax_recommended = r**2 * 1e4
        if xmax is None:
            xmax = max(1e5, xmax_recommended)
        elif xmax < xmax_recommended:
            import warnings
            warnings.warn(
                f"xmax={xmax:.0e} may be too small for r={r:.1f}; "
                f"recommend xmax >= r^2 * 1e4 = {xmax_recommended:.0e}")

        # Build x-grid up to xmax
        xf_tail = [x for x in [2e2, 5e2, 1e3, 2e3, 5e3, 1e4] if x < xmax]
        if xmax > 1e4:
            xf_tail.extend(np.logspace(np.log10(max(1e4, xf_tail[-1] if xf_tail else 1e4) * 2),
                                        np.log10(xmax), 10).tolist())
        xf_tail.append(xmax)
        xfList = np.concatenate([
            np.array([1e-2, 1e-1, 0.99]),
            np.linspace(1, 10, 20),
            np.arange(10.5, 100, 1),
            np.array(xf_tail),
        ])
        xfList = np.unique(xfList)

        sv_YXX_XX = model.sigmav2_YXX_to_XX()
        sv_YYX_YX = model.sigmav2_YYX_to_YX()

        state = {"x_FO": None, "x_switch_QSSA": None,
                 "x_switch_full": None, "x_converged": None}

        # --- Thermodynamics: full (Phase 1 & 2) ---

        def _thermo_full(x, lnxi, muX, muY):
            xi = np.exp(np.clip(lnxi, -20, 20))
            T  = mX / x;  TD = xi * T
            s  = cosmo.s_entropy(T);  H = cosmo.hubble(T)
            dlngSdlnT_val = cosmo.dlngSdlnT(T)
            g_tilde = 1.0 + dlngSdlnT_val / 3.0
            zX = mX / TD;  zY = mY / TD

            BX1 = cosmo.Bfac1(zX, mX);  BX2 = cosmo.Bfac2(zX, mX)
            dBX1 = cosmo.dBfac1dT(zX, mX)
            BY1 = cosmo.Bfac1(zY, mY);  BY2 = cosmo.Bfac2(zY, mY)
            dBY1 = cosmo.dBfac1dT(zY, mY)
            lamX = cosmo.lambda_i(zX);  lamY = cosmo.lambda_i(zY)

            lnYeqX = cosmo.ln_Yeq(zX, gX, mX, include_antiparticlesX, s)
            lnYeqY = cosmo.ln_Yeq(zY, gY, mY, include_antiparticlesY, s)
            lnYX = lnYeqX + muX;  lnYY = lnYeqY + muY
            YX = np.exp(np.clip(lnYX, -700, 700))
            YY = np.exp(np.clip(lnYY, -700, 700))
            nX = YX * s;  nY = YY * s

            rho_h = max(BX1 * nX + BY1 * nY, 0.0)
            Htot  = np.sqrt(H**2 + rho_h / (3.0 * Mpl**2))

            return dict(
                xi=xi, T=T, TD=TD, s=s, H=H, Htot=Htot,
                dlngSdlnT=dlngSdlnT_val, g_tilde=g_tilde, zX=zX, zY=zY,
                BX1=BX1, BX2=BX2, dBX1=dBX1,
                BY1=BY1, BY2=BY2, dBY1=dBY1,
                lnYeqX=lnYeqX, lnYeqY=lnYeqY, lnRXY=lnYeqX - lnYeqY,
                YX=YX, YY=YY, lnYX=lnYX, lnYY=lnYY, nX=nX, nY=nY,
                lamX=lamX, lamY=lamY,
                sv_YYY_YY=model.sigmav2_YYY_to_YY(TD),
                sv_XXYY=model.sigmav_XX_to_YY(TD),
            )

        # --- Thermodynamics: Y-in-eq (Phase 1.5) ---

        def _thermo_Yeq(x, lnxi, muX):
            xi = np.exp(np.clip(lnxi, -20, 20))
            T  = mX / x;  TD = xi * T
            s  = cosmo.s_entropy(T);  H = cosmo.hubble(T)
            dlngSdlnT_val = cosmo.dlngSdlnT(T)
            g_tilde = 1.0 + dlngSdlnT_val / 3.0
            zX = mX / TD;  zY = mY / TD

            BX1 = cosmo.Bfac1(zX, mX);  BX2 = cosmo.Bfac2(zX, mX)
            dBX1 = cosmo.dBfac1dT(zX, mX)
            BY1 = cosmo.Bfac1(zY, mY);  BY2 = cosmo.Bfac2(zY, mY)
            dBY1 = cosmo.dBfac1dT(zY, mY)
            lamX = cosmo.lambda_i(zX);  lamY = cosmo.lambda_i(zY)

            lnYeqX = cosmo.ln_Yeq(zX, gX, mX, include_antiparticlesX, s)
            lnYeqY = cosmo.ln_Yeq(zY, gY, mY, include_antiparticlesY, s)
            lnYX = lnYeqX + muX;  lnYY = lnYeqY
            YX = np.exp(np.clip(lnYX, -700, 700))
            YY = np.exp(np.clip(lnYY, -700, 700))
            nX = YX * s;  nY = YY * s

            rho_h = max(BX1 * nX + BY1 * nY, 0.0)
            Htot  = np.sqrt(H**2 + rho_h / (3.0 * Mpl**2))

            return dict(
                xi=xi, T=T, TD=TD, s=s, H=H, Htot=Htot,
                dlngSdlnT=dlngSdlnT_val, g_tilde=g_tilde, zX=zX, zY=zY,
                BX1=BX1, BX2=BX2, dBX1=dBX1,
                BY1=BY1, BY2=BY2, dBY1=dBY1,
                lnYeqX=lnYeqX, lnYeqY=lnYeqY,
                YX=YX, YY=YY, nX=nX, nY=nY,
                lamX=lamX, lamY=lamY,
            )

        # === Phase 1: Equilibrium ===

        def dlnxi_dx_equilibrium(x, lnxi):
            th = _thermo_full(x, lnxi, 0.0, 0.0)
            neqX, neqY = th['nX'], th['nY']
            Num = neqY * th['BY2'] + neqX * th['BX2']
            dneqX = cosmo.dneqdT_MB(th['zX'], gX, mX, include_antiparticlesX)
            dneqY = cosmo.dneqdT_MB(th['zY'], gY, mY, include_antiparticlesY)
            Den = (th['BX1'] * dneqX + th['dBX1'] * neqX
                   + th['BY1'] * dneqY + th['dBY1'] * neqY)
            if np.abs(Den) < 1e-300:
                Den = np.copysign(1e-300, Den)
            dlnxi = (1.0 / x) * (1.0 - (3.0 + th['dlngSdlnT']) * Num / (th['TD'] * Den))
            Gamma_ann = (th['s'] / th['Htot']
                         * (th['sv_XXYY'] / 2.0) * th['YX'])
            return dlnxi, Gamma_ann

        def BEQs_eq(t, y):
            dlnxi, _ = dlnxi_dx_equilibrium(t, y[0])
            return [np.clip(dlnxi, -10.0 / t, 10.0 / t)]

        # === Phase 1.5: QSSA (muY=0, evolve lnxi & muX) ===

        def compute_derivatives_QSSA(x, lnxi, muX):
            th = _thermo_Yeq(x, lnxi, muX)
            g_tilde = th['g_tilde']
            nX, nY = th['nX'], th['nY']

            YXeq = np.exp(np.clip(th['lnYeqX'], -700, 700))
            sv = model.sigmav_XX_to_YY(th['TD'])
            Gamma_coll = th['s'] * g_tilde / (th['Htot'] * x) * sv * YXeq

            A = -Gamma_coll * np.sinh(muX) + (th['lamX'] - 3.0 * g_tilde) / x
            B = -th['lamX']

            D_cal  = nX * th['dBX1'] + nY * th['dBY1']
            E_cal  = th['BX1'] * nX * th['lamX'] + th['BY1'] * nY * th['lamY']
            Dtilde = th['TD'] * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            Num_xi = nX * th['BX2'] + nY * th['BY2']
            C_xi = (1.0 / x) * (1.0 - 3.0 * g_tilde * Num_xi / Dtilde)
            D_xi = -th['BX1'] * nX / Dtilde

            denom = 1.0 - B * D_xi
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)

            dlnxi_dx = np.clip((C_xi + D_xi * A) / denom, -10.0 / x, 10.0 / x)
            dmuX_dx  = (A + B * C_xi) / denom

            if muX > 0.05 and state["x_FO"] is None:
                state["x_FO"] = x

            return dlnxi_dx, dmuX_dx

        def BEQs_QSSA(t, y):
            return compute_derivatives_QSSA(t, *y)

        def estimate_cannibal_rate(x, lnxi, muX):
            th = _thermo_Yeq(x, lnxi, muX)
            YX, YY = th['YX'], th['YY']
            sv_YYY = model.sigmav2_YYY_to_YY(th['TD'])
            S_3to2 = (YY**2 * sv_YYY
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            return th['s']**2 / th['Htot'] * S_3to2

        # === Phase 2: Full system ===

        def compute_derivatives_full(x, lnxi, muX, muY):
            th = _thermo_full(x, lnxi, muX, muY)
            g_tilde = th['g_tilde']
            YX, YY  = th['YX'], th['YY']
            nX, nY  = th['nX'], th['nY']

            Pfac = th['s'] * g_tilde / (th['Htot'] * x) * (th['sv_XXYY'] / 2.0)
            RXY2_YY2 = np.exp(np.clip(2.0 * th['lnRXY'] + 2.0 * th['lnYY'], -700, 700))

            coll_X = Pfac * (YX - RXY2_YY2 / YX) if YX > 1e-300 else 0.0
            A = -coll_X + (th['lamX'] - 3.0 * g_tilde) / x
            B = -th['lamX']

            coll_Y_2to2 = Pfac * (YX**2 - RXY2_YY2) / YY if YY > 1e-300 else 0.0
            S_3to2 = (YY**2 * th['sv_YYY_YY']
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            Gamma_can = th['s']**2 * g_tilde / (th['Htot'] * x) * S_3to2

            if muY > 50:
                one_m = 1.0
            elif muY < -50:
                one_m = -np.exp(-muY)
            else:
                one_m = 1.0 - np.exp(-muY)

            C = coll_Y_2to2 - Gamma_can * one_m + (th['lamY'] - 3.0 * g_tilde) / x
            D = -th['lamY']

            D_cal  = nX * th['dBX1'] + nY * th['dBY1']
            E_cal  = th['BX1'] * nX * th['lamX'] + th['BY1'] * nY * th['lamY']
            Dtilde = th['TD'] * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            rho_plus_P = nX * th['BX2'] + nY * th['BY2']
            E = (1.0 / x) * (1.0 - 3.0 * g_tilde * rho_plus_P / Dtilde)
            F = -th['BX1'] * nX / Dtilde
            G = -th['BY1'] * nY / Dtilde

            denom = 1.0 - F * B - G * D
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)

            dlnxi_dx = np.clip((E + F * A + G * C) / denom, -10.0 / x, 10.0 / x)
            dmuX_dx  = A + B * dlnxi_dx
            dmuY_dx  = C + D * dlnxi_dx

            if muX > 0.05 and state["x_FO"] is None:
                state["x_FO"] = x

            return dlnxi_dx, dmuX_dx, dmuY_dx

        def BEQs_full(t, y):
            return compute_derivatives_full(t, *y)

        # === Convergence check helper ===

        def _check_convergence(x0, lnxi_now, muX_now, g_tilde_now,
                               phase_label, dlnxi_num_val):
            nonlocal converged, lnYX_check_prev, x_check_prev

            if state["x_FO"] is None or x0 < 2.0 * state["x_FO"]:
                return False

            if convergence_mode == 'dlnxi_NR':
                dlnxi_NR = (1.0 - 2.0 * g_tilde_now) / x0
                rel_dev = np.abs(dlnxi_num_val - dlnxi_NR) * x0
                if verbose:
                    iterator.set_postfix({
                        'ph': phase_label, 'x': f'{x0:.0f}',
                        'muX': f'{muX_now:.1f}', 'dev': f'{rel_dev:.1e}'})
                if rel_dev < convergence_threshold:
                    converged = True
                    state["x_converged"] = x0
                    if verbose:
                        print(f"\n  Converged ({phase_label}, dlnxi_NR)"
                              f" at x = {x0:.1f}")
                    return True

            elif convergence_mode == 'dlnYX':
                xi_now = np.exp(np.clip(lnxi_now, -20, 20))
                T_now  = mX / x0
                zX_now = mX / (xi_now * T_now)
                lnYX_now = muX_now + cosmo.ln_Yeq(
                    zX_now, gX, mX, include_antiparticlesX,
                    cosmo.s_entropy(T_now))
                if (np.isfinite(lnYX_check_prev) and np.isfinite(lnYX_now)
                        and x0 > x_check_prev * 1.5):
                    dlnYX_dlnx = ((lnYX_now - lnYX_check_prev)
                                  / np.log(x0 / x_check_prev))
                    if verbose:
                        iterator.set_postfix({
                            'ph': phase_label, 'x': f'{x0:.0f}',
                            'muX': f'{muX_now:.1f}',
                            '|dlnYX|': f'{np.abs(dlnYX_dlnx):.1e}'})
                    if np.abs(dlnYX_dlnx) < convergence_threshold:
                        converged = True
                        state["x_converged"] = x0
                        if verbose:
                            print(f"\n  Converged ({phase_label}, dlnYX)"
                                  f" at x = {x0:.1f}")
                        return True
                    lnYX_check_prev = lnYX_now
                    x_check_prev = x0
            return False

        # === Initial conditions ===

        x0    = xmin
        lnxi0 = (1.0 / 3.0) * np.log(
            np.round(cosmo.gstarS(mX / xmin) / cosmo.gstarS(Tinf), 4)) + np.log(model.xi_infl)
        lnYX_check_prev = cosmo.ln_Yeq(
            mX / (mX / x0), gX, mX, include_antiparticlesX,
            cosmo.s_entropy(mX / x0))
        x_check_prev = x0

        x_all, lnxi_all = [x0], [lnxi0]
        muX_all, muY_all = [0.0], [0.0]

        phase = 1
        y0_eq   = [lnxi0]
        y0_QSSA = None
        y0_full = None

        if verbose:
            print("=" * 60)
            print("Boltzmann Solver -- QSSA Three-Phase")
            print("=" * 60)
            print(f"  mX = {mX:.2e} GeV, mY = {mY:.2e} GeV, r = {r:.2f}")
            print(f"  Phase 1->1.5: Gamma_ann < {Gamma_switch_QSSA:.0e}")
            print(f"  Phase 1.5->2: Gamma_can < {cannibal_switch_full:.0e}")
            print(f"  Convergence: {convergence_mode}, thr={convergence_threshold}")
            print("-" * 60)

        # === Main loop ===

        converged = False
        iterator = tqdm(range(len(xfList)), disable=not verbose, desc="Evolving")

        for j in iterator:
            xf = xfList[j]
            if x0 >= xf:
                continue

            if x0 < 1:
                xs = np.logspace(np.log10(x0 * 1.001), np.log10(xf * 0.999), n_points)
            else:
                xs = np.linspace(x0 * 1.001, xf * 0.999, n_points)

            # ----- PHASE 1 -----
            if phase == 1:
                sol = solve_ivp(BEQs_eq, (x0, xf), y0_eq, t_eval=xs,
                                rtol=rtol_value, atol=1e-12, method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (eq): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break
                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend([0.0] * len(sol.t))
                muY_all.extend([0.0] * len(sol.t))
                x0 = sol.t[-1];  y0_eq = [sol.y[0, -1]]

                _, Gamma_now = dlnxi_dx_equilibrium(x0, y0_eq[0])
                if Gamma_now < Gamma_switch_QSSA:
                    phase = 1.5
                    state["x_switch_QSSA"] = x0
                    y0_QSSA = [y0_eq[0], 0.0]
                    if verbose:
                        print(f"\n  -> Phase 1.5 (QSSA) at x = {x0:.2f}"
                              f" (Gamma_ann = {Gamma_now:.1e})")

            # ----- PHASE 1.5 -----
            elif phase == 1.5:
                sol = solve_ivp(BEQs_QSSA, (x0, xf), y0_QSSA, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (QSSA): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break
                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend(sol.y[1].tolist())
                muY_all.extend([0.0] * len(sol.t))
                x0 = sol.t[-1];  y0_QSSA = sol.y[:, -1].tolist()

                # Check Phase 2 trigger
                Gamma_can_now = estimate_cannibal_rate(x0, y0_QSSA[0], y0_QSSA[1])

                if Gamma_can_now < cannibal_switch_full:
                    phase = 2
                    state["x_switch_full"] = x0
                    y0_full = [y0_QSSA[0], y0_QSSA[1], 0.0]
                    if verbose:
                        print(f"\n  -> Phase 2 (full) at x = {x0:.2f}"
                              f" (Gamma_can={Gamma_can_now:.1e},"
                              f" muX={y0_QSSA[1]:.2f})")

                # Convergence check in QSSA phase
                if convergence_mode == 'dlnxi_NR':
                    dlnxi_num, _ = compute_derivatives_QSSA(
                        x0, y0_QSSA[0], y0_QSSA[1])
                    g_tilde_now = _thermo_Yeq(
                        x0, y0_QSSA[0], y0_QSSA[1])['g_tilde']
                    done = _check_convergence(
                        x0, y0_QSSA[0], y0_QSSA[1], g_tilde_now,
                        '1.5', dlnxi_num)
                else:
                    done = _check_convergence(
                        x0, y0_QSSA[0], y0_QSSA[1], None, '1.5', None)
                if done:
                    break

            # ----- PHASE 2 -----
            elif phase == 2:
                sol = solve_ivp(BEQs_full, (x0, xf), y0_full, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (full): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break
                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend(sol.y[1].tolist())
                muY_all.extend(sol.y[2].tolist())
                x0 = sol.t[-1];  y0_full = sol.y[:, -1].tolist()

                if convergence_mode == 'dlnxi_NR':
                    lnxi_now = sol.y[0, -1]
                    muX_now  = sol.y[1, -1]
                    muY_now  = sol.y[2, -1]
                    dlnxi_num, _, _ = compute_derivatives_full(
                        x0, lnxi_now, muX_now, muY_now)
                    g_tilde_now = _thermo_full(
                        x0, lnxi_now, muX_now, muY_now)['g_tilde']
                    done = _check_convergence(
                        x0, lnxi_now, muX_now, g_tilde_now,
                        '2', dlnxi_num)
                else:
                    muX_now = sol.y[1, -1]
                    done = _check_convergence(
                        x0, sol.y[0, -1], muX_now, None, '2', None)
                if done:
                    break

        # === Post-process ===

        x_arr    = np.array(x_all)
        lnxi_arr = np.array(lnxi_all)
        muX_arr  = np.array(muX_all)
        muY_arr  = np.array(muY_all)

        xi_arr = np.exp(np.clip(lnxi_arr, -20, 20))
        T_arr  = mX / x_arr
        TD_arr = xi_arr * T_arr
        zX_arr = mX / TD_arr;  zY_arr = mY / TD_arr

        lnYX_arr   = np.zeros_like(x_arr)
        lnYY_arr   = np.zeros_like(x_arr)
        lnYXeq_arr = np.zeros_like(x_arr)
        lnYYeq_arr = np.zeros_like(x_arr)
        for i in range(len(x_arr)):
            s_i = cosmo.s_entropy(T_arr[i])
            lnYXeq_arr[i] = cosmo.ln_Yeq(zX_arr[i], gX, mX, include_antiparticlesX, s_i)
            lnYYeq_arr[i] = cosmo.ln_Yeq(zY_arr[i], gY, mY, include_antiparticlesY, s_i)
            lnYX_arr[i] = muX_arr[i] + lnYXeq_arr[i]
            lnYY_arr[i] = muY_arr[i] + lnYYeq_arr[i]

        YX_arr   = np.exp(np.clip(lnYX_arr,   -700, 700))
        YY_arr   = np.exp(np.clip(lnYY_arr,   -700, 700))
        YXeq_arr = np.exp(np.clip(lnYXeq_arr, -700, 700))
        YYeq_arr = np.exp(np.clip(lnYYeq_arr, -700, 700))

        if verbose:
            print("-" * 60)
            print("RESULTS:")
            if state['x_switch_QSSA']:
                print(f"  x_switch (eq -> QSSA): {state['x_switch_QSSA']:.2f}")
            if state['x_switch_full']:
                print(f"  x_switch (QSSA -> full): {state['x_switch_full']:.2f}")
            else:
                print(f"  (QSSA carried to convergence, Phase 2 never entered)")
            if state['x_FO']:
                print(f"  x_FO (muX > 0.05): {state['x_FO']:.1f}")
            print(f"  Converged: {converged}"
                  + (f" at x = {state['x_converged']:.1f}" if converged else ""))
            print(f"  YX_relic = {YX_arr[-1]:.6e}")
            print("=" * 60)

        result = {
            'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
            'YXeq': YXeq_arr, 'YYeq': YYeq_arr,
            'mubar_X': muX_arr, 'mubar_Y': muY_arr,
            'YX_relic': YX_arr[-1], 'YY_final': YY_arr[-1], 'xi_final': xi_arr[-1],
            'x_FO': state['x_FO'],
            'x_switch_QSSA': state['x_switch_QSSA'],
            'x_switch_full': state['x_switch_full'],
            'x_converged': state['x_converged'], 'converged': converged,
        }
        if return_bg_ICs:
            result['bg_ICs'] = {
                'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
                'T': T_arr, 'TD': TD_arr,
            }
        return result

    # =============================================================
    #  Boltzmann solver — hybrid Y-variables -> mu-variables
    # =============================================================
    def solve_boltzmann_hybrid(
        self,
        xmin: float = 1e-3,
        n_points: int = 500,
        rtol_value: float = 1e-8,
        atol_value: float = 1e-15,
        Gamma_switch_Yield: float = 5e9,
        cannibal_switch_mu: float = 1e-2,
        delta_switch_mu: float = 0.05,
        delta_eq_threshold: float = 0.05,
        convergence_threshold: float = 1e-2,
        convergence_mode: str = "dlnxi_NR",
        Tinf: float = 1e14,
        return_bg_ICs: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Three-phase hybrid Boltzmann solver.

        Phase 1   : equilibrium (muX=muY=0, evolve ln xi only)
        Phase 1.5 : Y-variables (ln xi, YX, YY)
        Phase 2   : mu-variables (ln xi, muX, muY)
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        include_antiparticlesX = model.include_antiparticlesX
        include_antiparticlesY = model.include_antiparticlesY

        xfList = np.concatenate([
            np.array([1e-2, 1e-1, 0.99]),
            np.arange(1, 100, 1),
            np.array([2e2, 5e2, 1e3, 2e3, 5e3]),
        ])

        sv_YXX_XX = model.sigmav2_YXX_to_XX()
        sv_YYX_YX = model.sigmav2_YYX_to_YX()

        r = mX / mY
        state = {"x_FO": None, "x_switch_Yield": None,
                 "x_switch_mu": None, "x_converged": None,
                 "switch_reason": None}

        # === Shared thermodynamics ===

        def _thermo(x, lnxi):
            xi = np.exp(np.clip(lnxi, -20, 20))
            T  = mX / x;  TD = xi * T
            s  = cosmo.s_entropy(T);  H = cosmo.hubble(T)
            dlngSdlnT_val = cosmo.dlngSdlnT(T)
            g_tilde = 1.0 + dlngSdlnT_val / 3.0
            zX = mX / TD;  zY = mY / TD

            BX1 = cosmo.Bfac1(zX, mX);  BX2 = cosmo.Bfac2(zX, mX)
            dBX1 = cosmo.dBfac1dT(zX, mX)
            BY1 = cosmo.Bfac1(zY, mY);  BY2 = cosmo.Bfac2(zY, mY)
            dBY1 = cosmo.dBfac1dT(zY, mY)
            lamX = cosmo.lambda_i(zX);  lamY = cosmo.lambda_i(zY)

            neqX = cosmo.neq_MB(zX, gX, mX, include_antiparticlesX)
            neqY = cosmo.neq_MB(zY, gY, mY, include_antiparticlesY)
            YeqX = neqX / s if s > 0 else 0.0
            YeqY = neqY / s if s > 0 else 0.0

            lnYeqX = cosmo.ln_Yeq(zX, gX, mX, include_antiparticlesX, s)
            lnYeqY = cosmo.ln_Yeq(zY, gY, mY, include_antiparticlesY, s)

            sv_XXYY = model.sigmav_XX_to_YY(TD)
            sv_YYY  = model.sigmav2_YYY_to_YY(TD)
            return dict(
                xi=xi, T=T, TD=TD, s=s, H=H,
                dlngSdlnT=dlngSdlnT_val, g_tilde=g_tilde, zX=zX, zY=zY,
                BX1=BX1, BX2=BX2, dBX1=dBX1,
                BY1=BY1, BY2=BY2, dBY1=dBY1,
                lamX=lamX, lamY=lamY,
                neqX=neqX, neqY=neqY, YeqX=YeqX, YeqY=YeqY,
                lnYeqX=lnYeqX, lnYeqY=lnYeqY,
                sv_XXYY=sv_XXYY, sv_YYY=sv_YYY,
            )

        def _Htot(th, nX, nY):
            rho_h = max(th['BX1'] * nX + th['BY1'] * nY, 0.0)
            return np.sqrt(th['H']**2 + rho_h / (3.0 * Mpl**2))

        # === Phase 1: Equilibrium ===

        def dlnxi_dx_equilibrium(x, lnxi):
            th = _thermo(x, lnxi)
            neqX, neqY = th['neqX'], th['neqY']
            Htot = _Htot(th, neqX, neqY)
            Num = neqY * th['BY2'] + neqX * th['BX2']
            dneqX = cosmo.dneqdT_MB(th['zX'], gX, mX, include_antiparticlesX)
            dneqY = cosmo.dneqdT_MB(th['zY'], gY, mY, include_antiparticlesY)
            Den = (th['BX1'] * dneqX + th['dBX1'] * neqX
                   + th['BY1'] * dneqY + th['dBY1'] * neqY)
            if np.abs(Den) < 1e-300:
                Den = np.copysign(1e-300, Den)
            dlnxi = (1.0 / x) * (1.0 - (3.0 + th['dlngSdlnT']) * Num / (th['TD'] * Den))
            Gamma_ann = (th['s'] / Htot
                         * (th['sv_XXYY'] / 2.0) * th['YeqX'])
            return dlnxi, Gamma_ann

        def BEQs_eq(t, y):
            dlnxi, _ = dlnxi_dx_equilibrium(t, y[0])
            return [np.clip(dlnxi, -10.0 / t, 10.0 / t)]

        # === Phase 1.5: Y-variables ===

        def compute_derivatives_Yield(x, lnxi, YX, YY):
            th = _thermo(x, lnxi)
            g_tilde = th['g_tilde']
            s = th['s']
            nX = YX * s;  nY = YY * s
            Htot = _Htot(th, nX, nY)
            YeqX, YeqY = th['YeqX'], th['YeqY']

            PreFac = s * g_tilde / (Htot * x)

            if YeqY > 1e-300:
                ratio = YeqX / YeqY
            else:
                ratio = 0.0
            term1 = YX + ratio * YY
            term2 = YX - ratio * YY
            dYX_dx = -PreFac * (th['sv_XXYY'] / 2.0) * term1 * term2

            sv_YYY = th['sv_YYY']
            S_3to2 = (YY**2 * sv_YYY
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            PreFac2 = s**2 * g_tilde / (Htot * x)
            diffY = YY - YeqY
            dYY_dx = -dYX_dx - PreFac2 * diffY * S_3to2

            delta_X = np.abs(YX - YeqX) / YeqX if YeqX > 1e-300 else 1e10

            if delta_X < delta_eq_threshold:
                Num = nY * th['BY2'] + nX * th['BX2']
                dneqX = cosmo.dneqdT_MB(th['zX'], gX, mX, include_antiparticlesX)
                dneqY = cosmo.dneqdT_MB(th['zY'], gY, mY, include_antiparticlesY)
                Den = (th['BX1'] * dneqX + th['dBX1'] * nX
                       + th['BY1'] * dneqY + th['dBY1'] * nY)
            else:
                Num = (nY * th['TD'] + nX * th['TD']
                       + x * s * (th['BX1'] * dYX_dx + th['BY1'] * dYY_dx)
                       / (g_tilde * 3.0))
                Den = th['dBX1'] * nX + th['dBY1'] * nY

            if np.abs(Den) < 1e-300:
                Den = np.copysign(1e-300, Den)

            dlnxi_dx = (1.0 / x) * (1.0 - (3.0 + th['dlngSdlnT']) * Num
                                     / (th['TD'] * Den))
            dlnxi_dx = np.clip(dlnxi_dx, -10.0 / x, 10.0 / x)

            Gamma_can = s**2 / Htot * S_3to2

            return dlnxi_dx, dYX_dx, dYY_dx, delta_X, Gamma_can

        def BEQs_Yield(t, y):
            dlnxi, dYX, dYY, _, _ = compute_derivatives_Yield(t, y[0], y[1], y[2])
            return [dlnxi, dYX, dYY]

        # === Phase 2: mu-variables ===

        def compute_derivatives_mu(x, lnxi, muX, muY):
            th = _thermo(x, lnxi)
            g_tilde = th['g_tilde']
            s = th['s']

            lnYX = th['lnYeqX'] + muX;  lnYY = th['lnYeqY'] + muY
            YX = np.exp(np.clip(lnYX, -700, 700))
            YY = np.exp(np.clip(lnYY, -700, 700))
            nX = YX * s;  nY = YY * s
            Htot = _Htot(th, nX, nY)

            lnRXY = th['lnYeqX'] - th['lnYeqY']
            Pfac = s * g_tilde / (Htot * x) * (th['sv_XXYY'] / 2.0)
            RXY2_YY2 = np.exp(np.clip(2.0 * lnRXY + 2.0 * lnYY, -700, 700))

            coll_X = Pfac * (YX - RXY2_YY2 / YX) if YX > 1e-300 else 0.0
            A = -coll_X + (th['lamX'] - 3.0 * g_tilde) / x
            B = -th['lamX']

            coll_Y_2to2 = Pfac * (YX**2 - RXY2_YY2) / YY if YY > 1e-300 else 0.0
            S_3to2 = (YY**2 * th['sv_YYY']
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            Gamma_can = s**2 * g_tilde / (Htot * x) * S_3to2

            if muY > 50:
                one_m = 1.0
            elif muY < -50:
                one_m = -np.exp(-muY)
            else:
                one_m = 1.0 - np.exp(-muY)

            C = coll_Y_2to2 - Gamma_can * one_m + (th['lamY'] - 3.0 * g_tilde) / x
            D = -th['lamY']

            D_cal  = nX * th['dBX1'] + nY * th['dBY1']
            E_cal  = th['BX1'] * nX * th['lamX'] + th['BY1'] * nY * th['lamY']
            Dtilde = th['TD'] * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            rho_plus_P = nX * th['BX2'] + nY * th['BY2']
            E = (1.0 / x) * (1.0 - 3.0 * g_tilde * rho_plus_P / Dtilde)
            F = -th['BX1'] * nX / Dtilde
            G = -th['BY1'] * nY / Dtilde

            denom = 1.0 - F * B - G * D
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)

            dlnxi_dx = np.clip((E + F * A + G * C) / denom, -10.0 / x, 10.0 / x)
            dmuX_dx  = A + B * dlnxi_dx
            dmuY_dx  = C + D * dlnxi_dx

            if muX > 0.05 and state["x_FO"] is None:
                state["x_FO"] = x

            return dlnxi_dx, dmuX_dx, dmuY_dx

        def BEQs_mu(t, y):
            return compute_derivatives_mu(t, *y)

        # === Initial conditions ===

        x0    = xmin
        lnxi0 = (1.0 / 3.0) * np.log(
            np.round(cosmo.gstarS(mX / xmin) / cosmo.gstarS(Tinf), 4)) + np.log(model.xi_infl)

        T0 = mX / x0;  s0 = cosmo.s_entropy(T0)
        YX0 = cosmo.neq_MB(mX / T0, gX, mX, include_antiparticlesX) / s0
        YY0 = cosmo.neq_MB(mY / T0, gY, mY, include_antiparticlesY) / s0

        lnYX_check_prev = np.log(YX0) if YX0 > 0 else -700.0
        x_check_prev = x0

        x_all, lnxi_all = [x0], [lnxi0]
        muX_all, muY_all = [0.0], [0.0]
        YX_all, YY_all   = [YX0], [YY0]

        phase = 1
        y0_eq    = [lnxi0]
        y0_Yield = None
        y0_mu    = None

        if verbose:
            print("=" * 60)
            print("Boltzmann Solver -- Hybrid Y -> mu")
            print("=" * 60)
            print(f"  mX = {mX:.2e} GeV, mY = {mY:.2e} GeV, r = {r:.2f}")
            print(f"  Phase 1->1.5 (Yield): Gamma_ann < {Gamma_switch_Yield:.0e}")
            print(f"  Phase 1.5->2 (mu): Gamma_can < {cannibal_switch_mu:.0e}"
                  f" OR delta > {delta_switch_mu}")
            print(f"  Convergence: {convergence_mode}, thr={convergence_threshold}")
            print("-" * 60)

        # === Main loop ===

        converged = False
        iterator = tqdm(range(len(xfList)), disable=not verbose, desc="Evolving")

        for j in iterator:
            xf = xfList[j]
            if x0 >= xf:
                continue

            if x0 < 1:
                xs = np.logspace(np.log10(x0 * 1.001), np.log10(xf * 0.999), n_points)
            else:
                xs = np.linspace(x0 * 1.001, xf * 0.999, n_points)

            # ----- PHASE 1: Equilibrium -----
            if phase == 1:
                sol = solve_ivp(BEQs_eq, (x0, xf), y0_eq, t_eval=xs,
                                rtol=rtol_value, atol=1e-12, method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (eq): {sol.message}")
                    break
                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend([0.0] * len(sol.t))
                muY_all.extend([0.0] * len(sol.t))

                for xi_x in sol.t:
                    th_tmp = _thermo(xi_x, sol.y[0, np.searchsorted(sol.t, xi_x)])
                    YX_all.append(th_tmp['YeqX'])
                    YY_all.append(th_tmp['YeqY'])

                x0 = sol.t[-1];  y0_eq = [sol.y[0, -1]]

                _, Gamma_now = dlnxi_dx_equilibrium(x0, y0_eq[0])
                if Gamma_now < Gamma_switch_Yield:
                    phase = 1.5
                    state["x_switch_Yield"] = x0
                    th_sw = _thermo(x0, y0_eq[0])
                    y0_Yield = [y0_eq[0], th_sw['YeqX'], th_sw['YeqY']]
                    if verbose:
                        print(f"\n  -> Phase 1.5 (Yield) at x = {x0:.2f}"
                              f" (Gamma_ann = {Gamma_now:.1e})")

            # ----- PHASE 1.5: Y-variables -----
            elif phase == 1.5:
                sol = solve_ivp(BEQs_Yield, (x0, xf), y0_Yield, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (Yield): {sol.message}")
                    break

                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                YX_all.extend(sol.y[1].tolist())
                YY_all.extend(sol.y[2].tolist())
                muX_all.extend([0.0] * len(sol.t))
                muY_all.extend([0.0] * len(sol.t))
                x0 = sol.t[-1];  y0_Yield = sol.y[:, -1].tolist()

                _, _, _, delta_now, Gamma_can_now = compute_derivatives_Yield(
                    x0, y0_Yield[0], y0_Yield[1], y0_Yield[2])

                if verbose:
                    iterator.set_postfix({
                        'ph': '1.5', 'x': f'{x0:.1f}',
                        'delta': f'{delta_now:.1e}', 'G_c': f'{Gamma_can_now:.1e}'})

                if Gamma_can_now < cannibal_switch_mu or delta_now > delta_switch_mu:
                    phase = 2
                    state["x_switch_mu"] = x0
                    reason = (f"Gamma_can={Gamma_can_now:.1e}"
                              if Gamma_can_now < cannibal_switch_mu
                              else f"delta={delta_now:.1e}")
                    state["switch_reason"] = reason

                    lnxi_sw = y0_Yield[0]
                    YX_sw   = y0_Yield[1]
                    YY_sw   = y0_Yield[2]
                    th_sw   = _thermo(x0, lnxi_sw)
                    muX_sw  = np.log(max(YX_sw, 1e-300)) - th_sw['lnYeqX']
                    muY_sw  = np.log(max(YY_sw, 1e-300)) - th_sw['lnYeqY']

                    y0_mu = [lnxi_sw, muX_sw, muY_sw]
                    if verbose:
                        print(f"\n  -> Phase 2 (mu) at x = {x0:.2f}"
                              f" ({reason},"
                              f" muX={muX_sw:.3f}, muY={muY_sw:.3f})")

            # ----- PHASE 2: mu-variables -----
            elif phase == 2:
                sol = solve_ivp(BEQs_mu, (x0, xf), y0_mu, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (mu): {sol.message}")
                    break

                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend(sol.y[1].tolist())
                muY_all.extend(sol.y[2].tolist())

                for k in range(len(sol.t)):
                    th_k = _thermo(sol.t[k], sol.y[0, k])
                    YX_all.append(np.exp(np.clip(
                        th_k['lnYeqX'] + sol.y[1, k], -700, 700)))
                    YY_all.append(np.exp(np.clip(
                        th_k['lnYeqY'] + sol.y[2, k], -700, 700)))

                x0 = sol.t[-1];  y0_mu = sol.y[:, -1].tolist()

                # Convergence check
                if state["x_FO"] is not None and x0 > 2.0 * state["x_FO"]:
                    lnxi_now = sol.y[0, -1]
                    muX_now  = sol.y[1, -1]
                    muY_now  = sol.y[2, -1]

                    if convergence_mode == 'dlnxi_NR':
                        dlnxi_num, _, _ = compute_derivatives_mu(
                            x0, lnxi_now, muX_now, muY_now)
                        th_now = _thermo(x0, lnxi_now)
                        g_tilde_now = th_now['g_tilde']
                        dlnxi_NR = (1.0 - 2.0 * g_tilde_now) / x0
                        rel_dev = np.abs(dlnxi_num - dlnxi_NR) * x0
                        if verbose:
                            iterator.set_postfix({
                                'ph': '2', 'x': f'{x0:.0f}',
                                'muX': f'{muX_now:.1f}', 'dev': f'{rel_dev:.1e}'})
                        if rel_dev < convergence_threshold:
                            converged = True
                            state["x_converged"] = x0
                            if verbose:
                                print(f"\n  Converged (dlnxi_NR) at x = {x0:.1f}")
                            break

                    elif convergence_mode == 'dlnYX':
                        xi_now  = np.exp(np.clip(lnxi_now, -20, 20))
                        T_now   = mX / x0
                        zX_now  = mX / (xi_now * T_now)
                        lnYX_now = muX_now + cosmo.ln_Yeq(
                            zX_now, gX, mX, include_antiparticlesX,
                            cosmo.s_entropy(T_now))
                        if (np.isfinite(lnYX_check_prev) and np.isfinite(lnYX_now)
                                and x0 > x_check_prev * 1.5):
                            dlnYX_dlnx = ((lnYX_now - lnYX_check_prev)
                                          / np.log(x0 / x_check_prev))
                            if verbose:
                                iterator.set_postfix({
                                    'ph': '2', 'x': f'{x0:.0f}',
                                    '|dlnYX|': f'{np.abs(dlnYX_dlnx):.1e}'})
                            if np.abs(dlnYX_dlnx) < convergence_threshold:
                                converged = True
                                state["x_converged"] = x0
                                if verbose:
                                    print(f"\n  Converged (dlnYX) at x = {x0:.1f}")
                                break
                            lnYX_check_prev = lnYX_now
                            x_check_prev = x0

        # === Post-process ===

        x_arr    = np.array(x_all)
        lnxi_arr = np.array(lnxi_all)
        xi_arr   = np.exp(np.clip(lnxi_arr, -20, 20))
        T_arr    = mX / x_arr
        TD_arr   = xi_arr * T_arr
        YX_arr   = np.array(YX_all)
        YY_arr   = np.array(YY_all)
        muX_arr  = np.array(muX_all)
        muY_arr  = np.array(muY_all)

        zX_arr = mX / TD_arr;  zY_arr = mY / TD_arr
        YXeq_arr = np.zeros_like(x_arr)
        YYeq_arr = np.zeros_like(x_arr)
        for i in range(len(x_arr)):
            s_i = cosmo.s_entropy(T_arr[i])
            YXeq_arr[i] = np.exp(np.clip(
                cosmo.ln_Yeq(zX_arr[i], gX, mX, include_antiparticlesX, s_i), -700, 700))
            YYeq_arr[i] = np.exp(np.clip(
                cosmo.ln_Yeq(zY_arr[i], gY, mY, include_antiparticlesY, s_i), -700, 700))

        if verbose:
            print("-" * 60)
            print("RESULTS:")
            if state.get('x_switch_Yield'):
                print(f"  x_switch (eq -> Yield): {state['x_switch_Yield']:.2f}")
            if state.get('x_switch_mu'):
                print(f"  x_switch (Yield -> mu): {state['x_switch_mu']:.2f}"
                      f"  [{state['switch_reason']}]")
            else:
                print(f"  (Yield phase carried to end, Phase 2 never entered)")
            if state['x_FO']:
                print(f"  x_FO (muX > 0.05): {state['x_FO']:.1f}")
            print(f"  Converged: {converged}"
                  + (f" at x = {state['x_converged']:.1f}" if converged else ""))
            print(f"  YX_relic = {YX_arr[-1]:.6e}")
            print("=" * 60)

        result = {
            'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
            'YXeq': YXeq_arr, 'YYeq': YYeq_arr,
            'mubar_X': muX_arr, 'mubar_Y': muY_arr,
            'YX_relic': YX_arr[-1], 'YY_final': YY_arr[-1], 'xi_final': xi_arr[-1],
            'x_FO': state['x_FO'],
            'x_switch_Yield': state.get('x_switch_Yield'),
            'x_switch_mu': state.get('x_switch_mu'),
            'switch_reason': state.get('switch_reason'),
            'x_converged': state['x_converged'], 'converged': converged,
        }
        if return_bg_ICs:
            result['bg_ICs'] = {
                'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
                'T': T_arr, 'TD': TD_arr,
            }
        return result

    # =============================================================
    #  Boltzmann solver — Y in equilibrium
    # =============================================================
    def solve_boltzmann_Yeq(
        self,
        xmin: float = 1e-3,
        n_points: int = 500,
        rtol_value: float = 1e-8,
        atol_value: float = 1e-15,
        Gamma_switch_threshold: float = 1e7,
        convergence_threshold: float = 1e-2,
        Tinf: float = 1e14,
        return_bg_ICs: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Boltzmann solver assuming Y remains in thermal equilibrium.
        Solves for (ln xi, mubar_X) only -- two state variables.
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        include_antiparticlesX = model.include_antiparticlesX
        include_antiparticlesY = model.include_antiparticlesY

        xfList = np.concatenate([
            np.array([1e-2, 1e-1, 0.99]),
            np.arange(1, 100, 1),
            np.array([2e2, 5e2, 1e3, 2e3]),
        ])

        sv_XXYY = model.sigmav_XX_to_YY_swave()
        r = mX / mY
        state = {"x_FO": None, "x_switch": None, "x_converged": None}

        def _compute_thermo(x, lnxi, muX):
            xi = np.exp(np.clip(lnxi, -20, 20))
            T  = mX / x
            TD = xi * T
            s  = cosmo.s_entropy(T)
            H  = cosmo.hubble(T)
            dlngSdlnT_val = cosmo.dlngSdlnT(T)
            g_tilde = 1.0 + dlngSdlnT_val / 3.0

            zX = mX / TD;  zY = mY / TD

            BX1 = cosmo.Bfac1(zX, mX);  BX2 = cosmo.Bfac2(zX, mX);  dBX1 = cosmo.dBfac1dT(zX, mX)
            BY1 = cosmo.Bfac1(zY, mY);  BY2 = cosmo.Bfac2(zY, mY);  dBY1 = cosmo.dBfac1dT(zY, mY)

            lamX = cosmo.lambda_i(zX);  lamY = cosmo.lambda_i(zY)

            lnYeqX = cosmo.ln_Yeq(zX, gX, mX, include_antiparticlesX, s)
            lnYeqY = cosmo.ln_Yeq(zY, gY, mY, include_antiparticlesY, s)
            lnYX = lnYeqX + muX
            lnYY = lnYeqY

            YX = np.exp(np.clip(lnYX, -700, 700))
            YY = np.exp(np.clip(lnYY, -700, 700))

            nX = YX * s
            nY = YY * s

            rho_h = max(BX1 * nX + BY1 * nY, 0.0)
            Htot  = np.sqrt(H**2 + rho_h / (3.0 * Mpl**2))

            return dict(
                xi=xi, T=T, TD=TD, s=s, H=H, Htot=Htot,
                dlngSdlnT=dlngSdlnT_val, g_tilde=g_tilde, zX=zX, zY=zY,
                BX1=BX1, BX2=BX2, dBX1=dBX1,
                BY1=BY1, BY2=BY2, dBY1=dBY1,
                lnYeqX=lnYeqX, lnYeqY=lnYeqY,
                YX=YX, YY=YY, nX=nX, nY=nY,
                lamX=lamX, lamY=lamY,
            )

        # --- Phase 1 ---

        def dlnxi_dx_equilibrium(x, lnxi):
            th = _compute_thermo(x, lnxi, 0.0)
            Num = th['nY'] * th['BY2'] + th['nX'] * th['BX2']
            dneqX = cosmo.dneqdT_MB(th['zX'], gX, mX, include_antiparticlesX)
            dneqY = cosmo.dneqdT_MB(th['zY'], gY, mY, include_antiparticlesY)
            Den = (th['BX1'] * dneqX + th['dBX1'] * th['nX']
                   + th['BY1'] * dneqY + th['dBY1'] * th['nY'])
            if np.abs(Den) < 1e-300:
                Den = np.copysign(1e-300, Den)

            dlnxi = (1.0 / x) * (1.0 - (3.0 + th['dlngSdlnT']) * Num / (th['TD'] * Den))
            Gamma_ann = (th['s'] / th['Htot']
                         * sv_XXYY * th['YX'])
            return dlnxi, Gamma_ann

        def BEQs_eq(t, y):
            dlnxi, _ = dlnxi_dx_equilibrium(t, y[0])
            return [np.clip(dlnxi, -10.0 / t, 10.0 / t)]

        # --- Phase 2 ---

        def compute_derivatives_full(x, lnxi, muX):
            th = _compute_thermo(x, lnxi, muX)
            g_tilde = th['g_tilde']
            YX = th['YX']
            nX, nY = th['nX'], th['nY']

            YXeq = np.exp(np.clip(th['lnYeqX'], -700, 700))
            Gamma_coll = th['s'] * g_tilde / (th['Htot'] * x) * sv_XXYY * YXeq

            A = -Gamma_coll * np.sinh(muX) + (th['lamX'] - 3.0 * g_tilde) / x
            B = -th['lamX']

            D_cal  = nX * th['dBX1'] + nY * th['dBY1']
            E_cal  = th['BX1'] * nX * th['lamX'] + th['BY1'] * nY * th['lamY']
            Dtilde = th['TD'] * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            Num_xi = nX * th['BX2'] + nY * th['BY2']
            C = (1.0 / x) * (1.0 - 3.0 * g_tilde * Num_xi / Dtilde)
            D = -th['BX1'] * nX / Dtilde

            denom = 1.0 - B * D
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)

            dlnxi_dx = np.clip((C + D * A) / denom, -10.0 / x, 10.0 / x)
            dmuX_dx  = (A + B * C) / denom

            if muX > 0.05 and state["x_FO"] is None:
                state["x_FO"] = x

            return dlnxi_dx, dmuX_dx

        def BEQs_full(t, y):
            return compute_derivatives_full(t, *y)

        # --- ICs ---

        x0    = xmin
        lnxi0 = (1.0 / 3.0) * np.log(
            np.round(cosmo.gstarS(mX / xmin) / cosmo.gstarS(Tinf), 4)) + np.log(model.xi_infl)
        lnYX_check_prev = cosmo.ln_Yeq(
            mX / (mX / x0), gX, mX, include_antiparticlesX,
            cosmo.s_entropy(mX / x0))
        x_check_prev = x0

        x_all, lnxi_all, muX_all = [x0], [lnxi0], [0.0]
        in_equilibrium = True
        y0_eq   = [lnxi0]
        y0_full = None

        if verbose:
            print("=" * 60)
            print("Boltzmann Solver -- Y in Equilibrium")
            print("=" * 60)
            print(f"  mX = {mX:.2e} GeV, mY = {mY:.2e} GeV, r = {r:.2f}")
            print(f"  Phase 1->2 switch: Gamma_ann < {Gamma_switch_threshold:.0e}")
            print(f"  Convergence: |d ln YX / d ln x| < {convergence_threshold}")
            print("-" * 60)

        # --- Main loop ---

        converged = False
        iterator = tqdm(range(len(xfList)), disable=not verbose, desc="Evolving")

        for j in iterator:
            xf = xfList[j]
            if x0 >= xf:
                continue

            if x0 < 1:
                xs = np.logspace(np.log10(x0 * 1.001), np.log10(xf * 0.999), n_points)
            else:
                xs = np.linspace(x0 * 1.001, xf * 0.999, n_points)

            if in_equilibrium:
                sol = solve_ivp(BEQs_eq, (x0, xf), y0_eq, t_eval=xs,
                                rtol=rtol_value, atol=1e-12, method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (eq): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break

                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend([0.0] * len(sol.t))
                x0 = sol.t[-1];  y0_eq = [sol.y[0, -1]]

                _, Gamma_now = dlnxi_dx_equilibrium(x0, y0_eq[0])
                if Gamma_now < Gamma_switch_threshold:
                    in_equilibrium = False
                    state["x_switch"] = x0
                    y0_full = [y0_eq[0], 0.0]
                    if verbose:
                        print(f"\n  -> Full system at x = {x0:.2f}"
                              f" (Gamma_ann = {Gamma_now:.1e})")

            else:
                sol = solve_ivp(BEQs_full, (x0, xf), y0_full, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (full): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break

                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend(sol.y[1].tolist())
                x0 = sol.t[-1];  y0_full = sol.y[:, -1].tolist()

                if state["x_FO"] is not None and x0 > 2.0 * state["x_FO"]:
                    xi_now  = np.exp(np.clip(sol.y[0, -1], -20, 20))
                    muX_now = sol.y[1, -1]
                    T_now   = mX / x0
                    zX_now  = mX / (xi_now * T_now)
                    lnYX_now = muX_now + cosmo.ln_Yeq(
                        zX_now, gX, mX, include_antiparticlesX,
                        cosmo.s_entropy(T_now))

                    if (np.isfinite(lnYX_check_prev) and np.isfinite(lnYX_now)
                            and x0 > x_check_prev * 1.5):
                        dlnYX_dlnx = ((lnYX_now - lnYX_check_prev)
                                      / np.log(x0 / x_check_prev))
                        if verbose:
                            iterator.set_postfix({
                                'x': f'{x0:.0f}', 'muX': f'{muX_now:.1f}',
                                '|dlnYX/dlnx|': f'{np.abs(dlnYX_dlnx):.1e}'})
                        if np.abs(dlnYX_dlnx) < convergence_threshold:
                            converged = True
                            state["x_converged"] = x0
                            if verbose:
                                print(f"\n  Converged at x = {x0:.1f}")
                            break
                        lnYX_check_prev = lnYX_now
                        x_check_prev = x0

        # --- Post-process ---

        x_arr    = np.array(x_all)
        lnxi_arr = np.array(lnxi_all)
        muX_arr  = np.array(muX_all)

        xi_arr = np.exp(np.clip(lnxi_arr, -20, 20))
        T_arr  = mX / x_arr
        TD_arr = xi_arr * T_arr
        zX_arr = mX / TD_arr;  zY_arr = mY / TD_arr

        lnYX_arr   = np.zeros_like(x_arr)
        lnYXeq_arr = np.zeros_like(x_arr)
        lnYYeq_arr = np.zeros_like(x_arr)
        for i in range(len(x_arr)):
            s_i = cosmo.s_entropy(T_arr[i])
            lnYXeq_arr[i] = cosmo.ln_Yeq(zX_arr[i], gX, mX, include_antiparticlesX, s_i)
            lnYYeq_arr[i] = cosmo.ln_Yeq(zY_arr[i], gY, mY, include_antiparticlesY, s_i)
            lnYX_arr[i]   = muX_arr[i] + lnYXeq_arr[i]

        YX_arr   = np.exp(np.clip(lnYX_arr,   -700, 700))
        YXeq_arr = np.exp(np.clip(lnYXeq_arr, -700, 700))
        YYeq_arr = np.exp(np.clip(lnYYeq_arr, -700, 700))
        YY_arr   = YYeq_arr.copy()

        if verbose:
            print("-" * 60)
            print("RESULTS:")
            if state['x_switch']:
                print(f"  x_switch (eq -> full): {state['x_switch']:.2f}")
            if state['x_FO']:
                print(f"  x_FO (muX > 0.05): {state['x_FO']:.1f}")
            print(f"  Converged: {converged}"
                  + (f" at x = {state['x_converged']:.1f}" if converged else ""))
            print(f"  YX_relic = {YX_arr[-1]:.6e}")
            print("=" * 60)

        result = {
            'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
            'YXeq': YXeq_arr, 'YYeq': YYeq_arr,
            'mubar_X': muX_arr, 'mubar_Y': np.zeros_like(muX_arr),
            'YX_relic': YX_arr[-1], 'YY_final': YY_arr[-1], 'xi_final': xi_arr[-1],
            'x_FO': state['x_FO'], 'x_switch': state['x_switch'],
            'x_converged': state['x_converged'], 'converged': converged,
        }
        if return_bg_ICs:
            result['bg_ICs'] = {
                'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
                'T': T_arr, 'TD': TD_arr,
            }
        return result

    # =============================================================
    #  Background evolution solver (late-time Y decay)
    # =============================================================
    def solve_background(
        self,
        epsX: float,
        bg_ICs: Dict[str, np.ndarray],
        Nmax: float = 80.0,
        rtol_value: float = 1e-6,
        atol_value: float = 1e-18,
        stop_on_rhoY: bool = True,
        rhoY_threshold: float = 1e-12,
        verbose: bool = False,
    ) -> Dict[str, Any]:
        """
        Background evolution: late-time Y decay with entropy injection.

        Evolves three log-scaled variables in e-folds N = ln(a/a_0):
            y[0] = ln(T/T_0),  y[1] = ln(rhoY/rhoY_0),  y[2] = ln(rhoX/rhoX_0)
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY

        x0  = bg_ICs['x'][-1]
        xi0 = bg_ICs['xi'][-1]
        YX0 = bg_ICs['YX'][-1]
        YY0 = bg_ICs['YY'][-1]
        T0  = bg_ICs['T'][-1]
        TD0 = bg_ICs['TD'][-1]

        s0    = cosmo.s_entropy(T0)
        rhoX0 = mX * YX0 * s0
        rhoY0 = cosmo.Bfac1(mY / TD0, mY) * YY0 * s0

        sv_XXYY = model.sigmav_XX_to_YY_swave()
        GammaY  = model.decay_width_to_SM(epsX)

        def rho_r(T):
            return np.pi**2 / 30.0 * cosmo.gstar(T) * T**4

        def ode_system(N, y):
            T    = T0 * np.exp(np.clip(y[0], -50, 50))
            rhoY = rhoY0 * np.exp(np.clip(y[1], -700, 200))
            rhoX = rhoX0 * np.exp(np.clip(y[2], -700, 200))

            T    = max(T,    1e-30)
            rhoY = max(rhoY, 1e-300)
            rhoX = max(rhoX, 1e-300)

            rhor    = rho_r(T)
            rho_tot = rhor + rhoX + rhoY
            H       = np.sqrt(rho_tot / (3.0 * Mpl**2))
            s       = cosmo.s_entropy(T)

            dlngS_dlnT = cosmo.dlngSdlnT(T)
            g_tilde_val = 1.0 + dlngS_dlnT / 3.0

            dlnrhoY_dN = -(3.0 + GammaY / max(H, 1e-300))

            injection = rhoY * GammaY / (max(H, 1e-300) * 3.0 * s * T)
            dlnT_dN   = -1.0 / g_tilde_val * (1.0 - injection)

            dlnrhoX_dN = -3.0 - sv_XXYY * rhoX / (2.0 * mX * max(H, 1e-300))

            return [dlnT_dN, dlnrhoY_dN, dlnrhoX_dN]

        def rhoY_negligible(N, y):
            T    = T0 * np.exp(np.clip(y[0], -50, 50))
            rhoY = rhoY0 * np.exp(np.clip(y[1], -700, 200))
            rhoX = rhoX0 * np.exp(np.clip(y[2], -700, 200))
            rhor    = rho_r(max(T, 1e-30))
            rho_tot = rhor + max(rhoX, 1e-300) + max(rhoY, 1e-300)
            return rhoY / rho_tot - rhoY_threshold

        rhoY_negligible.terminal  = True
        rhoY_negligible.direction = -1

        y0_bg  = [0.0, 0.0, 0.0]
        events = [rhoY_negligible] if stop_on_rhoY else None

        sol = solve_ivp(
            ode_system, (0, Nmax), y0_bg,
            method='Radau', max_step=0.1,
            rtol=rtol_value, atol=atol_value,
            events=events,
        )

        if verbose:
            print(f"BG_Solver: N_final = {sol.t[-1]:.2f}, success = {sol.success}")
            if sol.t_events and len(sol.t_events[0]) > 0:
                print(f"  rhoY/rho_tot threshold hit at N = {sol.t_events[0][0]:.2f}")

        Ns      = sol.t
        Tsol    = T0 * np.exp(sol.y[0])
        rhoYsol = rhoY0 * np.exp(sol.y[1])
        rhoXsol = rhoX0 * np.exp(sol.y[2])
        rhoRsol = rho_r(Tsol)

        S_BG  = cosmo.s_entropy(Tsol) * np.exp(3.0 * Ns)
        S0_BG = S_BG[0]

        xi_BG = xi0 * np.exp(-2.0 * Ns - sol.y[0])

        if stop_on_rhoY:
            s_final  = cosmo.s_entropy(Tsol[-1])
            YX_final = rhoXsol[-1] / (mX * s_final)
            SRatio   = S_BG[-1] / S0_BG
            if verbose:
                print(f"  YX_final = {YX_final:.6e}")
                print(f"  SRatio   = {SRatio:.4f}")
            return {'YX_final': YX_final, 'SRatio': SRatio, 'T_final': Tsol[-1]}

        else:
            x_pre  = bg_ICs['x']
            T_pre  = bg_ICs['T']
            TD_pre = bg_ICs['TD']
            YX_pre = bg_ICs['YX']
            YY_pre = bg_ICs['YY']
            xi_pre = bg_ICs['xi']
            n_pre  = len(x_pre)

            rhoX_pre = np.array([
                cosmo.Bfac1(mX / TD_pre[i], mX) * YX_pre[i] * cosmo.s_entropy(T_pre[i])
                for i in range(n_pre)])
            rhoY_pre = np.array([
                cosmo.Bfac1(mY / TD_pre[i], mY) * YY_pre[i] * cosmo.s_entropy(T_pre[i])
                for i in range(n_pre)])
            rhoR_pre = np.pi**2 / 30.0 * cosmo.gstar(T_pre) * T_pre**4
            S_pre    = np.full(n_pre, S0_BG)

            x_BG = mX / Tsol[1:]

            x_out    = np.concatenate([x_pre,  x_BG])
            T_out    = mX / x_out
            rhoX_out = np.concatenate([rhoX_pre, rhoXsol[1:]])
            rhoY_out = np.concatenate([rhoY_pre, rhoYsol[1:]])
            rhoR_out = np.concatenate([rhoR_pre, rhoRsol[1:]])
            S_out    = np.concatenate([S_pre,  S_BG[1:]])
            xi_out   = np.concatenate([xi_pre, xi_BG[1:]])

            SoverSf = S_out / S_out[-1]
            SoverSi = S_out / S_out[0]

            return {
                'x':       x_out,
                'T':       T_out,
                'rhoX':    rhoX_out,
                'rhoY':    rhoY_out,
                'rhoR':    rhoR_out,
                'rhoTot':  rhoX_out + rhoY_out + rhoR_out,
                'SoverSf': SoverSf,
                'SoverSi': SoverSi,
                'xi':      xi_out,
            }

    # =============================================================
    #  Epsilon finder — correct DM relic abundance
    # =============================================================
    def find_epsilon_DMRelic(
        self,
        bg_ICs: Dict[str, np.ndarray],
        Och2_target: float = 0.12,
        log10eps_range: tuple = (-20, -6),
        root_rtol: float = 1e-3,
        bg_kw: Optional[Dict] = None, verbose=False
    ) -> float:
        if bg_kw is None:
            bg_kw = {}

        Y_target = Och2_target * rhoc / (s0_cosmo * self.model.mX)

        def f(log10eps):
            res = self.solve_background(
                epsX=10**log10eps, bg_ICs=bg_ICs,
                stop_on_rhoY=True, verbose=False, **bg_kw,
            )
            if verbose:
                print(res['YX_final']*self.model.mX)
            return np.log10(res['YX_final']) - np.log10(Y_target)

        log10eps = brentq(f, log10eps_range[0], log10eps_range[1], rtol=root_rtol)
        return 10**log10eps

    # =============================================================
    #  Epsilon finder — BBN constraint
    # =============================================================
    def find_epsilon_BBN(
        self,
        bg_ICs: Dict[str, np.ndarray],
        T_BBN: float = 10e-3,
        log10eps_range: tuple = (-20, -6),
        root_rtol: float = 1e-3,
        rhoY_threshold: float = 1e-3,
        bg_kw: Optional[Dict] = None,
    ) -> float:
        if bg_kw is None:
            bg_kw = {}

        def f(log10eps):
            res = self.solve_background(
                epsX=10**log10eps, bg_ICs=bg_ICs,
                stop_on_rhoY=True, verbose=False, rhoY_threshold=rhoY_threshold, **bg_kw,
            )
            return res['T_final'] - T_BBN

        log10eps = brentq(f, log10eps_range[0], log10eps_range[1], rtol=root_rtol)
        return 10**log10eps

    # =============================================================
    #  Epsilon finder — BBN constraint via DN_eff
    # =============================================================
    def find_epsilon_BBN_Neff(
        self,
        bg_ICs,
        T_eval=1e-3,
        DNeff_max=0.151,
        log10eps_range=(-13, -9),
        root_rtol=1e-3,
        bg_kw=None,
        verbose=False,
    ):
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY

        T0  = bg_ICs['T'][-1]
        TD0 = bg_ICs['TD'][-1]
        YX0 = bg_ICs['YX'][-1]
        YY0 = bg_ICs['YY'][-1]
        s0  = cosmo.s_entropy(T0)
        rhoX0 = mX * YX0 * s0
        rhoY0 = cosmo.Bfac1(mY / TD0, mY) * YY0 * s0
        sv    = model.sigmav_XX_to_YY_swave()

        rho_1nu = 7.0 * np.pi**2 / 120.0 * T_eval**4

        def _rhoY_at_Teval(epsX):
            GammaY = model.decay_width_to_SM(epsX)

            def rho_r(T):
                return np.pi**2 / 30.0 * cosmo.gstar(T) * T**4

            def ode(N, y):
                T    = T0 * np.exp(np.clip(y[0], -50, 50))
                rhoY = rhoY0 * np.exp(np.clip(y[1], -700, 200))
                rhoX = rhoX0 * np.exp(np.clip(y[2], -700, 200))
                T    = max(T, 1e-30)
                rhoY = max(rhoY, 1e-300)
                rhoX = max(rhoX, 1e-300)
                rhor = rho_r(T)
                H = np.sqrt((rhor + rhoX + rhoY) / (3.0 * Mpl**2))
                s = cosmo.s_entropy(T)
                g_tilde = 1.0 + cosmo.dlngSdlnT(T) / 3.0
                inj = rhoY * GammaY / (max(H, 1e-300) * 3.0 * s * T)
                return [
                    -1.0 / g_tilde * (1.0 - inj),
                    -(3.0 + GammaY / max(H, 1e-300)),
                    -3.0 - sv * rhoX / (2.0 * mX * max(H, 1e-300)),
                ]

            def T_hit(N, y):
                return T0 * np.exp(np.clip(y[0], -50, 50)) - T_eval
            T_hit.terminal  = True
            T_hit.direction = -1

            sol = solve_ivp(
                ode, (0, 80), [0.0, 0.0, 0.0],
                method='Radau', max_step=0.5,
                rtol=1e-6, atol=1e-18,
                events=[T_hit],
            )
            return rhoY0 * np.exp(sol.y[1, -1])

        def f(log10eps):
            rhoY_e = _rhoY_at_Teval(10**log10eps)
            DNeff  = rhoY_e / rho_1nu
            if verbose:
                print(f"  log10eps={log10eps:+.3f}  DNeff={DNeff:.4e}")
            return DNeff - DNeff_max

        log10eps = brentq(f, log10eps_range[0], log10eps_range[1], rtol=root_rtol)
        return 10**log10eps


# =================================================================
#  find_mXstar  (module-level function)
# =================================================================

def find_mXstar(
    cosmo: Cosmology,
    model_factory,
    *,
    r: float,
    # --- interpolation ---
    n_refine: int = 5,
    bracket_factor: tuple = (0.1, 5.0),
    epsX_ref: float = 1e-1,
    Och2_target: float = Och2,
    # --- epsX_min ---
    find_epsX_min: bool = True,
    delta_S: float = 0.001,
    log10eps_range: tuple = (-20, -4),
    # --- Boltzmann solver ---
    solver_method: str = "QSSA",
    solver_kw: Optional[Dict] = None,
    # --- background solver ---
    bg_kw: Optional[Dict] = None,
    verbose: bool = True,
    verbose_solver: bool = False,
) -> Dict[str, Any]:
    """
    Find mX* where freeze-out + prompt Y decay yields Och2_target,
    and optionally the minimum epsX for negligible entropy injection.

    Parameters
    ----------
    cosmo : Cosmology
    model_factory : callable(mX, mY) -> DarkSectorModel
        Factory that creates a model instance for given masses.
    r : float
        Mass ratio mX / mY.
    solver_method : str
        Which Boltzmann solver to use: "QSSA", "chempot", "hybrid", or "Yeq".
    solver_kw : dict
        Keyword arguments passed to the chosen solver method.
    verbose : bool
        Print find_mXstar progress.
    verbose_solver : bool
        Also print verbose output from each individual Boltzmann solve.
    """
    if solver_kw is None:
        solver_kw = {}
    _bg_kw = dict(stop_on_rhoY=True, verbose=False)
    if bg_kw:
        _bg_kw.update(bg_kw)

    # Estimate mX scale
    _model_tmp = model_factory(1.0, 1.0 / r)
    # For the estimate we need alphaX-like info; use sigmav_swave ~ pi * alpha^2 / mX^2
    # mX_est ~ alpha * sqrt(pi * mY_relic * Mpl)
    # We'll compute it from the model's s-wave cross section at unit mass
    sv1 = _model_tmp.sigmav_XX_to_YY_swave()  # pi * alpha^2 / 1^2 for VP
    alpha_eff = np.sqrt(sv1 / np.pi)
    mX_est = alpha_eff * np.sqrt(np.pi * mY_relic * Mpl)

    log10mX_lo = np.log10(bracket_factor[0] * mX_est)
    log10mX_hi = np.log10(bracket_factor[1] * mX_est)

    # Solver dispatch
    _skw = dict(return_bg_ICs=True, verbose=verbose_solver)
    _skw.update(solver_kw)

    def _solve(solver):
        if solver_method == "QSSA":
            return solver.solve_boltzmann_chempot_QSSA(**_skw)
        elif solver_method == "chempot":
            return solver.solve_boltzmann_chempot(**_skw)
        elif solver_method == "hybrid":
            return solver.solve_boltzmann_hybrid(**_skw)
        elif solver_method == "Yeq":
            return solver.solve_boltzmann_Yeq(**_skw)
        else:
            raise ValueError(f"Unknown solver_method: {solver_method}")

    def eval_yield(log10mX):
        mX = 10.0 ** log10mX
        mdl = model_factory(mX, mX / r)
        solver = BoltzmannSolver(cosmo, mdl)
        sol = _solve(solver)
        sol_bg = solver.solve_background(
            epsX=epsX_ref, bg_ICs=sol['bg_ICs'], **_bg_kw)
        return np.log10(sol_bg['YX_final'] * mX)

    def interp_mX(lm, lmYX, target):
        o = np.argsort(lmYX)
        return float(np.interp(target, lmYX[o], lm[o]))

    log10_target = np.log10(mY_relic)

    if verbose:
        print("=" * 60)
        print(f"find_mXstar  (r={r:.2f})")
        print(f"  mX_est = {mX_est:.2e} GeV  "
              f"bracket = [{10**log10mX_lo:.2e}, {10**log10mX_hi:.2e}]")
        print("-" * 60)

    # Pass 1: bracket endpoints
    lm = np.array([log10mX_lo, log10mX_hi])
    lmYX = np.array([eval_yield(lm[0]), eval_yield(lm[1])])

    if verbose:
        print(f"  lo: mX={10**lm[0]:.2e} -> mX*YX={10**lmYX[0]:.2e}")
        print(f"  hi: mX={10**lm[1]:.2e} -> mX*YX={10**lmYX[1]:.2e}")

    log10mX_guess = interp_mX(lm, lmYX, log10_target)
    if verbose:
        print(f"  guess: mX = {10**log10mX_guess:.4e}")

    # Pass 2: refine
    refine = np.linspace(log10mX_guess + np.log10(0.5),
                         log10mX_guess + np.log10(2.0), n_refine)
    for rv in refine:
        lm = np.append(lm, rv)
        lmYX = np.append(lmYX, eval_yield(rv))
        if verbose:
            print(f"  refine: mX={10**rv:.4e} -> mX*YX={10**lmYX[-1]:.4e}")

    log10mX_sol = interp_mX(lm, lmYX, log10_target)
    mXstar = 10.0 ** log10mX_sol

    o = np.argsort(lm)
    log10mYX_sol = np.interp(log10mX_sol, lm[o], lmYX[o])
    Och2_achieved = 10.0 ** log10mYX_sol * s0_cosmo / rhoc

    if verbose:
        print(f"  => mX* = {mXstar:.4e} GeV  mY* = {mXstar/r:.4e} GeV  "
              f"Och2 = {Och2_achieved:.6f}")

    result = {
        'mXstar': mXstar, 'mYstar': mXstar / r,
        'Och2_achieved': Och2_achieved, 'mX_est': mX_est,
        'grid_log10mX': lm, 'grid_log10mYX': lmYX,
    }

    if find_epsX_min:
        mdl_star = model_factory(mXstar, mXstar / r)
        solver_star = BoltzmannSolver(cosmo, mdl_star)
        sol_star = _solve(solver_star)
        S_plateau = solver_star.solve_background(
            epsX=epsX_ref, bg_ICs=sol_star['bg_ICs'], **_bg_kw)['SRatio']

        def obj(log10eps):
            S = solver_star.solve_background(
                epsX=10.0 ** log10eps,
                bg_ICs=sol_star['bg_ICs'], **_bg_kw)['SRatio']
            return S / S_plateau - 1.0 - delta_S

        log10eps_sol = brentq(
            obj, log10eps_range[0], log10eps_range[1], rtol=1e-4)
        result['epsX_min'] = 10.0 ** log10eps_sol
        result['SRatio_plateau'] = S_plateau

        if verbose:
            print(f"  => epsX_min = {result['epsX_min']:.4e}  "
                  f"(dS/S < {delta_S})  SRatio_plateau = {S_plateau:.4f}")

    if verbose:
        print("=" * 60)

    return result

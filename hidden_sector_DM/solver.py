"""
BoltzmannSolver — model-independent Boltzmann equation machinery.

All hidden-sector rates are read from a HiddenSectorModel instance.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Dict, Any
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import kv as Knu
from tqdm import tqdm

from .constants import Mpl, rhoc, s0_cosmo, Och2, mY_relic
from .cosmology import Cosmology
from .model import HiddenSectorModel


class BoltzmannSolver:
    """
    Model-independent Boltzmann solver for hidden sectors.

    Parameters
    ----------
    cosmo : Cosmology
        Shared SM-background instance.
    model : HiddenSectorModel
        Hidden-sector model providing masses, DOF, and rate methods.
    """

    def __init__(self, cosmo: Cosmology, model: HiddenSectorModel):
        self.cosmo = cosmo
        self.model = model

    # =============================================================
    #  Boltzmann solver — chemical-potential formulation (3-phase)
    # =============================================================
    def solve_boltzmann_chempot_3phase(
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

        Treats Y as stable through freeze-out (the secluded regime); the
        Y -> SM portal decay enters only at the post-processing step via
        `solve_background(..., correct_Y_decay=True)`.  When the portal
        coupling is large enough that Y decay is dynamically important
        during freeze-out (epsX >= ~0.5 * eps_therm), use
        `solve_boltzmann_joint` instead, which evolves the dark sector
        and the Y <-> SM coupling jointly.
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        include_antiparticlesX = model.include_antiparticlesX
        include_antiparticlesY = model.include_antiparticlesY

        # Symmetry factor: 1/2 for Dirac (X != Xbar), 1 for Majorana (X = Xbar)
        ann_sym = 0.5 if include_antiparticlesX else 1.0

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

        sv_YXX_XX = model.sigmav2_YXX_to_XX() * ann_sym   # 1/2 for Dirac (XX in initial state)
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
            Gamma_over_H_ann = (th['s'] / th['Htot']
                         * (th['sv_XXYY'] * ann_sym) * th['YX'])
            return dlnxi, Gamma_over_H_ann

        def BEQs_eq(t, y):
            dlnxi, _ = dlnxi_dx_equilibrium(t, y[0])
            return [np.clip(dlnxi, -10.0 / t, 10.0 / t)]

        def ann_switch_event(t, y):
            """Event: Gamma_over_H_ann - threshold = 0.  Terminal, direction=-1."""
            _, Gamma_over_H_ann = dlnxi_dx_equilibrium(t, y[0])
            return Gamma_over_H_ann - Gamma_switch_QSSA
        ann_switch_event.terminal  = True
        ann_switch_event.direction = -1

        # === Phase 1.5: QSSA (muY=0, evolve lnxi & muX) ===

        def compute_derivatives_QSSA(x, lnxi, muX):
            th = _thermo_Yeq(x, lnxi, muX)
            g_tilde = th['g_tilde']
            nX, nY = th['nX'], th['nY']

            YXeq = np.exp(np.clip(th['lnYeqX'], -700, 700))
            sv = model.sigmav_XX_to_YY(th['TD'])
            # 2*ann_sym matches the Phase-2 form Pfac*(YX - YXeq^2/YX)
            # = 2*ann_sym * sv * s * g_tilde * YXeq * sinh(muX) / (Htot*x).
            Gamma_coll = (2.0 * ann_sym * th['s'] * g_tilde
                          / (th['Htot'] * x) * sv * YXeq)

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

        def cannibal_switch_event(t, y):
            """Event: Gamma_over_H_can - threshold = 0.  Terminal, direction=-1."""
            return estimate_cannibal_rate(t, y[0], y[1]) - cannibal_switch_full
        cannibal_switch_event.terminal  = True
        cannibal_switch_event.direction = -1

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

            Pfac = th['s'] * g_tilde / (th['Htot'] * x) * (th['sv_XXYY'] * ann_sym)
            RXY2_YY2 = np.exp(np.clip(2.0 * th['lnRXY'] + 2.0 * th['lnYY'], -700, 700))

            coll_X = Pfac * (YX - RXY2_YY2 / YX) if YX > 1e-300 else 0.0
            A = -coll_X + (th['lamX'] - 3.0 * g_tilde) / x
            B = -th['lamX']

            coll_Y_2to2 = Pfac * (YX**2 - RXY2_YY2) / YY if YY > 1e-300 else 0.0
            S_3to2 = (YY**2 * th['sv_YYY_YY']
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            Gamma_over_H_can = th['s']**2 * g_tilde / (th['Htot'] * x) * S_3to2

            if muY > 50:
                one_m = 1.0
            elif muY < -50:
                one_m = -np.exp(-muY)
            else:
                one_m = 1.0 - np.exp(-muY)

            C = coll_Y_2to2 - Gamma_over_H_can * one_m + (th['lamY'] - 3.0 * g_tilde) / x
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

        # === Convergence events (continuous detection) ===

        if convergence_mode == 'dlnxi_NR':
            def _conv_event_QSSA(t, y):
                if state["x_FO"] is None or t < 2.0 * state["x_FO"]:
                    return 1.0
                dlnxi_num, _ = compute_derivatives_QSSA(t, y[0], y[1])
                g_t = _thermo_Yeq(t, y[0], y[1])['g_tilde']
                dlnxi_NR = (1.0 - 2.0 * g_t) / t
                return np.abs(dlnxi_num - dlnxi_NR) * t - convergence_threshold
            _conv_event_QSSA.terminal  = True
            _conv_event_QSSA.direction = -1

            def _conv_event_full(t, y):
                if state["x_FO"] is None or t < 2.0 * state["x_FO"]:
                    return 1.0
                dlnxi_num, _, _ = compute_derivatives_full(t, *y)
                g_t = _thermo_full(t, *y)['g_tilde']
                dlnxi_NR = (1.0 - 2.0 * g_t) / t
                return np.abs(dlnxi_num - dlnxi_NR) * t - convergence_threshold
            _conv_event_full.terminal  = True
            _conv_event_full.direction = -1
        else:
            _conv_event_QSSA = None
            _conv_event_full = None

        # === Initial conditions ===

        x0    = xmin
        lnxi0 = (1.0 / 3.0) * np.log(
            np.round(cosmo.gstarS(mX / xmin) / cosmo.gstarS(Tinf), 4)) + np.log(model.xi_ini)
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
            print(f"  Phase 1->1.5: Gamma_over_H_ann < {Gamma_switch_QSSA:.0e}")
            print(f"  Phase 1.5->2: Gamma_over_H_can < {cannibal_switch_full:.0e}")
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
                                rtol=rtol_value, atol=1e-12, method='Radau',
                                max_step=1, events=[ann_switch_event])
                if not sol.success:
                    if verbose:
                        print(f"  Warning (eq): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break
                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend([0.0] * len(sol.t))
                muY_all.extend([0.0] * len(sol.t))

                # Check if annihilation event fired (Phase 1.5 trigger)
                event_fired = (sol.t_events is not None
                               and len(sol.t_events) > 0
                               and len(sol.t_events[0]) > 0)

                if event_fired:
                    x_ev = sol.t_events[0][0]
                    y_ev = sol.y_events[0][0]
                    _, Gamma_now = dlnxi_dx_equilibrium(x_ev, y_ev[0])
                    phase = 1.5
                    x0 = x_ev
                    state["x_switch_QSSA"] = x0
                    y0_eq  = [y_ev[0]]
                    y0_QSSA = [y_ev[0], 0.0]
                    if verbose:
                        print(f"\n  -> Phase 1.5 (QSSA) at x = {x0:.2f}"
                              f" (Gamma_over_H_ann = {Gamma_now:.1e})")
                else:
                    x0 = sol.t[-1];  y0_eq = [sol.y[0, -1]]

            # ----- PHASE 1.5 -----
            elif phase == 1.5:
                events_15 = [cannibal_switch_event]
                if _conv_event_QSSA is not None:
                    events_15.append(_conv_event_QSSA)
                sol = solve_ivp(BEQs_QSSA, (x0, xf), y0_QSSA, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1,
                                events=events_15)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (QSSA): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break
                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend(sol.y[1].tolist())
                muY_all.extend([0.0] * len(sol.t))

                # Check if cannibal event fired (Phase 2 trigger)
                cannibal_fired = (sol.t_events is not None
                                  and len(sol.t_events) > 0
                                  and len(sol.t_events[0]) > 0)
                # Check if convergence event fired
                conv_fired = (_conv_event_QSSA is not None
                              and sol.t_events is not None
                              and len(sol.t_events) > 1
                              and len(sol.t_events[1]) > 0)

                if conv_fired and (not cannibal_fired
                                   or sol.t_events[1][0] <= sol.t_events[0][0]):
                    x_ev = sol.t_events[1][0]
                    converged = True
                    state["x_converged"] = x_ev
                    if verbose:
                        print(f"\n  Converged (1.5, dlnxi_NR, event)"
                              f" at x = {x_ev:.1f}")
                    break

                if cannibal_fired:
                    x_ev  = sol.t_events[0][0]
                    y_ev  = sol.y_events[0][0]
                    Gamma_over_H_can_now = estimate_cannibal_rate(x_ev, y_ev[0], y_ev[1])
                    phase = 2
                    x0 = x_ev
                    state["x_switch_full"] = x0
                    y0_QSSA = y_ev.tolist()
                    y0_full = [y0_QSSA[0], y0_QSSA[1], 0.0]
                    if verbose:
                        print(f"\n  -> Phase 2 (full) at x = {x0:.2f}"
                              f" (Gamma_over_H_can={Gamma_over_H_can_now:.1e},"
                              f" muX={y0_QSSA[1]:.2f})")
                else:
                    x0 = sol.t[-1];  y0_QSSA = sol.y[:, -1].tolist()

            # ----- PHASE 2 -----
            elif phase == 2:
                events_2 = [_conv_event_full] if _conv_event_full is not None else None
                sol = solve_ivp(BEQs_full, (x0, xf), y0_full, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1,
                                events=events_2)
                if not sol.success:
                    if verbose:
                        print(f"  Warning (full): x={x0:.2e}->{xf:.2e}: {sol.message}")
                    break

                # Check if convergence event fired
                conv_fired = (_conv_event_full is not None
                              and sol.t_events is not None
                              and len(sol.t_events) > 0
                              and len(sol.t_events[0]) > 0)

                if conv_fired:
                    x_ev = sol.t_events[0][0]
                    y_ev = sol.y_events[0][0]
                    # Append trajectory up to the event
                    x_all.extend(sol.t.tolist())
                    lnxi_all.extend(sol.y[0].tolist())
                    muX_all.extend(sol.y[1].tolist())
                    muY_all.extend(sol.y[2].tolist())
                    converged = True
                    state["x_converged"] = x_ev
                    if verbose:
                        print(f"\n  Converged (2, dlnxi_NR, event)"
                              f" at x = {x_ev:.1f}")
                    break

                x_all.extend(sol.t.tolist())
                lnxi_all.extend(sol.y[0].tolist())
                muX_all.extend(sol.y[1].tolist())
                muY_all.extend(sol.y[2].tolist())
                x0 = sol.t[-1];  y0_full = sol.y[:, -1].tolist()

                # Fallback discrete convergence check (for dlnYX mode)
                if convergence_mode == 'dlnYX':
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
    #  Joint dark + SM Boltzmann solver (Y decay + inverse decay, 4-phase)
    # =============================================================
    def solve_boltzmann_joint(
        self,
        epsX: float = None,
        use_inverse_decay: bool = True,
        *,
        eps_floor_ratio: float = None,
        xmin: float = 1e-3,
        xmax: float = 2e3,
        n_points: int = 500,
        rtol_value: float = 1e-8,
        atol_value: float = 1e-15,
        Tinf: float = 1e14,
        dlnYX_thresh: float = 1e-2,
        dlnxi_NR_thresh: float = 1e-2,
        Gamma_switch_QSSA: float = 5e9,
        cannibal_switch_full: float = 1.0,
        SMlock_thr: float = 1.0,
        GoH_exit_thresh: float = 1.0,
        min_eps_floor_ratio: float = 0.5,
        phase2_form: str = 'Y',
        return_bg_ICs: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """Joint hidden-sector + SM-bath Boltzmann solver (Y decay + inverse decay).

        Solves the dark sector (X, Y) Boltzmann equations *jointly* with the
        Y <-> SM portal coupling.  Both Y -> SM decay AND SM -> Y inverse
        decay are included in the Y number- and energy-drain terms when
        `use_inverse_decay=True`; the combination satisfies detailed balance
        against the SM bath.  With `use_inverse_decay=False`, only the
        forward decay piece is kept (legacy, breaks detailed balance below
        the thermalisation floor).

        Up to four phases are used:

          Phase 1     : full thermal equilibrium (mu_X = mu_Y = 0).
          Phase 1.5-D : cannibal QSSA (mu_Y = 0 pinned).
          Phase 1.5-S : SM-locked QSSA (xi = 1 and mu_Y = 0 pinned).
                        Entered when Gamma_Y * n_Y^eq(T_SM) / (H * n_Y) > 1.
          Phase 2     : full 3-ODE system in (lnxi, lnY_X, lnY_Y) [Y-form]
                        or (lnxi, mu_X, mu_Y) [mu-form, legacy].

        Use this method instead of `solve_boltzmann_chempot_3phase` when
        epsX is comparable to or above the thermalisation floor
        (epsX >= ~0.5 * eps_therm), where Y <-> SM transfer is dynamically
        important.  Below that floor (the secluded regime), use `_3phase`
        followed by `solve_background` with correct_Y_decay=True.  The
        post-processing step (`solve_background`) is unchanged for both
        workflows.

        Parameters
        ----------
        epsX : float, optional
            Absolute portal coupling.  Mutually exclusive with
            `eps_floor_ratio`; provide exactly one.
        eps_floor_ratio : float, optional (keyword-only)
            Portal coupling expressed as a ratio to the thermalisation floor
            at x_FO=20:  epsX = eps_floor_ratio * eps_therm, with
            eps_therm = sqrt(H(mX/20) / Gamma_Y(epsX=1)).
        use_inverse_decay : bool, default True
            If True, the Y number/energy drains include both Y -> SM decay
            and the SM -> Y inverse-decay piece, in the detailed-balance
            form -Gamma_Y (n_Y - n_Y^eq(T_SM)).  If False, only the forward
            decay -Gamma_Y n_Y is kept (legacy).
        SMlock_thr : float
            Threshold on the SM-lock ratio Gamma_Y * n_Y^eq(T_SM)/(H * n_Y);
            crossing it triggers Phase 1.5-S.
        min_eps_floor_ratio : float
            Validity floor.  Raises ValueError if
            epsX < min_eps_floor_ratio * eps_therm.  Phase 2 numerics
            become unreliable well below this.  Set to 0 to disable.
        phase2_form : {'Y', 'mu'}
            Phase 2 variables: lnY_X/lnY_Y (default, numerically robust at
            small alphaX) or mu_X/mu_Y (legacy).
        """
        # Resolve epsilon: caller provides exactly one of epsX or eps_floor_ratio.
        if (epsX is None) == (eps_floor_ratio is None):
            raise ValueError(
                "Provide exactly one of `epsX` (absolute portal coupling) or "
                "`eps_floor_ratio` (epsX / eps_therm at x_FO=20).")

        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        inc_aX = model.include_antiparticlesX
        inc_aY = model.include_antiparticlesY
        ann_sym = 0.5 if inc_aX else 1.0

        # Compute thermalisation floor (used for ratio mode and validity guard).
        GammaY_norm = model.decay_width_to_SM(1.0)
        if GammaY_norm > 0:
            eps_therm = np.sqrt(cosmo.hubble(mX / 20.0) / GammaY_norm)
        else:
            eps_therm = None

        if eps_floor_ratio is not None:
            if eps_therm is None:
                raise ValueError(
                    "Cannot use `eps_floor_ratio` when the model has zero "
                    "portal decay width (secluded mode); pass `epsX` directly.")
            epsX = eps_floor_ratio * eps_therm

        if (eps_therm is not None and min_eps_floor_ratio > 0
                and epsX < min_eps_floor_ratio * eps_therm):
            raise ValueError(
                f"epsX={epsX:.3e} is below "
                f"min_eps_floor_ratio*eps_therm="
                f"{min_eps_floor_ratio * eps_therm:.3e}; "
                f"the joint solver's Phase 2 is unreliable in this regime. "
                f"Use the secluded workflow (solve_boltzmann_chempot_3phase "
                f"+ solve_background with correct_Y_decay=True) instead.")

        GammaY = model.decay_width_to_SM(epsX)
        if verbose:
            print(f"      [solve_boltzmann_joint] epsX={epsX:.3e}  "
                  f"GammaY={GammaY:.3e} GeV  "
                  f"use_inverse_decay={use_inverse_decay}")

        sv_YXX_XX = model.sigmav2_YXX_to_XX() * ann_sym
        sv_YYX_YX = model.sigmav2_YYX_to_YX()

        # --- Thermodynamics helper ---
        def _thermo(x, lnxi, muX, muY):
            xi = np.exp(np.clip(lnxi, -20, 20))
            T  = mX / x
            TD = xi * T
            s  = cosmo.s_entropy(T)
            H  = cosmo.hubble(T)
            g_tilde = 1.0 + cosmo.dlngSdlnT(T) / 3.0
            zX = mX / TD;  zY = mY / TD
            zY_SM = mY / T
            BX1 = cosmo.Bfac1(zX, mX);  BX2 = cosmo.Bfac2(zX, mX)
            dBX1 = cosmo.dBfac1dT(zX, mX)
            BY1 = cosmo.Bfac1(zY, mY);  BY2 = cosmo.Bfac2(zY, mY)
            dBY1 = cosmo.dBfac1dT(zY, mY)
            lamX = cosmo.lambda_i(zX);  lamY = cosmo.lambda_i(zY)

            lnYeqX = cosmo.ln_Yeq(zX, gX, mX, inc_aX, s)
            lnYeqY = cosmo.ln_Yeq(zY, gY, mY, inc_aY, s)
            lnYX = lnYeqX + muX;  lnYY = lnYeqY + muY
            YX = np.exp(np.clip(lnYX, -700, 700))
            YY = np.exp(np.clip(lnYY, -700, 700))
            nX = YX * s;  nY = YY * s

            # n_Y^eq(T_SM); MB form, consistent with ln_Yeq
            pref = (gY * mY**2 * T / (2.0 * np.pi**2)) * (2.0 if inc_aY else 1.0)
            if zY_SM < 50.0:
                K2_SM = Knu(2, zY_SM)
            else:
                K2_SM = np.sqrt(np.pi / (2.0 * zY_SM)) * np.exp(-zY_SM) \
                        * (1.0 + 15.0 / (8.0 * zY_SM))
            nY_eq_SM = pref * K2_SM
            BY1_SM = cosmo.Bfac1(zY_SM, mY)

            rho_h = max(BX1 * nX + BY1 * nY, 0.0)
            Htot  = np.sqrt(H**2 + rho_h / (3.0 * Mpl**2))

            if zY < 50.0:
                K2_h = Knu(2, zY)
            else:
                K2_h = np.sqrt(np.pi / (2.0 * zY)) * np.exp(-zY) \
                       * (1.0 + 15.0 / (8.0 * zY))
            RY = (T / TD) * (K2_SM / max(K2_h, 1e-300))

            return dict(
                xi=xi, T=T, TD=TD, s=s, H=H, Htot=Htot,
                g_tilde=g_tilde, zX=zX, zY=zY, zY_SM=zY_SM,
                BX1=BX1, BX2=BX2, dBX1=dBX1,
                BY1=BY1, BY2=BY2, dBY1=dBY1,
                lnYeqX=lnYeqX, lnYeqY=lnYeqY, lnRXY=lnYeqX - lnYeqY,
                lnYX=lnYX, lnYY=lnYY, YX=YX, YY=YY, nX=nX, nY=nY,
                nY_eq_SM=nY_eq_SM, RY=RY,
                lamX=lamX, lamY=lamY, BY1_SM=BY1_SM,
                sv_YYY_YY=model.sigmav2_YYY_to_YY(TD),
                sv_XXYY=model.sigmav_XX_to_YY(TD),
            )

        def _dE_decay(BY1, nY, BY1_SM, nY_eq_SM):
            """Y -> SM energy drain per unit volume on the hidden sector.
            Inverse-decay form keeps detailed balance against the SM bath."""
            if use_inverse_decay:
                return GammaY * (BY1 * nY - BY1_SM * nY_eq_SM)
            return GammaY * BY1 * nY

        # === Phase 1: full equilibrium ===
        def dlnxi_dx_eq(x, lnxi):
            th = _thermo(x, lnxi, 0.0, 0.0)
            neqX, neqY = th['nX'], th['nY']
            Num = neqY * th['BY2'] + neqX * th['BX2']
            dneqX = cosmo.dneqdT_MB(th['zX'], gX, mX, inc_aX)
            dneqY = cosmo.dneqdT_MB(th['zY'], gY, mY, inc_aY)
            Den = (th['BX1'] * dneqX + th['dBX1'] * neqX
                   + th['BY1'] * dneqY + th['dBY1'] * neqY)
            if np.abs(Den) < 1e-300:
                Den = np.copysign(1e-300, Den)
            dlnxi = (1.0 / x) * (1.0 - (3.0 + cosmo.dlngSdlnT(th['T']))
                                 * Num / (th['TD'] * Den))
            Gamma_ann = (th['s'] / th['Htot']
                         * (th['sv_XXYY'] * ann_sym) * th['YX'])
            return dlnxi, Gamma_ann

        def rhs_eq(t, y):
            dlnxi, _ = dlnxi_dx_eq(t, y[0])
            return [np.clip(dlnxi, -10.0 / t, 10.0 / t)]

        def ev_ann(t, y):
            _, G = dlnxi_dx_eq(t, y[0])
            return G - Gamma_switch_QSSA
        ev_ann.terminal = True;  ev_ann.direction = -1

        # === Phase 1.5-D: QSSA (muY = 0) ===
        def rhs_QSSA(t, y):
            lnxi, muX = y
            th = _thermo(t, lnxi, muX, 0.0)
            g_tilde = th['g_tilde']
            nX, nY  = th['nX'], th['nY']

            YXeq = np.exp(np.clip(th['lnYeqX'], -700, 700))
            # 2*ann_sym matches the Phase-2 form Pfac*(YX - YXeq^2/YX)
            # = 2*ann_sym * sv * s * g_tilde * YXeq * sinh(muX) / (Htot*t).
            Pfac = (2.0 * ann_sym * th['s'] / th['Htot']
                    * th['sv_XXYY'] * (g_tilde / t))
            A = -Pfac * YXeq * np.sinh(muX) + (th['lamX'] - 3.0 * g_tilde) / t
            B = -th['lamX']

            D_cal = nX * th['dBX1'] + nY * th['dBY1']
            E_cal = th['BX1'] * nX * th['lamX'] + th['BY1'] * nY * th['lamY']
            Dtilde = th['TD'] * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            Num_xi = nX * th['BX2'] + nY * th['BY2']
            C_xi = (1.0 / t) * (1.0 - 3.0 * g_tilde * Num_xi / Dtilde)

            if GammaY > 0:
                dE = _dE_decay(th['BY1'], nY, th['BY1_SM'], th['nY_eq_SM'])
                C_xi -= g_tilde * dE / (t * th['Htot'] * Dtilde)
            D_xi = -th['BX1'] * nX / Dtilde

            denom = 1.0 - B * D_xi
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)

            dlnxi = np.clip((C_xi + D_xi * A) / denom, -10.0 / t, 10.0 / t)
            dmuX  = (A + B * C_xi) / denom
            return [dlnxi, dmuX]

        def cannibal_rate(t, y):
            lnxi, muX = y
            th = _thermo(t, lnxi, muX, 0.0)
            YX, YY = th['YX'], th['YY']
            S_3to2 = (YY**2 * th['sv_YYY_YY']
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            return th['s']**2 / th['Htot'] * S_3to2

        def ev_can(t, y):
            return cannibal_rate(t, y) - cannibal_switch_full
        ev_can.terminal = True;  ev_can.direction = -1

        # === SM-lock ratio ===
        def sm_lock_ratio(x, lnxi, muX, muY):
            if GammaY <= 0:
                return 0.0
            th = _thermo(x, lnxi, muX, muY)
            nY_safe = max(th['nY'], 1e-300)
            return GammaY * th['nY_eq_SM'] / (th['Htot'] * nY_safe)

        def ev_SMlock_in_P15(t, y):
            return sm_lock_ratio(t, y[0], y[1], 0.0) - SMlock_thr
        ev_SMlock_in_P15.terminal = True;  ev_SMlock_in_P15.direction = +1

        # === Phase 1.5-S: SM-locked, xi=1 pinned ===
        # Total g_*, g_*S (SM + hidden) are needed in rhs_SMlock.  Direct
        # evaluation is expensive (4 cosmo.species_thermo calls per rhs
        # eval, ×thousands of solver steps).  Cache: build cubic-spline
        # interpolators of ln(g_S_tot) and ln(g_rho_tot) vs ln(T) over a
        # grid sized to Phase 1.5-S's actual T range; spline derivative
        # gives g_tilde with no finite-difference noise.  Lazy build:
        # only pay the setup cost when Phase 1.5-S is actually entered.
        from scipy.interpolate import CubicSpline
        smlock_cache = {}

        def _build_gtot_splines(x_ev):
            """Construct g_S_tot, g_rho_tot, g_tilde_tot interpolators for
            T in the Phase 1.5-S range.  Called once on entry to '1.5S'."""
            T_min = mX / (xmax * 1.05)
            T_max = mX / (x_ev * 0.95)
            T_grid = np.logspace(np.log10(T_min), np.log10(T_max), 200)
            g_S_arr   = np.empty_like(T_grid)
            g_rho_arr = np.empty_like(T_grid)
            for i, T in enumerate(T_grid):
                rX = cosmo.species_thermo(mX, T, g=gX, statistics="fermion",
                                          include_antiparticles=inc_aX,
                                          what="all")
                rY = cosmo.species_thermo(mY, T, g=gY, statistics="boson",
                                          include_antiparticles=inc_aY,
                                          what="all")
                rho_h  = rX['rho'] + rY['rho']
                rhoP_h = rho_h + rX['P'] + rY['P']
                dg_rho = rho_h  / ((np.pi**2 / 30.0) * T**4)
                dg_S   = rhoP_h / ((2.0 * np.pi**2 / 45.0) * T**4)
                g_rho_arr[i] = cosmo.gstar(T)  + dg_rho
                g_S_arr[i]   = cosmo.gstarS(T) + dg_S
            lnT = np.log(T_grid)
            ln_gS_spl  = CubicSpline(lnT, np.log(g_S_arr))
            ln_gR_spl  = CubicSpline(lnT, np.log(g_rho_arr))
            smlock_cache['ln_gS']  = ln_gS_spl
            smlock_cache['ln_gR']  = ln_gR_spl
            smlock_cache['lnT_lo'] = lnT[0]
            smlock_cache['lnT_hi'] = lnT[-1]

        def rhs_SMlock(t, y):
            muX = y[0]
            T   = mX / t
            lnT = np.log(T)
            # Clamp lnT to the spline's interpolation range; CubicSpline
            # extrapolates by default but we keep it tight to avoid drift.
            lnT_clipped = min(max(lnT, smlock_cache['lnT_lo']),
                              smlock_cache['lnT_hi'])
            ln_gS = smlock_cache['ln_gS'](lnT_clipped)
            ln_gR = smlock_cache['ln_gR'](lnT_clipped)
            g_S_tot   = float(np.exp(ln_gS))
            g_rho_tot = float(np.exp(ln_gR))
            # g_tilde = 1 + (1/3) d ln g_S / d ln T  -- spline analytic derivative
            g_tilde_tot = 1.0 + float(smlock_cache['ln_gS'](lnT_clipped, 1)) / 3.0

            s_tot = (2.0 * np.pi**2 / 45.0) * g_S_tot * T**3
            H_tot = np.sqrt(np.pi**2 * g_rho_tot / 90.0) * T**2 / Mpl

            zX = mX / T
            lnYeqX_tot = cosmo.ln_Yeq(zX, gX, mX, inc_aX, s_tot)
            lnYX = lnYeqX_tot + muX
            YX   = float(np.exp(np.clip(lnYX, -700, 700)))

            sv = model.sigmav_XX_to_YY(T)
            Pfac = s_tot * g_tilde_tot / (H_tot * t) * (sv * ann_sym)
            YXeq = float(np.exp(np.clip(lnYeqX_tot, -700, 700)))
            coll_X = Pfac * (YX - YXeq**2 / YX) if YX > 1e-300 else 0.0

            lamX = cosmo.lambda_i(zX)
            dmuX = -coll_X + (lamX - 3.0 * g_tilde_tot) / t
            return [dmuX]

        # === Phase 2: full system (mu-form) ===
        def rhs_full(t, y):
            lnxi, muX, muY = y
            th = _thermo(t, lnxi, muX, muY)
            g_tilde = th['g_tilde']
            YX, YY  = th['YX'], th['YY']
            nX, nY  = th['nX'], th['nY']

            Pfac = th['s'] * g_tilde / (th['Htot'] * t) * (th['sv_XXYY'] * ann_sym)
            RXY2_YY2 = np.exp(np.clip(2.0 * th['lnRXY'] + 2.0 * th['lnYY'], -700, 700))

            coll_X = Pfac * (YX - RXY2_YY2 / YX) if YX > 1e-300 else 0.0
            A = -coll_X + (th['lamX'] - 3.0 * g_tilde) / t
            B = -th['lamX']

            coll_Y_2to2 = Pfac * (YX**2 - RXY2_YY2) / YY if YY > 1e-300 else 0.0
            S_3to2 = (YY**2 * th['sv_YYY_YY']
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            Gamma_can = th['s']**2 * g_tilde / (th['Htot'] * t) * S_3to2
            one_m = 1.0 if muY > 50 else (-np.exp(-muY) if muY < -50
                                          else 1.0 - np.exp(-muY))

            C = coll_Y_2to2 - Gamma_can * one_m + (th['lamY'] - 3.0 * g_tilde) / t
            decay_active = GammaY > 0
            if decay_active:
                if use_inverse_decay:
                    emu = np.exp(-min(muY, 50.0)) if muY > -50 else np.exp(50.0)
                    C -= GammaY * g_tilde / (t * th['Htot']) * (1.0 - th['RY'] * emu)
                else:
                    C -= GammaY * g_tilde / (t * th['Htot'])
            D = -th['lamY']

            D_cal  = nX * th['dBX1'] + nY * th['dBY1']
            E_cal  = th['BX1'] * nX * th['lamX'] + th['BY1'] * nY * th['lamY']
            Dtilde = th['TD'] * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            rho_plus_P = nX * th['BX2'] + nY * th['BY2']
            E = (1.0 / t) * (1.0 - 3.0 * g_tilde * rho_plus_P / Dtilde)
            if decay_active:
                dE = _dE_decay(th['BY1'], nY, th['BY1_SM'], th['nY_eq_SM'])
                E -= g_tilde * dE / (t * th['Htot'] * Dtilde)
            F = -th['BX1'] * nX / Dtilde
            G = -th['BY1'] * nY / Dtilde

            denom = 1.0 - F * B - G * D
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)

            dlnxi = np.clip((E + F * A + G * C) / denom, -10.0 / t, 10.0 / t)
            dmuX  = A + B * dlnxi
            dmuY  = C + D * dlnxi

            nonlocal x_FO
            if muX > 0.05 and x_FO is None:
                x_FO = t

            return [dlnxi, dmuX, dmuY]

        # === Phase 2 Y-form: variables (lnxi, lnYX, lnYY) ===
        def rhs_full_Y(t, y):
            lnxi, lnYX, lnYY = y
            th0 = _thermo(t, lnxi, 0.0, 0.0)
            g_tilde = th0['g_tilde']
            s, TD   = th0['s'], th0['TD']
            lamX, lamY = th0['lamX'], th0['lamY']
            BX1, BX2, dBX1 = th0['BX1'], th0['BX2'], th0['dBX1']
            BY1, BY2, dBY1 = th0['BY1'], th0['BY2'], th0['dBY1']
            BY1_SM, nY_eq_SM = th0['BY1_SM'], th0['nY_eq_SM']

            YX = float(np.exp(np.clip(lnYX, -700, 700)))
            YY = float(np.exp(np.clip(lnYY, -700, 700)))
            nX = YX * s;  nY = YY * s

            YYeq    = float(np.exp(np.clip(th0['lnYeqY'], -700, 700)))
            YYeq_SM = nY_eq_SM / s

            rho_h = max(BX1 * nX + BY1 * nY, 0.0)
            Htot  = np.sqrt(th0['H']**2 + rho_h / (3.0 * Mpl**2))

            Pfac = s * g_tilde / (Htot * t) * (th0['sv_XXYY'] * ann_sym)
            # Single-clipped reverse-rate exp to preserve tiny YXeq/YYeq ratio.
            rev = float(np.exp(np.clip(
                2.0 * th0['lnRXY'] + 2.0 * lnYY, -700, 700)))
            coll_X      = Pfac * (YX - rev / YX) if YX > 1e-300 else 0.0
            coll_Y_2to2 = Pfac * (YX ** 2 - rev) / YY if YY > 1e-300 else 0.0

            S_3to2 = (YY ** 2 * th0['sv_YYY_YY']
                      + YY * YX * sv_YYX_YX
                      + YX ** 2 * sv_YXX_XX)
            Gamma_can = s ** 2 * g_tilde / (Htot * t) * S_3to2
            ratio_D = YYeq / YY if YY > 1e-300 else 1e100
            ratio_D = min(ratio_D, 1e100)
            Gamma_can_term = Gamma_can * (1.0 - ratio_D)

            decay_active = GammaY > 0
            if decay_active:
                if use_inverse_decay:
                    ratio_SM = YYeq_SM / YY if YY > 1e-300 else 1e100
                    ratio_SM = min(ratio_SM, 1e100)
                    decay_term = GammaY * g_tilde / (t * Htot) * (1.0 - ratio_SM)
                else:
                    decay_term = GammaY * g_tilde / (t * Htot)
            else:
                decay_term = 0.0

            D_cal  = nX * dBX1 + nY * dBY1
            E_cal  = BX1 * nX * lamX + BY1 * nY * lamY
            Dtilde = TD * D_cal + E_cal
            if np.abs(Dtilde) < 1e-300:
                Dtilde = np.copysign(1e-300, Dtilde)

            rho_plus_P = nX * BX2 + nY * BY2
            E = (1.0 / t) * (1.0 - 3.0 * g_tilde * rho_plus_P / Dtilde)
            if decay_active:
                dE = _dE_decay(BY1, nY, BY1_SM, nY_eq_SM)
                E -= g_tilde * dE / (t * Htot * Dtilde)

            F = -BX1 * nX / Dtilde
            G = -BY1 * nY / Dtilde
            B = -lamX
            D = -lamY

            A = -coll_X + (lamX - 3.0 * g_tilde) / t
            C = (coll_Y_2to2 - Gamma_can_term - decay_term
                 + (lamY - 3.0 * g_tilde) / t)

            denom = 1.0 - F * B - G * D
            if np.abs(denom) < 1e-300:
                denom = np.copysign(1e-300, denom)
            dlnxi = np.clip((E + F * A + G * C) / denom, -10.0 / t, 10.0 / t)

            dlnYX = -coll_X
            dlnYY = coll_Y_2to2 - Gamma_can_term - decay_term

            muX_recon = lnYX - th0['lnYeqX']
            nonlocal x_FO
            if muX_recon > 0.05 and x_FO is None:
                x_FO = t

            return [dlnxi, dlnYX, dlnYY]

        # === Driver ===
        lnxi0 = (1.0 / 3.0) * np.log(
            np.round(cosmo.gstarS(mX / xmin) / cosmo.gstarS(Tinf), 4)
        ) + np.log(model.xi_ini)

        x_all = [xmin];  lnxi_all = [lnxi0]
        muX_all = [0.0]; muY_all = [0.0]
        phase = 1
        x0 = xmin
        y_eq = [lnxi0];  y_QSSA = None;  y_SMlock = None;  y_full = None
        x_FO = None
        x_switch_QSSA = None
        x_to_SMlock = None
        x_to_P2 = None
        x_converged = None
        converged = False

        xf_tail = list(np.arange(120, xmax + 1, 10))
        xfList = np.concatenate([
            np.array([1e-2, 1e-1, 0.99]),
            np.linspace(1, 10, 20),
            np.arange(10.5, 100, 1),
            np.array(xf_tail),
        ])
        xfList = np.unique(xfList[xfList > xmin])

        for xf in xfList:
            if x0 >= xf:
                continue
            xs = np.linspace(x0 * 1.001, xf * 0.999, n_points) if x0 > 1 \
                 else np.logspace(np.log10(x0 * 1.001),
                                  np.log10(xf * 0.999), n_points)

            if phase == 1:
                sol = solve_ivp(rhs_eq, (x0, xf), y_eq, t_eval=xs,
                                rtol=rtol_value, atol=1e-12, method='Radau',
                                max_step=1, events=[ev_ann])
                if not sol.success:
                    break
                x_all.extend(sol.t);  lnxi_all.extend(sol.y[0])
                muX_all.extend([0.0] * len(sol.t))
                muY_all.extend([0.0] * len(sol.t))

                fired = sol.t_events and len(sol.t_events[0]) > 0
                if fired:
                    x_ev = sol.t_events[0][0];  y_ev = sol.y_events[0][0]
                    x0 = x_ev
                    x_switch_QSSA = x_ev
                    ratio = sm_lock_ratio(x_ev, y_ev[0], 0.0, 0.0)
                    if use_inverse_decay and ratio > SMlock_thr:
                        phase = '1.5S'
                        y_SMlock = [0.0]
                        x_to_SMlock = x_ev
                        if not smlock_cache:
                            _build_gtot_splines(x_ev)
                        if verbose:
                            print(f"      Phase 1 -> 1.5S at x={x_ev:.2f}"
                                  f"  (SM-lock ratio = {ratio:.2e})")
                    else:
                        phase = 1.5
                        y_QSSA = [y_ev[0], 0.0]
                        if verbose:
                            print(f"      Phase 1 -> 1.5 at x={x_ev:.2f}"
                                  f"  (SM-lock ratio = {ratio:.2e})")
                else:
                    x0 = sol.t[-1];  y_eq = [sol.y[0, -1]]

            elif phase == 1.5:
                events_15 = [ev_can]
                if use_inverse_decay and GammaY > 0:
                    events_15.append(ev_SMlock_in_P15)
                sol = solve_ivp(rhs_QSSA, (x0, xf), y_QSSA, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1, events=events_15)
                if not sol.success:
                    break
                x_all.extend(sol.t)
                lnxi_all.extend(sol.y[0]);  muX_all.extend(sol.y[1])
                muY_all.extend([0.0] * len(sol.t))

                for muX_i in sol.y[1]:
                    if muX_i > 0.05 and x_FO is None:
                        x_FO = sol.t[np.argmax(sol.y[1] > 0.05)]
                        break

                can_fired = len(sol.t_events[0]) > 0
                smlock_fired = (use_inverse_decay and GammaY > 0
                                and len(sol.t_events) > 1
                                and len(sol.t_events[1]) > 0)

                if smlock_fired and (not can_fired
                        or sol.t_events[1][0] < sol.t_events[0][0]):
                    x_ev = sol.t_events[1][0];  y_ev = sol.y_events[1][0]
                    phase = '1.5S'
                    x0 = x_ev
                    y_SMlock = [y_ev[1]]
                    x_to_SMlock = x_ev
                    if not smlock_cache:
                        _build_gtot_splines(x_ev)
                    if verbose:
                        print(f"      Phase 1.5 -> 1.5S at x={x_ev:.2f}"
                              f"  muX={y_ev[1]:.3f}")
                elif can_fired:
                    x_ev = sol.t_events[0][0];  y_ev = sol.y_events[0][0]
                    phase = 2;  x0 = x_ev
                    x_to_P2 = x_ev
                    if phase2_form == 'Y':
                        th0 = _thermo(x_ev, y_ev[0], 0.0, 0.0)
                        y_full = [y_ev[0],
                                  y_ev[1] + th0['lnYeqX'],
                                  th0['lnYeqY']]
                    else:
                        y_full = [y_ev[0], y_ev[1], 0.0]
                    if verbose:
                        print(f"      Phase 1.5 -> 2 at x={x_ev:.2f}"
                              f"  lnxi={y_ev[0]:.3f}  muX={y_ev[1]:.3f}"
                              f"  (form='{phase2_form}')")
                else:
                    x0 = sol.t[-1];  y_QSSA = sol.y[:, -1].tolist()

            elif phase == '1.5S':
                # Phase 1.5-S is the same physics as solve_boltzmann_thermSM
                # (Lee-Weinberg with total g_*).  Mirror its xmax=500 so
                # the residual annihilation tail (Y drops ~30% from x_FO
                # to x~500) is captured.  Cap further integration here.
                PHASE_15S_XMAX = 500.0
                if x0 >= PHASE_15S_XMAX:
                    converged = True
                    x_converged = x0
                    if verbose:
                        print(f"      Phase 1.5S complete at x={x0:.1f}")
                    break
                xf_local = min(xf, PHASE_15S_XMAX)
                if x0 * 1.001 >= xf_local * 0.999:
                    continue
                xs_local = (np.linspace(x0 * 1.001, xf_local * 0.999, n_points)
                            if x0 > 1
                            else np.logspace(np.log10(x0 * 1.001),
                                             np.log10(xf_local * 0.999),
                                             n_points))

                # atol=1e-10 here: with the cubic-spline-based g_tilde
                # the FP-noise issue is mostly gone, but keeping the
                # looser tol is harmless and matches Phase 1.5-D.
                sol = solve_ivp(rhs_SMlock, (x0, xf_local), y_SMlock,
                                t_eval=xs_local,
                                rtol=rtol_value, atol=1e-10,
                                method='Radau', max_step=1)
                if not sol.success:
                    if verbose:
                        print(f"      Phase 1.5S FAILED at x={x0:.3f}: "
                              f"{sol.message}")
                    break
                x_all.extend(sol.t)
                lnxi_all.extend([0.0] * len(sol.t))
                muX_all.extend(sol.y[0])
                muY_all.extend([0.0] * len(sol.t))

                for muX_i in sol.y[0]:
                    if muX_i > 0.05 and x_FO is None:
                        x_FO = sol.t[np.argmax(sol.y[0] > 0.05)]
                        break

                x0 = sol.t[-1];  y_SMlock = [sol.y[0, -1]]

            elif phase == 2:
                rhs_P2 = rhs_full_Y if phase2_form == 'Y' else rhs_full

                MUX_EXIT_THRESH = 1.0
                if GammaY > 0:
                    if phase2_form == 'Y':
                        th_chk = _thermo(x0, y_full[0], 0.0, 0.0)
                        muX_chk = y_full[1] - th_chk['lnYeqX']
                        muY_chk = y_full[2] - th_chk['lnYeqY']
                    else:
                        muX_chk = y_full[1]
                        muY_chk = y_full[2]
                    th_start = _thermo(x0, y_full[0], muX_chk, muY_chk)
                    GoH_start   = GammaY / th_start['Htot']
                    rho_ratio_0 = (th_start['Htot'] / th_start['H']) ** 2 - 1.0
                    if (GoH_start >= GoH_exit_thresh or rho_ratio_0 >= 1.0) \
                            and muX_chk > MUX_EXIT_THRESH:
                        converged = True
                        x_converged = x0
                        if verbose:
                            reason = (f'Gamma/H>={GoH_exit_thresh:g}'
                                      if GoH_start >= GoH_exit_thresh
                                      else 'rho_h>=rho_SM')
                            print(f"      Phase 2 exit ({reason}, initial) "
                                  f"at x={x0:.2f}  muX={muX_chk:.2f}")
                        break

                def _state_thermo(t, y):
                    lnxi_now = y[0]
                    if phase2_form == 'Y':
                        th00 = _thermo(t, lnxi_now, 0.0, 0.0)
                        muX_now = y[1] - th00['lnYeqX']
                        muY_now = y[2] - th00['lnYeqY']
                    else:
                        muX_now, muY_now = y[1], y[2]
                    return _thermo(t, lnxi_now, muX_now, muY_now)

                def _muX(t, y):
                    if phase2_form == 'Y':
                        th00 = _thermo(t, y[0], 0.0, 0.0)
                        return y[1] - th00['lnYeqX']
                    return y[1]

                def ev_GoH(t, y):
                    if _muX(t, y) < MUX_EXIT_THRESH:
                        return -1.0
                    th_e = _state_thermo(t, y)
                    return GammaY / th_e['Htot'] - GoH_exit_thresh
                ev_GoH.terminal  = (GammaY > 0)
                ev_GoH.direction = +1

                def ev_rho(t, y):
                    if _muX(t, y) < MUX_EXIT_THRESH:
                        return -1.0
                    th_e = _state_thermo(t, y)
                    ratio = (th_e['Htot'] / th_e['H']) ** 2 - 1.0
                    return ratio - 1.0
                ev_rho.terminal  = True
                ev_rho.direction = +1

                def _conv_event_full(t, y):
                    if x_FO is None or t < 2.0 * x_FO:
                        return 1.0
                    rhs_now = rhs_full_Y if phase2_form == 'Y' else rhs_full
                    dlnxi_num = rhs_now(t, list(y))[0]
                    g_t = _thermo(t, y[0], 0.0, 0.0)['g_tilde']
                    dlnxi_NR = (1.0 - 2.0 * g_t) / t
                    return abs(dlnxi_num - dlnxi_NR) * t - dlnxi_NR_thresh
                _conv_event_full.terminal  = True
                _conv_event_full.direction = -1

                events_P2 = [ev_rho]
                if GammaY > 0:
                    events_P2.append(ev_GoH)
                events_P2.append(_conv_event_full)

                sol = solve_ivp(rhs_P2, (x0, xf), y_full, t_eval=xs,
                                rtol=rtol_value, atol=atol_value,
                                method='Radau', max_step=1, events=events_P2)
                if not sol.success:
                    break
                x_all.extend(sol.t)
                lnxi_all.extend(sol.y[0])
                if phase2_form == 'Y':
                    for i_pt, t_pt in enumerate(sol.t):
                        th0 = _thermo(t_pt, sol.y[0, i_pt], 0.0, 0.0)
                        muX_all.append(sol.y[1, i_pt] - th0['lnYeqX'])
                        muY_all.append(sol.y[2, i_pt] - th0['lnYeqY'])
                else:
                    muX_all.extend(sol.y[1])
                    muY_all.extend(sol.y[2])
                x0 = sol.t[-1];  y_full = sol.y[:, -1].tolist()

                if sol.t_events:
                    conv_idx = 2 if GammaY > 0 else 1
                    if len(sol.t_events[0]) > 0:
                        converged = True;  x_converged = x0
                        if verbose:
                            print(f"      Phase 2 exit (rho_h>rho_SM) "
                                  f"at x={x0:.2f}")
                        break
                    if (GammaY > 0 and len(sol.t_events) > 1
                            and len(sol.t_events[1]) > 0):
                        converged = True;  x_converged = x0
                        if verbose:
                            print(f"      Phase 2 exit (Gamma/H=1) "
                                  f"at x={x0:.2f}")
                        break
                    if (len(sol.t_events) > conv_idx
                            and len(sol.t_events[conv_idx]) > 0):
                        converged = True;  x_converged = x0
                        if verbose:
                            print(f"      Phase 2 exit (dlnxi_NR) "
                                  f"at x={x0:.2f}")
                        break

                if x_FO is not None and x0 > 2.0 * x_FO and len(x_all) >= 50:
                    x_prev = x_all[-50]
                    lnYX_prev = muX_all[-50] + cosmo.ln_Yeq(
                        mX / (np.exp(lnxi_all[-50]) * mX / x_prev),
                        gX, mX, inc_aX, cosmo.s_entropy(mX / x_prev))
                    lnYX_now = muX_all[-1] + cosmo.ln_Yeq(
                        mX / (np.exp(lnxi_all[-1]) * mX / x0),
                        gX, mX, inc_aX, cosmo.s_entropy(mX / x0))
                    if np.log(x0 / x_prev) > 0.5:
                        rate = abs((lnYX_now - lnYX_prev) / np.log(x0 / x_prev))
                        if rate < dlnYX_thresh:
                            converged = True;  x_converged = x0
                            if verbose:
                                print(f"      Phase 2 exit (dlnYX) "
                                      f"at x={x0:.2f}")
                            break

        # === Post-process ===
        x_arr   = np.array(x_all)
        xi_arr  = np.exp(np.clip(np.array(lnxi_all), -20, 20))
        muX_arr = np.array(muX_all)
        muY_arr = np.array(muY_all)
        T_arr   = mX / x_arr;  TD_arr = xi_arr * T_arr
        zX_arr  = mX / TD_arr;  zY_arr = mY / TD_arr

        lnYX_arr   = np.zeros_like(x_arr);  lnYY_arr   = np.zeros_like(x_arr)
        lnYXeq_arr = np.zeros_like(x_arr);  lnYYeq_arr = np.zeros_like(x_arr)
        for i in range(len(x_arr)):
            si = cosmo.s_entropy(T_arr[i])
            lnYXeq_arr[i] = cosmo.ln_Yeq(zX_arr[i], gX, mX, inc_aX, si)
            lnYYeq_arr[i] = cosmo.ln_Yeq(zY_arr[i], gY, mY, inc_aY, si)
            lnYX_arr[i] = muX_arr[i] + lnYXeq_arr[i]
            lnYY_arr[i] = muY_arr[i] + lnYYeq_arr[i]
        YX_arr   = np.exp(np.clip(lnYX_arr,   -700, 700))
        YY_arr   = np.exp(np.clip(lnYY_arr,   -700, 700))
        YXeq_arr = np.exp(np.clip(lnYXeq_arr, -700, 700))
        YYeq_arr = np.exp(np.clip(lnYYeq_arr, -700, 700))

        if verbose:
            print(f"      phase_final={phase}  x_final={x_arr[-1]:.2f}  "
                  f"muX_final={muX_arr[-1]:.3f}  "
                  f"muY_final={muY_arr[-1]:.3f}  xi_final={xi_arr[-1]:.3e}")

        result = {
            'x': x_arr, 'xi': xi_arr,
            'YX': YX_arr, 'YY': YY_arr,
            'YXeq': YXeq_arr, 'YYeq': YYeq_arr,
            'mubar_X': muX_arr, 'mubar_Y': muY_arr,
            'YX_relic': YX_arr[-1], 'YY_final': YY_arr[-1],
            'xi_final': xi_arr[-1],
            'x_FO': x_FO,
            'x_switch_QSSA': x_switch_QSSA,
            'x_switch_full': x_to_P2,
            'x_to_SMlock': x_to_SMlock,
            'x_to_P2': x_to_P2,
            'x_converged': x_converged,
            'converged': converged,
            'phase_final': phase,
            'GammaY': GammaY,
            'eps_therm': eps_therm,
        }
        if return_bg_ICs:
            result['bg_ICs'] = {
                'x': x_arr, 'xi': xi_arr, 'YX': YX_arr, 'YY': YY_arr,
                'T': T_arr, 'TD': TD_arr,
            }
        return result

    # =============================================================
    #  SM-thermalized (Lee-Weinberg) Boltzmann solver
    # =============================================================
    def solve_boltzmann_thermSM(
        self,
        *,
        xmin: float = 1.0,
        xmax: float = 500.0,
        n_points: int = 500,
        rtol_value: float = 1e-8,
        atol_value: float = 1e-15,
        include_hidden_in_gstar: bool = True,
        return_bg_ICs: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """SM-thermalized (Lee-Weinberg) relic-abundance solver.

        Solves the standard WIMP-like Lee-Weinberg equation in the limit
        where the dark sector (X, Y) is in full thermal equilibrium with
        the SM bath throughout (xi = T_h/T_SM = 1).  X annihilates to Y
        via 2 -> 2; Y stays in chemical equilibrium with the SM via the
        portal coupling (mu_Y = 0, Y_Y = Y_Y^eq(T_SM) implicit).  The
        ODE tracks only the comoving X yield Y_X.

        Total g_*, g_*S (SM + hidden) are used in s, H, g_tilde, matching
        Phase 1.5-S of `solve_boltzmann_joint` so the two agree at large
        eps_floor_ratio.  Use this method when the portal coupling is
        large enough to keep Y locked to the SM bath (eps >> eps_therm).

        Returns
        -------
        dict with keys
          x        : x-grid
          T        : SM temperature mX/x
          YX       : tracked yield (total Y_X = Y_X + Y_Xbar for Dirac)
          YXeq    : equilibrium yield at T (using total g_*S)
          YX_relic : YX[-1]
          x_FO     : first x where YX > 1.05 * YXeq (None if never)
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        inc_aX = model.include_antiparticlesX
        inc_aY = model.include_antiparticlesY
        ann_sym = 0.5 if inc_aX else 1.0

        def _hidden_gstar_contrib(T):
            rX = cosmo.species_thermo(mX, T, g=gX, statistics="fermion",
                                      include_antiparticles=inc_aX, what="all")
            rY = cosmo.species_thermo(mY, T, g=gY, statistics="boson",
                                      include_antiparticles=inc_aY, what="all")
            rho_h  = rX['rho'] + rY['rho']
            rhoP_h = rho_h + rX['P'] + rY['P']
            dg_rho = rho_h  / ((np.pi**2 / 30.0) * T**4)
            dg_S   = rhoP_h / ((2.0 * np.pi**2 / 45.0) * T**4)
            return dg_rho, dg_S

        def gstar_tot(T):
            if not include_hidden_in_gstar:
                return cosmo.gstar(T), cosmo.gstarS(T)
            dg_rho, dg_S = _hidden_gstar_contrib(T)
            return cosmo.gstar(T) + dg_rho, cosmo.gstarS(T) + dg_S

        def g_tilde_tot(T):
            # g_tilde = 1 + (1/3) d ln g_S / d ln T  (central log-step finite diff)
            d = 0.01
            _, gS_p = gstar_tot(T * np.exp(d))
            _, gS_m = gstar_tot(T * np.exp(-d))
            return 1.0 + (np.log(gS_p) - np.log(gS_m)) / (6.0 * d)

        def Yeq_at(T):
            _, gS = gstar_tot(T)
            s = (2.0 * np.pi**2 / 45.0) * gS * T**3
            zX = mX / T
            return float(np.exp(np.clip(
                cosmo.ln_Yeq(zX, gX, mX, inc_aX, s), -700, 700)))

        def rhs(x, y):
            T = mX / x
            gR, gS = gstar_tot(T)
            s = (2.0 * np.pi**2 / 45.0) * gS * T**3
            H = np.sqrt(np.pi**2 * gR / 90.0) * T**2 / Mpl
            gt = g_tilde_tot(T)
            Yeq = Yeq_at(T)
            sv = model.sigmav_XX_to_YY(T)
            Pfac = s * gt / (H * x) * sv * ann_sym
            Y = y[0]
            return [-Pfac * (Y * Y - Yeq * Yeq)]

        # Initial condition at xmin: equilibrium yield.
        Y0 = Yeq_at(mX / xmin)

        xs = np.logspace(np.log10(xmin * 1.001),
                         np.log10(xmax * 0.999), n_points)
        sol = solve_ivp(rhs, (xmin, xmax), [Y0], t_eval=xs,
                        rtol=rtol_value, atol=atol_value,
                        method='Radau', max_step=1.0)
        if not sol.success:
            raise RuntimeError(f"solve_boltzmann_thermSM failed: {sol.message}")

        x_arr = sol.t
        Y_arr = sol.y[0]
        T_arr = mX / x_arr
        Yeq_arr = np.array([Yeq_at(T) for T in T_arr])

        x_FO = None
        for i in range(len(x_arr)):
            if Y_arr[i] > 1.05 * Yeq_arr[i]:
                x_FO = float(x_arr[i])
                break

        if verbose:
            print(f"      [solve_boltzmann_thermSM] mX={mX:.3e} GeV  "
                  f"YX_relic={Y_arr[-1]:.4e}  x_FO={x_FO}")

        result = {
            'x':        x_arr,
            'T':        T_arr,
            'YX':       Y_arr,
            'YXeq':     Yeq_arr,
            'YX_relic': float(Y_arr[-1]),
            'x_FO':     x_FO,
        }
        if return_bg_ICs:
            # Y is in eq with SM at T_SM throughout; xi = 1 pinned.
            zY_arr = mY / T_arr
            YYeq_arr = np.array([
                float(np.exp(np.clip(
                    cosmo.ln_Yeq(zY_arr[i], gY, mY, inc_aY,
                                 (2.0 * np.pi**2 / 45.0)
                                 * gstar_tot(T_arr[i])[1] * T_arr[i]**3),
                    -700, 700)))
                for i in range(len(x_arr))])
            result['bg_ICs'] = {
                'x':  x_arr,
                'xi': np.ones_like(x_arr),
                'YX': Y_arr,
                'YY': YYeq_arr,
                'T':  T_arr,
                'TD': T_arr,
            }
        return result

    # =============================================================
    #  [Legacy] Boltzmann solver — chemical-potential formulation (2-phase)
    # =============================================================
    def solve_boltzmann_chempot_2phase(
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
        [Legacy] Two-phase Boltzmann solver using chemical-potential variables
        (ln xi, mubar_X, mubar_Y) with Y_i = Y_i^eq exp(mubar_i).

        Phase 1 (equilibrium): mubar = 0, evolve only ln xi
        Phase 2 (full):        evolve all three via decoupled A,B,C,D,E,F,G

        Note: Superseded by solve_boltzmann_chempot_3phase.  Retained for
        testing and comparison only.
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        include_antiparticlesX = model.include_antiparticlesX
        include_antiparticlesY = model.include_antiparticlesY

        ann_sym = 0.5 if include_antiparticlesX else 1.0

        # x-grid
        xfList = np.concatenate([
            np.array([1e-2, 1e-1, 0.99]),
            np.arange(1, 100, 1),
            np.array([2e2, 5e2, 1e3, 2e3, 5e3, 1e4]),
        ])

        # Temperature-independent cross sections
        sv_YXX_XX = model.sigmav2_YXX_to_XX() * ann_sym
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
            Gamma_over_H_ann = (th['s'] / th['Htot']
                         * (th['sv_XXYY'] * ann_sym) * th['YX'])
            return dlnxi, Gamma_over_H_ann

        def BEQs_eq(t, y):
            dlnxi, _ = dlnxi_dx_equilibrium(t, y[0])
            return [np.clip(dlnxi, -10.0 / t, 10.0 / t)]

        # --- Phase 2: full system ---

        def compute_derivatives_full(x, lnxi, muX, muY):
            th = _compute_thermo(x, lnxi, muX, muY)
            g_tilde = th['g_tilde']
            YX, YY  = th['YX'], th['YY']
            nX, nY  = th['nX'], th['nY']

            Pfac = th['s'] * g_tilde / (th['Htot'] * x) * (th['sv_XXYY'] * ann_sym)
            RXY2_YY2 = np.exp(np.clip(2.0 * th['lnRXY'] + 2.0 * th['lnYY'], -700, 700))

            coll_X = Pfac * (YX - RXY2_YY2 / YX) if YX > 1e-300 else 0.0
            A = -coll_X + (th['lamX'] - 3.0 * g_tilde) / x
            B = -th['lamX']

            coll_Y_2to2 = Pfac * (YX**2 - RXY2_YY2) / YY if YY > 1e-300 else 0.0

            S_3to2 = (YY**2 * th['sv_YYY_YY']
                      + YY * YX * sv_YYX_YX
                      + YX**2 * sv_YXX_XX)
            Gamma_over_H_can = th['s']**2 * g_tilde / (th['Htot'] * x) * S_3to2

            if muY > 50:
                one_m_expnmuY = 1.0
            elif muY < -50:
                one_m_expnmuY = -np.exp(-muY)
            else:
                one_m_expnmuY = 1.0 - np.exp(-muY)

            C = coll_Y_2to2 - Gamma_over_H_can * one_m_expnmuY + (th['lamY'] - 3.0 * g_tilde) / x
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
            np.round(cosmo.gstarS(mX / xmin) / cosmo.gstarS(Tinf), 4)) + np.log(model.xi_ini)
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
            print(f"  Phase 1->2 switch: Gamma_over_H_ann < {Gamma_switch_threshold:.0e}")
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
                              f" (Gamma_over_H_ann = {Gamma_now:.1e})")

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
        correct_Y_decay: bool = False,
        traj_has_decay: bool = False,
    ) -> Dict[str, Any]:
        """
        Background evolution: late-time Y decay with entropy injection.

        Evolves three log-scaled variables in e-folds N = ln(a/a_0):
            y[0] = ln(T/T_0),  y[1] = ln(rhoY/rhoY_0),  y[2] = ln(rhoX/rhoX_0)

        If correct_Y_decay is True, the initial conditions are corrected
        for the SM-entropy injection that occurred during the Boltzmann
        evolution.  The integrand is the net Hidden -> SM energy flux

            d(ln S)/dt = Gamma_Y / T_SM * [ B_1(m_Y/T_D)  * Y_Y_phys
                                         - B_1(m_Y/T_SM) * Y_Y^eq(T_SM) ]

        where B_1 = rho/n and the inverse-decay term vanishes identically
        when Y is in SM equilibrium.

        traj_has_decay toggles the interpretation of Y_Y in the trajectory:

          * traj_has_decay = False (default, legacy secluded workflow):
            the Boltzmann integration treated Y as stable, so
            Y_Y_trajectory is the counterfactual "no-decay" yield.  The
            physical Y yield is tracked via f_survive = exp(-int Gamma dt)
            and the forward term is weighted by it.  The handoff YY0 is
            multiplied by f_survive / S_ratio_corr.

          * traj_has_decay = True (prototype workflow, or any trajectory
            produced with Gamma_Y already in the Boltzmann ODE):
            Y_Y_trajectory is already the physical yield.  f_survive is
            set to 1 in the integrand and the handoff YY0 is just
            divided by S_ratio_corr.
        """
        cosmo = self.cosmo
        model = self.model
        mX, mY = model.mX, model.mY
        gX, gY = model.gX, model.gY
        include_antiparticlesX = model.include_antiparticlesX
        include_antiparticlesY = model.include_antiparticlesY
        ann_sym_bg = 0.5 if include_antiparticlesX else 1.0

        x0  = bg_ICs['x'][-1]
        xi0 = bg_ICs['xi'][-1]
        YX0 = bg_ICs['YX'][-1]
        YY0 = bg_ICs['YY'][-1]
        T0  = bg_ICs['T'][-1]
        TD0 = bg_ICs['TD'][-1]

        GammaY  = model.decay_width_to_SM(epsX)

        S_ratio_corr = 1.0  # Boltzmann-phase entropy correction (1 if disabled)
        f_survive    = 1.0  # Fraction of initial Y surviving in the "no-decay
                            # in Boltzmann" interpretation; left at 1 if
                            # traj_has_decay=True.
        # --- Correct for SM-entropy injection during Boltzmann evolution ---
        if correct_Y_decay and GammaY > 0:
            x_traj  = bg_ICs['x']
            T_traj  = bg_ICs['T']
            xi_traj = bg_ICs['xi']
            YX_traj = bg_ICs['YX']
            YY_traj = bg_ICs['YY']

            gY         = model.gY
            inc_anti_Y = model.include_antiparticlesY

            ln_S_corr = 0.0
            for i in range(1, len(x_traj)):
                x_mid   = 0.5 * (x_traj[i]  + x_traj[i-1])
                T_mid   = 0.5 * (T_traj[i]  + T_traj[i-1])
                xi_mid  = 0.5 * (xi_traj[i] + xi_traj[i-1])
                YY_mid  = 0.5 * (YY_traj[i] + YY_traj[i-1])
                TD_mid  = xi_mid * T_mid

                s_mid = cosmo.s_entropy(T_mid)
                H_sm  = cosmo.hubble(T_mid)
                rho_h = (mX * 0.5*(YX_traj[i]+YX_traj[i-1])
                         + mY * YY_mid) * s_mid
                Htot  = np.sqrt(H_sm**2 + rho_h / (3.0 * Mpl**2))
                g_t   = 1.0 + cosmo.dlngSdlnT(T_mid) / 3.0
                dx    = x_traj[i] - x_traj[i-1]
                dt    = dx / (x_mid * Htot * g_t)

                zY_D    = mY / max(TD_mid, 1e-300)
                zY_SM   = mY / max(T_mid,  1e-300)
                B1_D    = cosmo.Bfac1(zY_D,  mY)
                B1_SM   = cosmo.Bfac1(zY_SM, mY)
                lnYYeq  = cosmo.ln_Yeq(zY_SM, gY, mY, inc_anti_Y, s_mid)
                YYeq_SM = np.exp(np.clip(lnYYeq, -700, 700))

                # Physical Y yield used in the forward-decay term:
                #   * traj_has_decay=False: trajectory is counterfactual
                #     (no decay), so YY_phys = YY_traj * f_survive.
                #   * traj_has_decay=True : trajectory already includes
                #     decay, so YY_phys = YY_traj (f_survive stays at 1).
                YY_phys_mid = YY_mid if traj_has_decay else YY_mid * f_survive

                bracket = B1_D * YY_phys_mid - B1_SM * YYeq_SM
                ln_S_corr += GammaY * bracket / T_mid * dt

                if not traj_has_decay:
                    f_survive *= np.exp(-GammaY * dt)

            S_ratio_corr = float(np.exp(ln_S_corr))

            # Apply dilution + Y-survival corrections at the handoff.
            YX0 = YX0 / S_ratio_corr
            if traj_has_decay:
                YY0 = YY0 / S_ratio_corr
            else:
                YY0 = YY0 * f_survive / S_ratio_corr
            T0  = T0  * S_ratio_corr ** (1.0 / 3.0)

            if verbose:
                print(f"  Y-decay correction: S_ratio_corr={S_ratio_corr:.4f}, "
                      f"f_survive={f_survive:.4f}, "
                      f"traj_has_decay={traj_has_decay}")

        s0    = cosmo.s_entropy(T0)
        rhoX0 = mX * YX0 * s0
        rhoY0 = cosmo.Bfac1(mY / TD0, mY) * YY0 * s0

        sv_XXYY = model.sigmav_XX_to_YY_swave()

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

            # Reverse YY -> XX rate: matters when xi is still elevated and
            # Y is abundant.  Uses matter-cooling xi (T_h ~ a^-2) and
            # equilibrium densities at T_h.  Self-suppresses to ~0 for
            # xi <= 1 via the Boltzmann factor in (n_Xeq/n_Yeq)^2.
            xi_now = xi0 * np.exp(-2.0 * N - y[0])
            T_h = max(xi_now * T, 1e-30)
            n_Xeq_h = cosmo.neq_MB(mX / T_h, gX, mX, include_antiparticlesX)
            n_Yeq_h = cosmo.neq_MB(mY / T_h, gY, mY, include_antiparticlesY)
            reverse = 0.0
            if n_Yeq_h > 1e-300:
                n_X = rhoX / mX
                n_Y = rhoY / mY
                ratio_sq = (n_Xeq_h / n_Yeq_h) ** 2
                reverse = (sv_XXYY * ann_sym_bg * ratio_sq
                           * n_Y * n_Y / (n_X * max(H, 1e-300)))

            dlnrhoX_dN = (-3.0
                          - sv_XXYY * ann_sym_bg * rhoX / (mX * max(H, 1e-300))
                          + reverse)

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
            return {'YX_final': YX_final, 'SRatio': SRatio, 'T_final': Tsol[-1],
                    'S_ratio_corr': S_ratio_corr}

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
        bg_kw: Optional[Dict] = None, verbose=False,
        f_survive_floor: float = 1e-10,
        n_scan: int = 17,
    ) -> float:
        """
        Root-find epsX (Brent on log10) such that the background evolution
        from `bg_ICs` reproduces Omega_c h^2 = Och2_target.  Counterpart of
        `find_mX_DMRelic_joint` (which fixes epsX and root-finds mX).
        Caps log10eps_hi where Y has already decayed before bg starts
        (avoids a chaotic numerical regime); falls back to a sign-change
        scan over the uncapped range if brentq fails.
        """
        if bg_kw is None:
            bg_kw = {}

        Y_target = Och2_target * rhoc / (s0_cosmo * self.model.mX)

        # Dynamic cap on log10eps_hi: solve_background becomes numerically
        # chaotic when Y has fully decayed before the bg phase (f_survive
        # underflows to ~0, event never fires, Radau runs 80 e-folds of
        # pure dilution and returns garbage).  This regime is also
        # unphysical for the relic bound: if Y is gone, the answer is
        # saturated to YX0 and no root exists up there anyway.  So walk
        # log10eps_hi down until f_survive(upper) > f_survive_floor.
        t_total = self._bg_trajectory_time(bg_ICs)
        log10eps_lo, log10eps_hi = log10eps_range
        while log10eps_hi > log10eps_lo + 1.0:
            GammaY_hi = self.model.decay_width_to_SM(10 ** log10eps_hi)
            f_surv_hi = np.exp(-GammaY_hi * t_total) if GammaY_hi > 0 else 1.0
            if f_surv_hi > f_survive_floor:
                break
            log10eps_hi -= 0.5

        def f(log10eps):
            res = self.solve_background(
                epsX=10**log10eps, bg_ICs=bg_ICs,
                stop_on_rhoY=True, verbose=False, **bg_kw,
            )
            if verbose:
                print(res['YX_final']*self.model.mX)
            return np.log10(res['YX_final']) - np.log10(Y_target)

        try:
            log10eps = brentq(
                f, log10eps_lo, log10eps_hi, rtol=root_rtol)
        except ValueError:
            # Brentq on the capped bracket failed. Two failure modes:
            # (a) f is non-monotonic and both capped endpoints have
            #     the same sign (the dilution valley sits between two
            #     positive regions);
            # (b) the cap walked too far down and excluded the
            #     plateau region where f > 0, so the capped bracket
            #     is entirely on the dilution side (f < 0 throughout).
            # In both cases we scan the ORIGINAL (uncapped) range for
            # sign changes. For mX just above mXstar, the physical
            # root sits at large eps in the plateau regime, which the
            # cap excludes by construction. We pick the rightmost
            # sign change (largest eps -- physical "least dilution
            # needed" root).
            scan_l10 = np.linspace(
                log10eps_range[0], log10eps_range[1], n_scan)
            scan_f = np.array([f(l10) for l10 in scan_l10])
            sign_change = np.where(
                np.sign(scan_f[:-1]) * np.sign(scan_f[1:]) < 0)[0]
            if len(sign_change) == 0:
                raise ValueError(
                    f"find_epsilon_DMRelic: no root in "
                    f"log10eps_range={log10eps_range} for "
                    f"mX={self.model.mX:.4e} "
                    f"(scan f_min={scan_f.min():.3e}, "
                    f"f_max={scan_f.max():.3e})")
            i = int(sign_change[-1])
            log10eps = brentq(
                f, scan_l10[i], scan_l10[i + 1], rtol=root_rtol)
        return 10**log10eps

    def _bg_trajectory_time(self, bg_ICs: Dict[str, np.ndarray]) -> float:
        """Total cosmological time over the Boltzmann trajectory used by
        the Y-decay correction in solve_background.  Given this, the
        surviving Y fraction at the handoff point is just
        f_survive(epsX) = exp(-GammaY(epsX) * t_total)."""
        cosmo = self.cosmo
        mX, mY = self.model.mX, self.model.mY
        x_traj = bg_ICs['x']
        T_traj = bg_ICs['T']
        YX_traj = bg_ICs['YX']
        YY_traj = bg_ICs['YY']
        t_total = 0.0
        for i in range(1, len(x_traj)):
            x_mid = 0.5 * (x_traj[i] + x_traj[i-1])
            T_mid = 0.5 * (T_traj[i] + T_traj[i-1])
            YY_mid = 0.5 * (YY_traj[i] + YY_traj[i-1])
            s_mid = cosmo.s_entropy(T_mid)
            H_sm  = cosmo.hubble(T_mid)
            rho_h = (mX * 0.5*(YX_traj[i]+YX_traj[i-1])
                     + mY * YY_mid) * s_mid
            Htot  = np.sqrt(H_sm**2 + rho_h / (3.0 * Mpl**2))
            g_t   = 1.0 + cosmo.dlngSdlnT(T_mid) / 3.0
            dx    = x_traj[i] - x_traj[i-1]
            t_total += dx / (x_mid * Htot * g_t)
        return t_total

    # =============================================================
    #  Epsilon finder — BBN hadronic injection constraint
    # =============================================================
    def find_epsilon_BBN_hadronic(
        self,
        YY_frozen: float,
        log10eps_range: tuple = (-20, -6),
        root_rtol: float = 1e-3,
        verbose: bool = False,
    ) -> float:
        """
        Minimum portal coupling such that hadronic injection from Y decay
        satisfies the Kawasaki et al. (2017) BBN constraint on m_Y * Y_Y
        as a function of the Y lifetime tau = hbar / Gamma(Y -> SM).
        Bound is weighted by the hadronic branching fractions (uu vs bb).
        `YY_frozen` is the comoving Y yield well above BBN (post-cannibal
        freeze-out, pre-decay).
        """
        from .bbn import mY_bound_model

        model = self.model
        mY = model.mY
        branching_ratios = model.branching_ratios_to_SM()
        hbar_GeV_s = 6.582119514e-25   # hbar in GeV*s

        mY_times_YY = mY * YY_frozen

        def f(log10eps):
            epsX = 10.0 ** log10eps
            Gamma = model.decay_width_to_SM(epsX)
            if Gamma <= 0:
                return -1.0   # no decay -> certainly excluded
            tau = hbar_GeV_s / Gamma
            bound = mY_bound_model(mY, tau, branching_ratios)

            if verbose:
                print(f"  log10eps={log10eps:+.3f}  tau={tau:.3e} s  "
                      f"mY*YY={mY_times_YY:.3e}  bound={bound:.3e}")

            # f > 0 means allowed (bound > mY*YY)
            # f < 0 means excluded
            return bound - mY_times_YY

        log10eps = brentq(
            f, log10eps_range[0], log10eps_range[1], rtol=root_rtol,
        )
        return 10.0 ** log10eps

    # =============================================================
    #  Epsilon finder — thermalization floor (heuristic estimate)
    # =============================================================
    def find_epsilon_thermal_floor(self, x_FO: float = 20.0) -> float:
        """
        Heuristic ESTIMATE of the minimum portal coupling for the hidden
        sector to remain in thermal equilibrium with the SM down to
        freeze-out.  This is NOT a solver result -- it is the analytic
        criterion Gamma(Y -> SM) >= H(T_FO) at T_FO = mX / x_FO, giving
            eps_therm ~ sqrt(H_FO / decay_width_to_SM(eps=1))
        (Sec. 3.1.2 of arXiv:2509.08043).  Use as a seed/anchor only;
        for precise thermalization boundaries solve the joint Boltzmann
        evolution.
        """
        T_FO = self.model.mX / x_FO
        H_FO = self.cosmo.hubble(T_FO)

        Gamma_norm = self.model.decay_width_to_SM(1.0)
        if Gamma_norm <= 0:
            return np.inf

        return np.sqrt(H_FO / Gamma_norm)

    # =============================================================
    #  Epsilon finder — direct-detection constraint
    # =============================================================
    def find_epsilon_DD(
        self,
        constraint: str = 'LZ',
        log10eps_range: tuple = (-20, -1),
        root_rtol: float = 1e-6,
    ) -> float:
        """
        Portal coupling that saturates a direct-detection SI bound on
        the DM-nucleon cross section.

        `constraint` selects the limit function:
          'LZ'         — LZ 90% CL (9-10000 GeV).
          'XENONnT'    — XENONnT 2024 combined 90% CL.
          'nufloor_Xe' — Xenon SI neutrino floor (0.1-10000 GeV).
          callable     — custom sigma_limit(mX) -> cm^2.

        Returns np.inf if the model has no SI cross section (e.g. LiLjPortal),
        np.nan if mX falls outside the limit function's mX range.
        """
        # --- Resolve the limit function ---
        if callable(constraint):
            sigma_limit_fn = constraint
        elif constraint == 'LZ':
            from .direct_detection import sigma_SI_LZ
            sigma_limit_fn = sigma_SI_LZ
        elif constraint == 'XENONnT':
            from .direct_detection import sigma_SI_XENONnT
            sigma_limit_fn = sigma_SI_XENONnT
        elif constraint == 'nufloor_Xe':
            from .direct_detection import sigma_SI_nufloor_Xe
            sigma_limit_fn = sigma_SI_nufloor_Xe
        else:
            raise ValueError(
                f"Unknown constraint '{constraint}'. "
                "Use 'LZ', 'XENONnT', 'nufloor_Xe', or a callable.")

        mX = self.model.mX
        try:
            sigma_target = sigma_limit_fn(mX)
        except ValueError:
            return np.nan

        # Check if the model has a nonzero SI cross section
        eps_test = 10.0 ** (0.5 * (log10eps_range[0] + log10eps_range[1]))
        if self.model.sigma_SI(eps_test) == 0.0:
            return np.inf

        def f(log10eps):
            eps = 10.0 ** log10eps
            return np.log10(self.model.sigma_SI(eps)) - np.log10(sigma_target)

        log10eps = brentq(
            f, log10eps_range[0], log10eps_range[1], rtol=root_rtol,
        )
        return 10.0 ** log10eps


# =============================================================
#  find_mX_DMRelic_joint  (module-level function)
# =============================================================
def find_mX_DMRelic_joint(
    cosmo: Cosmology,
    model_factory,
    epsX: float,
    *,
    r: float,
    Och2_target: float = Och2,
    mX_guess: Optional[float] = None,
    bracket_factor: tuple = (0.25, 4.0),
    n_refine: int = 5,
    polish_tol: Optional[float] = None,
    polish_maxiter: int = 10,
    joint_kw: Optional[Dict] = None,
    bg_kw: Optional[Dict] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    For a fixed portal coupling epsX, find the DM mass mX that produces the
    target relic abundance Omega_c h^2 via the joint solver workflow:

        solve_boltzmann_joint(epsX=epsX, ...)
        --> solve_background(epsX=epsX, correct_Y_decay=True,
                             traj_has_decay=True, ...)

    Counterpart to `find_epsilon_DMRelic` (which fixes mX and root-finds
    eps).  Two-pass log-interpolation: pass-1 evaluates Och2 at the bracket
    endpoints; pass-2 refines with `n_refine` log-spaced points around the
    linearly-interpolated guess; final answer is the interpolation of
    log10(mX) vs log10(Och2) at the target.

    Also reports Gamma_Y/H at the joint-solver handoff x, so a caller can
    detect when epsX has dropped low enough that the joint workflow is no
    longer necessary (Gamma/H << 1 means Y decay is sub-Hubble at handoff
    and the secluded workflow is more appropriate).

    Parameters
    ----------
    cosmo : Cosmology
    model_factory : callable(mX, mY) -> HiddenSectorModel
    epsX : float
        Portal coupling, fixed.
    r : float
        Mass ratio mX / mY.
    mX_guess : float, optional
        Seed mass [GeV].  Default seeds from a thermalized sigma-v estimate.
    bracket_factor : (lo, hi)
        Initial bracket [mX_guess * lo, mX_guess * hi].
    n_refine : int
        Refinement points in pass 2.
    polish_tol : float, optional
        Relative tolerance |Och2/target - 1| above which a brentq polish
        is run on the existing grid bracket.  None (default) disables
        polishing and preserves legacy behavior.  Typical value: 5e-3.
    polish_maxiter : int
        Max brentq iterations (joint solves) during polish.
    joint_kw, bg_kw : dict, optional
        Extra kwargs forwarded to `solve_boltzmann_joint` and
        `solve_background`.  Defaults are the production-tested settings:
          joint_kw : min_eps_floor_ratio=0.0, phase2_form='Y',
                     return_bg_ICs=True
          bg_kw    : correct_Y_decay=True, traj_has_decay=True,
                     stop_on_rhoY=True

    Returns
    -------
    dict with keys
      mX, mY, epsX
      Och2_achieved : Omega_c h^2 at the interpolated mX (final-check eval)
      GoH_handoff   : Gamma_Y / H at the joint-solver handoff x
      mX_est        : seed used
      grid_log10mX, grid_log10Och2 : scan history
    """
    # --- defaults ---
    _joint_kw = dict(verbose=False, return_bg_ICs=True,
                     min_eps_floor_ratio=0.0, phase2_form='Y')
    if joint_kw:
        _joint_kw.update(joint_kw)
    _bg_kw = dict(correct_Y_decay=True, traj_has_decay=True,
                  stop_on_rhoY=True, verbose=False)
    if bg_kw:
        _bg_kw.update(bg_kw)

    # --- mX seed ---
    if mX_guess is None:
        _mtmp = model_factory(1.0, 1.0 / r)
        sv1 = _mtmp.sigmav_XX_to_YY_swave()
        alpha_eff = np.sqrt(sv1 / np.pi)
        mX_est = alpha_eff * np.sqrt(np.pi * mY_relic * Mpl)
    else:
        mX_est = mX_guess

    log10mX_lo = np.log10(bracket_factor[0] * mX_est)
    log10mX_hi = np.log10(bracket_factor[1] * mX_est)
    log10_target = np.log10(Och2_target)

    # mX -> diagnostic GoH at handoff (keyed by log10mX)
    GoH_at = {}

    def eval_lOch2(log10mX):
        mX = 10.0 ** log10mX
        mY = mX / r
        mdl = model_factory(mX, mY)
        solver = BoltzmannSolver(cosmo, mdl)
        sol = solver.solve_boltzmann_joint(epsX=epsX, **_joint_kw)
        # Phase 1.5-S: Y is in detailed balance with the SM throughout
        # (decay <-> inverse decay).  No net Y -> SM entropy injection;
        # solve_background's "free Y decay" assumption is wrong here.
        # Use the joint result directly.
        if sol['phase_final'] == '1.5S':
            Och2_val = float(sol['YX_relic']) * mX * s0_cosmo / rhoc
        else:
            sol_bg = solver.solve_background(
                epsX=epsX, bg_ICs=sol['bg_ICs'], **_bg_kw)
            Och2_val = float(sol_bg['YX_final']) * mX * s0_cosmo / rhoc
        x_handoff = float(sol['x'][-1])
        T_handoff = mX / x_handoff
        H_handoff = cosmo.hubble(T_handoff)
        GoH_at[log10mX] = float(sol['GammaY']) / H_handoff
        return np.log10(max(Och2_val, 1e-30))

    if verbose:
        print("=" * 60)
        print(f"find_mX_DMRelic_joint  "
              f"(epsX={epsX:.3e}, r={r:.2f}, target Och2={Och2_target:.4f})")
        print(f"  mX_est = {mX_est:.3e} GeV  "
              f"bracket = [{10**log10mX_lo:.3e}, {10**log10mX_hi:.3e}] GeV")
        print("-" * 60)

    # --- Pass 1: bracket endpoints ---
    lm = np.array([log10mX_lo, log10mX_hi])
    lOch2 = np.array([eval_lOch2(lm[0]), eval_lOch2(lm[1])])
    if verbose:
        for i in range(2):
            print(f"  pass-1: mX={10**lm[i]:.3e}  Och2={10**lOch2[i]:.4e}  "
                  f"GoH={GoH_at[lm[i]]:.2e}")

    if log10_target < lOch2.min() or log10_target > lOch2.max():
        raise ValueError(
            f"Target Och2={Och2_target:.3e} not bracketed by initial "
            f"endpoint Och2 in [{10**lOch2.min():.3e}, "
            f"{10**lOch2.max():.3e}]; widen `bracket_factor` or supply "
            f"`mX_guess`.")

    o = np.argsort(lOch2)
    log10mX_guess = float(np.interp(log10_target, lOch2[o], lm[o]))

    # --- Pass 2: refine ---
    refine = np.linspace(log10mX_guess + np.log10(0.5),
                         log10mX_guess + np.log10(2.0), n_refine)
    for rv in refine:
        lm = np.append(lm, rv)
        lOch2 = np.append(lOch2, eval_lOch2(rv))
        if verbose:
            print(f"  refine: mX={10**rv:.3e}  Och2={10**lOch2[-1]:.4e}  "
                  f"GoH={GoH_at[rv]:.2e}")

    o = np.argsort(lOch2)
    log10mX_sol = float(np.interp(log10_target, lOch2[o], lm[o]))
    mX_sol = 10.0 ** log10mX_sol

    # --- Final check at interpolated mX ---
    lOch2_check = eval_lOch2(log10mX_sol)
    Och2_achieved = 10.0 ** lOch2_check
    GoH_handoff = GoH_at[log10mX_sol]

    if verbose:
        print(f"  => mX = {mX_sol:.4e} GeV   Och2 = {Och2_achieved:.6f}   "
              f"GoH_handoff = {GoH_handoff:.2e}")

    # --- Optional brentq polish if interpolation residual exceeds tol ---
    if polish_tol is not None:
        residual = abs(Och2_achieved / Och2_target - 1.0)
        if residual > polish_tol:
            # Pool all evaluations (incl. final-check point) and find the
            # tightest sign-change bracket of (lOch2 - log10_target).
            all_lm = np.append(lm, log10mX_sol)
            all_lO = np.append(lOch2, lOch2_check)
            o = np.argsort(all_lm)
            lm_s, lO_s = all_lm[o], all_lO[o]
            f_s = lO_s - log10_target
            sgn = np.where(np.diff(np.sign(f_s)) != 0)[0]
            if sgn.size == 0:
                if verbose:
                    print(f"  polish: no sign change in grid; "
                          f"residual={residual:.2%} kept.")
            else:
                # Pick the bracket whose midpoint is closest to log10mX_sol
                mids = 0.5 * (lm_s[sgn] + lm_s[sgn + 1])
                k = sgn[np.argmin(np.abs(mids - log10mX_sol))]
                a, b = float(lm_s[k]), float(lm_s[k + 1])
                last = {'lm': log10mX_sol, 'lO': lOch2_check}

                def _f_polish(x):
                    last['lm'] = x
                    last['lO'] = eval_lOch2(x)
                    return last['lO'] - log10_target

                try:
                    log10mX_sol = float(brentq(
                        _f_polish, a, b,
                        xtol=1e-4, rtol=1e-4, maxiter=polish_maxiter))
                except Exception as exc:
                    if verbose:
                        print(f"  polish: brentq failed ({exc}); "
                              f"keeping interpolation result.")
                else:
                    mX_sol = 10.0 ** log10mX_sol
                    if abs(last['lm'] - log10mX_sol) > 1e-12:
                        lOch2_check = eval_lOch2(log10mX_sol)
                    else:
                        lOch2_check = last['lO']
                    Och2_achieved = 10.0 ** lOch2_check
                    GoH_handoff = GoH_at[log10mX_sol]
                    if verbose:
                        print(f"  polish: brentq -> mX={mX_sol:.4e}  "
                              f"Och2={Och2_achieved:.6f}  "
                              f"GoH_handoff={GoH_handoff:.2e}")

    if verbose:
        print("=" * 60)

    return {
        'mX': mX_sol,
        'mY': mX_sol / r,
        'epsX': epsX,
        'Och2_achieved': Och2_achieved,
        'GoH_handoff': GoH_handoff,
        'mX_est': mX_est,
        'grid_log10mX': lm,
        'grid_log10Och2': lOch2,
    }

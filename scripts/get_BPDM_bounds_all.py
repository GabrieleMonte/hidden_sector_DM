"""
Baryon Portal: relic-abundance + BBN bounds in (mX, eps) plane.

Mirrors get_vPDM_bounds_all.py.  Two-stage workflow:

  Stage 1 -- Joint piece (low mass / low eps):
    For each k in K_GRID, find mX such that the joint workflow gives
    Omega_c h^2 = 0.12.  Track GoH_handoff = Gamma_Y/H at the joint
    solver's handoff x; the joint sweep stops on the first row with
    GoH < GOH_THRESH.  Saved as a checkpoint .dat.

  Stage 2 -- Secluded piece (high mass / smaller eps):
    From mX_join_max upward, build a clustered grid (ZB-pole + Y -> q
    qbar decay thresholds + log tails) up to 10**log10mXmax.  Root-find
    eps via solve_boltzmann_chempot_3phase + find_epsilon_DMRelic.

Output:
  ../output/bP/relic_joint_part_aX{alphaX:.0e}_r{r}.dat   (Stage 1 only)
  ../output/bP/relic_joint_aX{alphaX:.0e}_r{r}.dat        (Stages 1 + 2)

Columns: mX  mY  epsX_relic  epsX_bbn  Och2  GoH_handoff  regime
where regime = 0 (joint), 1 (secluded).
"""

import sys
from pathlib import Path
import numpy as np
from tqdm import trange
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from hidden_sector_DM.HiddenSectorDM import (
    Cosmology, BaryonPortal, BoltzmannSolver, find_mX_DMRelic_joint,
)
from hidden_sector_DM.constants import (
    Mpl, Mc, Mb, Mt, mY_relic, s0_cosmo, rhoc,
    Och2 as OCH2_TARGET,
)


# ---- configuration ------------------------------------------------
RUN_SECLUDED = True

# k-grid: 30 points, log-spaced, densest near k=1.
K_GRID = np.unique(np.concatenate([
    np.logspace(-2, 0, 13),
    np.logspace(0, 1, 13),
    np.linspace(0.5, 2.0, 8),
]))

# joint -> secluded boundary.
GOH_THRESH = 0.05

# Brentq polish tolerance for find_mX_DMRelic_joint.
POLISH_TOL = 5e-3

# eps-search ranges for the secluded stage (BaryonPortal: g_B coupling
# multiplied by epsX, no resonance below 100 GeV mediator -> eps stays
# small in the relic line, same windows as VectorPortal).
EPS_RELIC_RANGE = (-18, -6)
EPS_BBN_RANGE   = (-20, -9)

# ZB pole: BaryonPortal kinetic-mixing pole at mY = m_ZB (default
# 1000 GeV).  No SM Z-pole here.
M_ZB = 1e4

# Y -> q qbar decay thresholds: BaryonPortal couples to quarks only.
DECAY_PARTNERS = (Mc, Mb, Mt)

cases = [
    (1e-3, 5.0, 4),
    (1e-3, 2.5, 4),
    (1e-2, 2.5, 4),
    (1e-2, 5.0, 4),
]


cosmo = Cosmology()


def find_mX_thermalized(factory, r):
    """Brent on log10(mX) using solve_boltzmann_thermSM, seeded from
    the simple sigmav_swave estimate."""
    sv1 = factory(1.0, 1.0 / r).sigmav_XX_to_YY_swave()
    alpha_eff = np.sqrt(sv1 / np.pi)
    mX_est = alpha_eff * np.sqrt(np.pi * mY_relic * Mpl)
    lo = max(mX_est / 5.0, 5.0)
    hi = mX_est * 5.0

    def f(log10mX):
        mX = 10.0 ** log10mX
        bs = BoltzmannSolver(cosmo, factory(mX, mX / r))
        sol = bs.solve_boltzmann_thermSM(verbose=False)
        Och2 = sol['YX_relic'] * mX * s0_cosmo / rhoc
        return np.log10(Och2 / OCH2_TARGET)
    return 10.0 ** brentq(f, np.log10(lo), np.log10(hi), rtol=1e-4)


FAST_PATH_K_MIN = 2.0
FAST_PATH_TOL   = 0.01


def _solve_at(factory, r, mX, eps):
    """Single joint-solver eval at fixed (mX, eps); returns
    (Och2, GoH_handoff, eps_bbn).  Phase 1.5-S: Y in detailed balance
    with SM (no bg post-step).  Phase 2: hand off to bg for late-time
    Y decay + dilution."""
    model = factory(mX, mX / r)
    solver = BoltzmannSolver(cosmo, model)
    sol = solver.solve_boltzmann_joint(
        epsX=eps, return_bg_ICs=True,
        min_eps_floor_ratio=0.0, phase2_form='Y',
        verbose=False)
    if sol['phase_final'] == '1.5S':
        Och2 = float(sol['YX_relic']) * mX * s0_cosmo / rhoc
    else:
        bg = solver.solve_background(
            epsX=eps, bg_ICs=sol['bg_ICs'],
            correct_Y_decay=True, traj_has_decay=True,
            stop_on_rhoY=True, verbose=False)
        Och2 = float(bg['YX_final']) * mX * s0_cosmo / rhoc
    YY_for_bbn = float(sol['YY_final'])
    x_handoff = float(sol['x'][-1])
    GoH = float(sol['GammaY']) / cosmo.hubble(mX / x_handoff)
    try:
        eps_bbn = solver.find_epsilon_BBN_hadronic(
            YY_for_bbn, log10eps_range=EPS_BBN_RANGE, root_rtol=1e-3)
    except Exception:
        eps_bbn = np.nan
    return Och2, GoH, eps_bbn


def joint_sweep(factory, r, eps_therm, mX_therm):
    """Sweep over K_GRID (high k to low k).  Fast path for k > 2 tries
    mX_therm first (deep thermal regime).  Otherwise fall through to
    find_mX_DMRelic_joint, seeded from previous k.  Stops on first row
    with GoH < GOH_THRESH (not recorded)."""
    rows = []
    print(f"  Joint sweep: {len(K_GRID)} k-values  "
          f"(fast-path: k>{FAST_PATH_K_MIN}, tol={FAST_PATH_TOL:.0%})",
          flush=True)
    mX_seed = mX_therm
    for k in sorted(K_GRID, reverse=True):
        eps = k * eps_therm

        if k > FAST_PATH_K_MIN:
            try:
                Och2_t, GoH_t, eps_bbn_t = _solve_at(
                    factory, r, mX_therm, eps)
            except Exception as exc:
                print(f"    k={k:>7.4f}  fast-path FAILED: {exc}",
                      flush=True)
                continue

            if abs(Och2_t / OCH2_TARGET - 1.0) < FAST_PATH_TOL:
                rows.append((mX_therm, mX_therm / r, eps,
                             eps_bbn_t, Och2_t, GoH_t, 0))
                print(f"    k={k:>7.4f}  [FAST]  mX={mX_therm:.3e}  "
                      f"Och2={Och2_t:.5f}  GoH={GoH_t:.2e}  "
                      f"eps_bbn={eps_bbn_t:.3e}", flush=True)
                continue

        try:
            res = find_mX_DMRelic_joint(
                cosmo, factory, eps, r=r,
                mX_guess=mX_seed, polish_tol=POLISH_TOL,
                verbose=False)
        except Exception as exc:
            print(f"    k={k:>7.4f}  search FAILED: {exc}", flush=True)
            continue
        mX = res['mX']
        mY = mX / r
        Och2 = res['Och2_achieved']
        GoH = res['GoH_handoff']
        try:
            _, _, eps_bbn = _solve_at(factory, r, mX, eps)
        except Exception:
            eps_bbn = np.nan
        if GoH < GOH_THRESH:
            print(f"    k={k:>7.4f}  [SCAN]  mX={mX:.3e}  GoH={GoH:.2e}"
                  f" < GOH_THRESH={GOH_THRESH}; not recorded. "
                  f"Stopping joint sweep, handing off to secluded.",
                  flush=True)
            break
        rows.append((mX, mY, eps, eps_bbn, Och2, GoH, 0))
        print(f"    k={k:>7.4f}  [SCAN]  mX={mX:.3e}  Och2={Och2:.5f}  "
              f"GoH={GoH:.2e}  eps_bbn={eps_bbn:.3e}", flush=True)
        mX_seed = mX
    return rows


def secluded_sweep(factory, r, mX_join_max, log10mXmax):
    """Clustered mX-grid: log-spaced near + ZB-pole + decay-thr + high
    log tail.  At each mX, find epsX_relic via
    solve_boltzmann_chempot_3phase + find_epsilon_DMRelic, plus eps_bbn."""
    mX_lo = mX_join_max * 1.06
    mX_hi = 10.0 ** log10mXmax
    if mX_lo >= mX_hi:
        print(f"  Secluded skipped (mX_join_max={mX_join_max:.3e} GeV "
              f">= mX_hi={mX_hi:.3e} GeV).", flush=True)
        return []

    mXs_near = np.logspace(np.log10(mX_lo),
                           np.log10(min(10 * mX_join_max, mX_hi)), 100)
    mX_res = r * M_ZB
    if mX_lo < mX_res * 1.3 and mX_res * 0.7 < mX_hi:
        mXs_res = np.linspace(max(mX_res * 0.8, mX_lo),
                              min(mX_res * 1.2, mX_hi), 30)
    else:
        mXs_res = np.array([])
    mXs_thr = []
    for mf in DECAY_PARTNERS:
        mX_thr = 2.0 * r * mf
        if mX_lo < mX_thr < mX_hi:
            mXs_thr.append(np.linspace(
                max(0.92 * mX_thr, mX_lo),
                min(1.25 * mX_thr, mX_hi), 15))
    mXs_thr = np.concatenate(mXs_thr) if mXs_thr else np.array([])
    mXs_high = np.logspace(np.log10(max(10.01 * mX_join_max, mX_lo)),
                           np.log10(mX_hi), 30)
    mXs = np.unique(np.concatenate([mXs_near, mXs_res, mXs_thr, mXs_high]))
    print(f"  Secluded sweep: {len(mXs)} mass points "
          f"({mXs[0]:.3e} -> {mXs[-1]:.3e} GeV)", flush=True)

    rows = []
    for i in trange(len(mXs)):
        mX = mXs[i]
        mY = mX / r
        model = factory(mX, mY)
        solver = BoltzmannSolver(cosmo, model)
        sol0 = solver.solve_boltzmann_chempot_3phase(
            return_bg_ICs=True, verbose=False,
            cannibal_switch_full=1, convergence_threshold=1e-2)
        try:
            eps_relic = solver.find_epsilon_DMRelic(
                sol0['bg_ICs'], log10eps_range=EPS_RELIC_RANGE,
                root_rtol=1e-4, bg_kw={'correct_Y_decay': True})
        except Exception:
            eps_relic = np.nan
        try:
            eps_bbn = solver.find_epsilon_BBN_hadronic(
                sol0['YY_final'], log10eps_range=EPS_BBN_RANGE,
                root_rtol=1e-3)
        except Exception:
            eps_bbn = np.nan
        rows.append((mX, mY, eps_relic, eps_bbn,
                     OCH2_TARGET, np.nan, 1))
        if (i + 1) % 20 == 0 or i == len(mXs) - 1:
            print(f"    secl [{i+1:>3}/{len(mXs)}]  mX={mX:.3e}  "
                  f"eps_relic={eps_relic:.3e}  eps_bbn={eps_bbn:.3e}",
                  flush=True)
        del model, solver
    return rows


# ---- main loop ----------------------------------------------------
for alphaX, r, log10mXmax in cases:
    print(f"\n=== BaryonPortal: alphaX={alphaX:.0e}, r={r}, "
          f"log10mXmax={log10mXmax} ===", flush=True)
    factory = lambda mX, mY, aX=alphaX: BaryonPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=aX,
        include_antiparticlesX=True, include_antiparticlesY=False)

    mX_therm = find_mX_thermalized(factory, r)
    bs0 = BoltzmannSolver(cosmo, factory(mX_therm, mX_therm / r))
    eps_therm = bs0.find_epsilon_thermal_floor(x_FO=20.0)
    print(f"  mX_therm = {mX_therm:.4e} GeV   "
          f"eps_therm = {eps_therm:.4e}", flush=True)

    out_part = (ROOT / "output" / "bP" /
                f"relic_joint_part_aX{alphaX:.0e}_r{r}.dat")
    out_full = (ROOT / "output" / "bP" /
                f"relic_joint_aX{alphaX:.0e}_r{r}.dat")
    out_part.parent.mkdir(parents=True, exist_ok=True)

    if out_part.exists():
        print(f"  Loading cached joint piece from {out_part}", flush=True)
        j_arr = np.loadtxt(out_part)
        if j_arr.ndim == 1:
            j_arr = j_arr[np.newaxis, :]
        j_rows = [tuple(row) for row in j_arr]
    else:
        j_rows = joint_sweep(factory, r, eps_therm, mX_therm)
        if j_rows:
            anchor = (mX_therm, mX_therm / r, 1e-1, np.nan,
                      OCH2_TARGET, np.nan, 0)
            j_arr = np.array([anchor] + list(j_rows))
            header_part = (
                f"alphaX={alphaX:.1e}  r={r}  mX_therm={mX_therm:.6e}  "
                f"eps_therm={eps_therm:.6e}  GoH_thresh={GOH_THRESH}\n"
                f"mX [GeV]    mY [GeV]    epsX_relic    epsX_bbn    "
                f"Och2    GoH_handoff   regime")
            np.savetxt(out_part, j_arr, header=header_part, fmt='%.8e')
            print(f"  Wrote joint-piece checkpoint: {out_part}", flush=True)

    if not j_rows:
        print("  No joint rows; skipping case.", flush=True)
        continue

    boundary = min(j_rows, key=lambda row: row[2])
    mX_join_max = boundary[0]
    print(f"  -> Joint regime ends at mX={mX_join_max:.4e} GeV  "
          f"(eps={boundary[2]:.3e}, GoH={boundary[5]:.3f})",
          flush=True)

    if not RUN_SECLUDED:
        print("  RUN_SECLUDED=False; stopping after joint piece.",
              flush=True)
        continue

    s_rows = secluded_sweep(factory, r, mX_join_max, log10mXmax)

    anchor = (mX_therm, mX_therm / r, 1e-1, np.nan,
              OCH2_TARGET, np.nan, 0)
    data = np.array([anchor] + list(j_rows) + s_rows)
    header = (
        f"alphaX={alphaX:.1e}  r={r}  mX_therm={mX_therm:.6e}  "
        f"eps_therm={eps_therm:.6e}  GoH_thresh={GOH_THRESH}  "
        f"mX_join_max={mX_join_max:.6e}\n"
        f"mX [GeV]    mY [GeV]    epsX_relic    epsX_bbn    "
        f"Och2    GoH_handoff   regime  (0=joint, 1=secluded)")
    np.savetxt(out_full, data, header=header, fmt='%.8e')
    print(f"  Wrote {out_full}", flush=True)

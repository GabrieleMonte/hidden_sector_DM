"""
Vector Portal: relic-abundance + BBN bounds in (mX, eps) plane.

Two-stage workflow (replaces legacy/get_vPDM_bounds_all.py):

  Stage 1 -- Joint piece (low mass / low eps):
    For each k in K_GRID, find mX such that the joint workflow gives
    Omega_c h^2 = 0.12.  Track GoH_handoff = Gamma_Y/H at the joint
    solver's handoff x; the boundary with the secluded regime is set
    by GoH_THRESH below.  Saved as a checkpoint .dat.

  Stage 2 -- Secluded piece (high mass / smaller eps):
    From mX_join_max upward, build the legacy clustered grid (Z-pole
    + Y -> f fbar decay thresholds + log tails) up to 10**log10mXmax.
    For each mX, root-find eps via solve_boltzmann_chempot_3phase +
    find_epsilon_DMRelic.

The two stages share the BBN-bound sweep (find_epsilon_BBN_hadronic).

Toggle RUN_SECLUDED below to control which stage(s) run.  When False
the script only does Stage 1 and writes the joint checkpoint.  When
True, it loads the cached checkpoint (or recomputes if missing) and
appends the secluded piece.

Output:
  ../output/vP/relic_joint_part_aX{alphaX:.0e}_r{r}.dat   (Stage 1 only)
  ../output/vP/relic_joint_aX{alphaX:.0e}_r{r}.dat        (Stages 1 + 2)

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
    Cosmology, VectorPortal, BoltzmannSolver, find_mX_DMRelic_joint,
)
from hidden_sector_DM.constants import (
    Mpl, Mmu, Mtau, Mc, Mb, Mt, mY_relic, s0_cosmo, rhoc,
    Och2 as OCH2_TARGET,
)


# ---- configuration ------------------------------------------------
RUN_SECLUDED = True    # flip to True to run Stage 2 after Stage 1

# k-grid: 30 points, log-spaced, densest near k=1, where k=eps/eps_therm.
K_GRID = np.unique(np.concatenate([
    np.logspace(-2, 0, 13),
    np.logspace(0, 1, 13),
    np.linspace(0.5, 2.0, 8),
]))

# joint -> secluded boundary: agreement between the two solvers degrades
# rapidly above this; below it Och2 agrees to ~1% (test_joint_secluded_threshold.py).
GOH_THRESH = 0.05

# Brentq polish tolerance: if the interp-based mX gives Och2 off-target by
# more than this (relative), find_mX_DMRelic_joint refines via brentq on
# the tight bracket from the existing grid.  Mostly bites the high-k
# (k > ~0.5) rows where the relic curve has the most curvature.
POLISH_TOL = 5e-3

# Final relic-abundance tolerance: if the search-converged Och2 still
# differs from 0.12 by more than this (relative), the point is skipped.
RELIC_TOL = 0.05

cases = [
    (1e-3, 5.0, 4),
    (1e-2, 5.0, 4),
    (1e-2, 2.5, 4),
    (1e-3, 2.5, 4),
]


cosmo = Cosmology()


def find_mX_thermalized(factory, r):
    """Brent on log10(mX) using solve_boltzmann_thermSM, seeded from
    the simple sigmav_swave estimate (same trick as find_mXstar /
    find_mX_DMRelic_joint)."""
    sv1 = factory(1.0, 1.0 / r).sigmav_XX_to_YY_swave()
    alpha_eff = np.sqrt(sv1 / np.pi)
    mX_est = alpha_eff * np.sqrt(np.pi * mY_relic * Mpl)
    lo = max(mX_est / 5.0, 5.0)   # avoid QCD region T_FO < ~0.5 GeV
    hi = mX_est * 5.0

    def f(log10mX):
        mX = 10.0 ** log10mX
        bs = BoltzmannSolver(cosmo, factory(mX, mX / r))
        sol = bs.solve_boltzmann_thermSM(verbose=False)
        Och2 = sol['YX_relic'] * mX * s0_cosmo / rhoc
        return np.log10(Och2 / OCH2_TARGET)
    return 10.0 ** brentq(f, np.log10(lo), np.log10(hi), rtol=1e-4)


FAST_PATH_K_MIN = 2.0    # only attempt the fast path for k > this
FAST_PATH_TOL   = 0.01   # accept mX_therm if |Och2/target - 1| < 1 %


def _solve_at(factory, r, mX, eps):
    """Single eval at fixed (mX, eps); returns (Och2, GoH_handoff, eps_bbn).
    Phase 1.5-S: Y in detailed balance with SM, so no bg post-step needed
    (joint result is the relic).  Phase 2: hand off to bg for late-time
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
    # BBN takes the *frozen* Y abundance (joint-solver handoff value,
    # before solve_background processes the decay).  solve_background
    # only tracks YX, so YY must come from the joint result.
    YY_for_bbn = float(sol['YY_final'])
    x_handoff = float(sol['x'][-1])
    GoH = float(sol['GammaY']) / cosmo.hubble(mX / x_handoff)
    try:
        eps_bbn = solver.find_epsilon_BBN_hadronic(
            YY_for_bbn, log10eps_range=(-20, -9), root_rtol=1e-3)
    except Exception:
        eps_bbn = np.nan
    return Och2, GoH, eps_bbn


def joint_sweep(factory, r, eps_therm, mX_therm):
    """Sweep over K_GRID (high k to low k).  For k > FAST_PATH_K_MIN, try
    mX_therm first; if Och2 is already within FAST_PATH_TOL of target,
    accept it (no search needed since deep in the thermalised regime).
    Otherwise fall through to find_mX_DMRelic_joint (seeded from the
    previous k's solution for tighter brackets / faster convergence)."""
    rows = []
    print(f"  Joint sweep: {len(K_GRID)} k-values  "
          f"(fast-path: k>{FAST_PATH_K_MIN}, tol={FAST_PATH_TOL:.0%})",
          flush=True)
    mX_seed = mX_therm   # rolling seed for the fallback path
    for k in sorted(K_GRID, reverse=True):  # high to low k
        eps = k * eps_therm

        # --- Fast path (high-k only): evaluate at mX_therm ---
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

        # --- Fallback: full mX search around the rolling seed ---
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
        if abs(Och2 / OCH2_TARGET - 1.0) > RELIC_TOL:
            print(f"    k={k:>7.4f}  [SCAN]  mX={mX:.3e}  Och2={Och2:.5f}"
                  f" off-target by >{RELIC_TOL:.0%}; skipping.",
                  flush=True)
            continue
        # BBN at the converged point (separate solve to extract YY_final).
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
        mX_seed = mX   # update seed for next (lower) k
    return rows


def secluded_sweep(factory, r, mX_join_max, log10mXmax):
    """Legacy clustered mX-grid: log-spaced near + Z-pole + decay-thr +
    high-mass log tail.  At each mX, find epsX_relic via
    solve_boltzmann_chempot_3phase + find_epsilon_DMRelic, plus epsX_bbn."""
    mX_lo = mX_join_max * 1.06
    mX_hi = 10.0 ** log10mXmax
    if mX_lo >= mX_hi:
        print(f"  Secluded skipped (mX_join_max={mX_join_max:.3e} GeV "
              f">= mX_hi={mX_hi:.3e} GeV).", flush=True)
        return []

    # Reduced point counts (target ~200 unique total).  Density is biased
    # toward the joint -> secluded boundary (mXs_near is log-spaced from
    # mX_lo = 1.06 * mX_join_max, so it crowds the low end naturally).
    mXs_near = np.logspace(np.log10(mX_lo),
                           np.log10(min(10 * mX_join_max, mX_hi)), 100)
    MZ = 91.1876
    mX_res = r * MZ
    if mX_lo < mX_res * 1.3 and mX_res * 0.7 < mX_hi:
        mXs_res = np.linspace(max(mX_res * 0.8, mX_lo),
                              min(mX_res * 1.2, mX_hi), 30)
    else:
        mXs_res = np.array([])
    mXs_thr = []
    for mf in (Mmu, Mtau, Mc, Mb, Mt):
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
                sol0['bg_ICs'], log10eps_range=(-18, -6),
                root_rtol=1e-4, bg_kw={'correct_Y_decay': True})
        except Exception:
            eps_relic = np.nan
        try:
            eps_bbn = solver.find_epsilon_BBN_hadronic(
                sol0['YY_final'], log10eps_range=(-20, -9),
                root_rtol=1e-3)
        except Exception:
            eps_bbn = np.nan
        rows.append((mX, mY, eps_relic, eps_bbn,
                     OCH2_TARGET, np.nan, 1))   # regime=1 (secluded)
        if (i + 1) % 20 == 0 or i == len(mXs) - 1:
            print(f"    secl [{i+1:>3}/{len(mXs)}]  mX={mX:.3e}  "
                  f"eps_relic={eps_relic:.3e}  eps_bbn={eps_bbn:.3e}",
                  flush=True)
        del model, solver
    return rows


# ---- main loop ----------------------------------------------------
for alphaX, r, log10mXmax in cases:
    print(f"\n=== VectorPortal: alphaX={alphaX:.0e}, r={r}, "
          f"log10mXmax={log10mXmax} ===", flush=True)
    factory = lambda mX, mY, aX=alphaX: VectorPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=aX,
        include_antiparticlesX=True, include_antiparticlesY=False)

    # --- anchor ---
    mX_therm = find_mX_thermalized(factory, r)
    bs0 = BoltzmannSolver(cosmo, factory(mX_therm, mX_therm / r))
    eps_therm = bs0.find_epsilon_thermal_floor(x_FO=20.0)
    print(f"  mX_therm = {mX_therm:.4e} GeV   "
          f"eps_therm = {eps_therm:.4e}", flush=True)

    out_part = (ROOT / "output" / "vP" /
                f"relic_joint_part_aX{alphaX:.0e}_r{r}.dat")
    out_full = (ROOT / "output" / "vP" /
                f"relic_joint_aX{alphaX:.0e}_r{r}.dat")
    out_part.parent.mkdir(parents=True, exist_ok=True)

    # --- Stage 1: joint piece (load checkpoint if present) ---
    if out_part.exists():
        print(f"  Loading cached joint piece from {out_part}", flush=True)
        j_arr = np.loadtxt(out_part)
        if j_arr.ndim == 1:
            j_arr = j_arr[np.newaxis, :]
        j_rows = [tuple(row) for row in j_arr]
    else:
        j_rows = joint_sweep(factory, r, eps_therm, mX_therm)
        if j_rows:
            # Preserve iteration order (high k first => eps descending);
            # joint_sweep already iterates in that order.  Prepend the
            # manual high-eps anchor (eps=0.1, mX=mX_therm, Och2=target).
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

    # --- joint -> secluded boundary ---
    # joint_sweep stops on the first GoH < GOH_THRESH (no record), so
    # every row in j_rows is already valid; the last (smallest eps) is
    # the boundary.
    boundary = min(j_rows, key=lambda row: row[2])
    mX_join_max = boundary[0]
    print(f"  -> Joint regime ends at mX={mX_join_max:.4e} GeV  "
          f"(eps={boundary[2]:.3e}, GoH={boundary[5]:.3f})",
          flush=True)

    if not RUN_SECLUDED:
        print("  RUN_SECLUDED=False; stopping after joint piece.",
              flush=True)
        continue

    # --- Stage 2: secluded piece ---
    s_rows = secluded_sweep(factory, r, mX_join_max, log10mXmax)

    # --- glue + save ---
    # Preserve computation order: joint rows in iteration order (high k
    # first => eps descending), then secluded rows by mX ascending (the
    # trange loop natural order).  No global sort.  A manual anchor row
    # at eps=0.1, mX=mX_therm, Och2=target is prepended -- this is the
    # high-eps saturation point (Y deeply locked to SM, no computation
    # needed).
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

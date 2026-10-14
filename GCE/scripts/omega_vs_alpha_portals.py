"""Do the BL / LiLj / Baryon grids need rebuilding, and did they ever get the
Phase-3 follow-up?  Two questions, one probe, at the campaign's own settings.

SECLUDED -- the shared driver's crossing extraction changed on 2026-10-13
(commit 7b75d867e3): the coupling is now read off the first adjacent pair of
solved samples straddling the target, not off a Pchip fit to the whole axis
sorted by Omega.  On the Higgs portal that moved 25% of nodes, by up to 500x,
because 48 of 60 probed Omega(lam) curves are not monotone and the global sort
fitted both branches together.  The vector-portal curves are believed well
behaved -- relic_alphaX_local matches relic_alphaX to <= 0.2% on complete axes
(GCE/docs/phase3_followup.md) -- but that was measured on the JOINT slices of
vP.  The secluded slices of BL, LiLj and B have never been looked at.  This
evaluates Omega(alpha_X) on the campaign axis and runs both extractions over
the same curves.

JOINT -- the follow-up lowered Sdot_switch_thresh from 0.05 to 1e-3 for the vP
grids.  BL / LiLj / B were built at the campaign's 0.05 and never revisited.
Lowering it is not automatically safer: at 1e-3 some alphas handed off early
without running Phase 3, and the legacy path then gave Omega 2-7x too low.
This records the exit phase at every alpha so the legacy-exit rate at 0.05 can
be counted instead of assumed.

Both parts read the campaign's own COUPLING axis and PHASE3_KW, so the answer
applies to the grids as they are stored, not to some other sampling.

    cosmo_env
    sbatch GCE/scripts/omega_vs_alpha_portals.sbatch

Writes GCE/output/diag/omega_vs_alpha_portals.npz.
"""

import sys
from pathlib import Path

import numpy as np
from tqdm import tqdm

import multiprocessing as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                    # noqa: E402

# ---- configuration ---------------------------------------------------

PORTALS_PROBED = ("BL", "LiLj", "B")

N_NODE = 16              # nodes per portal, spread evenly over the solved mesh
EPS_JOINT = (1e-9, 1e-10)

PROBE_TIMEOUT = 600      # s per solve; the campaign used 7200 for stiff nodes.
                         # A login-node smoke test measured 8-10 s at rv=0.12, so
                         # this only ever bites the stiff high-rv end.

MOVED_TOL = 0.05         # relative change in the crossing that counts as moved

# Exit-phase codes, the convention make_relic_grid_vP.py uses for phase_path.
PATH_FAILED = 0
PATH_LEGACY = 1
PATH_PHASE3 = 3
PATH_SMLOCK = 4

DIAG_DIR = base.OUT_DIR.parent / "diag"
OUT = DIAG_DIR / "omega_vs_alpha_portals.npz"


# ---- the pre-fix extraction, for comparison only ---------------------

def crossing_global(omega):
    """The crossing as make_relic_grid._crossing read it before 7b75d867e3.

    Copied verbatim from that commit so the comparison is against what
    actually built the stored grids, not a reconstruction of it.  Returns the
    coupling alone; the status codes are not what is being compared.
    """
    finite = np.isfinite(omega) & (omega > 0)
    if finite.sum() < 4:
        return np.nan

    om, ax = omega[finite], base.COUPLING[finite]
    if not (om.min() <= base.OCH2_TARGET <= om.max()):
        return np.nan

    order = np.argsort(om)
    log_om, log_ax = np.log(om[order]), np.log(ax[order])

    keep = np.concatenate(([True], np.diff(log_om) > 0.0))
    if keep.sum() < 2:
        return np.nan

    interp = base.PchipInterpolator(log_om[keep], log_ax[keep])
    crossing = float(np.exp(interp(np.log(base.OCH2_TARGET))))

    if not (ax.min() <= crossing <= ax.max()):
        return np.nan
    return crossing


# ---- node selection --------------------------------------------------

def pick_nodes(portal):
    """Flat indices of N_NODE solved nodes for `portal`, sampled evenly.

    Even sampling over the flat index spreads the picks across mX and rv
    rather than clustering them wherever the grid happens to be dense.
    """
    base._set_portal(portal, 1.0)
    path = base.OUT_DIR / f"relic_grid_{base.TAG}_eps_secluded.npz"
    if not path.exists():
        raise SystemExit(f"{path.name}: not found; nothing to compare against")

    with np.load(path, allow_pickle=True) as d:
        done = d["done_mask"]

    solved = np.flatnonzero((done & ~base.off_grid).ravel())
    if solved.size == 0:
        raise SystemExit(f"{path.name}: no solved node")
    return solved[np.linspace(0, solved.size - 1, min(N_NODE, solved.size)).astype(int)]


# ---- probes ----------------------------------------------------------

def probe_secluded(task):
    """(p, k, l, mX, mY) -> Omega_c h^2 for every eps in EPS_SECLUDED."""
    p, k, l, mX, mY = task
    base._set_portal(PORTALS_PROBED[p], 1.0)
    return p, k, l, base.omega_h2_secluded(mX, mY, base.COUPLING[l])


def probe_joint(task):
    """(p, k, l, e, mX, mY) -> (Omega_c h^2, exit-phase code) at one alpha.

    base.omega_h2_joint returns Omega alone, so its body is reproduced here
    with the exit phase kept.  The solver keywords, the Phase-3 keywords and
    the background fallback are the campaign's, so the Omega is the
    campaign's too.
    """
    p, k, l, e, mX, mY = task
    base._set_portal(PORTALS_PROBED[p], 1.0)
    eps = EPS_JOINT[e]
    alphaX = base.COUPLING[l]

    solver = base.BoltzmannSolver(base.cosmo, base._make_model(mX, mY, alphaX))
    try:
        joint = base.timed_solve(solver.solve_boltzmann_joint, epsX=eps,
                                 return_bg_ICs=True, min_eps_floor_ratio=0.0,
                                 GoH_exit_thresh=base.GOH_EXIT_THRESH,
                                 xmax=base.XMAX_JOINT, phase2_form="Y",
                                 verbose=False, **base.PHASE3_KW)
    except Exception:
        return p, k, l, e, np.nan, PATH_FAILED

    phase = joint["phase_final"]
    code = PATH_PHASE3 if phase == 3 else PATH_SMLOCK if phase == "1.5S" else PATH_LEGACY

    if phase in ("1.5S", 3):
        return p, k, l, e, joint["YX_relic"] * mX * base.s0_cosmo / base.rhoc, code

    try:
        yield_X = base.timed_solve(solver.solve_background, epsX=eps,
                                   bg_ICs=joint["bg_ICs"], correct_Y_decay=True,
                                   traj_has_decay=True, stop_on_rhoY=True,
                                   verbose=False)["YX_final"]
    except Exception:
        return p, k, l, e, np.nan, code
    return p, k, l, e, yield_X * mX * base.s0_cosmo / base.rhoc, code


# ---- reporting -------------------------------------------------------

def report(omega_sec, omega_joint, path_joint, portal_of):
    """Print the two answers the probe exists to give."""
    n_alpha = base.COUPLING.size

    print("\n==== SECLUDED: is the crossing extraction material? ====\n")
    print("portal   curves  non-monotone   both found   moved >5%   worst ratio")
    for p, name in enumerate(PORTALS_PROBED):
        sel = portal_of == p
        curves = omega_sec[sel].transpose(0, 2, 1).reshape(-1, n_alpha)

        nonmono = 0
        moved = 0
        both = 0
        worst = 1.0
        for om in curves:
            good = np.isfinite(om) & (om > 0)
            if good.sum() >= 2 and np.any(np.diff(om[good]) >= 0.0):
                nonmono += 1
            new = base._crossing(om)[0]
            old = crossing_global(om)
            if np.isfinite(new) and np.isfinite(old):
                both += 1
                ratio = new / old
                if abs(np.log10(ratio)) > np.log10(1.0 + MOVED_TOL):
                    moved += 1
                if abs(np.log10(ratio)) > abs(np.log10(worst)):
                    worst = ratio
        print(f"{name:8s} {len(curves):6d}  {nonmono:6d} ({100*nonmono/len(curves):4.0f}%)"
              f"  {both:10d}  {moved:6d} ({100*moved/max(both,1):4.0f}%)  {worst:12.4g}")

    print("\n==== JOINT: how often does Sdot=0.05 skip Phase 3? ====\n")
    print(f"PHASE3_KW = {base.PHASE3_KW}\n")
    print("portal   eps        solves   Phase3   SM-locked   LEGACY   failed")
    for p, name in enumerate(PORTALS_PROBED):
        sel = portal_of == p
        for e, eps in enumerate(EPS_JOINT):
            code = path_joint[sel, :, e].ravel()
            n = code.size
            counts = {c: int((code == c).sum()) for c in
                      (PATH_PHASE3, PATH_SMLOCK, PATH_LEGACY, PATH_FAILED)}
            print(f"{name:8s} {eps:.0e}  {n:7d}  {counts[PATH_PHASE3]:6d}"
                  f"  {counts[PATH_SMLOCK]:10d}  {counts[PATH_LEGACY]:6d}"
                  f"  {counts[PATH_FAILED]:6d}")


# ---- main ------------------------------------------------------------

def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else base.NCORES
    base.SOLVE_TIMEOUT = PROBE_TIMEOUT

    flat = []
    portal_of = []
    for p, name in enumerate(PORTALS_PROBED):
        idx = pick_nodes(name)
        flat.append(idx)
        portal_of.append(np.full(idx.size, p))
    flat = np.concatenate(flat)
    portal_of = np.concatenate(portal_of)

    i, j = np.unravel_index(flat, base.mX_grid.shape)
    mX, mY = base.mX_grid[i, j], base.mY_grid[i, j]
    rv = mY / mX

    n_node = flat.size
    n_alpha = base.COUPLING.size
    n_sec = len(base.EPS_SECLUDED)

    omega_sec = np.full((n_node, n_alpha, n_sec), np.nan)
    omega_joint = np.full((n_node, n_alpha, len(EPS_JOINT)), np.nan)
    path_joint = np.full((n_node, n_alpha, len(EPS_JOINT)), PATH_FAILED, np.int8)

    sec_tasks = [(portal_of[k], k, l, mX[k], mY[k])
                 for k in range(n_node) for l in range(n_alpha)]
    joint_tasks = [(portal_of[k], k, l, e, mX[k], mY[k])
                   for k in range(n_node) for l in range(n_alpha)
                   for e in range(len(EPS_JOINT))]

    print(f"{len(PORTALS_PROBED)} portals x {N_NODE} nodes = {n_node} nodes")
    print(f"alpha_X axis: {n_alpha} points, {base.COUPLING[0]:.3g} .. {base.COUPLING[-1]:.3g}")
    print(f"secluded {len(sec_tasks)} solves + joint {len(joint_tasks)} solves "
          f"on {ncores} workers\n")

    with mp.Pool(ncores) as pool:
        for p, k, l, om in tqdm(pool.imap_unordered(probe_secluded, sec_tasks, chunksize=1),
                                total=len(sec_tasks), desc="secluded", mininterval=10.0):
            omega_sec[k, l] = om

        for p, k, l, e, om, code in tqdm(pool.imap_unordered(probe_joint, joint_tasks,
                                                             chunksize=1),
                                         total=len(joint_tasks), desc="joint",
                                         mininterval=10.0):
            omega_joint[k, l, e] = om
            path_joint[k, l, e] = code

    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT, omega_sec=omega_sec, omega_joint=omega_joint,
                        path_joint=path_joint, coupling=base.COUPLING,
                        mX=mX, mY=mY, rv=rv, portal_of=portal_of,
                        portals=np.array(PORTALS_PROBED),
                        eps_secluded=np.array(base.EPS_SECLUDED),
                        eps_joint=np.array(EPS_JOINT))
    print(f"\nwrote {OUT}")

    report(omega_sec, omega_joint, path_joint, portal_of)


if __name__ == "__main__":
    main()

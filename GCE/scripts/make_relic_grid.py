"""Relic-abundance grid for any vector portal, on the shared (mX, rv) mesh.

    python GCE/scripts/make_relic_grid.py --portal BL --dry-run
    python GCE/scripts/make_relic_grid.py --portal LiLj --xi-ini 2.5 -n 128
    python GCE/scripts/make_relic_grid.py --portal vP --no-phase3     # legacy path

Generalises make_relic_grid_vP*.py, which this reproduces exactly for
--portal vP --no-phase3 (same mesh, same alpha_X axis, same solver keywords;
only the output file name differs).  The mass grid is deliberately identical
across portals so the slices can be compared node for node.

ALL FOUR VECTOR PORTALS SHARE THE HIDDEN SECTOR.  sigmav_XX_to_YY, the 3->2
cannibal rates and the freeze-out are byte-identical between VectorPortal,
BLPortal, LiLjPortal and BaryonPortal -- the classes differ only in
decay_width_to_SM, i.e. in Gamma_Y and in the branching ratios.  So the relic
grids differ ONLY through mediator decay: the entropy injected when Y decays,
and whether Phase 3 engages.  At a fixed epsX the widths sit within a factor
~3-14 of the vector portal, so the epsX axis below is in range for all of them.

The branching ratios, which the GCE spectra downstream do care about, are not
similar at all: LiLj(mu-tau) is 50% neutrinos and has no hadronic channel,
BaryonPortal is purely hadronic, BLPortal is both, VectorPortal is
charged-fermion only.

PHASE 3 IS ON BY DEFAULT HERE, unlike the restored make_relic_grid_vP*.py.
These are new grids with no legacy neighbours to match, and the 2026-10
campaign measured the conserved-S approximation biasing alpha_relic by
x1.07-1.20 at rv >= 0.465 (docs/phase3_rerun_campaign.md).  Pass --no-phase3
to reproduce the legacy path.

Outputs land in GCE/output/relic as relic_grid_<tag>_eps{...}.npz, where the
tag encodes portal, xi_ini and any portal parameter left off its default, so
two configurations can never collide in one file.
"""

import argparse
import inspect
import multiprocessing as mp
import shutil
import signal
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hidden_sector_DM.HiddenSectorDM import (
    BaryonPortal, BLPortal, BoltzmannSolver, Cosmology, LiLjPortal,
    VectorPortal, rhoc, s0_cosmo,
)

# ---- configuration ---------------------------------------------------

# ---- portals ---------------------------------------------------------
# name -> (class, portal-specific defaults).  The hidden-sector arguments
# (mX, mY, gX, gY, alphaX, include_antiparticles*, xi_ini) are common and are
# added by _make_model.  m_* defaults follow the model classes, except that
# LiLjPortal's m_ZLiLj default is 1 TeV in model.py while the other two are
# 10 TeV; it is raised here so every portal sits in the same decoupled-mediator
# limit and the slices stay comparable.  Override with --g-portal/--m-mediator.
PORTALS = {
    "vP":   (VectorPortal,  {}),
    "BL":   (BLPortal,      {"g_BL":   1.0, "m_ZBL":   1.0e4}),
    "LiLj": (LiLjPortal,    {"g_LiLj": 1.0, "m_ZLiLj": 1.0e4, "flavor": "mu-tau"}),
    "B":    (BaryonPortal,  {"g_B":    1.0, "m_ZB":    1.0e4}),
}
_G_KEY = {"BL": "g_BL", "LiLj": "g_LiLj", "B": "g_B"}       # --g-portal target
_M_KEY = {"BL": "m_ZBL", "LiLj": "m_ZLiLj", "B": "m_ZB"}    # --m-mediator target

# Set from the command line in main(), before any worker is forked.
PORTAL       = "vP"
PORTAL_KW    = {}
XI_INI       = 1.0                       # T_D/T_SM at Tinf
TAG          = "vP"                      # output-name stem; built by _set_portal
OCH2_TARGET  = 0.12
EPS_JOINT    = (1e-9, 1e-10)             # joint-solver epsilons
RUN_SECLUDED = True                      # also run the secluded slices
EPS_SECLUDED = (1e-11, 5e-12, 1e-12)     # one freeze-out feeds all three; rv=0.95 is
                                         # masked NaN in the output (solver does not
                                         # converge there for alpha >~ 0.05)

MX_LIM   = (15.0, 100.0)                 # GeV; matches notebooks/GCE_fit_vP
MY_LIM   = (5.0, 100.0)                  # GeV, capped from above by mX
RV_MAX   = 0.95                          # rv = mY/mX; XX -> YY shuts off at rv = 1
N_MX, N_RV, N_ALPHAX = 46, 67, 8         # coarse (mX, rv) mesh; alpha_X per node.
                                         # N_RV STAYS 67: the eps=1e-9 output is now
                                         # 67 wide, and 67 does not embed in 34, so
                                         # reverting would make load_state discard it
                                         # and re-solve all 1053 nodes.
ONLY_RV_COLUMNS = None                   # None = queue every rv column
ALPHAX_LIMS = {                          # per-output alpha_X axis; load_state checks it, so
    None:        (2e-5, 2e-2),           # changing an entry resets that output's file.
    1e-9:        (2e-5, 1.0),            # KEEP.  Was (2e-5, 2e-2); the Phase-3 campaign
                                         # moved crossings up to 1.6e-2 and the stored
                                         # files now carry THIS axis.  load_state compares
                                         # it exactly, so narrowing it again would reset
                                         # every eps=1e-9 output.
    1e-10:       (2e-5, 1.0),            # 1e-10 and secluded need alpha > 2e-2 at rv >= 0.87
    "secluded":  (2e-5, 1.0),            # (the old secluded file EXTRAPOLATED 78 such nodes);
}                                        # capped at 1 as for 1e-10. None = other joint eps

SOLVE_TIMEOUT = 7200                     # s; Phase-3 nodes were observed to
                                         # ~6200 s in the 2026-10 campaign
GOH_EXIT_THRESH = 1e99                   # converged mode: never exit on Gamma_Y/H (2026-09 audit:
                                         # gated handoffs bias Omega by 1.15-3.2x, smoothly in rv)
XMAX_JOINT = 5e3                         # converged to <0.3% vs 1.5e4 at both rv extremes
RETRY_RTOL    = 1e-6                     # one retry at this rtol_value after a timeout
SOLVER_KW = dict(cannibal_switch_full=1, convergence_threshold=1e-2,
                 rtol_value=1e-6)        # 1e-8 exceeds any timeout at rv=0.95 (>45 min);
                                         # 1e-6 solves those in ~1-6 min with identical
                                         # handoffs and Omega to 0.03% (2026-09 probe)


PHASE3_KW = dict(rho_switch_thresh=0.3,  # hand the joint solve to Phase 3 (the
                 Sdot_switch_thresh=0.05,  # in-solver entropy-injection era) when
                 enable_phase3=True,     # rho_D/rho_SM or dlnS/dN crosses these.
                 phase3_kw=dict(rtol_value=1e-6))  # Omega is rtol-stable over
                                         # 1e-6..1e-8 (validation D; campaign
                                         # probe P2 measured 1.00000 at x15.3
                                         # speedup).  --no-phase3 clears this.

NCORES = 4                               # default worker count; -n overrides
CHECKPOINT_EVERY = 8                     # completed nodes between saves
SEED_FROM_COARSE = True                  # reuse a coarser run's finished nodes
BACKUP_PREVIOUS  = False                 # history lives on cluster-relic-vP commits
SEC_PER_NODE = 3000.0                    # converged mode, rough; low-rv nodes are ~5x faster

# <sigma v>: GeV^-2 -> cm^3/s   (identical to GCE.spectrum.GEVM2_TO_CM3S)
GEVM2_TO_CM3S = 1.9732698e-14 ** 2 * 2.99792458e10

OUT_DIR = Path(__file__).resolve().parents[1] / "output" / "relic"

# ---- fixed inputs, built once --------------------------------------

cosmo = Cosmology()

log_mX = np.linspace(*np.log10(MX_LIM), N_MX)
log_rv = np.linspace(np.log10(MY_LIM[0] / MX_LIM[1]), np.log10(RV_MAX), N_RV)
mX_grid, rv_grid = np.meshgrid(10.0 ** log_mX, 10.0 ** log_rv, indexing="ij")
mY_grid = rv_grid * mX_grid
off_grid = mY_grid < MY_LIM[0]                       # below the mY floor
mX_grid[off_grid] = np.nan
mY_grid[off_grid] = np.nan
N_VALID = int(np.isfinite(mX_grid).sum())

def _set_alphax(eps=None):
    """Select the alpha_X axis for the output about to run.  Pool workers fork
    after this is called, so the module global is what node_* and relic_alphaX
    see; load_state and _save read the same one."""
    global ALPHAX
    ALPHAX = np.logspace(*np.log10(ALPHAX_LIMS.get(eps, ALPHAX_LIMS[None])),
                         N_ALPHAX)


_set_alphax()


# ---- Boltzmann solves ---------------------------------------------

def _raise_timeout(signum, frame):
    raise TimeoutError


signal.signal(signal.SIGALRM, _raise_timeout)


def _baseline_rtol(method):
    """The rtol_value `method` would use if we do not pass one."""
    try:
        return inspect.signature(method).parameters["rtol_value"].default
    except (ValueError, TypeError, KeyError):
        return None


def _log_nan(mX, mY, alphaX, eps, exc):
    """A failed node used to become NaN silently; the Phase-3 S_corr_max guard
    (ValueError) and timeouts must stay visible in the job log."""
    print(f"NaN node mX={mX:.3f} mY={mY:.3f} alphaX={alphaX:.3e} eps={eps:g}: "
          f"{type(exc).__name__}: {exc}", flush=True)


def timed_solve(method, **kwargs):
    """Call a BoltzmannSolver method, abandoning it after SOLVE_TIMEOUT seconds
    (the near-degenerate, low-mass corner can be pathologically stiff).

    A timed-out solve gets one more attempt at `rtol_value=RETRY_RTOL`, with a
    fresh SOLVE_TIMEOUT: the stiff corner usually crawls rather than diverges,
    and a looser tolerance is enough to get it over the line.  The retry is
    skipped when RETRY_RTOL is not actually looser than the tolerance already in
    force -- `solve_background` already defaults to 1e-6, so retrying it at the
    same tolerance would only burn another SOLVE_TIMEOUT for an identical run.
    """
    def _run(**kw):
        signal.alarm(SOLVE_TIMEOUT)
        try:
            return method(**kw)
        finally:
            signal.alarm(0)

    try:
        return _run(**kwargs)
    except TimeoutError:
        in_force = kwargs.get("rtol_value", _baseline_rtol(method))
        if in_force is None or RETRY_RTOL <= in_force:
            raise
        return _run(**{**kwargs, "rtol_value": RETRY_RTOL})


def _make_model(mX, mY, alphaX, antiY=False):
    """The portal object for one grid node, under the active configuration."""
    cls, _ = PORTALS[PORTAL]
    return cls(mX=mX, mY=mY, gX=2, gY=3, alphaX=alphaX, xi_ini=XI_INI,
               include_antiparticlesX=True, include_antiparticlesY=antiY,
               **PORTAL_KW)


def _set_portal(portal, xi_ini, g_portal=None, m_mediator=None, flavor=None):
    """Fix PORTAL / PORTAL_KW / XI_INI / TAG for this run.  The tag carries
    every setting that changes the numbers, so differing runs cannot share a
    file: portal, xi_ini, and any portal parameter left off its default."""
    global PORTAL, PORTAL_KW, XI_INI, TAG
    if portal not in PORTALS:
        raise SystemExit(f"unknown portal {portal!r}; choose from {sorted(PORTALS)}")
    PORTAL = portal
    PORTAL_KW = dict(PORTALS[portal][1])
    XI_INI = float(xi_ini)

    bits = [portal]
    if portal != "vP":
        if g_portal is not None:
            PORTAL_KW[_G_KEY[portal]] = float(g_portal)
        if m_mediator is not None:
            PORTAL_KW[_M_KEY[portal]] = float(m_mediator)
        if flavor is not None:
            if portal != "LiLj":
                raise SystemExit("--flavor applies to --portal LiLj only")
            PORTAL_KW["flavor"] = flavor
        defaults = PORTALS[portal][1]
        if "flavor" in PORTAL_KW:                       # always named, default or not
            bits.append(str(PORTAL_KW["flavor"]))
        for k in sorted(PORTAL_KW):
            if k != "flavor" and PORTAL_KW[k] != defaults[k]:
                bits.append(f"{k}_{PORTAL_KW[k]:g}")
    elif g_portal is not None or m_mediator is not None or flavor is not None:
        raise SystemExit("--g-portal/--m-mediator/--flavor do not apply to --portal vP")
    if XI_INI != 1.0:
        bits.append(f"xi_ini_{XI_INI:g}")
    TAG = "_".join(bits)


def omega_h2_secluded(mX, mY, alphaX):
    """Omega_c h^2 for each eps in EPS_SECLUDED, from one eps-independent
    secluded freeze-out with the background evolution re-run per eps.  A failed
    solve leaves NaN in that slot."""
    solver = BoltzmannSolver(cosmo, _make_model(mX, mY, alphaX))
    try:
        bg_ics = timed_solve(solver.solve_boltzmann_chempot_3phase,
                             return_bg_ICs=True, verbose=False,
                             **SOLVER_KW)["bg_ICs"]
    except Exception:
        return np.full(len(EPS_SECLUDED), np.nan)

    omega = np.full(len(EPS_SECLUDED), np.nan)
    for i, eps in enumerate(EPS_SECLUDED):
        try:
            yield_X = timed_solve(solver.solve_background, epsX=eps,
                                  bg_ICs=bg_ics, correct_Y_decay=True,
                                  stop_on_rhoY=True, verbose=False)["YX_final"]
            omega[i] = yield_X * mX * s0_cosmo / rhoc
        except Exception:
            pass
    return omega


def omega_h2_joint(mX, mY, alphaX, eps):
    """Omega_c h^2 from the joint solver.  The background step is skipped when
    the run ends SM-locked (phase 1.5S).  Returns NaN on a failed solve."""
    solver = BoltzmannSolver(cosmo, _make_model(mX, mY, alphaX))
    try:
        joint = timed_solve(solver.solve_boltzmann_joint, epsX=eps,
                            return_bg_ICs=True, min_eps_floor_ratio=0.0,
                            GoH_exit_thresh=GOH_EXIT_THRESH, xmax=XMAX_JOINT,
                            phase2_form="Y", verbose=False, **PHASE3_KW)
    except Exception as e:
        _log_nan(mX, mY, alphaX, eps, e)
        return np.nan

    # phase_final == 3: Phase 3 already applied decay + dilution in-solver and
    # YX_relic is Omega-ready, exactly like the SM-locked (1.5S) exit.
    if joint["phase_final"] in ("1.5S", 3):
        return joint["YX_relic"] * mX * s0_cosmo / rhoc
    try:
        yield_X = timed_solve(solver.solve_background, epsX=eps,
                              bg_ICs=joint["bg_ICs"], correct_Y_decay=True,
                              traj_has_decay=True, stop_on_rhoY=True,
                              verbose=False)["YX_final"]
    except Exception as e:
        _log_nan(mX, mY, alphaX, eps, e)
        return np.nan
    return yield_X * mX * s0_cosmo / rhoc


# ---- coupling at the relic target -------------------------------

def relic_alphaX(omega):
    """alpha_X giving Omega_c h^2 = OCH2_TARGET, from a cubic spline of
    log(alpha_X) against log(Omega h^2) over the sampled ALPHAX; NaN when the
    target is not bracketed by the sampled range (no extrapolation)."""
    finite = np.isfinite(omega) & (omega > 0)
    if finite.sum() < 4:
        return np.nan
    if not (omega[finite].min() <= OCH2_TARGET <= omega[finite].max()):
        return np.nan            # crossing outside the sampled range: report
                                 # honestly as absent instead of extrapolating
    order = np.argsort(omega[finite])                        # Omega h^2 ascending
    spline = CubicSpline(np.log(omega[finite][order]),
                         np.log(ALPHAX[finite][order]))
    return float(np.exp(spline(np.log(OCH2_TARGET))))


def sigmav_swave(mX, mY, alphaX):
    """s-wave <sigma v>(XX -> YY) at `alphaX`, in cm^3/s (NaN in -> NaN out)."""
    if not np.isfinite(alphaX):
        return np.nan
    return _make_model(mX, mY, alphaX).sigmav_XX_to_YY_swave() * GEVM2_TO_CM3S


# ---- one node per task (module level: Pool has to pickle these) ------

def node_secluded(task):
    """(i, j) -> alpha_X_relic and <sigma v> for every eps in EPS_SECLUDED."""
    i, j, mX, mY = task
    omega = np.array([omega_h2_secluded(mX, mY, a) for a in ALPHAX])   # (N_ALPHAX, 3)
    alpha = np.array([relic_alphaX(omega[:, e]) for e in range(len(EPS_SECLUDED))])
    return i, j, alpha, np.array([sigmav_swave(mX, mY, a) for a in alpha])


def node_joint(task):
    """(i, j) -> alpha_X_relic and <sigma v> at a single eps."""
    i, j, mX, mY, eps = task
    omega = np.array([omega_h2_joint(mX, mY, a, eps) for a in ALPHAX])
    alpha = relic_alphaX(omega)
    return i, j, alpha, sigmav_swave(mX, mY, alpha)


# ---- state: load / seed / save -------------------------------------

def _embed(old_log, new_log):
    """Index in `new_log` of each entry of `old_log`, or None when the old mesh
    is not a subset of the new one."""
    idx = []
    for v in old_log:
        k = int(np.argmin(np.abs(new_log - v)))
        if abs(new_log[k] - v) > 1e-9:
            return None
        idx.append(k)
    return np.array(idx)


def load_state(path, shape, eps, backup=True):
    """(alpha_relic, sigmav, done) on the CURRENT mesh, seeded from `path` when
    it holds the same alpha_X axis and eps on an embeddable coarser mesh.
    `done` is per node, shape (N_MX, N_RV); off-grid nodes start done."""
    alpha = np.full(shape, np.nan)
    sigmav = np.full(shape, np.nan)
    done = off_grid.copy()
    if not (path.exists() and SEED_FROM_COARSE):
        return alpha, sigmav, done

    with np.load(path) as saved:
        # A file must not mix Phase-3 and legacy nodes: they differ by up to
        # x1.2 in alpha_relic (docs/phase3_rerun_campaign.md), so appending
        # under the other setting would make one grid two physics.
        want = bool(PHASE3_KW)
        if "phase3" in saved:
            ok = bool(saved["phase3"]) == want
            was = "ON" if saved["phase3"] else "OFF"
        else:
            # No key: written before this driver existed, i.e. by
            # make_relic_grid_vP*.py.  Those are the campaign's own files and
            # their provenance cannot be read back, so only the legacy path may
            # extend them -- refusing here is what protects the committed grids.
            ok = not want
            was = "unrecorded (pre-dates this driver)"
        if not ok:
            raise SystemExit(
                f"\n{path.name} was built with Phase 3 {was}, but this run has it "
                f"{'ON' if want else 'OFF'}.\nAppending would put two different "
                f"physics in one grid.  Either pass "
                f"{'--no-phase3' if want else 'Phase 3 on (drop --no-phase3)'}, "
                f"or pick a --portal/--xi-ini whose tag differs so this run gets "
                f"its own file.")
        if not (np.array_equal(saved["alphaX"], ALPHAX)
                and np.array_equal(saved["eps"], np.atleast_1d(eps))):
            print(f"  {path.name}: alpha_X axis or eps changed -- starting fresh")
            return alpha, sigmav, done
        ix = _embed(saved["log_mX"], log_mX)
        jx = _embed(saved["log_rv"], log_rv)
        if ix is None or jx is None:
            print(f"  {path.name}: mesh does not embed in {N_MX}x{N_RV} "
                  f"-- starting fresh")
            return alpha, sigmav, done

        if "done_mask" in saved:                       # written by this version
            old_done = saved["done_mask"]
        else:                                          # 16x12 pass: column counter
            old_done = np.zeros((saved["log_mX"].size, saved["log_rv"].size), bool)
            old_done[:int(saved["done"])] = True

        sel = np.ix_(ix, jx)
        if alpha.ndim == 3:                            # leading eps axis
            alpha[np.ix_(np.arange(shape[0]), ix, jx)] = saved["alpha_relic"]
            sigmav[np.ix_(np.arange(shape[0]), ix, jx)] = saved["sigmav"]
        else:
            alpha[sel] = saved["alpha_relic"]
            sigmav[sel] = saved["sigmav"]
        done[sel] |= old_done

        n_old = int(saved["log_mX"].size * saved["log_rv"].size)
        bak = path.with_name(f"{path.stem}_{n_old}nodes_backup.npz")
        if backup and BACKUP_PREVIOUS and not bak.exists():
            shutil.copy2(path, bak)
            print(f"  {path.name}: previous run copied to {bak.name}")

    done |= off_grid                                    # never queue an off-grid node
    print(f"  {path.name}: reusing {int((done & ~off_grid).sum())} finished nodes")
    return alpha, sigmav, done


def _save(path, eps, alpha_relic, sigmav, done):
    np.savez(path, log_mX=log_mX, log_rv=log_rv, mX_grid=mX_grid, mY_grid=mY_grid,
             alphaX=ALPHAX, eps=np.atleast_1d(eps), alpha_relic=alpha_relic,
             sigmav=sigmav, done_mask=done, xi_ini=XI_INI,
             phase3=bool(PHASE3_KW), portal=PORTAL,
             n_done=int((done & ~off_grid).sum()))


# ---- grid runs ----------------------------------------------

def _pending(done):
    """(i, j, mX, mY) for every node still to solve.

    With ONLY_RV_COLUMNS set, nodes outside those rv columns are skipped but
    left NOT done, so they stay queued for a later run."""
    todo = ~done
    if ONLY_RV_COLUMNS is not None:
        keep = np.zeros(N_RV, bool)
        keep[list(ONLY_RV_COLUMNS)] = True
        todo = todo & keep[None, :]
    return [(i, j, float(mX_grid[i, j]), float(mY_grid[i, j]))
            for i, j in zip(*np.where(todo))]


def _run_pool(fn, tasks, ncores, desc, on_result, save):
    """Map `fn` over `tasks`, checkpointing every CHECKPOINT_EVERY results."""
    if not tasks:
        print(f"  {desc}: nothing to do")
        return
    n = 0
    with mp.Pool(ncores) as pool:
        for res in tqdm(pool.imap_unordered(fn, tasks, chunksize=1),
                        total=len(tasks), desc=desc):
            on_result(res)
            n += 1
            if n % CHECKPOINT_EVERY == 0:
                save()
    save()


def run_secluded(ncores):
    """alpha_X_relic and <sigma v> for every eps in EPS_SECLUDED, one output
    file with eps as the leading array axis."""
    path = OUT_DIR / f"relic_grid_{TAG}_eps_secluded.npz"
    eps = np.array(EPS_SECLUDED)
    _set_alphax("secluded")
    alpha_relic, sigmav, done = load_state(path, (len(eps), N_MX, N_RV), eps)

    def on_result(res):
        i, j, a, sv = res
        alpha_relic[:, i, j] = a
        sigmav[:, i, j] = sv
        done[i, j] = True

    _run_pool(node_secluded, _pending(done), ncores, "secluded",
              on_result, lambda: _save(path, eps, alpha_relic, sigmav, done))
    print(f"secluded: {np.isfinite(sigmav).sum()}/{len(eps) * N_VALID} solved "
          f"-> {path.name}")


def run_joint(eps, ncores):
    """alpha_X_relic and <sigma v> at a single eps via the joint solver."""
    path = OUT_DIR / f"relic_grid_{TAG}_eps{eps:.0e}.npz"
    _set_alphax(eps)
    alpha_relic, sigmav, done = load_state(path, (N_MX, N_RV), eps)

    def on_result(res):
        i, j, a, sv = res
        alpha_relic[i, j] = a
        sigmav[i, j] = sv
        done[i, j] = True

    tasks = [t + (eps,) for t in _pending(done)]
    _run_pool(node_joint, tasks, ncores, f"joint {eps:.0e}",
              on_result, lambda: _save(path, eps, alpha_relic, sigmav, done))
    print(f"joint {eps:.0e}: {np.isfinite(sigmav).sum()}/{N_VALID} solved "
          f"-> {path.name}")


def dry_run():
    """Grid, what each output would reuse, and the wall-clock estimate."""
    todo = 0
    outputs = [(OUT_DIR / f"relic_grid_{TAG}_eps{e:.0e}.npz", e, (N_MX, N_RV))
               for e in EPS_JOINT]
    if RUN_SECLUDED:
        outputs.insert(0, (OUT_DIR / f"relic_grid_{TAG}_eps_secluded.npz",
                           np.array(EPS_SECLUDED), (len(EPS_SECLUDED), N_MX, N_RV)))
    for path, eps, shape in outputs:
        _set_alphax(eps if np.isscalar(eps) else "secluded")
        alpha_relic, sigmav, done = load_state(path, shape, eps, backup=False)
        left = len(_pending(done))
        todo += left
        print(f"  {path.name}: {left} of {N_VALID} nodes still to solve"
              + ("" if ONLY_RV_COLUMNS is None
                 else f"  (queue restricted to rv columns {tuple(ONLY_RV_COLUMNS)}"
                      f" = rv {', '.join(f'{10 ** log_rv[k]:.4f}' for k in ONLY_RV_COLUMNS)})"))
    per_out = SEC_PER_NODE / 3.0
    print(f"\n{todo} node-solves left ~ {todo * per_out / 3600:.1f} h serial, "
          f"{todo * per_out / 3600 / NCORES:.1f} h on {NCORES} cores "
          f"(at {SEC_PER_NODE:.0f} s per node across all three outputs)")


def main():
    global PHASE3_KW
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--portal", default="vP", choices=sorted(PORTALS),
                    help="which vector portal (default vP)")
    ap.add_argument("--xi-ini", type=float, default=1.0,
                    help="T_D/T_SM at Tinf (default 1.0)")
    ap.add_argument("--g-portal", type=float, default=None,
                    help="portal gauge coupling g_BL / g_LiLj / g_B (default 1.0)")
    ap.add_argument("--m-mediator", type=float, default=None,
                    help="heavy mediator mass m_ZBL / m_ZLiLj / m_ZB in GeV "
                         "(default 1e4)")
    ap.add_argument("--flavor", default=None,
                    choices=["e-mu", "mu-tau", "e-tau"],
                    help="lepton pair for --portal LiLj (default mu-tau)")
    ap.add_argument("--no-phase3", action="store_true",
                    help="legacy path: clear PHASE3_KW, as the restored "
                         "make_relic_grid_vP*.py drivers run")
    ap.add_argument("-n", "--ncores", type=int, default=NCORES,
                    help=f"worker processes (default {NCORES})")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the grid, what is reusable and the ETA, then exit")
    args = ap.parse_args()

    # Must happen before any Pool is forked: the workers inherit these globals.
    _set_portal(args.portal, args.xi_ini, args.g_portal, args.m_mediator,
                args.flavor)
    if args.no_phase3:
        PHASE3_KW = {}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cls, _ = PORTALS[PORTAL]
    if PORTAL == "vP" and PHASE3_KW and not args.dry_run:
        print("NOTE: --portal vP writes relic_grid_vP*.npz, the files the 2026-10\n"
              "      campaign produced on the LEGACY path.  load_state will refuse\n"
              "      to mix the two; pass --no-phase3 to extend them.")
    print(f"portal: {cls.__name__}  xi_ini={XI_INI:g}  "
          f"{PORTAL_KW if PORTAL_KW else '(no extra parameters)'}  |  "
          f"Phase 3 {'ON' if PHASE3_KW else 'OFF'}  |  tag '{TAG}'")
    axes = ", ".join(
        ("default" if k is None else k if isinstance(k, str) else f"eps={k:g}")
        + f" [{v[0]:g}, {v[1]:g}]" for k, v in ALPHAX_LIMS.items())
    print(f"grid: {N_MX} x {N_RV} = {N_VALID} valid (mX, mY) nodes, "
          f"mX in [{MX_LIM[0]:g}, {MX_LIM[1]:g}] GeV, rv <= {RV_MAX}  |  "
          f"{N_ALPHAX} alpha_X, axis {axes}  |  out -> {OUT_DIR}")
    if args.dry_run:
        dry_run()
        return
    print(f"workers: {args.ncores}")
    if RUN_SECLUDED:
        run_secluded(args.ncores)
    for eps in EPS_JOINT:
        run_joint(eps, args.ncores)


if __name__ == "__main__":
    main()

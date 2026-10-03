"""
alpha_X that reproduces Omega_c h^2 = 0.12 for the Vector (Hypercharge) portal
cascade  X Xbar -> Y Y -> SM, and the matching s-wave <sigma v>(XX -> YY) in
cm^3/s, over the GCE (mX, mY) grid.

    cosmo_env
    cd hidden_sector_DM
    python GCE/scripts/make_relic_grid_vP_xi_ini_0.5.py --dry-run   # grid, reuse, ETA
    python GCE/scripts/make_relic_grid_vP_xi_ini_0.5.py -n 4        # the run

XI_INI VARIANT: identical to make_relic_grid_vP.py except that every
VectorPortal is built with xi_ini = 0.5 (the initial dark/SM temperature
ratio T_D/T_SM at Tinf = 1e14 GeV; the baseline grid is xi_ini = 1).  The
output filenames carry the same tag so this run CANNOT seed from, resume
into, or overwrite the xi_ini = 1 grids.

Three files land in GCE/output/relic/:

  relic_grid_vP_xi_ini_0.5_eps1e-09.npz      joint solver
  relic_grid_vP_xi_ini_0.5_eps1e-10.npz      joint solver
  relic_grid_vP_xi_ini_0.5_eps_secluded.npz  secluded solver, eps = 1e-11, 5e-12, 1e-12
                                  (one freeze-out per node, background re-run
                                  per eps -- eps is the leading array axis)

Which solver: for eps <= 1e-11 the secluded and joint workflows agree to <~4 %
across the grid; for eps >= 1e-9 they differ by up to ~2x at rv >~ 0.7, so the
two joint epsilons get the full joint solve.

Per (mX, mY, eps): Omega h^2 is evaluated on a log-spaced alpha_X axis and a
cubic spline of log(alpha_X) vs log(Omega h^2) is read off at the target.  No
chi^2 here -- the GCE overlay is a separate notebook step.

Grid resolution
---------------
N_MX x N_RV is 46 x 34 (1053 valid nodes), up from the 16 x 12 (128 nodes) of
the first pass.  The coarse mesh was not converged: reading it onto the fine
(mX, mY) triangle grid gave answers that differed by up to ~50 % for mY/mX >~
0.85 depending on which way it was interpolated, because a single 16x12 rv
interval (0.748 -> 0.980) spans a factor ~4.5 drop in <sigma v>.

Keep N_MX in {16, 31, 46, 61} and N_RV in {12, 23, 34, 45}: each such mesh
contains the 16 x 12 one, so a previous run's finished nodes are reused instead
of resolved (`SEED_FROM_COARSE`).  The previous file is copied to
`*_<n>nodes_backup.npz` before the first overwrite.

CAREFUL: those sets nest inside the 16 x 12 mesh, NOT inside each other.  Only
12 of the 34 rv columns land on the 45-point mesh, so N_RV 34 -> 45 fails
`_embed` and restarts all 1053 nodes from scratch.  To REFINE an existing 46 x
34 run, bisect: N_RV 34 -> 67 puts every old column on an even index (old j ->
new 2j, verified exact) and every new column at the log-midpoint of a pair, so
all 1053 finished nodes are reused.

Campaign mode
-------------
`ONLY_RV_COLUMNS` restricts a run to a few rv column indices, so a bisected
mesh can be filled one column at a time instead of solving all 33 new ones
(1518 nodes/output).  Nodes outside those columns are left untouched -- NOT
marked done -- so a later run picks them up with no mask clearing.  Set
`EPS_JOINT`/`RUN_SECLUDED` to match, and restore all three when the campaign
ends (as after the 2026-09-25 and 2026-09-29 campaigns).

2026-09-30 campaign (DONE, commit a03364657f): N_RV=67, ONLY_RV_COLUMNS=(65,).
Constants restored; N_RV deliberately left at 67.  BEWARE what that means for
the OTHER two outputs: they are still only 34 columns deep, so on the 67-mesh
--dry-run now honestly reports ~1030 unsolved nodes each.  That is not a bug
and nothing is lost -- their 34 columns seed fine (34 embeds in 67) -- but do
not launch a bare run without reading the dry-run first.  To fill more columns
for eps=1e-9 only, set ONLY_RV_COLUMNS and EPS_JOINT/RUN_SECLUDED as below.
Details of that campaign: rv = 0.9085, the
log-midpoint of the 0.8689 -> 0.9500 gap, eps = 1e-9 only.  Reason: a 3-node
probe at that rv found alpha_relic sits 7.7 % BELOW the log-log interpolation
the GCE notebook assumes across that gap (0.9231/0.9223/0.9226 at mX =
39.6/43.0/46.8 -- mX-independent to 0.1 %), and propagating that moves the
near-diagonal GCE lobe from dchi2 = 6.68 to 6.17, i.e. across the 2-sigma
threshold.  46 nodes x 8 alpha; solves measured at 113-129 s each.

Runs on `multiprocessing.Pool`, one node per task, `-n/--ncores` workers.
State is a per-node boolean mask rewritten every CHECKPOINT_EVERY completions,
so a killed run resumes from the nodes it had actually finished (tasks complete
out of order, so this is a mask, not a column counter).  A solve exceeding
SOLVE_TIMEOUT is retried once at RETRY_RTOL and only then abandoned (that node
stays NaN).
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
    BoltzmannSolver, Cosmology, VectorPortal, rhoc, s0_cosmo,
)

# ---- configuration ---------------------------------------------------

XI_INI       = 0.5                      # T_D/T_SM at Tinf; baseline grid is 1.0
TAG          = f"xi_ini_{XI_INI:g}"       # keeps these outputs off the xi_ini=1 files
OCH2_TARGET  = 0.12
EPS_JOINT    = (1e-9, 1e-10)             # full joint solve (see "Which solver" above)
RUN_SECLUDED = True                      # also refresh the secluded output
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
ONLY_RV_COLUMNS = None                   # restrict the run to these rv column indices
                                         # (None = all); see "Campaign mode" above.
ALPHAX_LIMS = {                          # per-output alpha_X axis; load_state checks it, so
    None:        (2e-5, 2e-2),           # changing an entry resets that output's file.
    1e-10:       (2e-5, 1.0),            # 1e-10 and secluded need alpha > 2e-2 at rv >= 0.87
    "secluded":  (2e-5, 1.0),            # (the old secluded file EXTRAPOLATED 78 such nodes);
}                                        # capped at 1 as for 1e-10. None = other joint eps

SOLVE_TIMEOUT = 2700                     # s; converged-mode solves run to xmax at high rv
GOH_EXIT_THRESH = 1e99                   # converged mode: never exit on Gamma_Y/H (2026-09 audit:
                                         # gated handoffs bias Omega by 1.15-3.2x, smoothly in rv)
XMAX_JOINT = 5e3                         # converged to <0.3% vs 1.5e4 at both rv extremes
RETRY_RTOL    = 1e-6                     # one retry at this rtol_value after a timeout
SOLVER_KW = dict(cannibal_switch_full=1, convergence_threshold=1e-2,
                 rtol_value=1e-6)        # 1e-8 exceeds any timeout at rv=0.95 (>45 min);
                                         # 1e-6 solves those in ~1-6 min with identical
                                         # handoffs and Omega to 0.03% (2026-09 probe)

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


def omega_h2_secluded(mX, mY, alphaX):
    """Omega_c h^2 for each eps in EPS_SECLUDED, from one eps-independent
    secluded freeze-out with the background evolution re-run per eps.  A failed
    solve leaves NaN in that slot."""
    solver = BoltzmannSolver(cosmo, VectorPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=alphaX, xi_ini=XI_INI,
        include_antiparticlesX=True, include_antiparticlesY=False))
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
    solver = BoltzmannSolver(cosmo, VectorPortal(
        mX=mX, mY=mY, gX=2, gY=3, alphaX=alphaX, xi_ini=XI_INI,
        include_antiparticlesX=True, include_antiparticlesY=False))
    try:
        joint = timed_solve(solver.solve_boltzmann_joint, epsX=eps,
                            return_bg_ICs=True, min_eps_floor_ratio=0.0,
                            GoH_exit_thresh=GOH_EXIT_THRESH, xmax=XMAX_JOINT,
                            phase2_form="Y", verbose=False)
    except Exception:
        return np.nan

    if joint["phase_final"] == "1.5S":
        return joint["YX_relic"] * mX * s0_cosmo / rhoc
    try:
        yield_X = timed_solve(solver.solve_background, epsX=eps,
                              bg_ICs=joint["bg_ICs"], correct_Y_decay=True,
                              traj_has_decay=True, stop_on_rhoY=True,
                              verbose=False)["YX_final"]
    except Exception:
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
    return VectorPortal(mX=mX, mY=mY, gX=2, gY=3, alphaX=alphaX, xi_ini=XI_INI,
                        include_antiparticlesX=True
                        ).sigmav_XX_to_YY_swave() * GEVM2_TO_CM3S


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
    path = OUT_DIR / f"relic_grid_vP_{TAG}_eps_secluded.npz"
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
    path = OUT_DIR / f"relic_grid_vP_{TAG}_eps{eps:.0e}.npz"
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
    outputs = [(OUT_DIR / f"relic_grid_vP_{TAG}_eps{e:.0e}.npz", e, (N_MX, N_RV))
               for e in EPS_JOINT]
    if RUN_SECLUDED:
        outputs.insert(0, (OUT_DIR / f"relic_grid_vP_{TAG}_eps_secluded.npz",
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
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("-n", "--ncores", type=int, default=NCORES,
                    help=f"worker processes (default {NCORES})")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the grid, what is reusable and the ETA, then exit")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
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

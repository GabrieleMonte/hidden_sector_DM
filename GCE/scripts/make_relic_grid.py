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
campaign measured the conserved-S approximation biasing coupling_relic by
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
from scipy.interpolate import PchipInterpolator
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
N_MX, N_RV, N_COUPLING = 46, 67, 8         # coarse (mX, rv) mesh; alpha_X per node.
                                         # N_RV STAYS 67: the eps=1e-9 output is now
                                         # 67 wide, and 67 does not embed in 34, so
                                         # reverting would make load_state discard it
                                         # and re-solve all 1053 nodes.
ONLY_RV_COLUMNS = None                   # None = queue every rv column

# Why a node has no coupling.  coupling_relic is NaN for all of these, which on
# its own cannot say whether the number is missing or whether no coupling
# exists -- opposite statements for a downstream fit.  Stored as `status`.
STATUS_UNRECORDED  = -1                  # written before status was recorded, or
                                         # never attempted (incl. off-grid nodes)
STATUS_OK          = 0                   # a crossing was found
STATUS_SOLVE_FAILED = 1                  # the solves failed where it mattered: too few
                                         # finite Omega samples to interpolate, or the
                                         # bracketing pair is not adjacent on the axis
                                         # because the samples between them failed
STATUS_NO_SOLUTION = 2                   # Omega never reaches the target on the
                                         # sampled axis; no coupling solves this
                                         # node.  Physics, not a shortfall.
STATUS_OFF_AXIS    = 3                   # bracketed, but the interpolant landed
                                         # outside the sampled axis

STATUS_NAMES = {STATUS_UNRECORDED: "unrecorded", STATUS_OK: "ok",
                STATUS_SOLVE_FAILED: "solve_failed",
                STATUS_NO_SOLUTION: "no_solution", STATUS_OFF_AXIS: "off_axis"}
COUPLING_LIMS = {                          # per-output alpha_X axis; load_state checks it, so
    None:        (2e-5, 2e-2),           # changing an entry resets that output's file.
    1e-9:        (2e-5, 1.0),            # KEEP.  Was (2e-5, 2e-2); the Phase-3 campaign
                                         # moved crossings up to 1.6e-2 and the stored
                                         # files now carry THIS axis.  load_state compares
                                         # it exactly, so narrowing it again would reset
                                         # every eps=1e-9 output.
    1e-10:       (2e-5, 1.0),            # 1e-10 and secluded need alpha > 2e-2 at rv >= 0.87
    "secluded":  (2e-5, 1.0),            # (the old secluded file EXTRAPOLATED 78 such nodes);
}                                        # capped at 1 as for 1e-10. None = other joint eps

COUPLING_SWEEP_STOP = True               # solve the axis downward and stop one point past
                                         # the crossing, leaving the rest NaN.  The axis
                                         # itself is unchanged, so stored files still load.
                                         # make_relic_grid_vP.py calls this ALPHAX_SWEEP_STOP
                                         # but sweeps upward, which trims the cheap end.

SOLVE_TIMEOUT = 600                      # s; nothing that SUCCEEDS comes near this.
                                         # Both 2026-10 probes capped every solve at
                                         # 600 s and timed out zero of 2052 solves,
                                         # across both solvers and all four portals;
                                         # and over 95 finished eps=1e-9/1e-10 nodes
                                         # the successful solves of a whole node
                                         # total 3.8 min at the median.  7200 s only
                                         # ever bought failures: 1862 of them in
                                         # 59773233 + 59784544, at 2 h each plus a
                                         # 2 h retry, which is the entire overrun.
GOH_EXIT_THRESH = 1e99                   # converged mode: never exit on Gamma_Y/H (2026-09 audit:
                                         # gated handoffs bias Omega by 1.15-3.2x, smoothly in rv)
XMAX_JOINT = 5e3                         # converged to <0.3% vs 1.5e4 at both rv extremes
RETRY_RTOL    = 1e-6                     # one retry at this rtol_value after a timeout
REQUIRE_CONVERGED = True                 # discard a solve whose own 'converged' flag is
                                         # False; see _NotConverged
ACCEPT_UNCONVERGED_PHASES = (2, 3)       # ...except a joint run that got this far.
                                         # solve_boltzmann_joint reports converged=False
                                         # for three unrelated outcomes and only one of
                                         # them is a failure:
                                         #
                                         #   phase 1 / 1.5  the phase-entry artifact --
                                         #     never integrated, muX pinned to 0.  Fixed
                                         #     at source in 1d99f58edb; still rejected.
                                         #   phase 2  completed Phase 2 and stopped at
                                         #     xmax because dlnS/dN never reached
                                         #     Sdot_switch_thresh.  That is the weak-
                                         #     injection regime where the conserved-S
                                         #     form is still valid and bg_ICs are good;
                                         #     omega_h2_joint's else branch finishes it.
                                         #     Measured 2026-10-22 at eps=1e-7 on 4 nodes
                                         #     against a forced Phase 3: Omega to 2.3-2.6%
                                         #     (0.6% on the coupling), control 0.14%.
                                         #   phase 3  Phase 3 ran and produced YX_relic,
                                         #     then exited without setting the flag.
                                         #     Measured at the 4 off-edge eps=1e-7 nodes
                                         #     that survived the top-up, each flanked on
                                         #     the axis by two CONVERGED samples: the
                                         #     rejected Omega lands on the power law its
                                         #     neighbours define (successive ratios
                                         #     7.04/7.26, 7.38/7.50, 6.65/6.57, 6.60/6.45
                                         #     against the 1.701^4 = 8.4 that Omega ~
                                         #     coupling^-4 gives per step), and the
                                         #     Phase-3-to-solve_background ratio is
                                         #     0.87-0.95 on them against 0.84-1.14 on the
                                         #     converged neighbours.  The flag carries no
                                         #     information about correctness here.  They
                                         #     also finish in 6 s, so there is no cost
                                         #     pathology either.
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
CHECKPOINT_EVERY = 1                     # completed nodes between saves; the stiff
                                         # eps=1e-9/1e-10 nodes cost ~8 h each, so 8
                                         # between saves leaves a shard writing nothing
                                         # for ~64 h and losing it all at the wall
RETRY_STATUSES = ()                      # status codes to re-queue from a stored grid
                                         # before solving; see _requeue.  Empty means
                                         # resume only the nodes never attempted.
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

def _set_coupling(eps=None):
    """Select the alpha_X axis for the output about to run.  Pool workers fork
    after this is called, so the module global is what node_* and relic_coupling
    see; load_state and _save read the same one."""
    global COUPLING
    COUPLING = np.logspace(*np.log10(COUPLING_LIMS.get(eps, COUPLING_LIMS[None])),
                         N_COUPLING)


_set_coupling()


# ---- Boltzmann solves ---------------------------------------------

class _NotConverged(RuntimeError):
    """A solve that ran to the end without meeting its convergence criterion.

    Both solvers report this themselves in result['converged'], and the value
    they return when it is False is not a freeze-out: if neither phase switch
    fires, X stays pinned to chemical equilibrium and the yield is Boltzmann-
    suppressed by orders of magnitude.  Treating those as valid is what put a
    flat, lam-independent shelf into the 2026-10 Higgs-portal grids and turned
    real crossings into "no solution".
    """


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
        freeze_out = timed_solve(solver.solve_boltzmann_chempot_3phase,
                                 return_bg_ICs=True, verbose=False,
                                 **SOLVER_KW)
    except Exception as e:
        _log_nan(mX, mY, alphaX, np.nan, e)
        return np.full(len(EPS_SECLUDED), np.nan)

    if REQUIRE_CONVERGED and not freeze_out["converged"]:
        _log_nan(mX, mY, alphaX, np.nan, _NotConverged("freeze-out"))
        return np.full(len(EPS_SECLUDED), np.nan)

    bg_ics = freeze_out["bg_ICs"]

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

    if (REQUIRE_CONVERGED and not joint["converged"]
            and joint["phase_final"] not in ACCEPT_UNCONVERGED_PHASES):
        _log_nan(mX, mY, alphaX, eps, _NotConverged("joint"))
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

def relic_coupling(omega):
    """The sampled coupling reproducing Omega_c h^2 = OCH2_TARGET, by
    shape-preserving interpolation of log(coupling) against log(Omega h^2)
    over COUPLING; NaN when the target is not bracketed (no extrapolation).

    The coupling is alpha_X for the vector portals and lam for the Higgs
    portal -- the name is from the former, which came first.

    Pchip, not a cubic spline.  A cubic spline is global: every knot is
    coupled through one linear system, so where Omega saturates at low
    coupling -- annihilation no longer setting the relic, a shelf carrying no
    crossing information -- its oscillation reaches the crossing at the far
    end of the axis.  Measured on the Higgs portal, that put 4-17% of nodes
    outside the sampled axis entirely, by up to 180 orders of magnitude.
    Pchip is local and shape-preserving, so on monotone data the result
    cannot leave the interval its neighbouring samples span.
    """
    return _crossing(omega)[0]


def _crossing(omega):
    """(coupling, status) -- relic_coupling's value and why it is what it is.

    Split out so the grid can record the reason without relic_coupling's own
    signature changing; see the STATUS_* constants.
    """
    finite = np.isfinite(omega) & (omega > 0)
    if finite.sum() < 2:
        return np.nan, STATUS_SOLVE_FAILED

    # The FIRST crossing, reading the coupling axis upward, and only the
    # samples around it.  Omega ~ coupling^-p should fall everywhere, but the
    # solved curves do not: 48 of 60 Higgs-portal probe curves rise somewhere
    # (2026-10-13).  Fitting the whole axis sorts those branches together, and
    # the interpolant can then return a coupling off the wrong one -- 10 of 30
    # crossings moved when the fit was localised, one of them by 17x, from
    # 8.9e-4 to the correct 1.5e-2.  A local window cannot make that mistake:
    # it never looks past the first crossing.  make_relic_grid_vP.py reaches
    # the same conclusion from the vector portals in relic_alphaX_local.
    # A crossing is an adjacent pair of SOLVED samples that straddles the
    # target downward.  Not simply the first sample below it: these curves can
    # open below the target on a branch that carries no crossing information,
    # and reading from there would report a coupling off that branch.
    idx = np.flatnonzero(finite)
    om_f = omega[idx]
    straddle = np.flatnonzero((om_f[:-1] >= OCH2_TARGET) & (om_f[1:] < OCH2_TARGET))
    if straddle.size == 0:
        return np.nan, STATUS_NO_SOLUTION   # never falls through the target
    c0 = int(straddle[0])

    lo, k = int(idx[c0]), int(idx[c0 + 1])  # the bracketing pair itself

    # ...and it has to be adjacent ON THE AXIS, not merely adjacent among the
    # samples that happen to be finite.  Both reasons a sample is NaN look the
    # same here, and they are not the same thing: _sweep SKIPS the low end
    # deliberately, always below the crossing, where bridging is harmless; a
    # FAILED solve can sit anywhere, including between the bracketing pair,
    # and then the crossing is read across a hole whose shape is unknown.
    # Measured by replaying each run's own NaN log against its crossings:
    # 56 of 2083 at eps=1e-9, 0 of 2083 at 1e-8, 384 of 1968 at 1e-7 -- and at
    # 1e-7 the bridged values came out ~2x below the same nodes solved with
    # the holes filled, with nothing in the file to show for it.  A log-log
    # read across one step of a near power law is accurate; across four steps
    # and four decades of Omega it is a guess.  Report the gap instead.
    if k != lo + 1:
        return np.nan, STATUS_SOLVE_FAILED

    win = idx[max(c0 - 1, 0):min(c0 + 3, idx.size)]
    om, ax = omega[win], COUPLING[win]

    if win.size >= 3 and lo >= win[0] and np.all(np.diff(om) < 0.0):
        order = np.argsort(om)                               # Omega h^2 ascending
        interp = PchipInterpolator(np.log(om[order]), np.log(ax[order]))
        crossing = float(np.exp(interp(np.log(OCH2_TARGET))))
    else:
        # Too few points, the bracket outside the window, or not monotone
        # through it: log-log linear on the bracketing pair, which is all the
        # crossing is actually determined by.
        t = np.log(OCH2_TARGET / omega[lo]) / np.log(omega[k] / omega[lo])
        crossing = float(np.exp(np.log(COUPLING[lo])
                                + t * np.log(COUPLING[k] / COUPLING[lo])))

    # It has to lie in the bracket it was read from.
    if not (COUPLING[lo] <= crossing <= COUPLING[k]):
        return np.nan, STATUS_OFF_AXIS
    return crossing, STATUS_OK


def sigmav_swave(mX, mY, alphaX):
    """s-wave <sigma v>(XX -> YY) at `alphaX`, in cm^3/s (NaN in -> NaN out)."""
    if not np.isfinite(alphaX):
        return np.nan
    return _make_model(mX, mY, alphaX).sigmav_XX_to_YY_swave() * GEVM2_TO_CM3S


# ---- one node per task (module level: Pool has to pickle these) ------

def _sweep(evaluate, ncol):
    """Omega_c h^2 over COUPLING, as a (N_COUPLING, ncol) array.

    `evaluate(coupling)` returns one Omega per column -- the secluded node
    solves every eps at once, the joint node a single one.

    With COUPLING_SWEEP_STOP the axis is solved from the LARGEST coupling
    downward and abandoned one point past the first Omega at or above
    OCH2_TARGET in every column, so the crossing is bracketed from above and
    the low end is never touched.  That end is where the joint solve is
    pathologically stiff: 91% of the 1862 timeouts in the 2026-10 campaign sat
    on the two lowest axis points, which cannot hold a crossing -- over 5257
    solved cells the smallest one found is 0.0022, and those points are 2e-5
    and 9.38e-5.  Skipped points stay NaN, which _crossing reads as "not
    sampled" rather than as a failed solve.

    Replayed against every complete curve on hand -- 96 joint and 144 secluded
    from the portals probe, 60 from the Higgs-portal lam probe -- this returns
    the identical crossing on all 257 of them, loses none, and leaves 29-45%
    of the axis unsolved.
    """
    omega = np.full((N_COUPLING, ncol), np.nan)
    if not COUPLING_SWEEP_STOP:
        for k, a in enumerate(COUPLING):
            omega[k] = evaluate(a)
        return omega

    bracketed = np.zeros(ncol, bool)
    for k in range(N_COUPLING - 1, -1, -1):
        omega[k] = evaluate(COUPLING[k])
        bracketed |= np.isfinite(omega[k]) & (omega[k] >= OCH2_TARGET)
        if bracketed.all():
            if k - 1 >= 0:                   # one more past the crossing
                omega[k - 1] = evaluate(COUPLING[k - 1])
            break
    return omega


def node_secluded(task):
    """(i, j) -> alpha_X_relic and <sigma v> for every eps in EPS_SECLUDED."""
    i, j, mX, mY = task
    omega = _sweep(lambda a: omega_h2_secluded(mX, mY, a),
                   len(EPS_SECLUDED))                                # (N_COUPLING, 3)
    found = [_crossing(omega[:, e]) for e in range(len(EPS_SECLUDED))]
    alpha = np.array([c for c, _ in found])
    status = np.array([s for _, s in found], np.int8)
    return i, j, alpha, np.array([sigmav_swave(mX, mY, a) for a in alpha]), status


def node_joint(task):
    """(i, j) -> alpha_X_relic and <sigma v> at a single eps."""
    i, j, mX, mY, eps = task
    omega = _sweep(lambda a: omega_h2_joint(mX, mY, a, eps), 1)[:, 0]
    alpha, status = _crossing(omega)
    return i, j, alpha, sigmav_swave(mX, mY, alpha), np.int8(status)


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
    """(coupling_relic, sigmav, done) on the CURRENT mesh, seeded from `path` when
    it holds the same alpha_X axis and eps on an embeddable coarser mesh.
    `done` is per node, shape (N_MX, N_RV); off-grid nodes start done."""
    alpha = np.full(shape, np.nan)
    sigmav = np.full(shape, np.nan)
    status = np.full(shape, STATUS_UNRECORDED, np.int8)
    done = off_grid.copy()
    if not (path.exists() and SEED_FROM_COARSE):
        return alpha, sigmav, status, done

    with np.load(path) as saved:
        # A file must not mix Phase-3 and legacy nodes: they differ by up to
        # x1.2 in coupling_relic (docs/phase3_rerun_campaign.md), so appending
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
        if not (np.array_equal(saved["coupling_axis"], COUPLING)
                and np.array_equal(saved["eps"], np.atleast_1d(eps))):
            print(f"  {path.name}: alpha_X axis or eps changed -- starting fresh")
            return alpha, sigmav, status, done
        ix = _embed(saved["log_mX"], log_mX)
        jx = _embed(saved["log_rv"], log_rv)
        if ix is None or jx is None:
            print(f"  {path.name}: mesh does not embed in {N_MX}x{N_RV} "
                  f"-- starting fresh")
            return alpha, sigmav, status, done

        if "done_mask" in saved:                       # written by this version
            old_done = saved["done_mask"]
        else:                                          # 16x12 pass: column counter
            old_done = np.zeros((saved["log_mX"].size, saved["log_rv"].size), bool)
            old_done[:int(saved["done"])] = True

        sel = np.ix_(ix, jx)
        if alpha.ndim == 3:                            # leading eps axis
            alpha[np.ix_(np.arange(shape[0]), ix, jx)] = saved["coupling_relic"]
            sigmav[np.ix_(np.arange(shape[0]), ix, jx)] = saved["sigmav"]
        else:
            alpha[sel] = saved["coupling_relic"]
            sigmav[sel] = saved["sigmav"]
        if "status" in saved:                          # absent in older grids,
            if alpha.ndim == 3:                        # which keep UNRECORDED
                status[np.ix_(np.arange(shape[0]), ix, jx)] = saved["status"]
            else:
                status[sel] = saved["status"]
        done[sel] |= old_done

        n_old = int(saved["log_mX"].size * saved["log_rv"].size)
        bak = path.with_name(f"{path.stem}_{n_old}nodes_backup.npz")
        if backup and BACKUP_PREVIOUS and not bak.exists():
            shutil.copy2(path, bak)
            print(f"  {path.name}: previous run copied to {bak.name}")

    done |= off_grid                                    # never queue an off-grid node
    print(f"  {path.name}: reusing {int((done & ~off_grid).sum())} finished nodes")
    return alpha, sigmav, status, done


def _save(path, eps, coupling_relic, sigmav, done, status=None):
    """Write one grid.  `status` says why each NaN in coupling_relic is NaN;
    omitting it records UNRECORDED rather than claiming something untrue."""
    if status is None:
        status = np.full(np.shape(coupling_relic), STATUS_UNRECORDED, np.int8)
    np.savez(path, log_mX=log_mX, log_rv=log_rv, mX_grid=mX_grid, mY_grid=mY_grid,
             coupling_axis=COUPLING, eps=np.atleast_1d(eps), coupling_relic=coupling_relic,
             sigmav=sigmav, done_mask=done, status=status, xi_ini=XI_INI,
             phase3=bool(PHASE3_KW), portal=PORTAL,
             n_done=int((done & ~off_grid).sum()))


# ---- grid runs ----------------------------------------------

def _requeue(done, status, codes):
    """Mark nodes carrying any of `codes` not-done, so _pending returns them.

    A top-up for nodes whose solves failed for a reason the next run changes
    -- a longer SOLVE_TIMEOUT, say.  Nodes that are off-grid, or that already
    carry a crossing, are never re-queued: `done` is only ever cleared where
    `status` says the node has no usable answer.  Returns how many.
    """
    if not codes:
        return 0
    hit = np.isin(status, list(codes))
    if hit.ndim > 2:                     # secluded: eps on the leading axis,
        hit = hit.any(axis=0)            # one `done` entry per (i, j)
    hit &= ~off_grid
    done[hit] = False
    return int(hit.sum())


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
    _set_coupling("secluded")
    coupling_relic, sigmav, status, done = load_state(path, (len(eps), N_MX, N_RV), eps)
    n_requeued = _requeue(done, status, RETRY_STATUSES)
    if n_requeued:
        print(f"  re-queued {n_requeued} node(s) with status "
              f"{[STATUS_NAMES[c] for c in RETRY_STATUSES]}")

    def on_result(res):
        i, j, a, sv, st = res
        coupling_relic[:, i, j] = a
        sigmav[:, i, j] = sv
        status[:, i, j] = st
        done[i, j] = True

    _run_pool(node_secluded, _pending(done), ncores, "secluded",
              on_result, lambda: _save(path, eps, coupling_relic, sigmav, done, status))
    print(f"secluded: {np.isfinite(sigmav).sum()}/{len(eps) * N_VALID} solved "
          f"-> {path.name}")


def run_joint(eps, ncores):
    """alpha_X_relic and <sigma v> at a single eps via the joint solver."""
    path = OUT_DIR / f"relic_grid_{TAG}_eps{eps:.0e}.npz"
    _set_coupling(eps)
    coupling_relic, sigmav, status, done = load_state(path, (N_MX, N_RV), eps)
    n_requeued = _requeue(done, status, RETRY_STATUSES)
    if n_requeued:
        print(f"  re-queued {n_requeued} node(s) with status "
              f"{[STATUS_NAMES[c] for c in RETRY_STATUSES]}")

    def on_result(res):
        i, j, a, sv, st = res
        coupling_relic[i, j] = a
        sigmav[i, j] = sv
        status[i, j] = st
        done[i, j] = True

    tasks = [t + (eps,) for t in _pending(done)]
    _run_pool(node_joint, tasks, ncores, f"joint {eps:.0e}",
              on_result, lambda: _save(path, eps, coupling_relic, sigmav, done, status))
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
        _set_coupling(eps if np.isscalar(eps) else "secluded")
        coupling_relic, sigmav, status, done = load_state(path, shape, eps, backup=False)
        _requeue(done, status, RETRY_STATUSES)   # or the estimate understates
        left = len(_pending(done))               # what the real run would solve
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
    ap.add_argument("--only", default=None,
                    choices=["secluded", "1e-9", "1e-10"],
                    help="run just one output.  The three write different "
                         "files, so they can be separate concurrent slices")
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
    if args.only is not None:
        global RUN_SECLUDED, EPS_JOINT
        RUN_SECLUDED = args.only == "secluded"
        EPS_JOINT = () if RUN_SECLUDED else (float(args.only),)

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
        + f" [{v[0]:g}, {v[1]:g}]" for k, v in COUPLING_LIMS.items())
    print(f"grid: {N_MX} x {N_RV} = {N_VALID} valid (mX, mY) nodes, "
          f"mX in [{MX_LIM[0]:g}, {MX_LIM[1]:g}] GeV, rv <= {RV_MAX}  |  "
          f"{N_COUPLING} alpha_X, axis {axes}  |  out -> {OUT_DIR}")
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

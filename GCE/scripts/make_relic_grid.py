"""Relic coupling and <sigma v> over the (mX, mY) mesh, for every portal in the paper.

    cosmo_env
    python GCE/scripts/make_relic_grid.py --portal bP --dry-run        # mesh, reuse, ETA
    python GCE/scripts/make_relic_grid.py --portal vP -n 128           # every output
    python GCE/scripts/make_relic_grid.py --portal vP --xi-ini 2.5 -n 128
    python GCE/scripts/make_relic_grid.py --portal hP --only 1e-7 -n 128

For one portal this finds the hidden-sector coupling that gives Omega_c h^2 = 0.12 --
alpha_X for the Z' portals (vP, bP, blP), lam for the Higgs portal (hP) -- and the
matching s-wave <sigma v>(XX -> YY), in three regimes of the portal coupling eps:

  thermal    eps >= eps_therm.  The hidden sector stays locked to the SM, so eps drops
             out (solve_boltzmann_thermSM).  One file, on the GCE fit grid itself.
  joint      eps near the thermalisation floor.  Y <-> SM decay and inverse decay act
             during freeze-out (solve_boltzmann_joint).  One file per eps.
  secluded   eps well below it.  One freeze-out with Y stable, then the late Y decay
             and its entropy dilution (solve_boltzmann_chempot_3phase, then
             solve_background).  One file for all of the portal's secluded eps.

Everything that differs between portals sits in the PORTALS table below; the mesh,
the solvers, the checkpointing and the crossing extraction are shared.

THE Z' PORTALS SHARE THE HIDDEN SECTOR.  sigmav_XX_to_YY, the 3->2 cannibal rates and
the freeze-out are byte-identical between VectorPortal, BLPortal and BaryonPortal --
the classes differ only in decay_width_to_SM, i.e. in Gamma_Y and in the branching
ratios.  So their joint and secluded grids differ ONLY through mediator decay: the
entropy injected when Y decays, and whether Phase 3 engages.  At a fixed eps the
widths sit within a factor ~3-14 of the vector portal, so one eps axis serves all
three.  The thermal relic never sees the decay at all, so it is solved with
VectorPortal for every Z' portal.

PHASE 3 IS ON BY DEFAULT.  The 2026-10 campaign measured the conserved-S
approximation biasing coupling_relic by x1.07-1.20 at rv >= 0.465.  Pass --no-phase3
to run the legacy path.

OUTPUTS, in GCE/output/relic/.  The tag is the portal, plus _xi_ini_<x> when
xi_ini != 1, so two configurations can never share a file.

    relic_grid_<tag>_thermal.npz        mX, mY [GeV] (the GCE fit grid), coupling_relic,
                                        sigmav_relic [cm^3/s]; the coarse mesh it is
                                        interpolated from as mX_coarse, mY_coarse,
                                        coupling_coarse (and *_coarse_ext for hP)
    relic_grid_<tag>_eps<eps>.npz       one joint eps on the coarse (mX, rv) mesh
    relic_grid_<tag>_eps_secluded.npz   every secluded eps, eps on the leading axis

The joint and secluded files carry log_mX, log_rv, mX_grid, mY_grid [GeV],
coupling_axis (the sampled couplings), eps, coupling_relic, sigmav [cm^3/s], status
(see STATUS_*), done_mask, n_done, xi_ini, phase3 and portal.  They double as
checkpoints: a rerun resumes from them, and an interrupted run is just relaunched.

PROVENANCE OF THE SHIPPED GRIDS.  The hP grids and every thermal grid are this
script's output.  The vP, bP and blP joint and secluded grids were solved on
2026-10-10..13 by earlier revisions of it (and, for vP, by the per-xi_ini drivers it
replaced): before the phase-entry fix of 2026-10-15 (1d99f58edb, solver.py and here),
the local crossing read, the downward sweep and the 600 s timeout.  Their mesh and
coupling axis are the ones below, so they load and resume here, but a rerun will not
reproduce them bit for bit.  Their status is unrecorded (-1).  The vP files carry no
`phase3` key because their joint nodes mix routes: Phase 3 everywhere at
xi_ini = 2.5 and at xi_ini = 1, eps = 1e-9; the legacy route below rv = 0.636 at
xi_ini = 1, eps = 1e-10 (under 1 % in <sigma v>), and below rv = 0.465 (eps = 1e-9)
and 0.636 (eps = 1e-10) at xi_ini = 0.5 (at most 4 %).

Their SECLUDED files predate the convergence check as well.  At high rv the secluded
freeze-out does not converge for alpha >~ 0.05, Omega(alpha) folds, and the
whole-axis spline those revisions used can return a wrong coupling inside the sampled
axis.  So they carry values in the rv = 0.95 row, and a minority of nodes at
rv >~ 0.6 differ from their neighbours by more than 20 % (most at rv >= 0.87).  This
script records such samples as solve_failed instead.
"""

import argparse
import inspect
import multiprocessing as mp
import shutil
import signal
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator, griddata
from scipy.optimize import brentq
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))   # sibling scripts

from hidden_sector_DM.HiddenSectorDM import (                  # noqa: E402
    BaryonPortal, BLPortal, BoltzmannSolver, Cosmology, HiggsPortal, VectorPortal,
    rhoc, s0_cosmo,
)

# ---- configuration ---------------------------------------------------

# ---- portals ---------------------------------------------------------
# Everything that differs between portals.  _set_portal copies one entry into the
# module globals below before any worker is forked, so the pool inherits it.
#
#   cls, kw         the model class and its portal-specific defaults.  The Z' gauge
#                   couplings default to 1, so eps is the product eps * g, and the
#                   heavy Z' sits at 10 TeV, in the decoupled limit.
#   gX, gY          DM and mediator degrees of freedom.
#   dirac_X         include_antiparticlesX: the 1/2 on XX -> YY and YXX -> XX.
#   coupling        the keyword the swept hidden-sector coupling is passed as.
#   coupling_lims   sampled coupling axis per output as (lo, hi), keyed by eps, by
#                   "secluded", or by None for any other eps.  load_state compares
#                   the axis exactly, so changing an entry resets that output's file.
#   n_coupling      points on that axis.
#   eps_joint       eps solved by the joint solver, one file each.
#   eps_secluded    eps served by one secluded freeze-out, one file together.
#   phase3_require  every joint sample must end in Phase 3; see PHASE3_REQUIRE.
#   extension       continue the mesh past mX = 100 GeV; see N_MX_EXT.
#   thermal         the thermal relic's model class, and the coupling range,
#                   starting guess and n in Omega ~ coupling^-n of its
#                   root-find; see thermal_coupling.
_ZPRIME = dict(
    gX=2, gY=3, dirac_X=True, coupling='alphaX',
    coupling_lims={
        None:       (2e-5, 2e-2),
        1e-9:       (2e-5, 1.0),   # was (2e-5, 2e-2); the Phase-3 campaign moved
                                   # crossings up to 1.6e-2
        1e-10:      (2e-5, 1.0),   # 1e-10 and secluded need alpha > 2e-2 at rv >= 0.87
        'secluded': (2e-5, 1.0),   # (the old secluded file EXTRAPOLATED 78 such
                                   # nodes); capped at the perturbative alpha_X = 1
    },
    n_coupling=8,
    eps_joint=(1e-9, 1e-10),
    eps_secluded=(1e-11, 5e-12, 1e-12),  # at rv = 0.95 the freeze-out does not
                                         # converge for alpha >~ 0.05; those samples
                                         # are recorded as solve_failed
    phase3_require=False,
    extension=False,
    # The thermal relic sees only the hidden sector, which the three share, so it is
    # solved with VectorPortal for all of them (checked: the same Omega h^2 to 4
    # digits); the shipped bP and blP thermal files are copies of vP's.
    thermal=dict(cls=VectorPortal, coupling_range=(1e-5, 1.0), guess=1e-2,
                 omega_power=2),      # sigmav ~ alpha_X^2
)

PORTALS = {
    'vP':  dict(_ZPRIME, cls=VectorPortal, kw={}),
    'bP':  dict(_ZPRIME, cls=BaryonPortal, kw={'g_B': 1.0, 'm_ZB': 1.0e4}),
    'blP': dict(_ZPRIME, cls=BLPortal, kw={'g_BL': 1.0, 'm_ZBL': 1.0e4}),

    # X is a Majorana fermion and phi a real scalar.  The swept coupling is lam, which
    # sets lam_s = lam_p, and eps is the Higgs mixing angle.  Gamma_phi comes from
    # HDECAY (hidden_sector_DM/phi_decay); every node has mY >= 5 GeV.
    'hP': dict(
        cls=HiggsPortal, kw={}, gX=2, gY=1, dirac_X=False, coupling='lam',
        # sigmav ~ lam^4, so Omega ~ lam^-4: the Z' axis (4.7 decades of coupling)
        # would span 18.8 decades of Omega, most of it where annihilation no longer
        # sets the relic.  Surveys of 30 nodes on the fixed solver put every crossing
        # at eps = 1e-10..1e-7 in [0.0048, 0.32], so (1e-3, 1) brackets them with a
        # factor >= 3 to spare, and 14 points give a factor 1.70 per step.  The same
        # axis serves every output.
        coupling_lims={None: (1e-3, 1.0)},
        n_coupling=14,
        # At the Z' secluded eps every phi decays at or after the onset of BBN
        # (tau = 2.8-37700 s), where the relic is not the binding constraint, so the
        # Higgs portal stops at 1e-10.  1e-9 and 1e-10 are served by the SECLUDED
        # route: the sectors have already decoupled there, its Omega matches the
        # joint route's to a median 0.997 over 16 nodes at eps = 1e-9 (0.07 % on lam),
        # and it costs 4-5 s a solve against ~750 s.
        eps_joint=(1e-7, 1e-8),
        eps_secluded=(1e-9, 1e-10),
        phase3_require=True,
        # The hP favoured region reaches the main grid's mX = 100 GeV edge along the
        # mY ~ mX diagonal with its contours still open; the vector portals' close
        # inside it.
        extension=True,
        thermal=dict(cls=HiggsPortal, coupling_range=(2e-5, 1.0), guess=0.1,
                     omega_power=4),  # sigmav ~ lam_s^2 lam_p^2 = lam^4
    ),
}

# The active portal, set from the command line by _set_portal before any worker is
# forked.  The defaults are vP's, so the module imports cleanly.
PORTAL       = 'vP'
XI_INI       = 1.0                       # T_D/T_SM at Tinf
TAG          = 'vP'                      # output-name stem; built by _set_portal
OCH2_TARGET  = 0.12
RUN_THERMAL  = True                      # also run the thermal output (xi_ini = 1 only)
RUN_SECLUDED = True                      # also run the secluded output
EPS_JOINT    = PORTALS['vP']['eps_joint']
EPS_SECLUDED = PORTALS['vP']['eps_secluded']
COUPLING_LIMS = PORTALS['vP']['coupling_lims']
N_COUPLING   = PORTALS['vP']['n_coupling']

# ---- the coarse (mX, rv) mesh ---------------------------------------
MX_LIM   = (15.0, 100.0)                 # GeV, the GCE fit grid's main block
MY_LIM   = (5.0, 100.0)                  # GeV, capped from above by mX
RV_MAX   = 0.95                          # top mY/mX row; XX -> YY shuts off at rv = 1
N_MX_MAIN, N_RV = 46, 67                 # columns uniform in log10 mX, rows in log10 rv
N_MX_EXT = 6                             # extension columns past mX = 100 GeV at the
                                         # main log-mX step: mX = 104.31 ... 128.78 GeV;
                                         # a seventh would land at 134.3, past 130
RV_MIN_EXT = 0.6                         # lowest mY/mX solved in the extension: for
                                         # mX >= 92 the GCE 2- and 3-sigma regions sit at
                                         # rv >= 0.70, and the spectra cache holds that
                                         # block over rv 0.6-1 only.  mY is NOT capped at
                                         # 100 GeV there; it reaches 122.34.
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

STATUS_NAMES = {STATUS_UNRECORDED: 'unrecorded', STATUS_OK: 'ok',
                STATUS_SOLVE_FAILED: 'solve_failed',
                STATUS_NO_SOLUTION: 'no_solution', STATUS_OFF_AXIS: 'off_axis'}

COUPLING_SWEEP_STOP = True               # solve the axis downward and stop one point past
                                         # the crossing, leaving the rest NaN.  The axis
                                         # itself is unchanged, so stored files still load.

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

PHASE3_REQUIRE = False                   # True: a joint solve must end in Phase 3 (or
                                         # the SM-locked 1.5S exit).  One that would take
                                         # the legacy solve_background step is re-solved
                                         # with the Sdot trigger lowered through
                                         # PHASE3_RETRY_SDOT; if none engages the legacy
                                         # value is kept and logged ("LEGACY coupling").
                                         # Set per portal by _set_portal from
                                         # PORTALS[...]["phase3_require"]: off for the Z'
                                         # portals, on for the Higgs portal, where
                                         # at eps=1e-7 freeze-out sits right at the
                                         # Phase-3 handoff, so which route a sample takes
                                         # is close to a coin flip between neighbouring
                                         # coupling or rv, and the two routes disagree.
                                         # Audited 2026-10-26 over 25 nodes per eps
                                         # at eps=1e-7, 24 of
                                         # 181 solved samples change route, 7 of 25 nodes
                                         # have one, but only 2 move the crossing by >1%
                                         # and 1 by >5% -- a sample matters only when it
                                         # is one of the two the crossing is read from.
                                         # Worst 5.74% at mX=80.994 rv=0.727, the node
                                         # the 2026-10-26 report names.  At eps=1e-8:
                                         # 4 of 168 samples, worst 2.74%.  The median is
                                         # 1.0000 at both, so this wobbles a contour
                                         # node-to-node rather than biasing it.
PHASE3_RETRY_SDOT = (1e-4, 1e-5)         # retry ladder for PHASE3_REQUIRE.  NOT lower:
                                         # at 1e-7 the joint solve can hand off at a
                                         # phase entry without the final-state check
                                         # passing, i.e. an early legacy handoff.

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

# <sigma v>: GeV^-2 -> cm^3/s   (identical to GCE.source.spectrum.GEVM2_TO_CM3S)
GEVM2_TO_CM3S = 1.9732698e-14**2*2.99792458e10

OUT_DIR = Path(__file__).resolve().parents[1] / 'output' / 'relic'

# ---- fixed inputs, built once --------------------------------------

cosmo = Cosmology()


def _set_mesh(extension):
    """Build the coarse (mX, rv) mesh into the module globals.

    Columns are uniform in log10 mX over MX_LIM, rows uniform in log10 rv = mY/mX from
    the mY floor at mX = 100 GeV up to RV_MAX.  Nodes below the mY floor are off-grid:
    mX_grid and mY_grid NaN there, and load_state marks them done so no sweep picks
    them up.

    With `extension` the columns continue at the same log-mX step for N_MX_EXT more,
    so the seam at mX = 100 GeV has neither a gap nor a repeated column, and those
    columns keep only the rows with rv >= RV_MIN_EXT; the rest are off-grid as well.

    The thermal relic below uses a mesh of its own whose top rows are placed by beta
    rather than by rv (see thermal_rv_rows).  These rows stay uniform in log10 rv:
    every joint and secluded grid was solved on them.

    Must run before load_state and before any worker is forked: _pending, load_state
    and _save read these globals, and the pool inherits them.
    """
    global log_mX, log_rv, mX_grid, mY_grid, off_grid, N_MX, N_VALID
    log_mX = np.linspace(*np.log10(MX_LIM), N_MX_MAIN)
    log_rv = np.linspace(np.log10(MY_LIM[0]/MX_LIM[1]), np.log10(RV_MAX), N_RV)
    if extension:
        du = log_mX[1] - log_mX[0]
        log_mX = np.concatenate([log_mX, log_mX[-1] + du*np.arange(1, N_MX_EXT + 1)])
    N_MX = log_mX.size

    mX_grid, rv_grid = np.meshgrid(10.0**log_mX, 10.0**log_rv, indexing='ij')
    mY_grid = rv_grid*mX_grid
    off_grid = mY_grid < MY_LIM[0]                   # below the mY floor
    if extension:
        off_grid[N_MX_MAIN:] |= log_rv[None, :] < np.log10(RV_MIN_EXT)
    mX_grid[off_grid] = np.nan
    mY_grid[off_grid] = np.nan
    N_VALID = int(np.isfinite(mX_grid).sum())


_set_mesh(extension=False)


def _set_coupling(eps=None):
    """Select the coupling axis for the output about to run.  Pool workers fork
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
        return inspect.signature(method).parameters['rtol_value'].default
    except (ValueError, TypeError, KeyError):
        return None


def _log_nan(mX, mY, coupling, eps, exc):
    """A failed node used to become NaN silently; the Phase-3 S_corr_max guard
    (ValueError) and timeouts must stay visible in the job log."""
    print(f'NaN node mX={mX:.3f} mY={mY:.3f} coupling={coupling:.3e} eps={eps:g}: '
          f'{type(exc).__name__}: {exc}', flush=True)


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
        in_force = kwargs.get('rtol_value', _baseline_rtol(method))
        if in_force is None or RETRY_RTOL <= in_force:
            raise
        return _run(**{**kwargs, 'rtol_value': RETRY_RTOL})


def _make_model(mX, mY, coupling, antiY=False):
    """The portal object for one grid node, under the active configuration.
    mX, mY in GeV; `coupling` is alpha_X for the Z' portals and lam for hP."""
    p = PORTALS[PORTAL]
    return p['cls'](mX=mX, mY=mY, gX=p['gX'], gY=p['gY'], xi_ini=XI_INI,
                    include_antiparticlesX=p['dirac_X'], include_antiparticlesY=antiY,
                    **{p['coupling']: coupling}, **p['kw'])


def _set_portal(portal, xi_ini):
    """Point the module at one portal: its model, coupling axis, eps sets, Phase-3
    requirement and mesh, and the output tag.  The tag carries every setting that
    changes the numbers -- the portal and xi_ini -- so differing runs cannot share
    a file.  Must run before any worker is forked."""
    global PORTAL, XI_INI, TAG, EPS_JOINT, EPS_SECLUDED, COUPLING_LIMS, N_COUPLING
    global PHASE3_REQUIRE
    if portal not in PORTALS:
        raise SystemExit(f'unknown portal {portal!r}; choose from {sorted(PORTALS)}')
    p = PORTALS[portal]
    PORTAL = portal
    XI_INI = float(xi_ini)
    TAG = portal if XI_INI == 1.0 else f'{portal}_xi_ini_{XI_INI:g}'

    EPS_JOINT = p['eps_joint']
    EPS_SECLUDED = p['eps_secluded']
    COUPLING_LIMS = p['coupling_lims']
    N_COUPLING = p['n_coupling']
    PHASE3_REQUIRE = p['phase3_require']
    _set_mesh(p['extension'])
    _set_coupling()


def omega_h2_secluded(mX, mY, coupling):
    """Omega_c h^2 for each eps in EPS_SECLUDED, from one eps-independent
    secluded freeze-out with the background evolution re-run per eps.  A failed
    solve leaves NaN in that slot."""
    solver = BoltzmannSolver(cosmo, _make_model(mX, mY, coupling))
    try:
        freeze_out = timed_solve(solver.solve_boltzmann_chempot_3phase,
                                 return_bg_ICs=True, verbose=False,
                                 **SOLVER_KW)
    except Exception as e:
        _log_nan(mX, mY, coupling, np.nan, e)
        return np.full(len(EPS_SECLUDED), np.nan)

    if REQUIRE_CONVERGED and not freeze_out['converged']:
        _log_nan(mX, mY, coupling, np.nan, _NotConverged('freeze-out'))
        return np.full(len(EPS_SECLUDED), np.nan)

    bg_ics = freeze_out['bg_ICs']

    omega = np.full(len(EPS_SECLUDED), np.nan)
    for i, eps in enumerate(EPS_SECLUDED):
        try:
            yield_X = timed_solve(solver.solve_background, epsX=eps,
                                  bg_ICs=bg_ics, correct_Y_decay=True,
                                  stop_on_rhoY=True, verbose=False)['YX_final']
            omega[i] = yield_X*mX*s0_cosmo/rhoc
        except Exception:
            pass
    return omega


def omega_h2_joint(mX, mY, coupling, eps):
    """Omega_c h^2 from the joint solver.  The background step is skipped when
    the run ends in Phase 3 or SM-locked (phase 1.5S).  With PHASE3_REQUIRE, a
    run that would take the legacy background step is first re-solved down
    PHASE3_RETRY_SDOT.  Returns NaN on a failed solve."""
    solver = BoltzmannSolver(cosmo, _make_model(mX, mY, coupling))
    tries = [PHASE3_KW]
    if PHASE3_REQUIRE:
        tries += [{**PHASE3_KW, 'Sdot_switch_thresh': s} for s in PHASE3_RETRY_SDOT]

    first = None                         # the PHASE3_KW run: the legacy fallback
    for kw in tries:
        try:
            joint = timed_solve(solver.solve_boltzmann_joint, epsX=eps,
                                return_bg_ICs=True, min_eps_floor_ratio=0.0,
                                GoH_exit_thresh=GOH_EXIT_THRESH, xmax=XMAX_JOINT,
                                phase2_form='Y', verbose=False, **kw)
        except Exception as e:
            _log_nan(mX, mY, coupling, eps, e)
            if first is None:            # the PHASE3_KW run itself failed
                return np.nan
            continue                     # a failed retry = did not engage

        if (REQUIRE_CONVERGED and not joint['converged']
                and joint['phase_final'] not in ACCEPT_UNCONVERGED_PHASES):
            if first is None:
                _log_nan(mX, mY, coupling, eps, _NotConverged('joint'))
                return np.nan
            continue

        # phase_final == 3: Phase 3 already applied decay + dilution in-solver
        # and YX_relic is Omega-ready, exactly like the SM-locked (1.5S) exit.
        if joint['phase_final'] in ('1.5S', 3):
            return joint['YX_relic']*mX*s0_cosmo/rhoc
        if first is None:
            first = joint

    if PHASE3_REQUIRE:
        print(f'LEGACY coupling mX={mX:.3f} mY={mY:.3f} coupling={coupling:.3e} '
              f'eps={eps:g}: Phase 3 never engaged (Sdot down to '
              f"{PHASE3_RETRY_SDOT[-1]:g}, phase_final={first['phase_final']}); "
              f'legacy value kept', flush=True)
    try:
        yield_X = timed_solve(solver.solve_background, epsX=eps,
                              bg_ICs=first['bg_ICs'], correct_Y_decay=True,
                              traj_has_decay=True, stop_on_rhoY=True,
                              verbose=False)['YX_final']
    except Exception as e:
        _log_nan(mX, mY, coupling, eps, e)
        return np.nan
    return yield_X*mX*s0_cosmo/rhoc


# ---- coupling at the relic target -------------------------------

def relic_coupling(omega):
    """The sampled coupling reproducing Omega_c h^2 = OCH2_TARGET, by
    shape-preserving interpolation of log(coupling) against log(Omega h^2)
    over COUPLING; NaN when the target is not bracketed (no extrapolation).

    The coupling is alpha_X for the Z' portals and lam for the Higgs portal.

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
    # it never looks past the first crossing.
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
        t = np.log(OCH2_TARGET/omega[lo])/np.log(omega[k]/omega[lo])
        crossing = float(np.exp(np.log(COUPLING[lo])
                                + t*np.log(COUPLING[k]/COUPLING[lo])))

    # It has to lie in the bracket it was read from.
    if not (COUPLING[lo] <= crossing <= COUPLING[k]):
        return np.nan, STATUS_OFF_AXIS
    return crossing, STATUS_OK


def sigmav_swave(mX, mY, coupling):
    """s-wave <sigma v>(XX -> YY) at `coupling`, in cm^3/s (NaN in -> NaN out)."""
    if not np.isfinite(coupling):
        return np.nan
    return _make_model(mX, mY, coupling).sigmav_XX_to_YY_swave()*GEVM2_TO_CM3S


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
    """(i, j) -> coupling_relic and <sigma v> for every eps in EPS_SECLUDED."""
    i, j, mX, mY = task
    omega = _sweep(lambda a: omega_h2_secluded(mX, mY, a),
                   len(EPS_SECLUDED))                                # (N_COUPLING, 3)
    found = [_crossing(omega[:, e]) for e in range(len(EPS_SECLUDED))]
    couplings = np.array([c for c, _ in found])
    status = np.array([s for _, s in found], np.int8)
    return i, j, couplings, np.array([sigmav_swave(mX, mY, a) for a in couplings]), status


def node_joint(task):
    """(i, j) -> coupling_relic and <sigma v> at a single eps."""
    i, j, mX, mY, eps = task
    omega = _sweep(lambda a: omega_h2_joint(mX, mY, a, eps), 1)[:, 0]
    coupling, status = _crossing(omega)
    return i, j, coupling, sigmav_swave(mX, mY, coupling), np.int8(status)


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
    it holds the same coupling axis and eps on an embeddable coarser mesh.
    `done` is per node, shape (N_MX, N_RV); off-grid nodes start done."""
    coupling_relic = np.full(shape, np.nan)
    sigmav = np.full(shape, np.nan)
    status = np.full(shape, STATUS_UNRECORDED, np.int8)
    done = off_grid.copy()
    if not (path.exists() and SEED_FROM_COARSE):
        return coupling_relic, sigmav, status, done

    with np.load(path) as saved:
        # A file must not mix Phase-3 and legacy nodes: they differ by up to
        # x1.2 in coupling_relic (the 2026-10 Phase-3 campaign), so appending
        # under the other setting would make one grid two physics.
        want = bool(PHASE3_KW)
        if 'phase3' in saved:
            ok = bool(saved['phase3']) == want
            was = 'ON' if saved['phase3'] else 'OFF'
        else:
            # No key: the shipped vP grids, whose joint nodes mix the legacy and
            # Phase-3 routes (see PROVENANCE above).  That cannot be read back
            # node by node, so only the legacy path may extend them -- refusing
            # here is what keeps one grid from becoming two physics.
            ok = not want
            was = 'unrecorded (the shipped vP grids mix routes)'
        if not ok:
            raise SystemExit(
                f'\n{path.name} was built with Phase 3 {was}, but this run has it '
                f"{'ON' if want else 'OFF'}.\nAppending would put two different "
                f'physics in one grid.  Either pass '
                f"{'--no-phase3' if want else 'Phase 3 on (drop --no-phase3)'}, "
                f'or pick a --portal/--xi-ini whose tag differs so this run gets '
                f'its own file.')
        if not (np.array_equal(saved['coupling_axis'], COUPLING)
                and np.array_equal(saved['eps'], np.atleast_1d(eps))):
            print(f'  {path.name}: coupling axis or eps changed -- starting fresh')
            return coupling_relic, sigmav, status, done
        ix = _embed(saved['log_mX'], log_mX)
        jx = _embed(saved['log_rv'], log_rv)
        if ix is None or jx is None:
            print(f'  {path.name}: mesh does not embed in {N_MX}x{N_RV} '
                  f'-- starting fresh')
            return coupling_relic, sigmav, status, done

        if 'done_mask' in saved:                       # written by this version
            old_done = saved['done_mask']
        else:                                          # 16x12 pass: column counter
            old_done = np.zeros((saved['log_mX'].size, saved['log_rv'].size), bool)
            old_done[:int(saved['done'])] = True

        sel = np.ix_(ix, jx)
        if coupling_relic.ndim == 3:                            # leading eps axis
            coupling_relic[np.ix_(np.arange(shape[0]), ix, jx)] = saved['coupling_relic']
            sigmav[np.ix_(np.arange(shape[0]), ix, jx)] = saved['sigmav']
        else:
            coupling_relic[sel] = saved['coupling_relic']
            sigmav[sel] = saved['sigmav']
        if 'status' in saved:                          # absent in older grids,
            if coupling_relic.ndim == 3:                        # which keep UNRECORDED
                status[np.ix_(np.arange(shape[0]), ix, jx)] = saved['status']
            else:
                status[sel] = saved['status']
        done[sel] |= old_done

        n_old = int(saved['log_mX'].size*saved['log_rv'].size)
        bak = path.with_name(f'{path.stem}_{n_old}nodes_backup.npz')
        if backup and BACKUP_PREVIOUS and not bak.exists():
            shutil.copy2(path, bak)
            print(f'  {path.name}: previous run copied to {bak.name}')

    done |= off_grid                                    # never queue an off-grid node
    print(f'  {path.name}: reusing {int((done & ~off_grid).sum())} finished nodes')
    return coupling_relic, sigmav, status, done


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
        print(f'  {desc}: nothing to do')
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
    """coupling_relic and <sigma v> for every eps in EPS_SECLUDED, one output
    file with eps as the leading array axis."""
    path = OUT_DIR / f'relic_grid_{TAG}_eps_secluded.npz'
    eps = np.array(EPS_SECLUDED)
    _set_coupling('secluded')
    coupling_relic, sigmav, status, done = load_state(path, (len(eps), N_MX, N_RV), eps)
    n_requeued = _requeue(done, status, RETRY_STATUSES)
    if n_requeued:
        print(f'  re-queued {n_requeued} node(s) with status '
              f'{[STATUS_NAMES[c] for c in RETRY_STATUSES]}')

    def on_result(res):
        i, j, a, sv, st = res
        coupling_relic[:, i, j] = a
        sigmav[:, i, j] = sv
        status[:, i, j] = st
        done[i, j] = True

    _run_pool(node_secluded, _pending(done), ncores, 'secluded',
              on_result, lambda: _save(path, eps, coupling_relic, sigmav, done, status))
    print(f'secluded: {np.isfinite(sigmav).sum()}/{len(eps) * N_VALID} solved '
          f'-> {path.name}')


def run_joint(eps, ncores):
    """coupling_relic and <sigma v> at a single eps via the joint solver."""
    path = OUT_DIR / f'relic_grid_{TAG}_eps{eps:.0e}.npz'
    _set_coupling(eps)
    coupling_relic, sigmav, status, done = load_state(path, (N_MX, N_RV), eps)
    n_requeued = _requeue(done, status, RETRY_STATUSES)
    if n_requeued:
        print(f'  re-queued {n_requeued} node(s) with status '
              f'{[STATUS_NAMES[c] for c in RETRY_STATUSES]}')

    def on_result(res):
        i, j, a, sv, st = res
        coupling_relic[i, j] = a
        sigmav[i, j] = sv
        status[i, j] = st
        done[i, j] = True

    tasks = [t + (eps,) for t in _pending(done)]
    _run_pool(node_joint, tasks, ncores, f'joint {eps:.0e}',
              on_result, lambda: _save(path, eps, coupling_relic, sigmav, done, status))
    print(f'joint {eps:.0e}: {np.isfinite(sigmav).sum()}/{N_VALID} solved '
          f'-> {path.name}')


# ---- the thermal relic ----------------------------------------------
# eps >= eps_therm: the hidden sector stays locked to the SM, eps drops out, and
# solve_boltzmann_thermSM integrates the DM yield alone.  Each node is a short
# root-find rather than a sampled axis, so the coupling is solved on a coarse
# (log10 mX, log10 rv) mesh of its own and log10 of it is interpolated linearly onto
# the GCE fit grid -- the (mX, mY) points make_channel_spectra defines -- where the
# s-wave <sigma v> is evaluated.  Not checkpointed: the whole output takes minutes.

THERMAL_RTOL = 3e-3                      # |Omega/0.12 - 1| accepted from a secant step
THERMAL_RV_MAX = 0.995                   # the solve fails at mY = mX
THERMAL_N_RV_BETA = 4                    # top rows placed by beta; see thermal_rv_rows
THERMAL_N_MX, THERMAL_N_RV = 46, 34      # main block's coarse mesh
THERMAL_N_MX_EXT = 8                     # the extension's, before one padding step on
THERMAL_N_RV_EXT = 8                     # each side in mX and one below in rv: it spans
                                         # 0.1139 dex in mX and 0.2218 in rv, so 8 x 8
                                         # keeps the main mesh's density


def _thermal_model(mX, mY, coupling):
    """The model the thermal relic is solved with; see PORTALS[...]["thermal"]."""
    p = PORTALS[PORTAL]
    return p['thermal']['cls'](mX=mX, mY=mY, gX=p['gX'], gY=p['gY'],
                               include_antiparticlesX=p['dirac_X'],
                               include_antiparticlesY=False, **{p['coupling']: coupling})


def _och2_thermal(mX, mY, coupling):
    """Omega_c h^2 with the hidden sector thermalised with the SM."""
    sol = BoltzmannSolver(cosmo, _thermal_model(mX, mY, coupling)).solve_boltzmann_thermSM(
        verbose=False, include_hidden_in_gstar=False)
    return sol['YX_relic']*mX*s0_cosmo/rhoc


def thermal_coupling(mX, mY):
    """The coupling giving Omega_c h^2 = OCH2_TARGET in the thermalised regime, NaN
    if it is outside the portal's coupling_range.  A guess rescaled by
    Omega ~ coupling^-n, one log-log secant step, then brentq as the fallback."""
    t = PORTALS[PORTAL]['thermal']
    a_lo, a_hi = t['coupling_range']
    a0 = t['guess']
    o0 = _och2_thermal(mX, mY, a0)
    a1 = float(np.clip(a0*(o0/OCH2_TARGET)**(1.0/t['omega_power']), a_lo, a_hi))
    o1 = _och2_thermal(mX, mY, a1)
    if abs(o1/OCH2_TARGET - 1.0) < THERMAL_RTOL:
        return a1

    la = (np.log(a1) + (np.log(OCH2_TARGET) - np.log(o1))
          *(np.log(a1) - np.log(a0))/(np.log(o1) - np.log(o0)))
    a2 = float(np.clip(np.exp(la), a_lo, a_hi))
    if abs(_och2_thermal(mX, mY, a2)/OCH2_TARGET - 1.0) < THERMAL_RTOL:
        return a2

    f = lambda l: np.log(_och2_thermal(mX, mY, np.exp(l))/OCH2_TARGET)
    for lo, hi in ((a2/3.0, a2*3.0), (a_lo, a_hi)):
        try:
            return float(np.exp(brentq(f, np.log(max(lo, a_lo)),
                                       np.log(min(hi, a_hi)), rtol=1e-3)))
        except ValueError:
            continue
    return np.nan


def thermal_rv_rows(lrv_lo, n_rv, rv_max=THERMAL_RV_MAX, n_beta=THERMAL_N_RV_BETA):
    """log10 of the thermal mesh's mY/mX rows, dense where the s-wave closes.

    A ladder of `n_rv` rows uniform in log10 rv from `lrv_lo` up to `rv_max`, whose
    top row is replaced by `n_beta` rows at a uniform step in log10 beta, with
    beta = sqrt(1 - rv^2).  The s-wave shuts off as beta -> 0, so the uniform ladder
    leaves one interval below `rv_max` that covers several fine-grid rows, and the
    <sigma v> interpolated across it came out wrong by up to 65 % (vP) or 9 % (hP)
    against an exact solve.  The same number of rows placed by beta holds both under
    4 %.  Interpolating in log beta with the rows left alone is worse for hP, whose
    lam barely moves with rv and turns over near `rv_max`: it is the row placement
    that matters, not the interpolation variable.
    """
    ladder = np.linspace(lrv_lo, np.log10(rv_max), n_rv)
    beta = np.sqrt(1.0 - (10.0**ladder[-2:])**2)
    log_beta = np.linspace(np.log10(beta[0]), np.log10(beta[1]), n_beta + 1)[1:]
    return np.concatenate([ladder[:-1], 0.5*np.log10(1.0 - (10.0**log_beta)**2)])


def _thermal_node(task):
    """One coarse node, for the pool: (i, j, mX, mY) -> (i, j, coupling, message).
    A failure is carried back rather than raised, so one bad node cannot kill the
    sweep."""
    i, j, mX, mY = task
    try:
        return i, j, thermal_coupling(mX, mY), ''
    except Exception as exc:
        return i, j, np.nan, f'  FAILED mX={mX:.2f} mY={mY:.2f}: {exc}'


def _thermal_block(mX_fine, mY_fine, lmx_c, lrv_c, keep, ncores, desc):
    """Solve the coarse mesh lmx_c x lrv_c at the nodes `keep(mX, mY)` accepts,
    interpolate log10 of the coupling onto the fine points (mX_fine, mY_fine) [GeV],
    and evaluate the s-wave <sigma v> there.

    Returns coupling and sigmav [cm^3/s] on the fine points, then the coarse mesh's
    mX, mY [GeV] and coupling.  The nodes are independent, so the sweep is a plain
    map and the numbers do not depend on `ncores`.
    """
    LMX_C, LRV_C = np.meshgrid(lmx_c, lrv_c, indexing='ij')
    tasks = []
    for i in range(LMX_C.shape[0]):
        for j in range(LMX_C.shape[1]):
            mX = 10.0**LMX_C[i, j]
            mY = mX*10.0**LRV_C[i, j]
            if keep(mX, mY):
                tasks.append((i, j, mX, mY))

    coupling_c = np.full(LMX_C.shape, np.nan)
    with mp.get_context('fork').Pool(ncores) as pool:
        for i, j, coupling, message in tqdm(pool.imap_unordered(_thermal_node, tasks,
                                                                chunksize=1),
                                            total=len(tasks), desc=desc):
            coupling_c[i, j] = coupling
            if message:
                print(message, flush=True)

    ok_c = np.isfinite(coupling_c)
    coupling = 10.0**griddata(
        np.column_stack([LMX_C[ok_c], LRV_C[ok_c]]), np.log10(coupling_c[ok_c]),
        np.column_stack([np.log10(mX_fine), np.log10(mY_fine/mX_fine)]), method='linear')

    sigmav = np.full(mX_fine.size, np.nan)
    for k in np.flatnonzero(np.isfinite(coupling)):
        m = _thermal_model(float(mX_fine[k]), float(mY_fine[k]), float(coupling[k]))
        sigmav[k] = m.sigmav_XX_to_YY_swave()*GEVM2_TO_CM3S
    # the s-wave closes at mY = mX: nothing reaches the GCE there, whatever the coupling
    sigmav[mY_fine >= mX_fine] = 0.0
    return coupling, sigmav, 10.0**LMX_C, 10.0**(LMX_C + LRV_C), coupling_c


def run_thermal(ncores):
    """Thermal coupling and <sigma v> on the GCE fit grid: the main block, plus the
    mX = 100-130 GeV extension for a portal that has one, appended after it in the
    order make_channel_spectra gives, so one file is the portal's whole fit grid."""
    # Imported here, not at the top: make_channel_spectra loads Pythia, which only
    # this output needs, and only for the grid definition.
    from make_channel_spectra import EXT_GRID, MAIN_GRID

    path = OUT_DIR / f'relic_grid_{TAG}_thermal.npz'

    # main block: the mesh spans the fit grid, padded one rv row below its mY floor
    mX, mY = MAIN_GRID
    lmx, lrv = np.log10(mX), np.log10(mY/mX)
    lmx_c = np.linspace(lmx.min(), lmx.max(), THERMAL_N_MX)
    lrv_c = thermal_rv_rows(lrv.min(), THERMAL_N_RV)
    rv_step = 10.0**(lrv_c[1] - lrv_c[0])
    keep = lambda m_x, m_y: MY_LIM[0]/rv_step <= m_y <= THERMAL_RV_MAX*m_x
    coupling, sigmav, mX_c, mY_c, coupling_c = _thermal_block(
        mX, mY, lmx_c, lrv_c, keep, ncores, f'thermal {PORTAL}')
    out = {'mX': mX, 'mY': mY, 'coupling_relic': coupling, 'sigmav_relic': sigmav,
           'mX_coarse': mX_c, 'mY_coarse': mY_c, 'coupling_coarse': coupling_c}

    if PORTALS[PORTAL]['extension']:
        # Padded by one step on either side in mX and one below in rv, so that every
        # target point falls inside the interpolation hull rather than on its edge
        # (the extension's masses are rounded to 2 decimals).  Its top row is rv = 1,
        # where <sigma v> = 0 whatever the coupling.
        mX_e, mY_e = EXT_GRID
        lmx, lrv = np.log10(mX_e), np.log10(mY_e/mX_e)
        dmx = (lmx.max() - lmx.min())/(THERMAL_N_MX_EXT - 1)
        drv = (lrv.max() - lrv.min())/(THERMAL_N_RV_EXT - 1)
        lmx_c = np.linspace(lmx.min() - dmx, lmx.max() + dmx, THERMAL_N_MX_EXT + 2)
        lrv_c = thermal_rv_rows(lrv.min() - drv, THERMAL_N_RV_EXT + 1)
        keep = lambda m_x, m_y: m_y <= THERMAL_RV_MAX*m_x
        coupling, sigmav, mX_c, mY_c, coupling_c = _thermal_block(
            mX_e, mY_e, lmx_c, lrv_c, keep, ncores, f'thermal {PORTAL} ext')
        for key, block in (('mX', mX_e), ('mY', mY_e), ('coupling_relic', coupling),
                           ('sigmav_relic', sigmav)):
            out[key] = np.concatenate([out[key], block])
        out.update(mX_coarse_ext=mX_c, mY_coarse_ext=mY_c, coupling_coarse_ext=coupling_c)

    np.savez(path, **out)
    print(f"thermal: {np.isfinite(out['sigmav_relic']).sum()}/{out['mX'].size} points "
          f'with a <sigma v> -> {path.name}')


# ---- sharding: one output split across concurrent jobs ------------------

def _shard_tag(tag, shard, nshards):
    """Output stem for one shard.  Shards write separate files so they can run
    concurrently; --merge folds them back into `tag`."""
    return f'{tag}_s{shard:02d}of{nshards:02d}'


def _apply_shard(shard, nshards):
    """Restrict this process to every nshards-th rv column, starting at `shard`.

    A stride rather than a block: the stiff nodes sit at high rv, so a block split
    would hand one shard most of the cost.  _pending leaves the other columns
    not-done, so each shard's file is the full mesh with only its own columns
    filled, which is what _merge relies on.
    """
    global ONLY_RV_COLUMNS, TAG
    ONLY_RV_COLUMNS = list(range(shard, N_RV, nshards))
    TAG = _shard_tag(TAG, shard, nshards)


def _merge(nshards, only):
    """Fold the shard files of output `only` back into its canonical grid.

    Each shard contributes only the nodes its own done_mask marks, so the shards never
    disagree about a node.  An existing canonical file is kept as the base, so nodes
    solved before sharding are not lost.
    """
    secluded = only == 'secluded'
    eps = np.array(EPS_SECLUDED) if secluded else float(only)
    stem = 'eps_secluded' if secluded else f'eps{float(only):.0e}'
    path = OUT_DIR / f'relic_grid_{TAG}_{stem}.npz'

    _set_coupling('secluded' if secluded else float(only))
    shape = (len(EPS_SECLUDED), N_MX, N_RV) if secluded else (N_MX, N_RV)
    coupling_relic, sigmav, status, done = load_state(path, shape, eps, backup=False)

    added = 0
    for k in range(nshards):
        sp = OUT_DIR / f'relic_grid_{_shard_tag(TAG, k, nshards)}_{stem}.npz'
        if not sp.exists():
            print(f'  {stem}: shard {k} missing ({sp.name}) -- skipped')
            continue
        with np.load(sp, allow_pickle=True) as d:
            take = d['done_mask'] & ~done
            coupling_relic[..., take] = d['coupling_relic'][..., take]
            sigmav[..., take] = d['sigmav'][..., take]
            status[..., take] = d['status'][..., take]
        done |= take
        added += int(take.sum())

    _save(path, eps, coupling_relic, sigmav, done, status)
    print(f'  {path.name}: merged {added} nodes, '
          f'{int((done & ~off_grid).sum())}/{N_VALID} done')


# ---- driver -----------------------------------------------------------

def dry_run():
    """Grid, what each output would reuse, and the wall-clock estimate."""
    if RUN_THERMAL:
        path = OUT_DIR / f'relic_grid_{TAG}_thermal.npz'
        print(f"  {path.name}: {'present' if path.exists() else 'missing'}; "
              f'not checkpointed, a run re-solves it in full (minutes)')
    todo = 0
    outputs = [(OUT_DIR / f'relic_grid_{TAG}_eps{e:.0e}.npz', e, (N_MX, N_RV))
               for e in EPS_JOINT]
    if RUN_SECLUDED:
        outputs.insert(0, (OUT_DIR / f'relic_grid_{TAG}_eps_secluded.npz',
                           np.array(EPS_SECLUDED), (len(EPS_SECLUDED), N_MX, N_RV)))
    for path, eps, shape in outputs:
        _set_coupling(eps if np.isscalar(eps) else 'secluded')
        coupling_relic, sigmav, status, done = load_state(path, shape, eps, backup=False)
        _requeue(done, status, RETRY_STATUSES)   # or the estimate understates
        left = len(_pending(done))               # what the real run would solve
        todo += left
        print(f'  {path.name}: {left} of {N_VALID} nodes still to solve'
              + ('' if ONLY_RV_COLUMNS is None
                 else f'  (queue restricted to rv columns {tuple(ONLY_RV_COLUMNS)}'
                      f" = rv {', '.join(f'{10 ** log_rv[k]:.4f}' for k in ONLY_RV_COLUMNS)})"))
    per_out = SEC_PER_NODE/3.0
    print(f'\n{todo} node-solves left ~ {todo * per_out / 3600:.1f} h serial, '
          f'{todo * per_out / 3600 / NCORES:.1f} h on {NCORES} cores '
          f'(at {SEC_PER_NODE:.0f} s per node across all three outputs)')


def main():
    global PHASE3_KW, RUN_THERMAL, RUN_SECLUDED, EPS_JOINT, SOLVE_TIMEOUT, RETRY_STATUSES
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--portal', default='vP', choices=list(PORTALS),
                    help='vP, bP, blP or hP (default vP)')
    ap.add_argument('--xi-ini', type=float, default=1.0,
                    help='T_D/T_SM at Tinf (default 1.0).  The thermal output has the '
                         'hidden sector at T_SM by construction and runs only at 1')
    ap.add_argument('--only', default=None, metavar='OUTPUT',
                    help="run one output: 'thermal', 'secluded', or one of the portal's "
                         'joint eps (1e-9, 1e-10 for vP/bP/blP; 1e-7, 1e-8 for hP).  '
                         'Each writes its own file, so they can run as concurrent jobs')
    ap.add_argument('--no-phase3', action='store_true',
                    help='legacy path: clear PHASE3_KW')
    ap.add_argument('-n', '--ncores', type=int, default=NCORES,
                    help=f'worker processes (default {NCORES})')
    # Run control, not configuration: like -n these say how to run, not what to
    # compute.  --retry-failed re-solves nodes a previous run could not finish;
    # --solve-timeout is what makes that worth doing.
    ap.add_argument('--nshards', type=int, default=1,
                    help='split the rv columns of one output across this many jobs')
    ap.add_argument('--shard', type=int, default=None,
                    help='which shard this process runs (0 .. nshards-1)')
    ap.add_argument('--merge', action='store_true',
                    help='fold the shard files back into the canonical grid and exit')
    ap.add_argument('--retry-failed', action='store_true',
                    help='also re-solve stored nodes whose status is solve_failed')
    ap.add_argument('--solve-timeout', type=int, default=None, metavar='S',
                    help=f'seconds per solve (default {SOLVE_TIMEOUT})')
    ap.add_argument('--dry-run', action='store_true',
                    help='report the mesh, what is reusable and the ETA, then exit')
    args = ap.parse_args()

    # Must happen before any Pool is forked: the workers inherit these globals.
    _set_portal(args.portal, args.xi_ini)

    only_eps = None
    if args.only not in (None, 'thermal', 'secluded'):
        try:
            only_eps = float(args.only)
        except ValueError:
            only_eps = None
        if only_eps not in EPS_JOINT:
            ap.error(f'--only for --portal {PORTAL} is thermal, secluded or one of '
                     f"{', '.join(f'{e:g}' for e in EPS_JOINT)}")
    if args.only == 'thermal' and XI_INI != 1.0:
        ap.error('the thermal output exists only at --xi-ini 1')
    if args.nshards < 1:
        ap.error('--nshards must be >= 1')
    if args.shard is not None and not 0 <= args.shard < args.nshards:
        ap.error(f'--shard must be in 0..{args.nshards - 1}')
    if (args.shard is not None or args.merge) and args.only in (None, 'thermal'):
        ap.error('--shard and --merge need --only secluded or a joint eps: shards of '
                 'different outputs are separate files')

    if args.no_phase3:
        PHASE3_KW = {}
    if args.retry_failed:
        RETRY_STATUSES = (STATUS_SOLVE_FAILED,)
    if args.solve_timeout is not None:
        SOLVE_TIMEOUT = args.solve_timeout
    if args.only is not None:
        RUN_THERMAL = args.only == 'thermal'
        RUN_SECLUDED = args.only == 'secluded'
        EPS_JOINT = () if only_eps is None else (only_eps,)
    RUN_THERMAL = RUN_THERMAL and XI_INI == 1.0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.merge:
        _merge(args.nshards, args.only)
        return
    if args.shard is not None:
        _apply_shard(args.shard, args.nshards)

    p = PORTALS[PORTAL]
    print(f"portal: {PORTAL} = {p['cls'].__name__}  gX={p['gX']} gY={p['gY']} "
          f"{'Dirac' if p['dirac_X'] else 'Majorana'} X  xi_ini={XI_INI:g}  |  "
          f"Phase 3 {'ON' if PHASE3_KW else 'OFF'}"
          f"{' (required)' if PHASE3_KW and PHASE3_REQUIRE else ''}  |  tag '{TAG}'")
    axes = ', '.join(
        ('default' if k is None else k if isinstance(k, str) else f'eps={k:g}')
        + f' [{v[0]:g}, {v[1]:g}]' for k, v in COUPLING_LIMS.items())
    print(f'mesh: {N_MX} x {N_RV} = {N_VALID} valid (mX, mY) nodes, mX in '
          f'[{10.0 ** log_mX[0]:.2f}, {10.0 ** log_mX[-1]:.2f}] GeV, rv <= {RV_MAX}  |  '
          f"{N_COUPLING} {p['coupling']} per node, axis {axes}  |  out -> {OUT_DIR}")
    outputs = ((['thermal'] if RUN_THERMAL else [])
               + ([f"secluded {', '.join(f'{e:g}' for e in EPS_SECLUDED)}"]
                  if RUN_SECLUDED else [])
               + [f'joint {e:g}' for e in EPS_JOINT])
    print(f"outputs: {' + '.join(outputs) if outputs else 'none'}")

    if args.dry_run:
        dry_run()
        return
    print(f'workers: {args.ncores}')
    if RUN_THERMAL:
        run_thermal(args.ncores)
    if RUN_SECLUDED:
        run_secluded(args.ncores)
    for eps in EPS_JOINT:
        run_joint(eps, args.ncores)


if __name__ == '__main__':
    main()

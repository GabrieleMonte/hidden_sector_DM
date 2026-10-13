"""Why the Higgs-portal secluded grid has no lam at most nodes: Omega(lam) itself.

The campaign leaves 81% of eps=1e-12 nodes without a crossing, yet it records
zero solve failures and only one crossing at the sampled axis floor.  So the
nodes are not failing and are not running off the end of a too-narrow axis.
This evaluates Omega_c h^2 on a much wider lam axis at a sample of those nodes
and asks which of two things is true:

  SATURATED  Omega flattens as lam falls, peaking below 0.12.  No lam solves
             the node and no widening of the axis would; the point is excluded.
  NARROW     Omega keeps rising as lam falls and crosses 0.12 below 2e-5.  The
             axis is simply too short and should be extended.

    cosmo_env
    sbatch GCE/scripts/omega_vs_lam_hP.sbatch

Writes GCE/output/diag/omega_vs_lam_hP.npz: omega (n_node, n_lam, n_eps),
with the probed lam axis and the (mX, mY, rv) of each node.
"""

import sys
from pathlib import Path

import numpy as np
from tqdm import tqdm

import multiprocessing as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                    # noqa: E402
import make_relic_grid_hP as hp                                   # noqa: E402

# ---- configuration ---------------------------------------------------

# Probed lam axis: 2.3 decades below the campaign floor of 2e-5, up to its
# ceiling of 1.  Omega ~ lam^-4, so this is ~28 decades of Omega -- far more
# than enough to see a crossing if one exists.
LAM_PROBE = np.logspace(-7.0, 0.0, 15)

N_FAILED = 16            # sampled nodes with no crossing at DIAG_EPS
N_CONTROL = 4            # sampled nodes that did get one, for comparison
DIAG_EPS = 1e-12         # the worst slice: 18.7% of nodes finite

PROBE_TIMEOUT = 600      # s; the campaign's 7200 only matters for stiff nodes
                         # and these recorded no failures at all.  Keeps the
                         # whole probe inside one debug-QOS allocation.

OMEGA_TARGET = 0.12      # Omega_c h^2 the grid solves for
FLAT_TOL = 0.05          # relative Omega change over the bottom decade that
                         # still counts as a saturated shelf

# The campaign checkpoints the canonical grid with np.savez, which rewrites it
# in place rather than renaming into position, so reading it mid-run can catch a
# half-written file.  The sbatch copies it here first and this prefers that copy.
DIAG_DIR = base.OUT_DIR.parent / "diag"
SNAPSHOT = DIAG_DIR / "relic_grid_hP_eps_secluded_snapshot.npz"
GRID = SNAPSHOT if SNAPSHOT.exists() else base.OUT_DIR / "relic_grid_hP_eps_secluded.npz"
OUT = DIAG_DIR / "omega_vs_lam_hP.npz"


# ---- node selection --------------------------------------------------

def pick_nodes():
    """(failed, control) index arrays into the (mX, rv) mesh, sampled evenly.

    Even sampling over the flat index spreads the picks across mX and rv
    instead of clustering them in whichever corner failed first.
    """
    with np.load(GRID, allow_pickle=True) as d:
        coupling, done, eps = d["coupling_relic"], d["done_mask"], d["eps"]

    e = int(np.argmin(np.abs(eps - DIAG_EPS)))
    solved = done & ~base.off_grid
    finite = np.isfinite(coupling[e])

    failed = np.flatnonzero((solved & ~finite).ravel())
    control = np.flatnonzero((solved & finite).ravel())
    if failed.size == 0:
        raise SystemExit(f"{GRID.name}: no solved node lacks a crossing at eps={DIAG_EPS:g}")

    take = lambda a, n: a[np.linspace(0, a.size - 1, min(n, a.size)).astype(int)]
    return take(failed, N_FAILED), take(control, N_CONTROL)


# ---- probe -----------------------------------------------------------

def probe(task):
    """(k, l, mX, mY) -> Omega_c h^2 for every eps in EPS_SECLUDED at one lam."""
    k, l, mX, mY = task
    return k, l, base.omega_h2_secluded(mX, mY, LAM_PROBE[l])


def verdict(omega_node, lam):
    """Classify one node's Omega(lam) curve as SATURATED, NARROW or UNRESOLVED.

    `omega_node` is (n_lam,) at a single eps, ordered by increasing lam.
    """
    good = np.isfinite(omega_node)
    if good.sum() < 3:
        return "UNRESOLVED", np.nan

    om, lm = omega_node[good], lam[good]
    peak = om.max()
    if peak >= OMEGA_TARGET:
        return "NARROW", peak

    bottom = lm <= 10.0 * lm.min()                  # the lowest decade probed
    if bottom.sum() >= 2:
        spread = np.ptp(om[bottom]) / max(om[bottom].max(), 1e-300)
        if spread <= FLAT_TOL:
            return "SATURATED", peak

    return "UNRESOLVED", peak


# ---- main ------------------------------------------------------------

def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else base.NCORES

    hp._configure(1.0)
    base.SOLVE_TIMEOUT = PROBE_TIMEOUT

    failed, control = pick_nodes()
    flat = np.concatenate([failed, control])
    kind = np.array(["failed"] * failed.size + ["control"] * control.size)

    i, j = np.unravel_index(flat, base.mX_grid.shape)
    mX, mY = base.mX_grid[i, j], base.mY_grid[i, j]
    rv = mY / mX

    n_eps = len(base.EPS_SECLUDED)
    omega = np.full((flat.size, LAM_PROBE.size, n_eps), np.nan)

    tasks = [(k, l, mX[k], mY[k]) for k in range(flat.size) for l in range(LAM_PROBE.size)]
    print(f"{flat.size} nodes ({failed.size} failed, {control.size} control) "
          f"x {LAM_PROBE.size} lam = {len(tasks)} solves on {ncores} workers")
    print(f"lam probed: {LAM_PROBE[0]:.3g} .. {LAM_PROBE[-1]:.3g}  "
          f"(campaign floor 2e-05)\n")

    with mp.Pool(ncores) as pool:
        for k, l, om in tqdm(pool.imap_unordered(probe, tasks, chunksize=1),
                             total=len(tasks), desc="omega", mininterval=10.0):
            omega[k, l] = om

    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT, omega=omega, lam=LAM_PROBE, mX=mX, mY=mY, rv=rv,
                        kind=kind, eps=np.array(base.EPS_SECLUDED),
                        diag_eps=DIAG_EPS, node_index=flat)

    # ---- report ----
    e = int(np.argmin(np.abs(np.array(base.EPS_SECLUDED) - DIAG_EPS)))
    print(f"\neps = {DIAG_EPS:g}, Omega target {OMEGA_TARGET}\n")
    print(f"  {'node':>5s} {'kind':>8s} {'mX':>8s} {'rv':>6s} "
          f"{'peak Omega':>11s} {'verdict':>11s}")

    tally = {}
    for k in range(flat.size):
        tag, peak = verdict(omega[k, :, e], LAM_PROBE)
        tally[tag] = tally.get(tag, 0) + 1
        print(f"  {flat[k]:5d} {kind[k]:>8s} {mX[k]:8.2f} {rv[k]:6.3f} "
              f"{peak:11.4g} {tag:>11s}")

    print("\n " + "  ".join(f"{t}: {n}" for t, n in sorted(tally.items())))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()

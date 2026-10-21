"""Where do Higgs-portal relic crossings actually sit, once the solver is fixed?

The phase-entry fix (1d99f58edb) changed Omega(lam) completely below
lam ~ 8e-3: what used to be a flat, lam-independent shelf orders of magnitude
BELOW the live branch is now a monotonically falling curve, and nodes the
grids recorded as NO_SOLUTION do have crossings.  Every axis choice in
make_relic_grid_hP.py was tuned against the artifact, so LAM_LIMS has to be
re-derived before the grids are rebuilt.

This samples a spread of (mX, rv) across the grid domain and reports, per eps,
the lam at which Omega_c h^2 crosses 0.12.  The answer sets the new axis: a
floor comfortably below the smallest crossing, a ceiling above the largest,
and no decades spent where nothing happens.

It covers the joint outputs too (eps = 1e-9, 1e-10).  The same fix was applied
to solve_boltzmann_joint but has not been shown to move any joint result, so
these columns are a check on whether those grids also need rebuilding.

    cosmo_env
    sbatch GCE/scripts/lam_crossing_survey_hP.sbatch

Writes GCE/output/diag/lam_crossing_survey_hP.npz.
"""

import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                    # noqa: E402
import make_relic_grid_hP as hp                                   # noqa: E402

# ---- configuration ---------------------------------------------------

# Probe axis.  Re-aimed 2026-10-20 for eps = 1e-8 and 1e-7, which sit ABOVE
# every eps solved so far.  The crossing climbs steadily with eps --
# 0.0043-0.019 at 1e-12, 0.0089-0.047 at 1e-11, 0.016-0.136 at 1e-10,
# 0.029-0.319 at 1e-9, so roughly lam ~ eps^0.4 at the top end -- which
# extrapolates to ~0.8 at 1e-8 and ~2 at 1e-7.  Both are at or over the
# production axis ceiling of 1.0, so the probe has to look past it; 3.98 is
# just inside perturbativity (sqrt(4 pi) = 3.54 is the usual marker, so
# anything found above ~3.5 is a result about the model, not the axis).
# The old axis was logspace(-4, 0, 13), aimed at 1e-9 and 1e-10.
LAM_PROBE = np.logspace(-2.0, 0.6, 14)
PROBE_TIMEOUT = 300      # s per solve; a healthy solve here takes 4-9 s
EPS_JOINT = (1e-8, 1e-7)

# The secluded column is skipped: base.EPS_SECLUDED is (1e-9, 1e-10) for this
# portal now and both are solved, so sweeping it again would only cost time.
INCLUDE_SECLUDED = False

N_MX, N_RV = 6, 5        # survey mesh over the grid's own domain
MX_LIMS = (16.0, 95.0)
RV_LIMS = (0.10, 0.93)

DIAG_DIR = base.OUT_DIR.parent / "diag"
OUT = DIAG_DIR / "lam_crossing_survey_hP.npz"


def survey_nodes():
    """(mX, rv) pairs spanning the grid domain, log-spaced in mX."""
    mX = np.logspace(*np.log10(MX_LIMS), N_MX)
    rv = np.linspace(*RV_LIMS, N_RV)
    return [(float(a), float(b)) for a in mX for b in rv]


def one_task(arg):
    """(k, kind) -> Omega over LAM_PROBE for one node and one output."""
    k, mX, mY, kind = arg
    t0 = time.perf_counter()
    ncol = len(base.EPS_SECLUDED) if kind == "secluded" else 1
    if kind == "secluded":
        omega = base._sweep(lambda a: base.omega_h2_secluded(mX, mY, a), ncol)
    else:
        eps = float(kind)
        omega = base._sweep(lambda a: base.omega_h2_joint(mX, mY, a, eps), 1)
    return k, kind, omega, time.perf_counter() - t0


def main():
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else base.NCORES

    hp._configure(1.0)
    base.SOLVE_TIMEOUT = PROBE_TIMEOUT
    base.COUPLING = LAM_PROBE
    base.N_COUPLING = LAM_PROBE.size
    base.COUPLING_SWEEP_STOP = True          # stop one point past the crossing

    nodes = survey_nodes()
    kinds = (["secluded"] if INCLUDE_SECLUDED else []) \
        + [f"{e:g}" for e in EPS_JOINT]
    tasks = [(k, mX, rv * mX, kind)
             for k, (mX, rv) in enumerate(nodes) for kind in kinds]

    print(f"{len(nodes)} nodes x {len(kinds)} outputs = {len(tasks)} sweeps")
    print(f"lam probe: {LAM_PROBE.size} points "
          f"[{LAM_PROBE[0]:.3g}, {LAM_PROBE[-1]:.3g}], {ncores} workers",
          flush=True)

    n_sec = len(base.EPS_SECLUDED)
    omega_sec = np.full((len(nodes), LAM_PROBE.size, n_sec), np.nan)
    omega_joint = np.full((len(nodes), LAM_PROBE.size, len(EPS_JOINT)), np.nan)

    with mp.Pool(min(ncores, len(tasks))) as pool:
        for k, kind, omega, dt in pool.imap_unordered(one_task, tasks):
            if kind == "secluded":
                omega_sec[k] = omega
            else:
                omega_joint[k, :, EPS_JOINT.index(float(kind))] = omega[:, 0]

    # ---- crossings -------------------------------------------------
    labels, columns = [], []
    if INCLUDE_SECLUDED:
        labels += [f"secluded eps={e:g}" for e in base.EPS_SECLUDED]
        columns += [omega_sec[:, :, e] for e in range(n_sec)]
    labels += [f"joint    eps={e:g}" for e in EPS_JOINT]
    columns += [omega_joint[:, :, e] for e in range(len(EPS_JOINT))]

    print(f"\n{'output':>22} {'found':>7} {'min lam':>11} {'max lam':>11}")
    allx = []
    for lab, col in zip(labels, columns):
        xs = []
        for row in col:
            a, s = base._crossing(row)
            if s == base.STATUS_OK:
                xs.append(a)
        allx += xs
        lo = f"{min(xs):11.4g}" if xs else f"{'-':>11}"
        hi = f"{max(xs):11.4g}" if xs else f"{'-':>11}"
        print(f"{lab:>22} {len(xs):3d}/{len(col):<3d} {lo} {hi}")

    if allx:
        print(f"\nALL crossings: {len(allx)} in "
              f"[{min(allx):.4g}, {max(allx):.4g}]")
        print(f"suggested LAM_LIMS floor (min/3): {min(allx) / 3:.2g}")

    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, omega_secluded=omega_sec, omega_joint=omega_joint,
             lam=LAM_PROBE, eps_secluded=np.array(base.EPS_SECLUDED),
             eps_joint=np.array(EPS_JOINT),
             mX=np.array([n[0] for n in nodes]),
             rv=np.array([n[1] for n in nodes]))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()

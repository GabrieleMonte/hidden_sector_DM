"""lam that reproduces Omega_c h^2 = 0.12 for the Higgs-portal cascade
X Xbar -> phi phi -> SM, and the matching s-wave <sigma v>, over the GCE
(mX, mY) grid.

    cosmo_env
    python GCE/scripts/make_relic_grid_hP.py --dry-run   # grid, reuse, ETA
    python GCE/scripts/make_relic_grid_hP.py --only 1e-9 -n 128

Three files land in GCE/output/relic/, in the format the vector-portal grids
use: relic_grid_hP_eps1e-09.npz, relic_grid_hP_eps1e-10.npz and
relic_grid_hP_eps_secluded.npz (eps = 1e-11, 5e-12, 1e-12 on the leading axis).

The mesh, solver settings, checkpointing and crossing extraction are
make_relic_grid.py's and are used unchanged.  This script supplies only what
the Higgs portal does differently: the model, and the sampled coupling axis.

MODEL
-----
X is a Majorana fermion and phi a real scalar, so gX = 2, gY = 1 and X is its
own antiparticle -- not the vector portals' gY = 3 and Dirac X.  The swept
coupling is `lam` (which sets lam_s = lam_p), not alphaX, and epsX here is the
Higgs mixing angle sin(theta) rather than a kinetic mixing.  Gamma_phi comes
from HDECAY via hidden_sector_DM/phi_decay; every solved node has mY >= 5 GeV,
so only HDECAY's regime is reached, never the Winkler or pQCD ones.

AXIS
----
sigmav ~ lam_s^2 lam_p^2 = lam^4, so Omega ~ lam^-4 where the vector portals
give Omega ~ alpha^-2.  Their (2e-5, 1.0) axis is 4.7 decades of coupling,
which is 9.4 decades of Omega for them but 18.8 decades here: its lower half
falls where annihilation no longer sets the relic, Omega saturates, and the
cubic spline through that flat shelf oscillates (it returned lam = 0.0 at the
rv = 0.950 corner).  LAM_LIMS keeps every sample in the live region instead.
Measured crossings sit near lam ~ 0.05-0.3.

The secluded outputs leave most nodes without a lam -- 74 / 59 / 19 % finite at
eps = 1e-11 / 5e-12 / 1e-12 -- and that is the physics, not a short axis.  Do
not widen LAM_LIMS to chase it.  omega_vs_lam_hP.py probed 1e-7 to 1, which is
2.3 decades below the floor and ~28 decades of Omega, at 16 nodes carrying no
crossing at eps = 1e-12: all 16 peak below 0.12 (highest 0.094), and below
lam ~ 1e-5 Omega is flat to the last bit, because X stops annihilating at all
and the yield can no longer rise.  phi at sin(theta) = 1e-12 then dilutes that
ceiling by many orders.  Across 60 (node, eps) cases, "grid found no lam"
and "Omega never reaches 0.12" agree 58 times.  Crossings survive only near
rv -> 1 (controls 0.727-0.950, no-crossing nodes 0.146-0.695).
"""

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                    # noqa: E402

from hidden_sector_DM.HiddenSectorDM import HiggsPortal           # noqa: E402

# ---- configuration ---------------------------------------------------

# Sampled lam axis per output, as (lo, hi); lam is dimensionless.  load_state
# records the axis and compares it exactly, so changing an entry resets that
# output's file.  See AXIS above for why this is not the vector portals'.
LAM_LIMS = {
    None:       (2e-5, 1.0),
    1e-9:       (2e-5, 1.0),
    1e-10:      (2e-5, 1.0),
    "secluded": (2e-5, 1.0),
}

GX_MAJORANA = 2          # Majorana fermion DM: 2 spin states, no antiparticle
GY_REAL_SCALAR = 1       # real scalar mediator: one degree of freedom

TAG_STEM = "hP"          # output-name stem; xi_ini != 1 is appended


def _make_model(mX, mY, lam, antiY=False):
    """The Higgs-portal model at one grid node.

    mX, mY in GeV.  `lam` is the swept dimensionless coupling, setting
    lam_s = lam_p.  The signature matches make_relic_grid._make_model so the
    shared solve functions call it unchanged.
    """
    return HiggsPortal(mX=mX, mY=mY, gX=GX_MAJORANA, gY=GY_REAL_SCALAR,
                       lam=lam, xi_ini=base.XI_INI,
                       include_antiparticlesX=False,
                       include_antiparticlesY=antiY)


def _configure(xi_ini):
    """Point the shared driver at the Higgs portal.

    Must run before any worker is forked: the pool inherits these globals,
    exactly as make_relic_grid._set_portal relies on.
    """
    base.XI_INI = float(xi_ini)
    base._make_model = _make_model
    base.COUPLING_LIMS = LAM_LIMS
    base.PORTAL = TAG_STEM
    base.PORTAL_KW = {}
    base.TAG = TAG_STEM if base.XI_INI == 1.0 else f"{TAG_STEM}_xi_ini_{base.XI_INI:g}"


def _shard_tag(tag, shard, nshards):
    """Output stem for one shard.  Shards write separate files so they can run
    concurrently; `--merge` folds them back into `tag`."""
    return f"{tag}_s{shard:02d}of{nshards:02d}"


def _apply_shard(shard, nshards):
    """Restrict this process to every nshards-th rv column, starting at `shard`.

    A stride rather than a block: the stiff nodes sit at high rv, so a block
    split would hand one shard most of the cost.  _pending leaves the other
    columns not-done, so each shard's file is the full grid with only its own
    columns filled, which is what _merge relies on.
    """
    base.ONLY_RV_COLUMNS = list(range(shard, base.N_RV, nshards))
    base.TAG = _shard_tag(base.TAG, shard, nshards)


def _merge(nshards, outputs):
    """Fold the shard files for `outputs` back into the canonical grids.

    Each shard contributes only the nodes its own done_mask marks, so the
    shards never disagree about a node.  An existing canonical file is kept as
    the base, so nodes solved before sharding are not lost.
    """
    for only in outputs:
        base.RUN_SECLUDED = only == "secluded"
        base.EPS_JOINT = () if base.RUN_SECLUDED else (float(only),)
        eps = np.array(base.EPS_SECLUDED) if base.RUN_SECLUDED else float(only)
        stem = "eps_secluded" if base.RUN_SECLUDED else f"eps{float(only):.0e}"
        path = base.OUT_DIR / f"relic_grid_{TAG_STEM}_{stem}.npz"

        base._set_coupling("secluded" if base.RUN_SECLUDED else float(only))
        shape = ((len(base.EPS_SECLUDED), base.N_MX, base.N_RV) if base.RUN_SECLUDED
                 else (base.N_MX, base.N_RV))
        coupling_relic, sigmav, status, done = base.load_state(path, shape, eps,
                                                               backup=False)

        added = 0
        for k in range(nshards):
            sp = base.OUT_DIR / f"relic_grid_{_shard_tag(TAG_STEM, k, nshards)}_{stem}.npz"
            if not sp.exists():
                print(f"  {stem}: shard {k} missing ({sp.name}) -- skipped")
                continue
            with np.load(sp, allow_pickle=True) as d:
                take = d["done_mask"] & ~done
                coupling_relic[..., take] = d["coupling_relic"][..., take]
                sigmav[..., take] = d["sigmav"][..., take]
                if "status" in d:                      # shards predating status
                    status[..., take] = d["status"][..., take]
            done |= take
            added += int(take.sum())

        base._save(path, eps, coupling_relic, sigmav, done, status)
        print(f"  {path.name}: merged {added} nodes, "
              f"{int((done & ~base.off_grid).sum())}/{base.N_VALID} done")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--xi-ini", type=float, default=1.0,
                    help="T_D/T_SM at Tinf (default 1.0)")
    ap.add_argument("--only", default=None,
                    choices=["secluded", "1e-9", "1e-10"],
                    help="run just one output.  The three write different "
                         "files, so they can be separate concurrent slices")
    ap.add_argument("--no-phase3", action="store_true",
                    help="legacy path: clear PHASE3_KW")
    ap.add_argument("--nshards", type=int, default=1,
                    help="split the rv columns across this many concurrent jobs")
    ap.add_argument("--shard", type=int, default=None,
                    help="which shard this process runs (0 .. nshards-1)")
    ap.add_argument("--merge", action="store_true",
                    help="fold the shard files back into the canonical grids "
                         "and exit; needs --nshards and --only")
    ap.add_argument("-n", "--ncores", type=int, default=base.NCORES,
                    help=f"worker processes (default {base.NCORES})")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the grid, what is reusable and the ETA, then exit")
    args = ap.parse_args()

    if args.nshards < 1:
        ap.error("--nshards must be >= 1")
    if args.shard is not None and not 0 <= args.shard < args.nshards:
        ap.error(f"--shard must be in 0..{args.nshards - 1}")
    if args.shard is not None and args.only is None:
        ap.error("--shard needs --only: shards of different outputs are "
                 "separate files and must not share one run")

    _configure(args.xi_ini)

    if args.no_phase3:
        base.PHASE3_KW = {}
    if args.only is not None:
        base.RUN_SECLUDED = args.only == "secluded"
        base.EPS_JOINT = () if base.RUN_SECLUDED else (float(args.only),)

    base.OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.merge:
        if args.only is None:
            ap.error("--merge needs --only")
        _merge(args.nshards, [args.only])
        return

    if args.shard is not None:
        _apply_shard(args.shard, args.nshards)

    axes = ", ".join(
        ("default" if k is None else k if isinstance(k, str) else f"eps={k:g}")
        + f" [{v[0]:g}, {v[1]:g}]" for k, v in LAM_LIMS.items())
    print(f"portal: HiggsPortal  xi_ini={base.XI_INI:g}  "
          f"gX={GX_MAJORANA} gY={GY_REAL_SCALAR} Majorana X  |  "
          f"Phase 3 {'ON' if base.PHASE3_KW else 'OFF'}  |  tag '{base.TAG}'")
    print(f"grid: {base.N_MX} x {base.N_RV} = {base.N_VALID} valid (mX, mY) nodes, "
          f"mX in [{base.MX_LIM[0]:g}, {base.MX_LIM[1]:g}] GeV, "
          f"rv <= {base.RV_MAX}  |  {base.N_COUPLING} lam, axis {axes}  |  "
          f"out -> {base.OUT_DIR}")

    if args.dry_run:
        base.dry_run()
        return

    print(f"workers: {args.ncores}")
    if base.RUN_SECLUDED:
        base.run_secluded(args.ncores)
    for eps in base.EPS_JOINT:
        base.run_joint(eps, args.ncores)


if __name__ == "__main__":
    main()

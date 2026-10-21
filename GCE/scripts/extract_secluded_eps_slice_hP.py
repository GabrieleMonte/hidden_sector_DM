"""Write one eps of the Higgs-portal secluded grid out in the joint layout.

eps = 1e-10 is a JOINT eps -- GCE_fits.ipynb lists it in EPS_JOINT and reads
it from relic_grid_hP_eps1e-10.npz, a file with no leading eps axis -- but it
is SOLVED by the secluded route, because at that coupling the sectors have
long decoupled and the joint route only costs more.  Measured 2026-10-20 over
all 2083 nodes at eps = 1e-9, where both routes run: lam_secluded/lam_joint
has median 0.9988, is inside 5% on 98.7% of nodes and inside 10% on 99.95%,
with one disagreeing node at the rv = 0.950 edge.  The cost is not marginal --
24 minutes on 16 nodes for both eps, against a two-hour joint run that
reached 167 of 2083 nodes.  See make_relic_grid_hP.EPS_SECLUDED.

So the solve writes relic_grid_hP_eps_secluded.npz with eps = (1e-9, 1e-10)
on the leading axis, and this lifts one eps out of it into the single-eps
layout the notebook expects.  Every other field is copied unchanged; only
`eps`, `coupling_relic`, `sigmav` and `status` carry the eps axis.

    cosmo_env
    python GCE/scripts/extract_secluded_eps_slice_hP.py 1e-10
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_relic_grid as base                                    # noqa: E402
import make_relic_grid_hP as hp                                   # noqa: E402

# Fields that carry the leading eps axis in the secluded file; everything
# else is shared across eps and copied through untouched.
PER_EPS = ("coupling_relic", "sigmav", "status")

SRC = base.OUT_DIR / f"relic_grid_{hp.TAG_STEM}_eps_secluded.npz"


def extract(eps):
    """Write relic_grid_hP_eps{eps:.0e}.npz from the secluded grid's `eps`."""
    if not SRC.exists():
        raise SystemExit(f"missing {SRC}")
    src = np.load(SRC, allow_pickle=True)

    hit = np.flatnonzero(src["eps"] == eps)      # exact: every eps is a
    if hit.size != 1:                            # decade, atol would alias
        raise SystemExit(f"eps={eps:g} not in {src['eps']} (found {hit.size})")
    k = int(hit[0])

    out = {name: src[name] for name in src.files if name not in PER_EPS}
    out["eps"] = np.array([eps])
    for name in PER_EPS:
        out[name] = src[name][k]

    dst = base.OUT_DIR / f"relic_grid_{hp.TAG_STEM}_eps{eps:.0e}.npz"
    np.savez(dst, **out)

    ok = out["status"] == 0
    lam = out["coupling_relic"][ok]
    print(f"{dst.name}: {int(ok.sum())} of {int(out['n_done'])} nodes OK, "
          f"lam {lam.min():.5g} - {lam.max():.5g}")
    return dst


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__.strip().splitlines()[-1].strip())
    extract(float(sys.argv[1]))

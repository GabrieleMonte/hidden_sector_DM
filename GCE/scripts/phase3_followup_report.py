"""Phase-3 follow-up (docs/phase3_followup.md): what the run changed.  Read-only.

    python GCE/scripts/phase3_followup_report.py              # vs the pre-invalidation commit
    python GCE/scripts/phase3_followup_report.py --ref <sha>  # vs any commit

For each of the four re-solved grids, compares the working-tree file with the
same file at --ref (default: the commit recorded in the invalidation marker,
else HEAD), over the re-solved block only:
  - nodes still queued / NaN in the block (both should be ~0);
  - the spread of alpha_new / alpha_old;
  - at mX = 34.86, alpha_new / alpha_old against the shift the 2026-10-12
    diagnosis predicts, sqrt(Omega_forced/0.12) (Omega ~ alpha^-2 locally);
  - isolated nodes (> 20 % off the mean of their 4 mesh neighbours) before/after.
"""
import argparse
import io
import re
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "output" / "relic"
MARKER = OUT_DIR / "phase3_followup_invalidated.txt"

#        file                                    cols [lo, hi)
FILES = [("relic_grid_vP_xi_ini_2.5_eps1e-09.npz", (0, 50)),
         ("relic_grid_vP_xi_ini_2.5_eps1e-10.npz", (0, 57)),
         ("relic_grid_vP_eps1e-09.npz",            (0, 67)),
         ("relic_grid_vP_eps1e-10.npz",            (57, 67))]

# Omega_forced/0.12 at the stored crossing, mX index 20 (phase3_followup_preflight F2
# plus the rest of the 2026-10-12 bias map), keyed by file then rv column.
PREDICTED = {
    "relic_grid_vP_xi_ini_2.5_eps1e-09.npz": {26: 1.128, 34: 1.046, 42: 0.915, 48: 0.833},
    "relic_grid_vP_xi_ini_2.5_eps1e-10.npz": {26: 1.403, 34: 1.342, 42: 1.287, 48: 1.252,
                                              52: 1.233},
    "relic_grid_vP_eps1e-09.npz":            {26: 0.918, 34: 0.927, 42: 0.995, 48: 1.041,
                                              52: 1.000, 57: 1.001, 62: 0.996, 66: 1.003},
    "relic_grid_vP_eps1e-10.npz":            {57: 0.998, 62: 0.999, 66: 1.122},
}
I_REF = 20                                     # mX = 34.86


def _ref_default():
    if MARKER.exists():
        m = re.search(r"at commit ([0-9a-f]{7,40})", MARKER.read_text())
        if m:
            return m.group(1)
    return "HEAD"


def _load_at(ref, name):
    blob = subprocess.run(["git", "show", f"{ref}:GCE/output/relic/{name}"],
                          cwd=HERE, capture_output=True, check=True).stdout
    with np.load(io.BytesIO(blob)) as z:
        return {k: z[k] for k in z.files}


def _isolated(la, lo, hi):
    n = 0
    for i in range(1, la.shape[0] - 1):
        for j in range(max(lo, 1), min(hi, la.shape[1] - 1)):
            nb = np.array([la[i - 1, j], la[i + 1, j], la[i, j - 1], la[i, j + 1]])
            if np.all(np.isfinite(nb)) and np.isfinite(la[i, j]) \
                    and abs(la[i, j] - nb.mean()) > np.log(1.2):
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="Phase-3 follow-up: new vs old crossings")
    ap.add_argument("--ref", default=None, help="commit holding the pre-follow-up grids")
    ref = ap.parse_args().ref or _ref_default()
    print(f"comparing GCE/output/relic against commit {ref}\n")

    for name, (lo, hi) in FILES:
        new = dict(np.load(OUT_DIR / name))
        old = _load_at(ref, name)
        valid = np.isfinite(new["mX_grid"])
        blk = np.zeros(valid.shape, bool)
        blk[:, lo:hi] = True
        blk &= valid
        an, ao = new["coupling_relic"], old["coupling_relic"]
        queued = int((~new["done_mask"] & blk).sum())
        nan = int((~np.isfinite(an) & blk).sum())
        r = (an / ao)[blk & np.isfinite(an) & np.isfinite(ao)]
        print(f"{name}  cols {lo}-{hi - 1}, {int(blk.sum())} nodes")
        print(f"  still queued {queued}, NaN {nan}")
        if r.size:
            q = np.percentile(r, [0, 5, 50, 95, 100])
            print(f"  alpha_new/alpha_old: min {q[0]:.3f}  5% {q[1]:.3f}  median {q[2]:.3f}"
                  f"  95% {q[3]:.3f}  max {q[4]:.3f}")
        rv = 10.0 ** new["log_rv"]
        for j, om in PREDICTED[name].items():
            got = an[I_REF, j] / ao[I_REF, j]
            print(f"  mX=34.86 rv={rv[j]:.3f}: alpha ratio {got:.4f}  "
                  f"predicted {np.sqrt(om):.4f}  (<sigma v> x{got ** 2:.3f})")
        print(f"  isolated nodes in the block: before {_isolated(np.log(ao), lo, hi)}, "
              f"after {_isolated(np.log(an), lo, hi)}\n")


if __name__ == "__main__":
    main()

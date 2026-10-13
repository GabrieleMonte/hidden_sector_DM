"""Phase-3 follow-up ROUND 2 (docs/phase3_followup.md): what the run changed.
Read-only.

    python GCE/scripts/phase3_followup_report.py              # vs the round-1 result
    python GCE/scripts/phase3_followup_report.py --ref <sha>  # vs any commit

For the two re-solved xi=1 grids, over the re-solved block only:
  - nodes still queued / NaN (both should be ~0);
  - the path each alpha took (phase_path): Phase 3 at the 1e-3 trigger, Phase 3
    via the PHASE3_REQUIRE retry ladder, SM-locked, legacy (never engaged),
    failed, or skipped past the crossing -- and every node with an alpha
    legacy or failed, listed;
  - isolated nodes (> 20 % off the mean of their 4 mesh neighbours): round 1
    left 44 (eps=1e-9) and 14 (eps=1e-10);
  - alpha_new / alpha_round1 (default --ref: the commit in the round-2 marker);
  - alpha_new / alpha_campaign over the 2026-10 campaign's columns (commit in
    the round-1 marker): there the campaign values were right except at its own
    glitch nodes, so this should sit at ~1 with a short tail;
  - at mX = 34.86, alpha_new / alpha_campaign against sqrt(Omega_forced/0.12)
    from the 2026-10-12 diagnosis.
"""
import argparse
import io
import re
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "output" / "relic"
MARKER1 = OUT_DIR / "phase3_followup_invalidated.txt"      # round 1: pre-follow-up grids
MARKER2 = OUT_DIR / "phase3_followup2_invalidated.txt"     # round 2: round-1 grids

#        file                          re-solved cols [lo, hi)  campaign cols [lo, hi)
FILES = [("relic_grid_vP_eps1e-09.npz", (0, 67),               (50, 67)),
         ("relic_grid_vP_eps1e-10.npz", (57, 67),              (57, 67))]

# Omega_forced/0.12 at the pre-follow-up crossing, mX index 20 (2026-10-12 bias map).
# eps=1e-10 col 66 (rv=0.950) does NOT reproduce on the 8-point alpha axis: Omega
# alpha^2 has a ~20 % bump between the sampled alphas there.  rv=0.95 is above the
# notebook's eps=1e-10 ceiling (0.9), so it is masked downstream anyway.
PREDICTED = {
    "relic_grid_vP_eps1e-09.npz": {26: 0.918, 34: 0.927, 42: 0.995, 48: 1.041,
                                   52: 1.000, 57: 1.001, 62: 0.996, 66: 1.003},
    "relic_grid_vP_eps1e-10.npz": {57: 0.998, 62: 0.999, 66: 1.122},
}
P3, SMLOCK, LEGACY, FAILED, SKIPPED = 3.0, 4.0, 1.0, 0.0, -1.0   # == the drivers' PATH_*
I_REF = 20                                                       # mX = 34.86


def _marker_commit(marker, fallback):
    if marker.exists():
        m = re.search(r"at commit ([0-9a-f]{7,40})", marker.read_text())
        if m:
            return m.group(1)
    return fallback


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


def _quantiles(r):
    q = np.percentile(r, [0, 5, 50, 95, 100])
    return (f"min {q[0]:.3f}  5% {q[1]:.3f}  median {q[2]:.3f}  95% {q[3]:.3f}  "
            f"max {q[4]:.3f}")


def main():
    ap = argparse.ArgumentParser(description="Phase-3 follow-up round 2: new vs old crossings")
    ap.add_argument("--ref", default=None, help="commit holding the round-1 grids")
    ref1 = ap.parse_args().ref or _marker_commit(MARKER2, "HEAD")
    ref0 = _marker_commit(MARKER1, None)
    print(f"round-1 grids from commit {ref1}; campaign grids from commit {ref0}\n")

    for name, (lo, hi), (clo, chi) in FILES:
        new = dict(np.load(OUT_DIR / name))
        r1 = _load_at(ref1, name)
        valid = np.isfinite(new["mX_grid"])
        blk = np.zeros(valid.shape, bool)
        blk[:, lo:hi] = True
        blk &= valid
        an, a1 = new["coupling_relic"], r1["coupling_relic"]
        print(f"{name}  cols {lo}-{hi - 1}, {int(blk.sum())} nodes")
        print(f"  still queued {int((~new['done_mask'] & blk).sum())}, "
              f"NaN {int((~np.isfinite(an) & blk).sum())}")

        if "phase_path" in new:
            pp, sd = new["phase_path"][blk], new["p3_sdot"][blk]      # (nodes, alphas)
            print(f"  alpha paths ({int(np.isfinite(pp).sum())} on record): "
                  f"{int(((pp == P3) & (sd == 1e-3)).sum())} Phase 3 at Sdot 1e-3, "
                  f"{int(((pp == P3) & (sd < 1e-3)).sum())} via the retry ladder, "
                  f"{int((pp == SMLOCK).sum())} SM-locked, {int((pp == LEGACY).sum())} legacy, "
                  f"{int((pp == FAILED).sum())} failed, {int((pp == SKIPPED).sum())} skipped")
            odd = [(i, j) for i, j in zip(*np.where(blk))
                   if np.any(np.isin(new["phase_path"][i, j], (LEGACY, FAILED)))]
            print(f"  nodes with an alpha legacy or failed: {len(odd)}"
                  + "".join(f"\n    i={i} j={j} mX={10 ** new['log_mX'][i]:.2f} "
                            f"rv={10 ** new['log_rv'][j]:.3f} "
                            f"paths={new['phase_path'][i, j].tolist()}" for i, j in odd[:20])
                  + ("\n    ..." if len(odd) > 20 else ""))
        else:
            print("  no phase_path in this file -- not written by a round-2 driver")

        print(f"  isolated nodes in the block: round 1 {_isolated(np.log(a1), lo, hi)}, "
              f"now {_isolated(np.log(an), lo, hi)}")
        r = (an / a1)[blk & np.isfinite(an) & np.isfinite(a1)]
        if r.size:
            print(f"  alpha_new/alpha_round1: {_quantiles(r)};  "
                  f"{int((np.abs(r - 1) > 0.05).sum())} nodes moved > 5 %")

        if ref0:
            a0 = _load_at(ref0, name)["coupling_relic"]
            cb = np.zeros(valid.shape, bool)
            cb[:, clo:chi] = True
            cb &= valid
            r = (an / a0)[cb & np.isfinite(an) & np.isfinite(a0)]
            if r.size:
                print(f"  alpha_new/alpha_campaign, campaign cols {clo}-{chi - 1}: "
                      f"{_quantiles(r)};  {int((np.abs(r - 1) > 0.05).sum())} nodes > 5 % "
                      f"(the campaign's own glitch nodes are expected here)")
            rv = 10.0 ** new["log_rv"]
            for j, om in PREDICTED[name].items():
                got = an[I_REF, j] / a0[I_REF, j]
                print(f"  mX=34.86 rv={rv[j]:.3f}: alpha_new/alpha_pre-follow-up {got:.4f}  "
                      f"predicted {np.sqrt(om):.4f}")
        print()


if __name__ == "__main__":
    main()

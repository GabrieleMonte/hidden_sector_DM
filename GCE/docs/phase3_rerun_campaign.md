# Phase-3 rerun campaign — v3 (turnkey)

> **Executed** (job 59610830, `813de606bc`). Its rv cut left Phase-3-off values
> below it, wrong by up to 40 % at xi_ini = 2.5. **Follow-up:
> `phase3_followup.md`** (`sbatch GCE/scripts/phase3_followup.sbatch`).

v2 resolved the five on-cluster review findings; v3 ships the campaign as
code so nothing is edited by hand on the cluster. Solver: this branch at
`94a43d37b1` or later. Secluded slices are NOT rerun (no in-era injection).

## Quickstart (on the cluster, after pulling this branch)

```bash
# 1. Stage 1 — preflight, interactive/debug node, ~1 h on >= 30 cores
python GCE/scripts/phase3_preflight.py          # prints PASS/FAIL + GATES

# 2. Stage 2 — only after the GATES line says GO
python GCE/scripts/phase3_invalidate.py         # clears rerun columns, swaps 1e-9 axis

# 3. Stage 3 — pass A (eps=1e-9), one slice per job, xi=2.5 first
sbatch -q regular -t 24:00:00 GCE/scripts/relic_phase3.sbatch make_relic_grid_vP_xi_ini_2.5.py
sbatch -q regular -t 24:00:00 GCE/scripts/relic_phase3.sbatch make_relic_grid_vP.py
sbatch -q regular -t 24:00:00 GCE/scripts/relic_phase3.sbatch make_relic_grid_vP_xi_ini_0.5.py
# each job prints --dry-run first: the queue MUST be 782 nodes (460 in pass B)

# 4. Stage 4 — pass B (eps=1e-10), ONLY if preflight P3 said GO: in all three
#    drivers flip the two PASS-B constants (EPS_JOINT, ONLY_RV_COLUMNS, both
#    commented in place), then repeat step 3.
```

Everything below is the why; the quickstart is the what.

## What was committed on this branch (no hand edits left)

- The three drivers `make_relic_grid_vP*.py` carry the campaign config:
  `PHASE3_KW` (rho 0.3 / Sdot 0.05 switches, `enable_phase3=True`,
  `phase3_kw=dict(rtol_value=1e-6)`), the `phase_final in ("1.5S", 3)`
  Omega-ready return, NaN-node logging (`_log_nan` — the Phase-3
  `S_corr_max=2` ValueError is no longer swallowed silently),
  `ALPHAX_LIMS[1e-9] = (2e-5, 1.0)`, `SOLVE_TIMEOUT = 7200`, and PASS-A
  `EPS_JOINT` / `ONLY_RV_COLUMNS` / `RUN_SECLUDED=False`. Every campaign
  constant is tagged `CAMPAIGN` in a comment with its restore value.
  NOTE: the base driver was also refreshed — the copy previously on this
  branch was still frozen in its 2026-09-30 single-column campaign config.
- `phase3_preflight.py` — Stage-1 probes P1-P4 (see below), GATES summary.
- `phase3_invalidate.py` — Stage-2 state surgery, idempotent.
- `relic_phase3.sbatch` — one slice per job, `--requeue` + checkpoint-resume
  (a worst-case grid node is ~32 h, so 24 h slices simply requeue).

## Blocks being recomputed

A. eps = 1e-9, rv >= 0.465 (cols 50-66), xi_ini = 1 / 0.5 / 2.5: 782 nodes/slice.
   Why: B35 legacy conserved-S bias (+9/+43/+47 % in Omega at xi=1,
   rv = 0.556/0.869/0.909); at xi=2.5 the rv > 0.61 values are seam garbage
   outright (x3.55 off at the validated anchor).
B. eps = 1e-10, rv >= 0.636 (cols 57-66), all three slices: 460 nodes/slice.
   Nothing is NaN there and the stored alpha axis already reaches 1.0 —
   "untrusted" means the RV_CEILING[1e-10] = 0.61 ordering-test bound in
   notebooks/GCE_fit_vP_xi_test.ipynb (and the muX-gate audit's rv >= 0.73).
   Runs ONLY if preflight P3 shows Phase 3 engages there.

## Status of the hold in 94a43d37b1

Owner decision 2026-10-09: proceed. The hold wanted a scope map and an
independent review before re-deriving any grid; the xi=1/0.5 rerun doubles as
that test: B35 predicts alpha_relic shifts of x1.04 / x1.19 / x1.21 at
rv = 0.556 / 0.869 / 0.909 (xi=1, local slope -2). Landing there confirms the
legacy bias; landing at x1.00 falsifies it. Same compute either way.

## The alpha_X axis

Corrected crossings move UP (x1.04-1.21 at xi=1, x3.55 at the xi=2.5 corner)
and the old ALPHAX_LIMS[None] = (2e-5, 2e-2) leaves only x1.3-1.9 headroom in
the block; relic_alphaX never extrapolates, so production would have written
NaN at exactly the nodes the campaign pays for. The drivers now pin
`1e-9: (2e-5, 1.0)` (the axis eps=1e-10 and secluded already use; crossing-
interpolation error is negligible — measured |d2 lnOmega / d lnalpha2| ~ 0.001).
Raising the limit normally resets the file (load_state compares the saved
axis exactly), so phase3_invalidate.py ALSO rewrites the stored axis: the
kept low-rv crossings are axis-independent results, and after the swap
load_state resumes with them intact. KEEP the 1e-9 entry after the campaign.

## The preflight probes

- P1 (gates everything): alpha_relic at the two validated xi=2.5 anchors,
  end-to-end through the edited driver — 8.38e-3 (rv=0.950) and 6.40e-3
  (rv=0.9085) at mX=34.86, to ~2 %.
- P2 (cost): Phase-3 rtol 1e-6 vs 1e-8 on the rv=0.950 anchor node; Omega
  must agree to 0.5 % (validation D said it does); the timing ratio is the
  real per-node cost lever.
- P3 (gates pass B): eps=1e-10 at rv=0.9085, xi_ini = 1 and 2.5, three alphas
  spanning the crossing — if phase_final never comes back 3, Phase 3 would
  reproduce the legacy numbers there and pass B is dropped.
- P4 (gates pass A): the worst stored corner (xi=0.5, mX=71.37, rv=0.95,
  alpha_relic = 1.58e-2) must stay bracketed below 1.0 on the new axis.

## File-shape caveats

relic_grid_vP_eps1e-10.npz is 34 columns wide: load_state embeds it into the
67 mesh, so 230 of pass B's 460 xi=1 nodes are new midpoints, not recomputes,
and after the rerun the file is 67 wide with empty low-rv odd columns until
Stage 5. Nothing the notebook used is lost (all 34 old columns land on even
indices); griddata interpolates across the holes exactly as it did across the
coarser mesh. relic_grid_vP_eps1e-09.npz is itself only 2098/3082 done.

- **Stage 5 (after the campaign)**: backfill the xi=1 holes on the LEGACY
  path (restore the campaign constants first, so the backfilled nodes match
  their trusted neighbours): 662 nodes at rv < 0.465 in eps1e-09 (odd 67-mesh
  columns the 34->67 bisection never filled); optionally the empty midpoint
  columns below col 57 that the embed leaves in eps1e-10.

## Resources

~30k solves x 700-1500 s = ~45-100 node-hours at 256 workers/node against
~2,200 left on m3166: budget is fine, a redo is not — hence the preflight
gates and one-slice-at-a-time launches.

## Sanity anchors (production must reproduce to ~1-2 %)

mX = 34.86, xi_ini = 2.5, eps = 1e-9: alpha_relic = 8.38e-3 at rv = 0.950,
6.40e-3 at rv = 0.9085 (Phase-3 validated; local slope d lnOmega/d lnalpha = -2.0).

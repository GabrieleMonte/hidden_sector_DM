# Phase-3 follow-up (2026-10-12)

The 2026-10 campaign (`phase3_rerun_campaign.md`) re-solved only rv >= 0.465
(eps=1e-9) and rv >= 0.636 (eps=1e-10). Below that cut the grids still hold
Phase-3-OFF (legacy) crossings, and at xi_ini = 2.5 those are wrong by up to
40 % in <sigma v>. This job re-solves them, plus the xi=1 nodes whose alpha
axis straddled the trigger.

## Run it

```bash
sbatch GCE/scripts/phase3_followup.sbatch      # one job, 3 nodes, requeue-safe
```

The job runs the preflight, then the invalidation, then four slices, then a report.
Rehearse first with `PHASE3_FOLLOWUP_DRYRUN=1 bash GCE/scripts/phase3_followup.sbatch`.
When it ends, commit `GCE/output/relic` (including
`phase3_followup_invalidated.txt`) on this branch.

## What is re-solved (Phase 3 forced: Sdot_switch_thresh 1e-3)

| slice | rv columns | nodes | why |
|---|---|---|---|
| xi=2.5, eps=1e-9 | 0-49 (rv < 0.465) | 1301 | legacy below the cut: -17 % .. +13 % |
| xi=2.5, eps=1e-10 | 0-56 (rv < 0.636) | 1623 | legacy below the cut: +23 % .. +40 % |
| xi=1, eps=1e-9 | 0-66 (all) | 2083 | below the cut: -8 % .. +4 %; above it, 18 glitch nodes (e.g. rv=0.760, 0.7x their neighbours) |
| xi=1, eps=1e-10 | 57-66 (rv >= 0.636) | 460 | mixed path at rv=0.95: +12 % |

The percentages are the change in <sigma v>_relic: Omega_forced/0.12 at the
stored crossing, measured at mX=34.86 (rv 0.159-0.950). They are not
re-solved: xi=0.5 (at most 4 %, inside the trigger's design tolerance),
xi=1 eps=1e-10 below col 57 (under 1 %), and every secluded grid.

## Why the trigger is fine and forcing it is still right

Phase 3 fires when dlnS/dN >= 0.05. The legacy error grows with the entropy
injected by Y decays, which is ~xi^3, so at xi=2.5 the trigger fires at every
rv. Re-solving the stored crossings with the UNCHANGED campaign trigger engages
Phase 3 and agrees with the forced trigger to < 0.2 %. The large errors exist
only because the campaign's cut left those columns on the legacy path.

The trigger only does damage where it fires for some alphas of a node and not
others. The crossing then splines across two different solvers: those are the
xi=1 glitch nodes, plus the eps=1e-10 rv=0.95 column. Forcing (Sdot 1e-3) puts
every alpha on one path, and it reproduces the campaign wherever the campaign
engaged. That is why it is used for all four slices.

Full-axis check (8 alphas, both triggers, 2026-10-12):
- At xi=2.5 (eps=1e-9 rv=0.326; eps=1e-10 rv=0.05 and 0.426) the campaign
  trigger engages at all 8 alphas and agrees with forced to <= 0.2 %.
- At xi=1, eps=1e-9, rv=0.326 it leaves 2 of 8 alphas on the legacy path,
  including alpha=2.06e-3 next to the crossing: <sigma v> off by 5.5 %. So
  below the cut at xi=1 the campaign setting would itself make mixed-path
  nodes.
- At the glitch node mX=24.9, rv=0.760 (xi=1, eps=1e-9, inside the GCE 1-sigma
  region) the one legacy alpha puts the crossing at 0.77x the forced one:
  <sigma v> 41 % too low.
- Forced engaged at all 56 alphas tested, with no failures, at 30-95 s mean
  per solve (max 265 s).

## Gates

- `phase3_followup_preflight.py` (all through the generated drivers):
  - F1: the validated anchor alpha_relic = 8.386e-3 (xi=2.5, eps=1e-9,
    mX=34.86, rv=0.950) to 2 %.
  - F2: seven stored-crossing Omega/0.12 values reproduce the diagnosis
    to 1 %.
  - F3: the glitch node at alpha=2.065e-3 comes back on the Phase-3 branch
    (Omega 0.308), not the legacy one (0.181).
- `phase3_followup_invalidate.py` is run-once. It checks shape, alpha axis and
  eps first, then writes a marker, so a requeued job resumes instead of
  re-clearing finished nodes.
- Each slice aborts if its driver's `--dry-run` says "starting fresh" or queues
  more nodes than its block.

## Cost

5467 nodes x 8 alphas, about 20-200 s per solve. Expect ~4-6 h of wall clock
on 3 nodes, about 15-25 node-hours.

## After the run

`phase3_followup_report.py` (Stage 4 of the job) compares every file with the
pre-invalidation commit. Expect:
- 0 queued nodes;
- at mX=34.86, alpha_new/alpha_old ~ sqrt(Omega_forced/0.12), e.g. x1.18 at
  xi=2.5 eps=1e-10 rv=0.159;
- the xi=1 eps=1e-9 isolated-node count dropping from 18 toward 0.

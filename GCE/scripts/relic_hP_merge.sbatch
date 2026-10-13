#!/bin/bash
#SBATCH -J relic-hP-merge
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_merge_%x_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/relic_hP_merge_%x_%j.err
#SBATCH -C cpu
#SBATCH -A m3166
#SBATCH -N 1
#SBATCH -q debug
#SBATCH -t 00:20:00
#
# Fold one output's shard files into its canonical grid.  Submitted with
# --dependency=afterany on the shard job so it runs whatever happened to it:
# a merge that only runs on success is a merge that does not run, which is how
# 59695521 and 59713507 both ended with finished shards and no grid.
#
#   sbatch --dependency=afterany:<jobid> -J hP-secluded \
#          --export=ALL,SLICE=secluded GCE/scripts/relic_hP_merge.sbatch
#
# Each shard contributes only the nodes its own done_mask marks, so a partly
# finished shard folds in what it has and the rest stays queued for a resubmit.

module load python
source "$HOME/cosmo_env/bin/activate"
cd /global/u1/g/gab97/secluded_DM
: "${SLICE:?set SLICE=secluded|1e-9|1e-10 via --export}"

python GCE/scripts/make_relic_grid_hP.py --only "$SLICE" --nshards 16 --merge

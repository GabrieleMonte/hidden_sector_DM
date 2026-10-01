#!/bin/bash
#SBATCH -J chanspec_finish
#SBATCH -o /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/chanspec_finish_%j.out
#SBATCH -e /global/u1/g/gab97/secluded_DM/GCE/scripts/logs/chanspec_finish_%j.err
#SBATCH -q xfer
#SBATCH -A m3166
#SBATCH -t 12:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=montefalcone@wisc.edu
#
# Runs on a data transfer node (no compute allocation charged, $HOME mounted
# natively rather than over DVS).  Waits for every compute rank to leave the
# queue, reports how much of the grid landed, then copies the N400000 set from
# $PSCRATCH back into the repo.  The N200000 files already in the repo are a
# different filename set, so nothing is overwritten.

SCRATCH_DIR=$PSCRATCH/channel_spectra
REPO_DIR=/global/u1/g/gab97/secluded_DM/GCE/channel_spectra
MANIFEST=$PSCRATCH/channel_spectra_expected.txt

echo "waiting for compute ranks to finish ($(date))"
# A job whose dependency can never be met sits in PD as DependencyNeverSatisfied
# forever rather than being reaped, so it must not count as work still to come.
while squeue -u "$USER" -h -n chanspec -t PD,R,CG -O JobID,Reason 2>/dev/null \
        | grep -qv DependencyNeverSatisfied; do
    sleep 120
done
echo "compute queue clear ($(date))"

have=$(ls "$SCRATCH_DIR" | grep -c N400000)
want=$(wc -l < "$MANIFEST")
echo "grid completeness: $have / $want"

# Both sides must use the same collation or comm reports phantom differences:
# the manifest is written in byte order, so force LC_ALL=C on each.
ls "$SCRATCH_DIR" | grep N400000 | LC_ALL=C sort > /tmp/have.$$
LC_ALL=C sort "$MANIFEST" > /tmp/want.$$
comm -23 /tmp/want.$$ /tmp/have.$$ > "$PSCRATCH/channel_spectra_missing.txt"
missing=$(wc -l < "$PSCRATCH/channel_spectra_missing.txt")
echo "missing: $missing  (listed in \$PSCRATCH/channel_spectra_missing.txt)"
rm -f /tmp/have.$$ /tmp/want.$$

echo "copying to repo ($(date))"
rsync -a --include='*N400000*' --exclude='*' "$SCRATCH_DIR/" "$REPO_DIR/"
echo "in repo now: $(ls "$REPO_DIR" | grep -c N400000) N400000 files"
echo "done ($(date))"

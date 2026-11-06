#!/usr/bin/env bash
# Step 3 (alternative) of docs/new-cohort.md: no vendor target BED? Infer one for ONE assay from the
# coverage of a few of its BAMs. Plain bash - wrap it in your scheduler (it wants ~8 CPUs, ~16 GB, a few hours;
# e.g. `sbatch --wrap "THREADS=8 bash infer_targets.sh ..."` on SLURM).
#
#   bash infer_targets.sh <bam_list.txt> <ref.fai> <out_dir> [min_depth=20] [min_frac=0.5]
#
#   bam_list.txt  one BAM path per line, all from the SAME assay and panel version (8-12 is plenty;
#                 mix tumours and normals). Each BAM needs its .bai next to it.
#   ref.fai       index of the reference the BAMs were aligned to (contig names must match the BAMs).
#   out_dir       gets inferred_targets.bed + inferred_targets.summary.tsv.
#   min_depth     a base counts as covered in a sample at >= this depth (default 20x).
#   min_frac      a base is kept if covered in at least this fraction of samples (default 0.5).
#
# Check the summary before using the BED: consensus_Mb should be close to the vendor's stated panel size and
# every sample should have most of its covered bp inside the consensus. The summary holds only aggregate
# numbers (no intervals, no sample names); the BED itself derives from the BAMs - treat it like the data.
# b37-named output for a chr-named --ref? Run bed_to_ref_names.sh on it afterwards.
set -euo pipefail
LIST=${1:?bam_list.txt}; FAI=${2:?ref.fai}; OUT=${3:?out_dir}; MIN=${4:-20}; FRAC=${5:-0.5}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MOSDEPTH=${MOSDEPTH:-mosdepth}
THREADS=${THREADS:-4}

mkdir -p "$OUT/mosdepth"
n=0
while IFS= read -r bam; do
    [[ -z $bam || $bam == \#* ]] && continue
    n=$((n + 1))
    # Index-based prefix: no sample IDs in output names.
    MOSDEPTH_Q0=LOW MOSDEPTH_Q1=HIGH "$MOSDEPTH" -t "$THREADS" -n -x --mapq 20 \
        --quantize "0:${MIN}:" "$OUT/mosdepth/s${n}" "$bam"
    echo "sample $n done"
done < "$LIST"
(( n >= 2 )) || { echo "need at least 2 BAMs in $LIST" >&2; exit 1; }

python3 -E "$HERE/consensus_targets.py" --out "$OUT" --fai "$FAI" --min-frac "$FRAC" \
    "$OUT"/mosdepth/s*.quantized.bed.gz
echo "BED: $OUT/inferred_targets.bed   summary: $OUT/inferred_targets.summary.tsv"

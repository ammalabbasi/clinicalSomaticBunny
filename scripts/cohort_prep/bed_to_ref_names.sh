#!/usr/bin/env bash
# Step 3 of docs/new-cohort.md. Make a target BED use the same contig names as --ref, which --targets requires
# (the branch stops at startup otherwise: "--targets uses contigs that are not in --ref").
# b37-style names (1, X, MT) become chr1, chrX, chrM when that name is in REF.fai, and the reverse for a
# no-chr reference. Contigs absent from REF.fai are dropped and counted. Padding and merging are left to the
# pipeline (--target_padding). Prints counts only.
#
#   bash bed_to_ref_names.sh <targets.bed[.gz]> <ref.fasta.fai> <out.bed>
set -euo pipefail
BED=${1:?usage: bed_to_ref_names.sh targets.bed ref.fai out.bed}
FAI=${2:?ref.fai}
OUT=${3:?out.bed}
[[ -s $FAI ]] || { echo "ERROR: $FAI missing"; exit 1; }
rm -f "$OUT.dropped_contigs"
cat_bed() { if [[ $BED == *.gz ]]; then gzip -dc "$BED"; else cat "$BED"; fi; }

cat_bed | grep -vE '^(track|browser|#)' | awk -v OFS='\t' 'NF>=3{print $1,$2,$3}' \
  | awk -v OFS='\t' -v drop="$OUT.dropped_contigs" '
      NR==FNR { ok[$1]=1; next }
      { c=$1
        if (!(c in ok)) {
          if (c=="MT") c2="chrM"; else if (c=="chrM") c2="MT"; else if (c ~ /^chr/) c2=substr(c,4); else c2="chr" c
          if (c2 in ok) c=c2
        }
        if (c in ok) { $1=c; print } else { print $1 > drop } }' "$FAI" - \
  | sort -k1,1 -k2,2n > "$OUT"

n_in=$(cat_bed | grep -vcE '^(track|browser|#)' || true)
n_drop=0; [[ -f $OUT.dropped_contigs ]] && n_drop=$(wc -l < "$OUT.dropped_contigs")
echo "input intervals: $n_in   written: $(wc -l < "$OUT")   dropped (contig not in ref): $n_drop"
echo "bp: $(awk '{s+=$3-$2} END{print s+0}' "$OUT")   contigs: $(cut -f1 "$OUT" | sort -u | tr '\n' ' ')"
[[ -s $OUT ]] || { echo "ERROR: no interval matched $FAI (wrong build?)"; exit 1; }

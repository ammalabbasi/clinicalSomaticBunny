#!/usr/bin/env bash
# Step 2 of docs/new-cohort.md: your BAMs are b37-style (contigs 1, 2, MT) but the pipeline database is
# chr-named? Make a chr-named copy of the BAMs' OWN reference FASTA for --ref. Sequence is unchanged;
# only names change: 1..22,X,Y -> chr1..chrY, MT -> chrM; decoy/alt contigs (GL*, hs37d5, ...) keep their names.
#
#   SRC=/path/to/the/bams/reference.fasta OUT=/path/to/outdir bash make_chr_reference.sh
#   optional: SAMTOOLS=/path/to/samtools GATK=/path/to/gatk
#
# NOTE: gatk CreateSequenceDictionary needs Java 17+ on PATH (many cluster nodes default to Java 8:
#   export JAVA_HOME=/path/to/jdk17; export PATH="$JAVA_HOME/bin:$PATH").
# Use the FASTA the BAMs were aligned to (check_bams_aggregate.sh prints its UR), not a different build:
# GATK verifies contig ORDER as well as names and lengths.
set -euo pipefail
SRC=${SRC:?SRC=/path/to/reference.fasta}
OUT=${OUT:?OUT=/path/to/outdir}
SAMTOOLS=${SAMTOOLS:-samtools}
GATK=${GATK:-gatk}
FA=$OUT/$(basename "${SRC%.fa*}").chr.fasta
mkdir -p "$OUT"
awk '/^>/{ n=substr($1,2); if (n=="MT") n="chrM"; else if (n ~ /^([0-9]+|X|Y)$/) n="chr" n; print ">" n; next } {print}' "$SRC" > "$FA"
"$SAMTOOLS" faidx "$FA"
"$GATK" CreateSequenceDictionary -R "$FA" -O "${FA%.fasta}.dict"
awk -v OFS='\t' '{n=$1; if (n=="MT") n="chrM"; else if (n ~ /^([0-9]+|X|Y)$/) n="chr" n; print $1, n}' "$SRC.fai" > "$OUT/b37_to_chr.map" 2>/dev/null || true
echo "wrote $FA (+ .fai, .dict)"
echo "run the pipeline with: --ref $FA --ref_fai $FA.fai --ref_dict ${FA%.fasta}.dict --contig_rename auto"

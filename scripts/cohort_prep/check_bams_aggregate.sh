#!/usr/bin/env bash
# Step 1 of docs/new-cohort.md: what are these BAMs? RUN THIS AS THE DATA OWNER.
# Reads BAM *headers* only and prints aggregate facts per row - never read groups or sample names.
#
#   bash check_bams_aggregate.sh manifest.tsv            (header: patient sample status bam, tab-separated)
#   SAMTOOLS=/path/to/samtools REF_FAI=/path/to/ref.fai bash check_bams_aggregate.sh manifest.tsv
#
# Answers, per BAM: reference build (chr1 length: GRCh38=248956422, GRCh37=249250621), chr-prefix style,
# whether the contig set matches your --ref (REF_FAI), and which programs touched the BAM
# (dup-marking / BQSR / UMI consensus), which decides the entry point and whether --umi applies.
set -euo pipefail
MANIFEST=${1:?usage: check_bams_aggregate.sh manifest.tsv}
SAMTOOLS=${SAMTOOLS:-samtools}
REF_FAI=${REF_FAI:-}

i=0
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r patient sample status bam; do
    i=$((i+1))
    echo "== row $i ($status)"
    if [[ ! -s $bam ]]; then echo "  MISSING bam"; continue; fi
    bai=""; for c in "$bam.bai" "${bam%.bam}.bai"; do [[ -s $c ]] && bai=$c && break; done
    echo "  size: $(du -h --apparent-size "$bam" | cut -f1)   bai: ${bai:+present}${bai:-MISSING}"
    hdr=$("$SAMTOOLS" view -H "$bam")
    nsq=$(grep -c '^@SQ' <<<"$hdr" || true)
    first=$(grep -m1 '^@SQ' <<<"$hdr" | grep -o 'SN:[^[:space:]]*' | cut -d: -f2)
    chr1len=$(grep -P '^@SQ\tSN:(chr)?1\t' <<<"$hdr" | grep -o 'LN:[0-9]*' | cut -d: -f2)
    echo "  @SQ: $nsq contigs, first=$first, chr1 length=$chr1len  (GRCh38=248956422, GRCh37=249250621)"
    if [[ -n $REF_FAI ]]; then
        if diff -q <(grep '^@SQ' <<<"$hdr" | grep -o 'SN:[^[:space:]]*' | cut -d: -f2 | sort) \
                   <(cut -f1 "$REF_FAI" | sort) >/dev/null; then
            echo "  contigs MATCH $REF_FAI"
        else
            echo "  contigs DIFFER from $REF_FAI -> pass the BAMs' own FASTA as --ref (see docs/new-cohort.md step 2)"
        fi
    fi
    echo "  UR (reference named in header): $(grep -m1 -o 'UR:[^[:space:]]*' <<<"$hdr" || echo none)"
    echo "  @PG programs: $(grep '^@PG' <<<"$hdr" | grep -o 'PN:[^[:space:]]*' | cut -d: -f2 | sort -u | tr '\n' ' ')"
    echo "    (MarkDuplicates/sambamba/samblaster = dup-marked; ApplyBQSR/PrintReads = BQSR'd; fgbio = UMI consensus)"
done

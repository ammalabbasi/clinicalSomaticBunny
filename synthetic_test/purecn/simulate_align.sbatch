#!/usr/bin/env bash
#SBATCH --job-name=sb-purecn-sim
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --open-mode=append
#SBATCH --output=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/alignment_benchmark/somaticbunny_synthetic/agents/purecn/job-sim-%j.log
#
# Synthetic copy-number test data for PureCN/CNVkit (public GRCh38 + public gnomAD sites + simulated reads only):
# gen_cn_synth.py -> FASTQ -> bwa-mem2 (full GRCh38.d1.vd1 index) -> sorted, indexed BAMs in $D/bam.
set -euo pipefail
D=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/alignment_benchmark/somaticbunny_synthetic/agents/purecn
# sbatch runs a spooled copy, so the generator's directory is passed in (default: this repo's synthetic_test/purecn)
CODE=${PCN_CODE:-$SLURM_SUBMIT_DIR}
DB=/tscc/projects/ps-lalexandrov/shared/EVC_nextflow
REF=$DB/GRCh38_ref/GRCh38.d1.vd1.fa          # NOT bwa_mem_2/ (truncated index)
GNOMAD=$DB/Databases/GRCh38/af-only-gnomad.hg38_no_alt.vcf.gz
T=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/envs/bwa_env/bin
BCF=/tscc/projects/ps-lalexandrov/amabbasi/somaticbunny_envs/purecn_inspect/bin/bcftools
export PATH=/usr/local/bin:/usr/bin:/bin
mkdir -p "$D/bam"
cd "$D"
echo "== $(date -Is) simulate"
if [[ ! -s truth/samples.tsv ]]; then
    /tscc/projects/ps-lalexandrov/amabbasi/anaconda3/bin/python3 -E "$CODE/gen_cn_synth.py" "$REF" "$GNOMAD" "$BCF" "$T/samtools" "$D"
fi
echo "== $(date -Is) align"
for r1 in fastq/*_1.fq; do
    id=$(basename "$r1" _1.fq)
    [[ -s bam/$id.bam.bai ]] && continue
    "$T/bwa-mem2" mem -t 16 -R "@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA\tLB:${id}" "$REF" "$r1" "${r1%_1.fq}_2.fq" \
        | "$T/samtools" sort -@ 4 -m 2G -o "bam/$id.bam" -
    "$T/samtools" index "bam/$id.bam"
    "$T/samtools" flagstat "bam/$id.bam" | head -5
done
echo "== $(date -Is) done"
ls -l bam

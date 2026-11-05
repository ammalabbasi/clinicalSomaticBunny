#!/usr/bin/env bash
#SBATCH --job-name=sb-purecn-test
#SBATCH --account=ddp302
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --open-mode=append
#SBATCH --output=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/alignment_benchmark/somaticbunny_synthetic/agents/purecn/job-test-%j.log
#
# PureCN + CNVkit on the SYNTHETIC copy-number data of simulate_align.sbatch (public GRCh38, simulated reads).
#   run A (results_A): PCN01 tumour + matched normal (paired), PCN07 = the same tumour BAM as a tumour-only patient,
#                      PCN02..PCN06 normal-only pool rows. Builds and publishes the PureCN normal DB (6 normals).
#   run B (results_B): PCN08 = the same tumour alone, sheet purity 55 %, --purecn_normaldb results_A/PureCN/normaldb.
# Then evaluate_purecn.py scores purity, ploidy, integer CN, LOH and CNVkit log2 direction against truth/.
# usage: sbatch --export=ALL,PIPE=<pipeline checkout> run_purecn_synthetic.sbatch   (default PIPE: the submit dir's repo)
set -euo pipefail
if [[ -z ${SB_CLEAN_ENV:-} ]]; then
    keep=(SB_CLEAN_ENV=1 HOME="$HOME" USER="$USER" LOGNAME="${LOGNAME:-$USER}" SHELL=/bin/bash LANG=C.UTF-8
          PIPE="${PIPE:-}" SLURM_SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}" RUNS="${RUNS:-A B}")
    while IFS= read -r v; do keep+=("$v"); done < <(env | grep -E '^SLURM_[A-Z0-9_]*=')
    exec env -i "${keep[@]}" /bin/bash "$0" "$@"
fi

D=/tscc/lustre/restricted/alexandrov-ddn/users/amabbasi/alignment_benchmark/somaticbunny_synthetic/agents/purecn
PIPE=${PIPE:-$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)}
CODE=$PIPE/synthetic_test/purecn
DB=/tscc/projects/ps-lalexandrov/shared/EVC_nextflow
REF=$DB/GRCh38_ref/GRCh38.d1.vd1.fa
cd "$D"
export PATH=/usr/local/bin:/usr/bin:/bin
B=$D/bam
cat > samples_A.csv <<EOF
patient,sample,status,bam,bai
PCN01,N,normal,$B/PCN01_N_normal.bam,$B/PCN01_N_normal.bam.bai
PCN01,T,tumor,$B/PCN01_T_tumor.bam,$B/PCN01_T_tumor.bam.bai
PCN07,T,tumor,$B/PCN01_T_tumor.bam,$B/PCN01_T_tumor.bam.bai
PCN02,N,normal,$B/PCN02_N_normal.bam,$B/PCN02_N_normal.bam.bai
PCN03,N,normal,$B/PCN03_N_normal.bam,$B/PCN03_N_normal.bam.bai
PCN04,N,normal,$B/PCN04_N_normal.bam,$B/PCN04_N_normal.bam.bai
PCN05,N,normal,$B/PCN05_N_normal.bam,$B/PCN05_N_normal.bam.bai
PCN06,N,normal,$B/PCN06_N_normal.bam,$B/PCN06_N_normal.bam.bai
EOF
cat > samples_B.csv <<EOF
patient,sample,status,bam,bai,purity
PCN08,T,tumor,$B/PCN01_T_tumor.bam,$B/PCN01_T_tumor.bam.bai,55
EOF
md5sum bam/*.bam > inputs.md5

echo "== $(date -Is) SomaticBunny @ $(cd "$PIPE" && /usr/bin/git rev-parse --short HEAD) ($PIPE)"
source /tscc/projects/ps-lalexandrov/amabbasi/anaconda3/etc/profile.d/conda.sh
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL
export JAVA_HOME=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/envs/nf-env/lib/jvm
export PATH="$PATH:$JAVA_HOME/bin"
export NXF_VER=25.04.8
export TMPDIR=$D/tmp NXF_CONDA_CACHEDIR=/tscc/projects/ps-lalexandrov/amabbasi/somaticbunny_conda_cache_v2
export NXF_OPTS="-Djava.io.tmpdir=$TMPDIR -Xms1g -Xmx6g"
mkdir -p "$TMPDIR"
NF=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/envs/nf-env/bin/nextflow
common=(-profile conda -work-dir "$D/work" -resume
        --type panel --targets "$D/truth/targets.bed" --tool purecn,cnvkit --purecn_gc_norm false
        --genome GRCh38 --first_step variant_calling --database_path "$DB"
        --ref "$REF" --ref_fai "$REF.fai" --ref_dict "$DB/GRCh38_ref/GRCh38.d1.vd1.dict"
        --bed "$DB/GRCh38_ref/hg38_chr.bed.gz" --bed_tbi "$DB/GRCh38_ref/hg38_chr.bed.gz.tbi"
        --max_cpus 16 --max_memory 60.GB --max_time 8.h)

for run in $RUNS; do
    echo "== $(date -Is) run $run"
    extra=()
    [[ $run == B ]] && extra=(--purecn_normaldb "$D/results_A/PureCN/normaldb")
    "$NF" run "$PIPE/main.nf" "${common[@]}" --sample "$D/samples_$run.csv" --outdir "$D/results_$run" "${extra[@]}" \
        -with-trace "$D/trace_$run.txt" -with-report "$D/report_$run.html" || echo "RUN $run FAILED (exit $?)"
done

echo "== $(date -Is) input BAMs unchanged?"
md5sum -c inputs.md5

echo "== $(date -Is) evaluate"
P3=/tscc/projects/ps-lalexandrov/amabbasi/anaconda3/bin/python3
"$P3" -E "$CODE/evaluate_purecn.py" --truth "$D/truth" --out "$D/eval_purecn.tsv" \
    --sample "paired=$D/results_A/PureCN/PCN01_T" \
    --sample "tumor_only=$D/results_A/PureCN/PCN07_T" \
    --sample "tumor_only_prebuilt_db=$D/results_B/PureCN/PCN08_T" \
    --cnvkit "pooled_paired=$D/results_A/CNVkit/PCN01_T_calls.cns" \
    --cnvkit "pooled_tumor_only=$D/results_A/CNVkit/PCN07_T_calls.cns" \
    --cnvkit "flat_tumor_only=$D/results_B/CNVkit/PCN08_T_calls.cns"
echo "== $(date -Is) done"

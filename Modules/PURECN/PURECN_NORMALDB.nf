nextflow.enable.dsl=2

// Pool of normals for PureCN (NormalDB.R): coverage denoising (tangent normal), interval weights and filtering.
// Built once per run from all normals' GC-normalised coverages (needs >= 2). Published as
// ${params.purecn_dir}/normaldb/ together with the intervals it was built on, so a later run (e.g. a tumour-only
// cohort of the same assay) can reuse it with --purecn_normaldb <that dir>.
process PURECN_NORMALDB {
    scratch true
    label 'process_low'
    conda "${params.purecn_env}"
    publishDir("${params.purecn_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(coverages, stageAs: 'input/cov/*')
    path(intervals, stageAs: 'input/intervals/*')

    output:
    path("normaldb"), emit: normaldb

    script:
    def assay = params.purecn_assay ? "--assay ${params.purecn_assay}" : ""
    """
    set -euo pipefail
    unset R_LIBS R_LIBS_USER R_LIBS_SITE
    X=\$(Rscript --vanilla -e 'cat(system.file("extdata", package = "PureCN"))')
    mkdir -p normaldb
    ls -1 \$PWD/input/cov/*.txt.gz > normaldb/normal_coverages.list
    n=\$(wc -l < normaldb/normal_coverages.list)
    if [ "\$n" -lt 2 ]; then
        echo "ERROR: PureCN NormalDB.R needs at least 2 normal coverages, got \$n" >&2
        exit 1
    fi

    Rscript --vanilla \$X/NormalDB.R \\
        --out-dir normaldb \\
        --coverage-files normaldb/normal_coverages.list \\
        --genome ${params.purecn_genome} \\
        ${assay} --force

    ls normaldb/normalDB*.rds > /dev/null
    cp ${intervals} normaldb/intervals.txt
    # keep only file names in the published list (work paths are meaningless after the run)
    sed -i 's#.*/##' normaldb/normal_coverages.list
    echo "PureCN normal DB from \$n normals"
    """
}

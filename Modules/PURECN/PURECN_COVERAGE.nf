nextflow.enable.dsl=2

// Per-BAM target coverage and GC normalisation (Coverage.R). Runs on every BAM of the run: normals (paired and
// normal-only pool rows, which feed the normal DB) and tumours. Published so later runs can build a larger DB.
// --purecn_gc_norm false skips the GC loess (Coverage.R --skip-gc-norm) and emits the raw coverage. Default true:
// real data has GC bias. The synthetic test turns it off because simulated coverage has no GC signal (the loess then
// raised the per-target log-ratio SD from 0.13 to 1.4).
process PURECN_COVERAGE {
    tag "${meta.patient}_${meta.sample}_${meta.status}"
    scratch true
    label 'process_low'
    conda "${params.purecn_env}"
    publishDir("${params.purecn_dir}/coverage", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*')
    path(intervals, stageAs: 'input/intervals/*')

    output:
    tuple val(meta), path("${meta.patient}_${meta.sample}_${meta.status}_coverage${params.purecn_gc_norm ? '_loess' : ''}.txt.gz"), emit: coverage
    path("${meta.patient}_${meta.sample}_${meta.status}_coverage.txt.gz"), optional: true
    path("${meta.patient}_${meta.sample}_${meta.status}_coverage_loess*.{png,txt}"), optional: true

    script:
    def id = "${meta.patient}_${meta.sample}_${meta.status}"
    """
    set -euo pipefail
    unset R_LIBS R_LIBS_USER R_LIBS_SITE
    X=\$(Rscript --vanilla -e 'cat(system.file("extdata", package = "PureCN"))')
    mkdir -p cov

    # Coverage.R names its output after the BAM's real path (it resolves symlinks), so rename afterwards.
    Rscript --vanilla \$X/Coverage.R \\
        --out-dir cov \\
        --bam ${bam} \\
        --bai ${bai} \\
        --intervals ${intervals} \\
        ${params.purecn_gc_norm ? '' : '--skip-gc-norm'} \\
        --force

    shopt -s nullglob
    raw=(cov/*_coverage.txt.gz)
    loess=(cov/*_coverage_loess.txt.gz)
    want=${params.purecn_gc_norm ? 1 : 0}
    if [ \${#raw[@]} -ne 1 ] || [ \${#loess[@]} -ne \$want ]; then
        echo "ERROR: expected one coverage file and \$want GC-normalised file(s) from Coverage.R" >&2
        ls -l cov >&2
        exit 1
    fi
    mv "\${raw[0]}" ${id}_coverage.txt.gz
    if [ \$want -eq 1 ]; then mv "\${loess[0]}" ${id}_coverage_loess.txt.gz; fi
    for f in cov/*_coverage_loess.png cov/*_coverage_loess_qc.txt; do
        mv "\$f" "${id}_coverage_loess\${f##*_coverage_loess}"
    done
    """
}

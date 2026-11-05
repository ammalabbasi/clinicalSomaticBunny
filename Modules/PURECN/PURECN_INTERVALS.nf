nextflow.enable.dsl=2

// PureCN interval file from the capture BED (unpadded baits): GC content, gene symbols (TxDb for hg38/hg19),
// optional mappability and off-target bins. Built once per run; a prebuilt --purecn_normaldb brings its own.
process PURECN_INTERVALS {
    scratch true
    label 'process_low'
    conda "${params.purecn_env}"
    publishDir("${params.purecn_dir}/intervals", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(targets, stageAs: 'input/*')
    path(mappability, stageAs: 'input/mappability/*')

    output:
    path("intervals.txt"), emit: intervals
    path("intervals_optimized.bed")

    script:
    def map_arg = mappability ? "--mappability ${mappability}" : ""
    def off_target = params.purecn_off_target ? "--off-target" : ""
    """
    set -euo pipefail
    unset R_LIBS R_LIBS_USER R_LIBS_SITE
    X=\$(Rscript --vanilla -e 'cat(system.file("extdata", package = "PureCN"))')

    Rscript --vanilla \$X/IntervalFile.R \\
        --in-file ${targets} \\
        --fasta ${params.ref} \\
        --genome ${params.purecn_genome} \\
        --out-file intervals.txt \\
        --export intervals_optimized.bed \\
        ${map_arg} ${off_target} --force

    test -s intervals.txt
    echo "PureCN intervals: \$((\$(wc -l < intervals.txt) - 1))"
    """
}

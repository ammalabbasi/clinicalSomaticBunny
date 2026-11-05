nextflow.enable.dsl=2

// One `somalier relate` across every BAM of the run, then the identity checker
// (Modules/SOMALIER/somalier_check.py): same-patient samples with low relatedness (swap), different patients with
// high relatedness (duplicate / mislabel), and inferred vs declared sex when the sheet has a `sex` column.
process SOMALIER_RELATE {
    conda "${params.somalier_env}"
    scratch true
    label 'process_low'
    publishDir("${params.somalier_dir}", mode: 'copy', pattern: "somalier.*")
    publishDir("${params.somalier_dir}", mode: 'copy', pattern: "identity_sex.tsv")
    publishDir("${params.report_dir}", mode: 'copy', pattern: "identity_flags.tsv")
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(digests, stageAs: 'extract/*')
    path(sample_sheet, stageAs: 'input/*')

    output:
    path("somalier.pairs.tsv"), emit: pairs
    path("somalier.samples.tsv"), emit: samples
    path("somalier.html"), optional: true, emit: html
    path("identity_flags.tsv"), emit: flags
    path("identity_sex.tsv"), emit: sex

    script:
    """
    set -euo pipefail
    somalier relate ${params.somalier_relate_args} -o somalier extract/*.somalier
    python3 ${projectDir}/Modules/SOMALIER/somalier_check.py \\
        --pairs somalier.pairs.tsv --samples somalier.samples.tsv --sheet ${sample_sheet} \\
        --swap-max ${params.somalier_swap_max_relatedness} \\
        --dup-min ${params.somalier_dup_min_relatedness} \\
        --min-sites ${params.somalier_min_sites} \\
        --min-x-sites ${params.somalier_min_x_sites} \\
        --flags identity_flags.tsv --sex-out identity_sex.tsv
    """
}

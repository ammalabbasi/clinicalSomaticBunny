nextflow.enable.dsl=2

include { msiStatusCmd } from './msi_status'

// Tumour-only MSI: msisensor-pro pro against a baseline of normals from the same assay.
process MSI_TUMOR_ONLY {
    conda "${params.msi_env}"
    scratch true
    label 'process_low'
    publishDir("${params.msi_dir}/tumor_only/${map.patient}_${map.tumor_meta.sample}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)
    path(baseline, stageAs: 'input/*')

    output:
    tuple val(map), path("${prefix}.msi_status.tsv"), emit: status
    tuple val(map), path("${prefix}*"), emit: MSI_out

    script:
    prefix = "${map.patient}_${map.tumor_meta.sample}"
    """
    set -euo pipefail
    msisensor-pro pro -d ${baseline} -t ${map.tumor} -g ${params.ref} \\
        -c ${params.msi_min_coverage} -b ${task.cpus} -o ${prefix}
    ${msiStatusCmd(prefix, "${prefix}.msi_status.tsv", map.patient, map.tumor_meta.sample, 'tumor_only', params.msi_threshold_tumor_only)}
    """
}

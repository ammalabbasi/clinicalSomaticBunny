nextflow.enable.dsl=2

include { msiStatusCmd } from './msi_status'

// Paired tumour-normal MSI: msisensor-pro msi (per-site distribution comparison between tumour and its normal).
process MSI_PAIRED {
    conda "${params.msi_env}"
    scratch true
    label 'process_low'
    publishDir("${params.msi_dir}/paired/${map.patient}_${map.tumor_meta.sample}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)
    path(sites, stageAs: 'input/*')

    output:
    tuple val(map), path("${prefix}.msi_status.tsv"), emit: status
    tuple val(map), path("${prefix}*"), emit: MSI_out

    script:
    prefix = "${map.patient}_${map.tumor_meta.sample}"
    """
    set -euo pipefail
    msisensor-pro msi -d ${sites} -n ${map.normal} -t ${map.tumor} -g ${params.ref} \\
        -c ${params.msi_min_coverage} -b ${task.cpus} -o ${prefix}
    ${msiStatusCmd(prefix, "${prefix}.msi_status.tsv", map.patient, map.tumor_meta.sample, 'paired', params.msi_threshold_paired)}
    """
}

nextflow.enable.dsl=2

// Per-normal site distributions (msisensor-pro pro) that MSI_BASELINE pools into a tumour-only baseline.
// Normal-only sheet rows (pool normals) are included.
process MSI_PRO_NORMAL {
    conda "${params.msi_env}"
    scratch true
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(id), val(bam)
    path(sites, stageAs: 'input/*')

    output:
    path("${id}_all"), emit: all

    script:
    """
    set -euo pipefail
    msisensor-pro pro -d ${sites} -t ${bam} -g ${params.ref} -c ${params.msi_min_coverage} -b ${task.cpus} -o ${id}
    test -s ${id}_all
    """
}

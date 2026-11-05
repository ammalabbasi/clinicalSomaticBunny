nextflow.enable.dsl=2

// Genotype-like allele counts at the sites for one BAM (tumour, normal, or normal-only pool row). The digest is
// named <patient>_<sample>_<status> (SOMALIER_SAMPLE_NAME), not by the BAM's @RG SM, so vendor BAMs that reuse an SM
// (or have none) still give one unique, sheet-traceable name per row.
process SOMALIER_EXTRACT {
    conda "${params.somalier_env}"
    scratch true
    label 'process_low'
    publishDir("${params.somalier_dir}/extract", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*')
    path(sites)

    output:
    path("${meta.patient}_${meta.sample}_${meta.status}.somalier"), emit: somalier

    script:
    def id = "${meta.patient}_${meta.sample}_${meta.status}"
    """
    set -euo pipefail
    export SOMALIER_SAMPLE_NAME="${id}"
    somalier extract -d extracted --sites ${sites} -f ${params.ref} ${bam}
    mv "extracted/${id}.somalier" "${id}.somalier"
    """
}

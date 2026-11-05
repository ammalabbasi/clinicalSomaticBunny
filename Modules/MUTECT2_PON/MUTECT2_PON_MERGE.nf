nextflow.enable.dsl=2

// Cohort panel of normals, step 2: one bgzipped, indexed VCF per normal (GenomicsDBImport takes one VCF per sample).
process MUTECT2_PON_MERGE {
    scratch true
    conda "${params.gatk_env}"
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 3

    input:
    // Staged under input/ so no shard name can collide with the merged output.
    tuple val(id), path(vcfs, stageAs: 'input/*')

    output:
    tuple val(id), path("${id}.pon_normal.vcf.gz"), path("${id}.pon_normal.vcf.gz.tbi"), emit: normal_vcf

    script:
    def shards = vcfs instanceof List ? vcfs : [vcfs]
    def inputs = shards.collect { "-I ${it}" }.join(' ')
    """
    ${params.gatk} --java-options "-Xmx4g" MergeVcfs ${inputs} -O ${id}.pon_normal.vcf.gz
    """
}

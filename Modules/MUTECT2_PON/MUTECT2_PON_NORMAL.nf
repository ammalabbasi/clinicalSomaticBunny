nextflow.enable.dsl=2

// Cohort panel of normals, step 1: Mutect2 in tumour-only mode on one normal BAM over one interval shard
// (GATK PON recipe: --max-mnp-distance 0 so GenomicsDBImport can merge the calls). No germline resource
// and no PON here: the point is to see every recurrent site the assay produces in normals.
process MUTECT2_PON_NORMAL {
    scratch true
    conda "${params.gatk_env}"
    label 'MUTECT2_CALLING'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)
    each chunk
    path interval_dir

    output:
    tuple val("${patient}_${meta.sample}"), path("*.pon_normal.vcf"), emit: shard_vcf

    script:
    def xmx = task.memory ? Math.max(2, (task.memory.toGiga() * 0.8) as int) : 8
    """
    ${params.gatk} --java-options "-Xmx${xmx}g" Mutect2 \\
        -R ${params.ref} \\
        -I ${bam} \\
        --max-mnp-distance 0 \\
        --native-pair-hmm-threads ${task.cpus} \\
        -L ${interval_dir}/${chunk}-scattered.interval_list \\
        -O ${patient}_${meta.sample}_${chunk}.pon_normal.vcf
    """
}

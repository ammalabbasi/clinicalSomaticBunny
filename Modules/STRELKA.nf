nextflow.enable.dsl=2

process STRELKA {
    conda "${params.strelka_env}"
    scratch true
    label 'STRELKA'
    publishDir("${params.strelka_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)
    tuple path(call_regions), path(call_regions_tbi)

    output:
    path("*.somatic_indels.vcf.gz"), emit: STRELKA_indels
    path("*.somatic_indels.vcf.gz.tbi"), emit: STRELKA_indels_tbi
    path("*.somatic_snvs.vcf.gz"), emit: STRELKA_snvs
    path("*.somatic_snvs.vcf.gz.tbi"), emit: STRELKA_snvs_tbi
    tuple val(map),path("*.somatic_snvs.vcf.gz"),path("*.somatic_indels.vcf.gz"), emit: STRELKA_out


    script:
    // Clinical: --callRegions = padded --targets when given, else the callable-regions --bed (as upstream);
    // --exome for exome and panel (no depth filters).
    def exome_flag = params.targeted ? "--exome" : ""
    def prefix     = "${map.patient}_${map.tumor_meta.sample}"
    """
    configureStrelkaSomaticWorkflow.py --referenceFasta ${params.ref} --normalBam ${map.normal} --tumorBam ${map.tumor} \\
        --callRegions ${call_regions} ${exome_flag} --runDir strelka

    python2 strelka/runWorkflow.py -m local -j ${task.cpus}

    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi
    """
}

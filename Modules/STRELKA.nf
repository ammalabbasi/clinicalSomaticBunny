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
    // genome: callable regions (--bed). exome/panel: --exome (no depth filters), restricted to --targets when given.
    def exome_flag = params.targeted ? "--exome" : ""
    def regions = (params.targeted && !params.targets) ? "" : "--callRegions ${call_regions}"
    """
    configureStrelkaSomaticWorkflow.py --referenceFasta ${params.ref} --normalBam ${map.normal} --tumorBam ${map.tumor} --runDir strelka ${exome_flag} ${regions}
    python2 strelka/runWorkflow.py -m local -j ${task.cpus}
    mv strelka/results/variants/somatic.indels.vcf.gz     ${map.patient}_${map.tumor_meta.sample}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${map.patient}_${map.tumor_meta.sample}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${map.patient}_${map.tumor_meta.sample}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${map.patient}_${map.tumor_meta.sample}.somatic_snvs.vcf.gz.tbi
    """
}

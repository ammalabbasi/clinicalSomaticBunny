nextflow.enable.dsl=2

process MOSDEPTH {
    conda "${params.mosdepth_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.mosdepth_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3
        
    input:
    val(map)

    output:
    tuple val(map.patient), path("*txt"), emit: coverage
    tuple val(map.patient), path("*txt"), emit: MOSDEPTH_out
    tuple val(map.patient), path("*regions.bed.gz"), emit:regions_bed

    script:
    def by = params.targeted ? "--by ${params.mosdepth_bed}" : "-n --fast-mode --by 500"
    def normal_cmd = map.normal ? "mosdepth -t ${task.cpus} ${by} ${map.patient}_${map.tumor_meta.sample}_normal ${map.normal}" : ""
    """
    mosdepth -t ${task.cpus} ${by} ${map.patient}_${map.tumor_meta.sample}_tumor ${map.tumor}
    ${normal_cmd}
    """
}

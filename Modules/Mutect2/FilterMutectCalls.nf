nextflow.enable.dsl=2

process FilterMutectCalls {
    scratch true
    label 'process_low'
    conda "${params.gatk_env}"
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy = 'retry'
    maxRetries 3

    input:
    tuple val(map),
          path(unfiltered_vcf),
          path(contamination_table),
          path(segments_table),
          path(read_orientation_model_tar),
          path(merged_stats)

    output:
    tuple val(map), path("*vcf"), path("*idx"), path("*stats"), emit: MUTECT2_final_out
    tuple val(map), path("*vcf"), emit: Mutect2_out

    script:
    """
    CONT_ARGS=""
    if [ -s ${contamination_table} ]; then CONT_ARGS="--contamination-table ${contamination_table}"; fi
    if [ -s ${segments_table} ]; then CONT_ARGS="\$CONT_ARGS --tumor-segmentation ${segments_table}"; fi
    ${params.gatk} FilterMutectCalls \
    -R ${params.ref} \
    -V ${unfiltered_vcf} \
    \$CONT_ARGS \
    --ob-priors ${read_orientation_model_tar} \
    -O ${map.patient}_${map.tumor_meta.sample}_mutect2_filtered.vcf \
    --stats ${merged_stats} \
    --filtering-stats ${map.patient}_${map.tumor_meta.sample}_mutect2_filtered.stats
    """
}

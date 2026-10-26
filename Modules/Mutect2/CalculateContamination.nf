nextflow.enable.dsl=2

process CalculateContamination {
    scratch true
    label 'process_low'
    conda "${params.gatk_env}"
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    tuple val(map), path("*contamination.table"), emit: MUTECT2_contamination_table
    tuple val(map), path("*segments.table"), emit: MUTECT2_segments_table
    path("*contamination.FALLBACK.txt"), optional: true

    script:
    def matched = map.normal ? "-matched ${map.normal}" : ""
    def prefix = "${map.patient}_${map.tumor_meta.sample}"
    // Panels and tumour-only samples can have too few common germline SNPs for an estimate. For those only,
    // write empty tables plus a published FALLBACK marker; FilterMutectCalls then runs without contamination
    // instead of the whole run stopping. Paired WGS/WES keep the old behaviour (the task fails).
    // Signals (exit >= 128, e.g. an OOM kill) always fail the task so errorStrategy can retry it.
    def allow_fallback = (!map.normal || params.type == 'panel') ? 1 : 0
    """
    rc=0
    ${params.gatk} CalculateContamination -I ${map.tumor} ${matched} -O ${prefix}_contamination.table --tumor-segmentation ${prefix}_segments.table || rc=\$?
    if [ \$rc -ge 128 ] || { [ \$rc -ne 0 ] && [ ${allow_fallback} -eq 0 ]; }; then exit \$rc; fi
    if [ \$rc -ne 0 ]; then
        echo "WARNING: CalculateContamination failed for ${prefix}; filtering without a contamination estimate" >&2
        echo "CalculateContamination exit \$rc: ${prefix} filtered without contamination/segmentation tables" > ${prefix}_contamination.FALLBACK.txt
        : > ${prefix}_contamination.table
        : > ${prefix}_segments.table
    fi
    """
}

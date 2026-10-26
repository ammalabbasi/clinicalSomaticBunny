nextflow.enable.dsl=2

// 20 Mutect2/GetPileupSummaries shards over the padded targets instead of whole chromosomes.
// Same layout and 1-based names as GENERATE_INTERVALS, so the 20-chunk contract downstream holds.
process GENERATE_TARGET_INTERVALS {
    conda "${params.gatk_env}"
    scratch true
    label 'process_low'
    publishDir("${params.intervals_dir}", mode: 'copy')

    input:
    path(padded_bed)
    path(ref)
    path(ref_fai)
    path(ref_dict)

    output:
    path "interval_list", emit: interval_dir

    script:
    """
    mkdir -p interval_list
    # INTERVAL_SUBDIVISION (default) splits intervals as needed, so even a small panel gives 20 shards.
    ${params.gatk} SplitIntervals \\
        -R ${ref} \\
        -L ${padded_bed} \\
        --scatter-count 20 \\
        -O interval_list

    cd interval_list
    for f in *-scattered.interval_list; do
        idx=\$(echo "\${f}" | grep -oP '^\\d+')
        num=\$((10#\${idx} + 1))
        mv "\${f}" "\${num}-scattered.interval_list"
    done
    n=\$(ls *-scattered.interval_list | wc -l)
    if [ "\$n" -ne 20 ]; then
        echo "ERROR: expected 20 interval shards over --targets, got \$n" >&2
        exit 1
    fi
    """
}

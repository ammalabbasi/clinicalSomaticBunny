nextflow.enable.dsl=2

// Group reads into molecular families by mapping position + UMI (fgbio GroupReadsByUmi); sets the MI tag.
// The family-size histogram is the file MultiQC's fgbio module reads (detected by its header).
process UMI_GROUP {
    conda "${params.umi_env}"
    scratch true
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 3
    publishDir("${params.umi_dir}", mode: 'copy', pattern: '*.{txt,tsv}')

    input:
    tuple val(meta), path(bam, stageAs: 'input/*')

    output:
    tuple val(meta), path("*_umi_grouped.bam"), emit: grouped
    path("*.fgbio_group_family_size_histogram.txt"), emit: family_sizes
    path("*.fgbio_group_metrics.txt"), emit: metrics

    script:
    def id = "${meta.patient}_${meta.sample}_${meta.status}"
    """
    set -o pipefail
    fgbio -Xmx${Math.max(2, task.memory.toGiga() - 2)}g -Djava.io.tmpdir=. --compression 1 GroupReadsByUmi \\
        --input ${bam} --output ${id}_umi_grouped.bam \\
        --strategy ${params.umi_group_strategy_resolved} --edits ${params.umi_edits} \\
        --min-map-q ${params.umi_min_map_q} --raw-tag RX --assign-tag MI \\
        --family-size-histogram ${id}.fgbio_group_family_size_histogram.txt \\
        --grouping-metrics ${id}.fgbio_group_metrics.txt
    """
}

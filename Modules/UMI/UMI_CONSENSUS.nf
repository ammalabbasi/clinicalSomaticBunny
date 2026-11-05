nextflow.enable.dsl=2

// One consensus read pair per molecular family (unmapped output; UMI_ALIGN 'consensus' re-aligns and filters it).
//   simplex: CallMolecularConsensusReads --min-reads <first value of --umi_min_reads>
//   duplex : CallDuplexConsensusReads (needs GroupReadsByUmi --strategy paired and a UMI on both ends) and
//            CollectDuplexSeqMetrics on the grouped reads.
// The consensus read group keeps ID <patient> and SM <patient>_<sample>_<status>, like BWA_MEM's.
process UMI_CONSENSUS {
    conda "${params.umi_env}"
    scratch true
    label 'process_medium'
    errorStrategy 'retry'
    maxRetries 3
    publishDir("${params.umi_dir}", mode: 'copy', pattern: '*.{txt,pdf}')

    input:
    tuple val(meta), path(bam, stageAs: 'input/*')

    output:
    tuple val(meta), path("*_umi_consensus_unmapped.bam"), emit: ubam
    path("*.fgbio_duplex.*"), emit: duplex_metrics, optional: true

    script:
    def id = "${meta.patient}_${meta.sample}_${meta.status}"
    def mem = Math.max(2, task.memory.toGiga() - 2)
    def min_reads = params.umi_min_reads.toString().trim().split(/\s+/)
    def common = "--input ${bam} --output ${id}_umi_consensus_unmapped.bam --read-group-id ${meta.patient} " +
                 "--read-name-prefix ${id} --threads ${task.cpus}"
    if (params.umi_duplex)
        """
        set -o pipefail
        fgbio -Xmx${mem}g -Djava.io.tmpdir=. --compression 1 CallDuplexConsensusReads ${common} \\
            --min-reads ${min_reads.join(' ')} --error-rate-pre-umi 45 --error-rate-post-umi 40
        fgbio -Xmx${mem}g -Djava.io.tmpdir=. CollectDuplexSeqMetrics --input ${bam} --output ${id}.fgbio_duplex
        """
    else
        """
        set -o pipefail
        fgbio -Xmx${mem}g -Djava.io.tmpdir=. --compression 1 CallMolecularConsensusReads ${common} \\
            --min-reads ${min_reads[0]} --error-rate-pre-umi 45 --error-rate-post-umi 40
        """
}

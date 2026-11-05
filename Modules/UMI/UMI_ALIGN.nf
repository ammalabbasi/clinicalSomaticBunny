nextflow.enable.dsl=2

// Align an unmapped BAM with bwa-mem2 and merge the alignments back onto the unmapped records with fgbio
// ZipperBams, which keeps RX/MI and the per-base consensus tags (reversed for reverse-strand reads).
//   stage 'raw'       : raw reads -> template-coordinate order, the input order GroupReadsByUmi wants
//   stage 'consensus' : consensus reads -> FilterConsensusReads (needs mapped reads + --ref) -> coordinate-sorted,
//                       indexed BAM; this replaces MarkDuplicates' output (duplicates are already collapsed).
process UMI_ALIGN {
    conda "${params.umi_env}"
    scratch true
    label 'BWA_MEM'
    errorStrategy 'retry'
    maxRetries 3
    publishDir("${params.umi_dir}", mode: 'copy', pattern: '*.{txt,tsv}')
    if (params.publish.tokenize(',').contains('markdup_bam')) {
        publishDir("${params.umi_dir}", mode: 'copy', pattern: '*_umi_consensus.bam*')
    }

    input:
    tuple val(meta), path(ubam, stageAs: 'input/*')
    val(stage)

    output:
    tuple val(meta), path("*_umi_raw_mapped.bam"), emit: raw, optional: true
    tuple val(meta.patient), val(meta), path("*_umi_consensus.bam"), path("*_umi_consensus.bam.bai"), emit: pair_mutect, optional: true
    tuple val(meta), path("*_umi_consensus.bam"), path("*_umi_consensus.bam.bai"), emit: mkdup_bam, optional: true
    path("*.umi_consensus_filter_counts.tsv"), emit: counts, optional: true
    val(meta), emit: cleanup_trigger

    script:
    def id = "${meta.patient}_${meta.sample}_${meta.status}"
    def sort_threads = Math.max(1, task.cpus.intdiv(4))
    def bwa_threads = Math.max(1, task.cpus - sort_threads)
    def batch_bases = 10000000L * task.cpus
    def zipper = "fgbio -Xmx4g -Djava.io.tmpdir=. --compression 0 ZipperBams --unmapped ${ubam} --ref ${params.ref} --tags-to-reverse Consensus --tags-to-revcomp Consensus"
    def align = "samtools fastq -@ 2 ${ubam} | bwa-mem2 mem -t ${bwa_threads} -p -K ${batch_bases} -Y ${params.ref} - | ${zipper}"
    if (stage == 'raw')
        """
        set -o pipefail
        ${align} \\
          | samtools sort -@ ${sort_threads} --template-coordinate -T ${id}_tc -o ${id}_umi_raw_mapped.bam -
        """
    else
        """
        set -o pipefail
        ${align} \\
          | fgbio -Xmx8g -Djava.io.tmpdir=. --compression 0 FilterConsensusReads \\
                --input /dev/stdin --output /dev/stdout --ref ${params.ref} \\
                --min-reads ${params.umi_min_reads} \\
                --min-base-quality ${params.umi_min_base_quality} \\
                --max-base-error-rate ${params.umi_max_base_error_rate} \\
                --max-no-calls ${params.umi_max_no_call_fraction} \\
          | samtools sort -@ ${sort_threads} -T ${id}_co -o ${id}_umi_consensus.bam -
        samtools index ${id}_umi_consensus.bam
        printf 'sample\\tconsensus_reads_in\\tconsensus_reads_kept\\n' > ${id}.umi_consensus_filter_counts.tsv
        printf '%s\\t%s\\t%s\\n' ${id} "\$(samtools view -c ${ubam})" "\$(samtools view -c -F 0x900 ${id}_umi_consensus.bam)" \\
            >> ${id}.umi_consensus_filter_counts.tsv
        """
}

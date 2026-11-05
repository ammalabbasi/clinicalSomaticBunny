nextflow.enable.dsl=2

// FASTQ(s) -> unmapped BAM with the UMI in the RX tag (fgbio FastqToBam).
// UMI source (params.umi_source, resolved in main.nf):
//   inline : UMI bases inside R1/R2, cut out by --umi_read_structure (e.g. "3M2S+T 3M2S+T")
//   fastq  : a separate UMI read (sheet column fastq_umi), structure default "+T +T +M"
//   header : UMI is the last ':' field of the read name (Illumina BCL Convert style); fgbio moves it to RX
// Read group matches BWA_MEM: ID <patient>, SM <patient>_<sample>_<status>, PL ILLUMINA.
process UMI_FASTQ_TO_BAM {
    conda "${params.umi_env}"
    scratch true
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(meta), val(fastq_umi)

    output:
    tuple val(meta), path("*_unmapped.bam"), emit: ubam

    script:
    def id = "${meta.patient}_${meta.sample}_${meta.status}"
    def inputs = [meta.fastq_1, meta.fastq_2] + (params.umi_source == 'fastq' ? [fastq_umi] : [])
    def header_opt = params.umi_source == 'header' ? '--extract-umis-from-read-names true' : ''
    """
    set -o pipefail
    fgbio -Xmx${Math.max(2, task.memory.toGiga() - 2)}g -Djava.io.tmpdir=. --compression 1 FastqToBam \\
        --input ${inputs.join(' ')} \\
        --read-structures ${params.umi_read_structures_resolved} \\
        --umi-tag RX ${header_opt} \\
        --sample ${id} --library ${id} \\
        --read-group-id ${meta.patient} --platform ILLUMINA \\
        --output ${id}_unmapped.bam
    """
}

nextflow.enable.dsl=2

// Tumour mutational burden for one tumour: callable territory from mosdepth --quantize on the tumour (and the
// matched normal for pairs), restricted to --tmb_regions; eligible consensus variants counted by tmb.py.
// mode: 'paired' | 'tumor_only'. normal_bam/normal_bai/germline may be [] (not staged).
process TMB {
    conda "${params.tmb_env}"
    scratch true
    label 'process_low'
    publishDir("${params.tmb_dir}", mode: 'copy', pattern: "${patient_sample}.tmb.*")
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient_sample), val(mode), val(source),
          path(vcf, stageAs: 'input/vcf/*'),
          path(tumor_bam, stageAs: 'input/tumor/*'), path(tumor_bai, stageAs: 'input/tumor/*'),
          path(normal_bam, stageAs: 'input/normal/*'), path(normal_bai, stageAs: 'input/normal/*'),
          path(germline, stageAs: 'input/germline/*')
    path(regions, stageAs: 'input/regions/*')

    output:
    tuple val(patient_sample), path("${patient_sample}.tmb.tsv"), emit: tmb_tsv
    tuple val(patient_sample), path("${patient_sample}.tmb.variants.tsv"), path("${patient_sample}.tmb.callable.bed"), emit: details

    script:
    def paired = (mode == 'paired')
    def flags = "-n -F 3844 --mapq ${params.tmb_min_mapq}"
    def germ_arg = (!paired && germline) ? "--germline_calls ${germline} --min_somatic_posterior ${params.tmb_min_somatic_posterior}" : ""
    def coding_arg = params.tmb_coding_only ? "--coding_only" : ""
    """
    set -euo pipefail
    # mosdepth/htslib find the index next to the BAM; vendor indexes may be named X.bai instead of X.bam.bai
    [[ -e "${tumor_bam}.bai" ]] || ln -s "\$(basename ${tumor_bai})" "${tumor_bam}.bai"
    ${paired ? "[[ -e \"${normal_bam}.bai\" ]] || ln -s \"\$(basename ${normal_bai})\" \"${normal_bam}.bai\"" : ""}

    MOSDEPTH_Q0=LOW MOSDEPTH_Q1=CALLABLE \\
        mosdepth -t ${task.cpus} ${flags} --quantize 0:${params.tmb_min_depth}: tmbq_tumor ${tumor_bam}
    ${paired ? "MOSDEPTH_Q0=LOW MOSDEPTH_Q1=CALLABLE mosdepth -t ${task.cpus} ${flags} --quantize 0:${params.tmb_min_normal_depth}: tmbq_normal ${normal_bam}" : ""}

    python3 -E ${projectDir}/Modules/TMB/tmb.py \\
        --sample ${patient_sample} --mode ${mode} --source ${source} \\
        --vcf ${vcf} --tumor_bam ${tumor_bam} --tumor_bai ${tumor_bai} \\
        --regions ${regions} \\
        --tumor_quantized tmbq_tumor.quantized.bed.gz \\
        ${paired ? "--normal_quantized tmbq_normal.quantized.bed.gz" : ""} \\
        --min_depth ${params.tmb_min_depth} --min_normal_depth ${params.tmb_min_normal_depth} \\
        --min_vaf ${params.tmb_min_vaf} --min_mapq ${params.tmb_min_mapq} --min_baseq ${params.tmb_min_baseq} \\
        ${coding_arg} ${germ_arg} \\
        --out ${patient_sample}.tmb
    """
}

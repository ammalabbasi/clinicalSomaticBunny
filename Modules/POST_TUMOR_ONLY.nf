nextflow.enable.dsl=2

// Consensus for tumour-only samples. Strelka2 somatic and MuSE need a matched normal, so only Mutect2
// (tumour-only, PON + gnomAD) and SAGE (tumour-only) run. Keep variants that PASS in both, and soft-flag
// likely germline sites by population frequency (Mutect2 POPAF = -log10 AF; POPAF < 3 means AF > 0.001).
process POST_TUMOR_ONLY {
    conda "${params.bcftools_env}"
    scratch true
    label 'process_low'
    publishDir("${params.post_dir}/tumor_only/${patient_sample}", mode: 'copy')

    input:
    tuple val(patient_sample), path(mutect2_vcf), path(sage_vcf)

    output:
    tuple val(patient_sample), path("${patient_sample}_tumor_only_consensus.vcf.gz"), path("${patient_sample}_tumor_only_consensus.vcf.gz.tbi"), emit: consensus
    tuple val(patient_sample), path("${patient_sample}_tumor_only_consensus.PASS.vcf.gz"), path("${patient_sample}_tumor_only_consensus.PASS.vcf.gz.tbi"), emit: consensus_pass

    script:
    """
    bcftools view -f PASS ${mutect2_vcf} | bcftools norm -m -any -f ${params.ref} -Oz -o mutect2.pass.vcf.gz
    bcftools index -t mutect2.pass.vcf.gz
    bcftools view -f PASS ${sage_vcf} | bcftools norm -m -any -f ${params.ref} -Oz -o sage.pass.vcf.gz
    bcftools index -t sage.pass.vcf.gz

    # Records called by both; Mutect2's record (with POPAF) is kept.
    bcftools isec -n=2 -w1 mutect2.pass.vcf.gz sage.pass.vcf.gz -Oz -o both.vcf.gz

    bcftools filter -m + -s germline_popaf -e 'INFO/POPAF<3' both.vcf.gz -Oz -o ${patient_sample}_tumor_only_consensus.vcf.gz
    bcftools index -t ${patient_sample}_tumor_only_consensus.vcf.gz
    bcftools view -f PASS ${patient_sample}_tumor_only_consensus.vcf.gz -Oz -o ${patient_sample}_tumor_only_consensus.PASS.vcf.gz
    bcftools index -t ${patient_sample}_tumor_only_consensus.PASS.vcf.gz
    """
}

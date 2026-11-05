nextflow.enable.dsl=2

// Panel-specific somalier sites (--somalier_sites targets): common biallelic gnomAD SNPs
// (params.somalier_min_af <= AF <= params.somalier_max_af) from the Mutect2 germline resource, inside the padded,
// merged targets. Built once per run (one assay per run). The public somalier release sites sit mostly in exomes and
// can miss a small panel almost entirely; this set uses every common SNP the panel actually covers.
process SOMALIER_SITES {
    conda "${params.somalier_env}"
    scratch true
    label 'process_low'
    publishDir("${params.somalier_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(padded_bed, stageAs: 'input/*')
    tuple path(germline), path(germline_tbi)

    output:
    path("somalier_sites.targets.vcf.gz"), emit: sites

    script:
    """
    set -euo pipefail
    bcftools view -R ${padded_bed} -v snps -m2 -M2 -f .,PASS \\
        -i 'INFO/AF>=${params.somalier_min_af} && INFO/AF<=${params.somalier_max_af}' ${germline} \\
      | bcftools annotate -x '^INFO/AF' \\
      | bcftools sort -Oz -o somalier_sites.targets.vcf.gz
    # -R can emit a site twice when padded targets overlap after merging per contig; keep one record.
    bcftools norm -d exact somalier_sites.targets.vcf.gz -Oz -o sites.dedup.vcf.gz
    mv sites.dedup.vcf.gz somalier_sites.targets.vcf.gz
    bcftools index -t somalier_sites.targets.vcf.gz
    n=\$(bcftools view -H somalier_sites.targets.vcf.gz | wc -l)
    echo "somalier targets sites: \$n common SNPs (AF ${params.somalier_min_af}-${params.somalier_max_af}) in the padded targets"
    if [ "\$n" -lt ${params.somalier_min_sites} ]; then
        echo "WARNING: only \$n sites; somalier relatedness on this panel will be noisy (min ${params.somalier_min_sites})" >&2
    fi
    if [ "\$n" -eq 0 ]; then
        echo "ERROR: no common SNPs in the targets; use --somalier_sites release or a sites VCF" >&2
        exit 1
    fi
    """
}

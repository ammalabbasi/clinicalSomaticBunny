nextflow.enable.dsl=2

// Common germline SNP sites for PureCN's B-allele frequencies: biallelic SNVs with population AF >= --purecn_min_pop_af
// in the public gnomAD resource that Mutect2 already uses (params.mutect2_germline), inside the padded targets.
// Built once per run. Outputs:
//   purecn_sites.vcf.gz         sites + INFO/AF (annotation source for DB / POP_AF)
//   purecn_sites.alleles.tsv.gz CHROM POS REF,ALT for `bcftools call -C alleles`
//   purecn_regions.bed          padded, merged targets (mpileup -R)
process PURECN_SITES {
    scratch true
    label 'process_low'
    conda "${params.purecn_env}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(targets, stageAs: 'input/*')
    path(germline_vcf)
    path(germline_tbi)

    output:
    tuple path("purecn_sites.vcf.gz"), path("purecn_sites.vcf.gz.tbi"), path("purecn_sites.alleles.tsv.gz"),
          path("purecn_sites.alleles.tsv.gz.tbi"), path("purecn_regions.bed"), emit: sites

    script:
    """
    set -euo pipefail
    awk -v OFS='\\t' -v pad=${params.target_padding} '!/^(track|browser|#)/ && NF>=3 {
            s=\$2-pad; if (s<0) s=0; print \$1, s, \$3+pad }' ${targets} \\
      | sort -k1,1 -k2,2n \\
      | awk -v OFS='\\t' 'c==\$1 && \$2<=e { if (\$3>e) e=\$3; next }
                         { if (c!="") print c, s, e; c=\$1; s=\$2; e=\$3 }
                         END { if (c!="") print c, s, e }' > purecn_regions.bed

    bcftools view -R purecn_regions.bed -v snps -m2 -M2 -i 'INFO/AF>=${params.purecn_min_pop_af}' ${germline_vcf} -Ou \\
      | bcftools annotate -x ID,QUAL,FILTER,^INFO/AF -Ou \\
      | bcftools sort -Oz -o purecn_sites.vcf.gz
    bcftools index -t purecn_sites.vcf.gz
    n=\$(bcftools index -n purecn_sites.vcf.gz)
    if [ "\$n" -eq 0 ]; then
        echo "ERROR: no gnomAD SNP with AF >= ${params.purecn_min_pop_af} inside the targets." >&2
        echo "       Do the contig names of ${germline_vcf} match the targets (chr1 vs 1)?" >&2
        exit 1
    fi
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' purecn_sites.vcf.gz | bgzip -c > purecn_sites.alleles.tsv.gz
    tabix -s1 -b2 -e2 purecn_sites.alleles.tsv.gz
    echo "PureCN SNP sites: \$n in \$(wc -l < purecn_regions.bed) padded target regions"
    """
}

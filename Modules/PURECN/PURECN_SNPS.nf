nextflow.enable.dsl=2

// Allele counts for PureCN, one VCF per tumour, WITHOUT touching the Mutect2 somatic calls (they are only read):
//   sites  = common gnomAD SNPs in the padded targets (PURECN_SITES) + this tumour's Mutect2 PASS biallelic SNVs
//   counts = bcftools mpileup (FORMAT/AD) on the tumour [and its matched normal], reduced to ref/alt counts of the
//            given alleles by resources/purecn_counts.awk (no genotyper); sites with no alt read in any sample dropped
//   INFO   = DB (site in the gnomAD set, i.e. population AF >= --purecn_min_pop_af) and POP_AF (gnomAD AF, 0-1 scale)
//            -> PureCN's germline prior; paired runs also get SOMATIC on the Mutect2 PASS sites (the matched normal
//            already vetted them), so PureCN uses its SOMATIC-based priors and drops non-het normal germline sites.
//            Tumour-only runs get no SOMATIC flag: PureCN then classifies every non-DB variant itself.
// Sample columns: <patient>_<sample> (tumour; equals PureCN --sampleid) [and <patient>_<sample>_matched_normal].
process PURECN_SNPS {
    tag "${map.patient}_${map.tumor_meta.sample}"
    scratch true
    label 'process_low'
    conda "${params.purecn_env}"
    publishDir("${params.purecn_dir}/${map.patient}_${map.tumor_meta.sample}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(map), path(mutect2_vcf, stageAs: 'input/m2/*'),
          path(tumor_bam, stageAs: 'input/t/*'), path(tumor_bai, stageAs: 'input/t/*'),
          path(normal_bam, stageAs: 'input/n/*'), path(normal_bai, stageAs: 'input/n/*')
    tuple path(sites_vcf), path(sites_tbi), path(alleles_tsv), path(alleles_tbi), path(regions_bed)

    output:
    tuple val(map), path("${map.patient}_${map.tumor_meta.sample}_purecn_input.vcf.gz"),
          path("${map.patient}_${map.tumor_meta.sample}_purecn_input.vcf.gz.tbi"), emit: vcf

    script:
    def id = "${map.patient}_${map.tumor_meta.sample}"
    def paired = normal_bam ? true : false
    def bams = paired ? "${tumor_bam}\n${normal_bam}" : "${tumor_bam}"
    def names = paired ? "${id}\n${id}_matched_normal" : "${id}"
    """
    set -euo pipefail

    # Mutect2 PASS biallelic SNVs of this tumour (read only)
    bcftools view -f PASS -v snps -m2 -M2 ${mutect2_vcf} -Ou \\
      | bcftools annotate -x ID,INFO,FORMAT -Ou \\
      | bcftools view -G -Ou \\
      | bcftools sort -Oz -o somatic_sites.vcf.gz
    bcftools index -t somatic_sites.vcf.gz
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' somatic_sites.vcf.gz > somatic.alleles.tsv
    echo "Mutect2 PASS SNVs: \$(wc -l < somatic.alleles.tsv)"

    # Sites to count: gnomAD sites not at a somatic position, plus the somatic SNVs (their own alleles)
    { bgzip -dc ${alleles_tsv} | awk -F'\\t' 'NR==FNR { s[\$1":"\$2]=1; next } !((\$1":"\$2) in s)' somatic.alleles.tsv -
      cat somatic.alleles.tsv; } | sort -k1,1 -k2,2n -u > alleles.tsv
    bgzip -f alleles.tsv
    tabix -f -s1 -b2 -e2 alleles.tsv.gz
    bgzip -dc alleles.tsv.gz | cut -f1,2 > positions.tsv
    # mpileup regions: padded targets plus 1-bp windows at the somatic sites (merged, so no record is emitted twice)
    { cat ${regions_bed}; awk -v OFS='\\t' '{ print \$1, \$2-1, \$2 }' somatic.alleles.tsv; } \\
      | sort -k1,1 -k2,2n \\
      | awk -v OFS='\\t' 'c==\$1 && \$2<=e { if (\$3>e) e=\$3; next }
                         { if (c!="") print c, s, e; c=\$1; s=\$2; e=\$3 }
                         END { if (c!="") print c, s, e }' > regions.bed

    printf '${bams}\\n' > bams.txt
    printf '${names}\\n' > names.txt
    bcftools mpileup -f ${params.ref} -b bams.txt --ignore-RG \\
        -R regions.bed -T positions.tsv \\
        -a FORMAT/AD -I -q ${params.purecn_min_mapq} -Q ${params.purecn_min_bq} -d 100000 -Ov -o pileup.vcf
    bgzip -dc alleles.tsv.gz > alleles.txt
    awk -f ${moduleDir}/resources/purecn_counts.awk alleles.txt pileup.vcf > counts.raw.vcf
    rm -f pileup.vcf alleles.txt
    bcftools reheader -s names.txt counts.raw.vcf \\
      | bcftools view -Oz -o counts.vcf.gz
    bcftools index -t counts.vcf.gz

    printf '##INFO=<ID=POP_AF,Number=A,Type=Float,Description="gnomAD population allele frequency (from ${params.mutect2_germline})">\\n' > popaf.hdr
    bcftools annotate -a ${sites_vcf} -c INFO/POP_AF:=INFO/AF -h popaf.hdr \\
        -m '+DB' counts.vcf.gz -Oz -o db.vcf.gz
    bcftools index -t db.vcf.gz

    if [ "${paired}" = "true" ]; then
        bcftools annotate -a somatic_sites.vcf.gz -m '+SOMATIC' db.vcf.gz -Oz -o ${id}_purecn_input.vcf.gz
    else
        mv db.vcf.gz ${id}_purecn_input.vcf.gz
    fi
    bcftools index -f -t ${id}_purecn_input.vcf.gz
    echo "PureCN input: \$(bcftools index -n ${id}_purecn_input.vcf.gz) sites," \\
         "\$(bcftools view -H -i 'INFO/DB=1' ${id}_purecn_input.vcf.gz | wc -l) DB-flagged"
    """
}

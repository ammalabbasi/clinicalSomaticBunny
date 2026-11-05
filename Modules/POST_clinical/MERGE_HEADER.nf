nextflow.enable.dsl=2

// Clinical SomaticBunny copy of Modules/MERGE_HEADER.nf (original kept unchanged).
// 1. Writes the original standard VCF exactly as before (same commands, same file name).
// 2. Also writes <id>_<scenario>_standard_PASSed.fixed.vcf.gz(.tbi): a valid sites-only VCF made from it
//    (Modules/POST_clinical/fix_standard_vcf.awk, then bcftools sort/index/view). Its header records the exact
//    commands used: the original merge, the fix, and the bcftools commands (bcftools also adds its own
//    ##bcftools_*Command lines).
process MERGE_HEADER {
    conda "${params.POST_rename_env}"
    scratch true
    label 'process_low'
    publishDir("${params.post_dir}/standard_vcf/${map[0]}", mode: 'copy', pattern: "*_standard_PASSed.vcf")
    publishDir("${params.post_dir}/standard_vcf/${map[0]}", mode: 'copy', pattern: "*_standard_PASSed.fixed.vcf.gz*")

    input:
    tuple val(map), val(scenario)

    output:
    path("*_standard_PASSed.vcf"), emit: final_standard_vcf
    tuple path("*_standard_PASSed.fixed.vcf.gz"), path("*_standard_PASSed.fixed.vcf.gz.tbi"), emit: fixed_standard_vcf

    script:
    def orig = "${map[0]}_${scenario}_standard_PASSed.vcf"
    def fixed = "${map[0]}_${scenario}_standard_PASSed.fixed.vcf.gz"
    def merge_cmd = "bash merge_header.sh ${map[1].name} ${map[3].name} ${map[2].name} ${map[4].name} > merged_header.vcf; cat merged_header.vcf ${map[5].name} > ${orig}"
    def fix_cmd = "awk -v provfile=provenance.txt -f fix_standard_vcf.awk ${orig} > fixed.unsorted.vcf"
    def sort_cmd = "bcftools sort -Oz -o ${fixed} fixed.unsorted.vcf"
    def index_cmd = "bcftools index -t ${fixed}"
    def check_cmd = "bcftools view ${fixed} > /dev/null"
    """
    # 1. original standard VCF, unchanged
    bash ${projectDir}/merge_header.sh \\
        ${map[1]} \\
        ${map[3]} \\
        ${map[2]} \\
        ${map[4]} > merged_header.vcf

    cat merged_header.vcf ${map[5]} > ${orig}

    # 2. fixed, valid VCF with the exact commands in its header
    cat > provenance.txt <<'PROV'
##SomaticBunny_source=${orig}
##SomaticBunny_originalMergeCommand=${merge_cmd}
##SomaticBunny_fixCommand=${fix_cmd}
##SomaticBunny_bcftoolsSortCommand=${sort_cmd}
##SomaticBunny_bcftoolsIndexCommand=${index_cmd}
##SomaticBunny_bcftoolsValidateCommand=${check_cmd}
PROV
    awk -v provfile=provenance.txt -f ${projectDir}/Modules/POST_clinical/fix_standard_vcf.awk ${orig} > fixed.unsorted.vcf
    ${sort_cmd}
    ${index_cmd}
    ${check_cmd}
    echo "fixed VCF records: \$(bcftools view -H ${fixed} | wc -l); original records: \$(grep -vc '^#' ${orig} || true)"
    """
}

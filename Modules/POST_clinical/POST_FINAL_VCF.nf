nextflow.enable.dsl=2

// Clinical SomaticBunny: for each tumour and consensus scenario, write the upstream POST_SCENARIOS final VCF as a
// sorted, indexed <id>_final.vcf.gz whose header records the exact commands that produced it (the consensus
// command from POST_SCENARIOS and every bcftools command used here; bcftools also adds its own ##bcftools_* lines),
// and a consensus_status.tsv row (n_snv, n_indel, OK | NO_CONSENSUS_MUTATIONS). The upstream plain final VCF is
// left unchanged.
process POST_FINAL_VCF {
    conda "${params.POST_env}"
    scratch true
    label 'process_low'
    tag "${patient_sample} ${scenario}"
    publishDir("${params.post_dir}/${patient_sample}/${scenario}", mode: 'copy',
               pattern: "${patient_sample}_*{final.vcf.gz,final.vcf.gz.tbi,consensus_status.tsv}")

    input:
    tuple val(patient_sample), val(scenario), path(final_vcf, stageAs: 'input/*'), path(provenance)

    output:
    tuple val(patient_sample), val(scenario), path("${patient_sample}_final.vcf.gz"), path("${patient_sample}_final.vcf.gz.tbi"), emit: vcf
    path("${patient_sample}_${scenario}_consensus_status.tsv"), emit: status

    script:
    def out = "${patient_sample}_final.vcf.gz"
    def ann = "bcftools annotate -h prov.txt -Ov -o with_prov.vcf ${final_vcf.name}"
    def srt = "bcftools sort -Oz -o ${out} with_prov.vcf"
    def idx = "bcftools index -t ${out}"
    def chk = "bcftools view ${out} > /dev/null"
    """
    {
        grep '^##' ${provenance} || true
        echo "##SomaticBunny_scenario=${scenario}"
        echo "##SomaticBunny_source=${scenario}/${final_vcf.name} (POST_SCENARIOS)"
        echo "##SomaticBunny_bcftoolsAnnotateCommand=${ann}"
        echo "##SomaticBunny_bcftoolsSortCommand=${srt}"
        echo "##SomaticBunny_bcftoolsIndexCommand=${idx}"
        echo "##SomaticBunny_bcftoolsValidateCommand=${chk}"
    } > prov.txt

    bcftools annotate -h prov.txt -Ov -o with_prov.vcf ${final_vcf}
    ${srt}
    ${idx}
    ${chk}

    n_snv=\$(bcftools view -H -v snps ${out} | wc -l)
    n_indel=\$(bcftools view -H -v indels ${out} | wc -l)
    status=OK
    if [ "\$n_snv" -eq 0 ] && [ "\$n_indel" -eq 0 ]; then status=NO_CONSENSUS_MUTATIONS; fi
    printf 'patient_sample\\tscenario\\tn_snv\\tn_indel\\tstatus\\n%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${patient_sample}" "${scenario}" "\$n_snv" "\$n_indel" "\$status" > ${patient_sample}_${scenario}_consensus_status.tsv
    echo "${patient_sample} ${scenario}: \$n_snv SNVs, \$n_indel indels (\$status)"
    """
}

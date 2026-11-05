nextflow.enable.dsl=2

// Cohort panel of normals, step 3: GenomicsDBImport over the padded targets (or the Mutect2 interval shards
// without --targets), then CreateSomaticPanelOfNormals with the germline resource, so common germline sites
// are kept out of the panel. Default --min-sample-count 2: a site enters the PON when >= 2 normals carry it.
process MUTECT2_PON_CREATE {
    scratch true
    conda "${params.gatk_env}"
    label 'process_medium'
    publishDir("${params.mutect2_cohort_pon_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(normal_vcfs, stageAs: 'input/*')
    // padded targets BED, or the directory of *-scattered.interval_list shards
    path(intervals, stageAs: 'intervals/*')

    output:
    tuple path("${params.mutect2_cohort_pon_name}.pon.vcf.gz"), path("${params.mutect2_cohort_pon_name}.pon.vcf.gz.tbi"), emit: pon
    path("${params.mutect2_cohort_pon_name}.pon.normals.txt"), emit: normals

    script:
    def name = params.mutect2_cohort_pon_name
    def files = normal_vcfs instanceof List ? normal_vcfs : [normal_vcfs]
    def vcfs = files.findAll { it.name.endsWith('.vcf.gz') }
    if (vcfs.size() < 2) {
        error "MUTECT2_PON_CREATE needs at least 2 normal VCFs, got ${vcfs.size()}"
    }
    def xmx = task.memory ? Math.max(4, (task.memory.toGiga() * 0.8) as int) : 8
    """
    : > intervals.args
    if [ -d ${intervals} ]; then
        for f in ${intervals}/*-scattered.interval_list; do printf -- '-L\n%s\n' "\$f" >> intervals.args; done
    else
        printf -- '-L\n%s\n' "${intervals}" >> intervals.args
    fi
    if [ ! -s intervals.args ]; then
        echo "ERROR: no intervals for GenomicsDBImport in ${intervals}" >&2
        exit 1
    fi
    cp intervals.args import.args
    for v in ${vcfs.join(' ')}; do printf -- '-V\n%s\n' "\$v" >> import.args; done
    mkdir -p tmp
    ${params.gatk} --java-options "-Xmx${xmx}g -Djava.io.tmpdir=\$PWD/tmp" GenomicsDBImport \\
        -R ${params.ref} \\
        --arguments_file import.args \\
        --merge-input-intervals \\
        --genomicsdb-shared-posixfs-optimizations true \\
        --tmp-dir tmp \\
        --genomicsdb-workspace-path pon_db

    ${params.gatk} --java-options "-Xmx${xmx}g -Djava.io.tmpdir=\$PWD/tmp" CreateSomaticPanelOfNormals \\
        -R ${params.ref} \\
        --germline-resource ${params.mutect2_germline} \\
        -V gendb://pon_db \\
        -O ${name}.pon.vcf.gz

    for v in ${vcfs.join(' ')}; do basename "\$v" .pon_normal.vcf.gz; done > ${name}.pon.normals.txt
    echo "Mutect2 cohort PON ${name}: \$(wc -l < ${name}.pon.normals.txt) normals, \$(zcat ${name}.pon.vcf.gz | grep -vc '^#') sites"
    """
}

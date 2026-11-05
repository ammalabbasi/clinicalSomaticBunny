nextflow.enable.dsl=2

include { PURECN_INTERVALS } from './PURECN_INTERVALS'
include { PURECN_SITES } from './PURECN_SITES'
include { PURECN_COVERAGE } from './PURECN_COVERAGE'
include { PURECN_NORMALDB } from './PURECN_NORMALDB'
include { PURECN_SNPS } from './PURECN_SNPS'
include { PURECN } from './PURECN'

// PureCN for paired and tumour-only tumours (--tool purecn). Called from both calling blocks of main.nf.
//   samples: [patient, meta, bam, bai] for every BAM of the run (normals incl. normal-only rows, and tumours)
//   mutect2: FilterMutectCalls.out.Mutect2_out, [map, vcf]; map.normal is null for tumour-only patients
workflow PURECN_WF {
    take:
    samples
    mutect2

    main:
    if (params.purecn_normaldb) {
        intervals_ch = Channel.value(file("${params.purecn_normaldb}/intervals.txt", checkIfExists: true))
        normaldb_ch = Channel.value(file(params.purecn_normaldb, type: 'dir', checkIfExists: true))
    } else {
        def mappability = params.purecn_mappability ? file(params.purecn_mappability, checkIfExists: true) : []
        intervals_ch = PURECN_INTERVALS(file(params.purecn_targets), mappability).intervals
    }

    PURECN_SITES(file(params.purecn_targets), file(params.mutect2_germline), file("${params.mutect2_germline}.tbi"))

    PURECN_COVERAGE(samples.map { patient, meta, bam, bai -> [meta, bam, bai] }, intervals_ch)
    coverage = PURECN_COVERAGE.out.coverage

    if (!params.purecn_normaldb) {
        if (params.purecn_build_normaldb) {
            normaldb_ch = PURECN_NORMALDB(
                coverage.filter { meta, cov -> meta.status == 'normal' }.map { meta, cov -> cov }.collect(),
                intervals_ch).normaldb
        } else {
            // a single normal: no DB; its paired tumour is normalised against the matched normal
            normaldb_ch = Channel.value([])
        }
    }

    tumor_bams = samples.filter { it[1].status == 'tumor' }
        .map { patient, meta, bam, bai -> ["${patient}_${meta.sample}", bam, bai] }
    normal_bams = samples.filter { it[1].status == 'normal' }
        .map { patient, meta, bam, bai -> [patient, bam, bai] }
    tumor_cov = coverage.filter { meta, cov -> meta.status == 'tumor' }
        .map { meta, cov -> ["${meta.patient}_${meta.sample}", cov] }
    normal_cov = coverage.filter { meta, cov -> meta.status == 'normal' }
        .map { meta, cov -> [meta.patient, cov] }

    m2 = mutect2.map { map, vcf -> ["${map.patient}_${map.tumor_meta.sample}", map, vcf] }
    snps_paired = m2.filter { it[1].normal }
        .join(tumor_bams)
        .map { key, map, vcf, tbam, tbai -> [map.patient, map, vcf, tbam, tbai] }
        .combine(normal_bams, by: 0)
        .map { patient, map, vcf, tbam, tbai, nbam, nbai -> [map, vcf, tbam, tbai, nbam, nbai] }
    snps_tonly = m2.filter { !it[1].normal }
        .join(tumor_bams)
        .map { key, map, vcf, tbam, tbai -> [map, vcf, tbam, tbai, [], []] }

    PURECN_SNPS(snps_paired.mix(snps_tonly), PURECN_SITES.out.sites)

    with_tcov = PURECN_SNPS.out.vcf
        .map { map, vcf, tbi -> ["${map.patient}_${map.tumor_meta.sample}", map, vcf, tbi] }
        .join(tumor_cov)
    purecn_paired = with_tcov.filter { it[1].normal }
        .map { key, map, vcf, tbi, tcov -> [map.patient, map, vcf, tbi, tcov] }
        .combine(normal_cov, by: 0)
        .map { patient, map, vcf, tbi, tcov, ncov -> [map, vcf, tbi, tcov, ncov] }
    purecn_tonly = with_tcov.filter { !it[1].normal }
        .map { key, map, vcf, tbi, tcov -> [map, vcf, tbi, tcov, []] }

    PURECN(purecn_paired.mix(purecn_tonly), normaldb_ch, intervals_ch)

    emit:
    purity_ploidy = PURECN.out.purity_ploidy
    segments = PURECN.out.segments
}

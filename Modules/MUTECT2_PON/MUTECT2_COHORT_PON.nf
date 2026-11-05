nextflow.enable.dsl=2

include { MUTECT2_PON_NORMAL } from './MUTECT2_PON_NORMAL'
include { MUTECT2_PON_MERGE } from './MUTECT2_PON_MERGE'
include { MUTECT2_PON_CREATE } from './MUTECT2_PON_CREATE'

// --mutect2_cohort_pon build: a process-matched panel of normals from every normal BAM of the run (paired
// normals and normal-only pool rows), over the same Mutect2 shards the tumours are called on.
// Emits a VALUE channel [pon.vcf.gz, pon.vcf.gz.tbi], so every MUTECT2_CALLING task waits for it.
workflow MUTECT2_COHORT_PON {
    take:
    normal          // [patient, meta, bam, bai], status normal
    chunk           // 1..20
    interval_dir    // Mutect2 shard directory
    pon_intervals   // padded targets BED (--targets) or the shard directory

    main:
    MUTECT2_PON_NORMAL(normal, chunk, interval_dir)
    MUTECT2_PON_MERGE(MUTECT2_PON_NORMAL.out.shard_vcf.groupTuple(by: 0, size: 20))
    MUTECT2_PON_CREATE(
        MUTECT2_PON_MERGE.out.normal_vcf.map { id, vcf, tbi -> [vcf, tbi] }.flatten().collect(),
        pon_intervals
    )
    pon = MUTECT2_PON_CREATE.out.pon.collect()
    pon.subscribe { files ->
        log.info "Mutect2: tumours are called against the cohort PON built in this run: " +
                 "${params.mutect2_cohort_pon_dir}/${files[0].name}"
    }

    emit:
    pon
}

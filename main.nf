// Author: Ting Yang, George Wu
// Lab: The Alexandrov Lab @ UCSD
// Date: 2026.9
// Version: 4.0

nextflow.enable.dsl=2

params.database_path = null
params.sample = "sample.csv"
params.genome = ""
params.tool = ""
params.cleanup = false 
params.type = "genome"
// Clinical-data options (feature/clinical-bam)
params.outdir = null              // results root; default $projectDir/RESULTS (previous behaviour)
params.targets = null             // capture BED for --type exome|panel: Mutect2 shards, Strelka regions, mosdepth, CNVkit, BQSR
params.target_padding = 100       // bp added on each side of --targets
params.contig_rename = "auto"     // auto: rename b37 <-> chr contigs in BAM headers when lengths match --ref; none: never
params.cnvkit_reference = "auto"  // auto: pooled from the run's normals, flat if the sheet has none; or pooled | flat
// >>> msi
params.msi_threshold_paired = 10      // MSI-H when % unstable sites >= this, paired (msisensor-pro msi; MSK-IMPACT MSIsensor cut-off)
params.msi_threshold_tumor_only = 20  // MSI-H cut-off for tumour-only (msisensor-pro pro vs baseline); calibrate per assay
params.msi_min_sites = 50             // fewer sites evaluated -> qc LOW_SITES (MSIsensor-pro paper: >= 50 sites)
params.msi_min_coverage = 20          // per-site read depth for a site to count (msisensor-pro -c; 20 for WES/panels)
params.msi_sites = null               // prebuilt msisensor-pro scan list (skips MSI_SCAN)
params.msi_baseline = null            // prebuilt tumour-only baseline (MSI/baseline/msi_baseline.tsv of an earlier run, same assay)
params.msi_baseline_min_normals = 10  // baseline: sites must be covered in this many normals (capped at the number of normals)
// <<< msi
// >>> mutect2_pon
params.mutect2_cohort_pon = null          // build: Mutect2 PON from this run's normal BAMs, used for its own tumours
params.mutect2_cohort_pon_name = "cohort" // published as Mutect2_PON/<name>.pon.vcf.gz(.tbi)
params.mutect2_cohort_pon_warn_below = 10 // warn when the PON is built from fewer normals (error below 2)
// --mutect2_pon <vcf.gz> (declared further down, default = database PON) replaces the database PON everywhere.
// <<< mutect2_pon
// >>> umi
// UMI-aware duplicates (fgbio molecular consensus) for --first_step mapping. Off by default. When on, the
// consensus BAM replaces BWA_MEM + MarkDuplicates and goes straight to BQSR and the callers.
// Exactly one UMI source: inline (--umi_read_structure with M segments), a UMI FASTQ (--umi_fastq_column,
// sheet column fastq_umi) or the read name (--umi_from_header).
params.umi = false
params.umi_read_structure = null        // fgbio read structures, one per FASTQ, e.g. "3M2S+T 3M2S+T"
params.umi_fastq_column = false         // UMI read in sheet column fastq_umi (default structures "+T +T +M")
params.umi_from_header = false          // UMI = last ':' field of the read name (default structures "+T +T")
params.umi_group_strategy = null        // identity | edit | adjacency | paired; default adjacency, paired with --umi_duplex
params.umi_edits = 1                    // UMI mismatches allowed when grouping
params.umi_min_map_q = 1                // reads below this MAPQ are not grouped
params.umi_min_reads = 1                // reads per family for a consensus; duplex accepts up to 3 values, e.g. "2 1 1"
params.umi_duplex = false               // CallDuplexConsensusReads (needs a UMI on both ends of the molecule)
params.umi_min_base_quality = 30        // FilterConsensusReads: consensus bases below this become N
params.umi_max_base_error_rate = 0.2    // FilterConsensusReads: mask bases where more reads disagree with consensus
params.umi_max_no_call_fraction = 0.2   // FilterConsensusReads: drop consensus reads with more Ns than this
if (params.umi) {
    if (params.first_step != 'mapping') {
        error "--umi needs --first_step mapping: UMIs are read from FASTQs (processed BAMs are already deduplicated)."
    }
    def seg = /(\d+|\+)[TBMSC]/
    def structs = params.umi_read_structure ? params.umi_read_structure.toString().trim().split(/\s+/) as List : []
    def bad = structs.findAll { !(it ==~ /(${seg})+/) }
    if (bad) error "--umi_read_structure: not an fgbio read structure: ${bad.join(' ')} (e.g. 3M2S+T 3M2S+T)"
    def nM = { String s -> (s =~ /(\d+|\+)M/).size() }
    def inline = structs && structs.any { nM(it) > 0 } && !params.umi_fastq_column
    def sources = [inline ? 'inline' : null, params.umi_fastq_column ? 'fastq' : null, params.umi_from_header ? 'header' : null].findAll()
    if (sources.size() != 1) {
        error "--umi needs exactly one UMI source: --umi_read_structure with M segments (inline UMI), " +
              "--umi_fastq_column (sheet column fastq_umi) or --umi_from_header. Got: ${sources ?: 'none'}"
    }
    def src = sources[0]
    if (src == 'fastq') {
        structs = structs ?: ['+T', '+T', '+M']
        if (structs.size() != 3 || nM(structs[2]) == 0) {
            error "--umi_fastq_column: give three read structures (R1 R2 UMI-read) with M in the third, e.g. '+T +T +M'"
        }
        def hdr = file(params.sample).exists() ? file(params.sample).readLines()[0].split(',').collect { it.trim() } : []
        if (!('fastq_umi' in hdr)) error "--umi_fastq_column: sample sheet ${params.sample} has no fastq_umi column"
    } else if (src == 'header') {
        structs = structs ?: ['+T', '+T']
        if (structs.size() != 2 || structs.any { nM(it) > 0 }) {
            error "--umi_from_header: --umi_read_structure may only trim (two structures, no M), e.g. '5S+T 5S+T'"
        }
    } else if (structs.size() != 2) {
        error "--umi_read_structure needs one structure per FASTQ (two), e.g. '3M2S+T 3M2S+T'"
    }
    def strategy = params.umi_group_strategy ?: (params.umi_duplex ? 'paired' : 'adjacency')
    if (!(strategy in ['identity', 'edit', 'adjacency', 'paired'])) {
        error "--umi_group_strategy must be identity, edit, adjacency or paired"
    }
    def minReads = params.umi_min_reads.toString().trim().split(/\s+/) as List
    if (!minReads.every { it ==~ /\d+/ } || minReads.size() > (params.umi_duplex ? 3 : 1)) {
        error "--umi_min_reads: one integer (simplex) or up to three for --umi_duplex, e.g. '2 1 1'"
    }
    if (params.umi_duplex) {
        if (strategy != 'paired') error "--umi_duplex needs --umi_group_strategy paired (the default with --umi_duplex)"
        if (src == 'inline' && !structs.every { nM(it) > 0 }) {
            error "--umi_duplex: duplex consensus needs a UMI on both reads (M in both read structures)"
        }
        if (src == 'fastq' && structs.sum { nM(it) } < 2) {
            error "--umi_duplex: the UMI read structures carry one UMI; duplex needs two (e.g. '+T +T 8M8M')"
        }
        if (src == 'header') log.warn "--umi_duplex with --umi_from_header: read names must carry two UMIs (e.g. ACGT+TTGA)"
    } else if (strategy == 'paired') {
        log.warn "--umi_group_strategy paired groups both strands of a molecule; consensus is still simplex (add --umi_duplex for duplex)"
    }
    params.umi_source = src
    params.umi_read_structures_resolved = structs.join(' ')
    params.umi_group_strategy_resolved = strategy
    log.info "UMI consensus: source ${src}, read structures '${structs.join(' ')}', ${strategy} grouping " +
             "(edits ${params.umi_edits}), ${params.umi_duplex ? 'duplex' : 'simplex'} consensus, min reads ${minReads.join(' ')}; MarkDuplicates is skipped"
}
// <<< umi
params.results_dir = params.outdir ?: "$projectDir/RESULTS"
// >>> tmb
// Tumour mutational burden (--tool tmb): PASS consensus SNVs+indels with tumour VAF >= tmb_min_vaf inside the
// callable territory (tmb_regions with tumour depth >= tmb_min_depth and, for pairs, normal depth >= tmb_min_normal_depth).
params.tmb_regions = null                 // BED to count over; default the cleaned --targets (pass a coding-exon BED for coding TMB)
params.tmb_min_depth = 50                 // tumour depth for a base to be callable
params.tmb_min_normal_depth = 20          // matched-normal depth (pairs only)
params.tmb_min_vaf = 0.05                 // tumour VAF, recomputed from the tumour BAM
params.tmb_min_mapq = 20                  // reads counted for depth and VAF
params.tmb_min_baseq = 20                 // bases counted for VAF
params.tmb_scenario = 4                   // paired consensus: POST scenario 4, 5 or 6 (see README)
params.tmb_coding_only = false            // count only coding consequences; needs a CSQ/ANN-annotated VCF, else ignored with a warning
params.tmb_germline_calls = null          // tumour-only: directory of PureCN <patient>_<sample>_variants.csv (optional)
params.tmb_min_somatic_posterior = 0.5    // tumour-only + PureCN: drop variants with POSTERIOR.SOMATIC below this
params.tmb_dir = "${params.results_dir}/TMB"
params.tmb_env = "$projectDir/yml/tmb.yml"
if (!(params.tmb_scenario.toString() in ['4', '5', '6'])) {
    error "--tmb_scenario must be 4, 5 or 6 (scenario 3 output is not PASS-gated)."
}
// <<< tmb
params.log_file = params.outdir ? "${params.outdir}/pipeline.log" : "$projectDir/pipeline.log"
// >>> vep
// --tool vep: Ensembl VEP (offline cache) + vcf2maf on the paired and tumour-only consensus VCFs (see README)
params.vep_cache = null            // VEP cache root (--dir), containing homo_sapiens/<version>_<assembly>/
params.vep_cache_version = 112     // must match the ensembl-vep release in yml/vep.yml
params.vep_assembly = null         // default: --genome (GRCh38 | GRCh37)
params.vep_fasta = null            // default: --ref (same contig names as the VCFs)
params.vep_scenario = "scenario_5" // paired consensus to annotate: scenario_3 | scenario_4 | scenario_5 | scenario_6
params.vep_env = "$projectDir/yml/vep.yml"
params.vep_dir = "${params.results_dir}/VEP"
params.vep_maf_dir = "${params.results_dir}/MAF"
if (params.tool && params.tool.toString().split(',').collect { it.trim() }.contains('vep')) {
    def vepAssembly = params.vep_assembly ?: params.genome
    if (!(vepAssembly in ['GRCh38', 'GRCh37'])) {
        error "--tool vep annotates human calls only: --vep_assembly/--genome must be GRCh38 or GRCh37 (got '${vepAssembly}')."
    }
    if (!(params.vep_scenario in ['scenario_3', 'scenario_4', 'scenario_5', 'scenario_6'])) {
        error "--vep_scenario must be scenario_3, scenario_4, scenario_5 or scenario_6 (the MERGE_HEADER standard VCFs)."
    }
    if (!params.vep_cache) {
        error "--tool vep needs --vep_cache <dir> holding homo_sapiens/${params.vep_cache_version}_${vepAssembly}/ (an indexed Ensembl VEP cache)."
    }
    def vepCacheDir = file("${params.vep_cache}/homo_sapiens/${params.vep_cache_version}_${vepAssembly}")
    if (!vepCacheDir.isDirectory()) {
        error "VEP cache not found: ${vepCacheDir} (check --vep_cache, --vep_cache_version, --vep_assembly)."
    }
}
// <<< vep
new File(params.log_file).absoluteFile.parentFile.mkdirs()

// Initialize the log file with a header
new File(params.log_file).text = """
==============================================
Pipeline Run - ${new Date()}
==============================================

"""
// Validate database being set
if (!params.database_path) {
    error "ERROR: --database_path is required. Please provide the path to your SomaticBunny_database directory.\n" +
          "Usage: nextflow run main.nf --database_path /path/to/SomaticBunny_database"
}

params.genomes = [
    'GRCh38': [
        ref: "${params.database_path}/GRCh38_ref/GRCh38.d1.vd1.fa",
        bed: "${params.database_path}/GRCh38_ref/hg38_chr.bed.gz",
        conpair_marker: "${params.database_path}/Databases/GRCh38/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed",
        conpair_marker_txt: "${params.database_path}/Databases/GRCh38/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt",
        mosdepth_bed: "${params.database_path}/Databases/GRCh38/GRCh38_exome.bed",
        muse_vcf: "${params.database_path}/Databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        recal_knownsite1: "${params.database_path}/Databases/GRCh38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf",
        recal_knownsite2: "${params.database_path}/Databases/GRCh38/Homo_sapiens_assembly38.known_indels.vcf.gz",
        recal_interval_wes: "${params.database_path}/Databases/GRCh38/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38_canonical.targets.interval_list",
        sage: [
            sage_ref_dir: "${params.database_path}/SAGE/v5_34/ref/38",
            sage_hotspots: "${params.database_path}/SAGE/v5_34/ref/38/variants/KnownHotspots.somatic.38.vcf.gz",
            sage_panel_bed: "${params.database_path}/SAGE/v5_34/ref/38/variants/ActionableCodingPanel.38.bed.gz",
            sage_high_confidence_bed: "${params.database_path}/SAGE/v5_34/ref/38/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz",
            sage_ref_genome_version: "38" 
        ],
        ASCAT: [
            ascat_ref_genome_version: "hg38",
            ascat_allele_prefix: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/Alleles/G1000_alleles_hg38_chr",
            ascat_loci_prefix: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/Loci/G1000_loci_hg38_chr",
            ascat_GCcontentfile: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/GC_Correction/GC_G1000_hg38.txt",
            ascat_replictimingfile: "${params.database_path}/Databases/GRCh38/ASCAT/WGS/hg38/RT_Correction/RT_G1000_hg38.txt",
			ascat_allele_prefix_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/Alleles/G1000_alleles_hg38_chr",
			ascat_loci_prefix_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/Loci/G1000_loci_hg38_chr",
            ascat_GCcontentfile_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/GC_Correction/GC_G1000_hg38.txt",
            ascat_replictimingfile_wes: "${params.database_path}/Databases/GRCh38/ASCAT/WES/hg38/RT_Correction/RT_G1000_hg38.txt"
        ],
        database_subdir: "GRCh38",
        mutect2_pon: "MuTect2.PON.5210.vcf.gz",
        mutect2_pon_wes: "MuTect2.PON.5210.vcf.gz",
        mutect2_germline: "af-only-gnomad.hg38_no_alt.vcf.gz",
        mutect2_interval_dir: "GRCh38_interval_list_20",
        delly_excl: "${params.database_path}/Databases/GRCh38/Delly/human.hg38.excl.tsv",
        tools: ["fastqc", "bwa_mem", "mkdup", "recalibrate", "sage", "strelka", "muse", "mutect2", "ascat", "delly", "cnvkit", "mosdepth", "conpair", "manta"]
    ],
    'GRCh37': [
        ref: "${params.database_path}/GRCh37_ref/GRCh37.fa", 
        bed: "${params.database_path}/GRCh37_ref/GRCh37_chr.bed.gz",
        conpair_marker: "${params.database_path}/Databases/GRCh37/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed",
        conpair_marker_txt: "${params.database_path}/Databases/GRCh37/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt",
        mosdepth_bed: "${params.database_path}/Databases/GRCh37/GRCh37_exome.bed",
        muse_vcf: "${params.database_path}/Databases/GRCh37/Homo_sapiens_assembly19.dbsnp.vcf",
        recal_knownsite1: "${params.database_path}/Databases/GRCh37/Homo_sapiens_assembly19.dbsnp.vcf.gz",
        recal_knownsite2: "${params.database_path}/Databases/GRCh37/Homo_sapiens_assembly19.known_indels.vcf.gz",
        recal_interval_wes: "${params.database_path}/Databases/GRCh37/intervals_Broad.human.exome.b37.interval_list",
        sage: [
            sage_ref_dir: "${params.database_path}/SAGE/v5_34/ref/37",
            sage_hotspots: "${params.database_path}/SAGE/v5_34/ref/37/dna/variants/KnownHotspots.somatic.37.vcf.gz",
            sage_panel_bed: "${params.database_path}/SAGE/v5_34/ref/37/dna/variants/ActionableCodingPanel.37.bed.gz",
            sage_high_confidence_bed: "${params.database_path}/SAGE/v5_34/ref/37/dna/variants/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz",
            sage_ref_genome_version: "37" 
        ],
        ASCAT: [
            ascat_ref_genome_version: "hg19",
            ascat_allele_prefix: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/Alleles/G1000_alleles_hg19_chr",
            ascat_loci_prefix: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/Loci/G1000_loci_hg19_chr",
            ascat_GCcontentfile: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/GC_Correction/GC_G1000_hg19.txt",
            ascat_replictimingfile: "${params.database_path}/Databases/GRCh37/ASCAT/WGS/hg19/RT_Correction/RT_G1000_hg19.txt",
            ascat_allele_prefix_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/Alleles/G1000_alleles_hg19_chr",
            ascat_loci_prefix_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/Loci/G1000_loci_hg19_chr",
            ascat_GCcontentfile_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/GC_Correction/GC_G1000_hg19.txt",
            ascat_replictimingfile_wes: "${params.database_path}/Databases/GRCh37/ASCAT/WES/hg19/RT_Correction/RT_G1000_hg19.txt"
        ],
        database_subdir: "GRCh37",
        mutect2_pon: "Mutect2-WGS-panel-b37.vcf",
        mutect2_pon_wes: "Mutect2-exome-panel_b37.vcf",
        mutect2_germline: "af-only-gnomad.raw.sites.grch37.vcf.gz",
        mutect2_interval_dir: "GRCh37_interval_list_20",
        delly_excl: "${params.database_path}/Databases/GRCh37/Delly/human.hg19.excl.tsv",
        tools: ["fastqc", "bwa_mem", "mkdup", "recalibrate", "sage", "strelka", "muse", "mutect2", "ascat", "delly", "cnvkit", "mosdepth", "conpair", "manta"]
    ],
    'mm39': [
        ref: "${params.database_path}/mm39_ref/mm39.fa",
        bed: "${params.database_path}/mm39_ref/mm39_chr.bed.gz",
        database_subdir: "mm39",
        mosdepth_bed: "${params.database_path}/Databases/mm39/mm39_exome.bed",
        recal_interval_wes: "${params.database_path}/Databases/mm39/mm39_exome.interval_list",
        muse_vcf: "${params.database_path}/Databases/mm39/af_only_mgp_mm39_unique.vcf.gz",
        mutect2_pon: "PoN.mm39.vcf.gz",
        mutect2_pon_wes: "PoN.mm39.vcf.gz",
        mutect2_germline: "af-only-mgp.mm39.vcf.gz",
        mutect2_interval_dir: "mm39_interval_list_20",
        tools: ["fastqc", "bwa_mem", "mkdup", "strelka", "muse", "mutect2", "cnvkit", "mosdepth", "manta"]
    ],
    'RN7': [
        ref: "${params.database_path}/RN7_ref/rn7.fa",
        bed: "${params.database_path}/RN7_ref/rn7_chr.bed.gz",
        database_subdir: "RN7",
        mosdepth_bed: "${params.database_path}/Databases/RN7/rn7_exome.bed",
        recal_interval_wes: "${params.database_path}/Databases/RN7/rn7_exome.interval_list",
        muse_vcf: "",
        mutect2_pon: "",
        mutect2_pon_wes: "",
        mutect2_germline: "af-only_rn7.vcf.gz",
        mutect2_interval_dir: "RN7_interval_list_20",
        tools: ["fastqc", "bwa_mem", "mkdup", "strelka", "muse", "mutect2", "cnvkit", "mosdepth", "manta"]
    ]
]

// Custom ref/bed is REQUIRED when starting from markdup, recalibration, or variant_calling
if (params.first_step in ['markdup', 'recalibration', 'variant_calling']) {
    def missing = []
    if (!params.ref)      missing << '--ref'
    if (!params.ref_fai)  missing << '--ref_fai'
    if (!params.ref_dict) missing << '--ref_dict'
    if (!params.bed)      missing << '--bed'
    if (!params.bed_tbi)  missing << '--bed_tbi'
    if (missing) {
        error "When --first_step is '${params.first_step}', you must provide all custom reference files: " +
              "--ref, --ref_fai, --ref_dict, --bed, --bed_tbi. Missing: ${missing.join(', ')}"
    }
}

// If any custom ref/bed provided, validate they come as a complete set
if (params.ref || params.bed || params.bed_tbi || params.ref_fai || params.ref_dict) {
    def missing = []
    if (!params.ref)      missing << '--ref'
    if (!params.ref_fai)  missing << '--ref_fai'
    if (!params.ref_dict) missing << '--ref_dict'
    if (!params.bed)      missing << '--bed'
    if (!params.bed_tbi)  missing << '--bed_tbi'
    if (missing) {
        error "When providing custom reference files, all five must be specified: " +
              "--ref, --ref_fai, --ref_dict, --bed, --bed_tbi. Missing: ${missing.join(', ')}"
    }
    if (!file(params.ref).exists())      error "Custom ref file not found: ${params.ref}"
    if (!file(params.ref_fai).exists())  error "Custom ref_fai file not found: ${params.ref_fai}"
    if (!file(params.ref_dict).exists()) error "Custom ref_dict file not found: ${params.ref_dict}"
    if (!file(params.bed).exists())      error "Custom bed file not found: ${params.bed}"
    if (!file(params.bed_tbi).exists())  error "Custom bed_tbi file not found: ${params.bed_tbi}"

    log.info """
    ==============================================
    Using CUSTOM reference files:
      ref      : ${params.ref}
      ref_fai  : ${params.ref_fai}
      ref_dict : ${params.ref_dict}
      bed      : ${params.bed}
      bed_tbi  : ${params.bed_tbi}
    ==============================================
    """
}


// Validate genome selection
if (!params.genomes.containsKey(params.genome)) {
    error "Invalid genome: ${params.genome}. Available options are: ${params.genomes.keySet().join(', ')}"
}

// Validate analysis type. "panel" is targeted like "exome" but requires --targets and skips genome-wide CN/SV callers.
if (!(params.type in ['genome', 'exome', 'panel'])) {
    error "Invalid --type '${params.type}'. Use genome, exome or panel."
}
params.targeted = params.type in ['exome', 'panel']
if (params.type == 'panel' && !params.targets) {
    error "--type panel needs --targets <capture BED> (same contig names as --ref)."
}
if (params.targets && !params.targeted) {
    error "--targets is only used with --type exome or --type panel."
}
if (params.targets && !file(params.targets).exists()) {
    error "--targets file not found: ${params.targets}"
}
if (params.type == 'panel') {
    def tools = params.tool ? params.tool.toString().split(',').collect { it.trim() } : []
    def unsuitable = tools.findAll { it in ['ascat', 'delly'] }
    if (unsuitable) {
        error "--tool ${unsuitable.join(',')} is not supported with --type panel: ASCAT needs genome-wide allele frequencies and " +
              "Delly needs genome-wide coverage. Use cnvkit (and manta only if the panel baits fusion introns)."
    }
}
if (!(params.contig_rename in ['auto', 'none'])) {
    error "--contig_rename must be auto or none."
}
if (!(params.cnvkit_reference in ['auto', 'pooled', 'flat'])) {
    error "--cnvkit_reference must be auto, pooled or flat."
}

// Check if user inputs custom ref or not
params.custom_ref = (params.ref && params.ref != params.genomes[params.genome].ref)

if (params.custom_ref) {
    log.warn """
    WARNING: You are using a custom reference genome.
    Ensure your reference uses the same contig naming convention 
    (e.g., chr1 vs 1) as the pipeline's database files for genome '${params.genome}'.
    Mismatched contig names will cause tools like Mutect2, MuSE, BQSR, 
    SAGE, Conpair, and Delly to fail.
    """
}

// Set genome-specific parameters
params.ref = params.ref ?: params.genomes[params.genome].ref
params.ref_fai = params.ref_fai ?: "${params.ref}.fai"
params.ref_dict = params.ref_dict ?: params.ref.replaceAll(/\.fa(sta)?$/, '.dict')
params.bed = params.bed ?: params.genomes[params.genome].bed
params.bed_tbi = params.bed_tbi ?: "${params.bed}.tbi"
params.conpair_marker = params.conpair_marker ?: params.genomes[params.genome].conpair_marker
params.conpair_marker_txt = params.conpair_marker_txt ?: params.genomes[params.genome].conpair_marker_txt

// Clean --targets once (gunzip, drop track/browser/# lines, keep 3 columns) and check its contigs against --ref,
// so mosdepth, CNVkit and BQSR read the same clean BED that PREPARE_TARGETS pads for Mutect2 and Strelka.
def targetsClean = null
if (params.targets) {
    def faiContigs = file(params.ref_fai).readLines().collect { it.split('\t')[0] } as Set
    def src = file(params.targets)
    def rawLines = src.name.endsWith('.gz') ?
        new java.util.zip.GZIPInputStream(src.newInputStream()).newReader().readLines() : src.readLines()
    def rows = rawLines.findAll { it.trim() && !(it ==~ /^(track|browser|#).*/) }
        .collect { it.split('\t') }
        .findAll { it.size() >= 3 }
    if (!rows) {
        error "--targets has no BED intervals: ${params.targets}"
    }
    def missing = rows.collect { it[0] }.unique().findAll { !(it in faiContigs) }
    if (missing) {
        error "--targets uses contigs that are not in --ref: ${missing.take(5).join(' ')} " +
              "(chr prefix mismatch? the BED must use the same contig names as the reference)"
    }
    def out = new File("${params.results_dir}/custom_intervals/targets.clean.bed")
    out.parentFile.mkdirs()
    // Rewrite only when the content changes: a new timestamp on every launch invalidated -resume for every step
    // downstream of the targets.
    def cleanText = rows.collect { "${it[0]}\t${it[1]}\t${it[2]}" }.join('\n') + '\n'
    if (!out.exists() || out.text != cleanText) {
        out.text = cleanText
    }
    targetsClean = out.absolutePath
}
params.targets_clean = targetsClean
params.mosdepth_bed = params.mosdepth_bed ?: (params.targets_clean ?: params.genomes[params.genome].mosdepth_bed)
// >>> tmb
if (params.tool && params.tool.toString().split(',').collect { it.trim() }.contains('tmb')) {
    if (!(params.tmb_regions ?: params.targets_clean)) {
        error "--tool tmb needs --tmb_regions <bed> (or --targets): TMB is mutations per callable Mb of a defined territory."
    }
    if (params.tmb_regions && !file(params.tmb_regions).exists()) {
        error "--tmb_regions file not found: ${params.tmb_regions}"
    }
    if (params.tmb_germline_calls && !file(params.tmb_germline_calls).isDirectory()) {
        error "--tmb_germline_calls must be a directory holding PureCN <patient>_<sample>_variants.csv files"
    }
    if (!(params.genome in ['GRCh38', 'GRCh37'])) {
        error "--tool tmb counts the POST consensus, which runs only for --genome GRCh38 or GRCh37."
    }
}
// <<< tmb
params.muse_vcf = params.muse_vcf ?: params.genomes[params.genome].muse_vcf
params.recal_knownsite1 = params.recal_knownsite1 ?: params.genomes[params.genome].recal_knownsite1
params.recal_knownsite2 = params.recal_knownsite2 ?: params.genomes[params.genome].recal_knownsite2
params.recal_interval_wes = params.recal_interval_wes ?: (params.targets_clean ?: params.genomes[params.genome].recal_interval_wes)
params.delly_excl = params.delly_excl ?: params.genomes[params.genome].delly_excl

// Set SAGE-specific parameters
params.SAGE_ref_dir = params.SAGE_ref_dir ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_ref_dir : null)
params.SAGE_hotspots = params.SAGE_hotspots ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_hotspots : null)
params.SAGE_panel_bed = params.SAGE_panel_bed ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_panel_bed : null)
params.SAGE_high_confidence_bed = params.SAGE_high_confidence_bed ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_high_confidence_bed : null)
params.SAGE_ref_genome_version = params.SAGE_ref_genome_version ?: (params.genomes[params.genome].containsKey('sage') ? params.genomes[params.genome].sage.sage_ref_genome_version : null)

// Set ASCAT-specific parameters
params.ascat_ref_genome_version = params.ascat_ref_genome_version ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_ref_genome_version : null)
params.ascat_allele_prefix = params.ascat_allele_prefix ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_allele_prefix : null)
params.ascat_loci_prefix = params.ascat_loci_prefix ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_loci_prefix : null)
params.ascat_GCcontentfile = params.ascat_GCcontentfile ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_GCcontentfile : null)
params.ascat_replictimingfile = params.ascat_replictimingfile ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_replictimingfile : null)

params.ascat_allele_prefix_wes = params.ascat_allele_prefix_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_allele_prefix_wes : null)
params.ascat_loci_prefix_wes = params.ascat_loci_prefix_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_loci_prefix_wes : null)
params.ascat_GCcontentfile_wes = params.ascat_GCcontentfile_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_GCcontentfile_wes : null)
params.ascat_replictimingfile_wes = params.ascat_replictimingfile_wes ?: (params.genomes[params.genome].containsKey('ASCAT') ? params.genomes[params.genome].ASCAT.ascat_replictimingfile_wes : null)


params.database_dir = "${params.database_path}/Databases"
params.genome_database_dir = "${params.database_dir}/${params.genomes[params.genome].database_subdir}"

// Genome-specific file paths for Mutect2
params.mutect2_pon = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_pon}"
params.mutect2_interval_dir = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_interval_dir}"
params.mutect2_germline = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_germline}"
// >>> mutect2_pon
// Which PON MUTECT2_CALLING gets (paired and tumour-only): the database one, --mutect2_pon <vcf> (a command-line
// value wins over the assignment above), or the cohort PON built in this run. It is passed to the process as a
// staged [vcf, index] input, so the cohort PON is an ordinary channel dependency.
params.mutect2_cohort_pon_dir = "${params.results_dir}/Mutect2_PON"
def mutect2PonDefault = "${params.genome_database_dir}/${params.genomes[params.genome].mutect2_pon}".toString()
def mutect2PonUser = params.mutect2_pon && params.mutect2_pon.toString() != mutect2PonDefault
if (params.mutect2_cohort_pon && params.mutect2_cohort_pon != 'build') {
    error "--mutect2_cohort_pon accepts only 'build' (got '${params.mutect2_cohort_pon}'). " +
          "To use a PON built earlier, pass it with --mutect2_pon <pon.vcf.gz>."
}
if (params.mutect2_cohort_pon && mutect2PonUser) {
    error "Use either --mutect2_cohort_pon build or --mutect2_pon <vcf>, not both."
}
def mutect2PonFiles = []
if (!params.mutect2_cohort_pon) {
    def ponPath = params.mutect2_pon ? params.mutect2_pon.toString() : ''
    if (ponPath && !ponPath.endsWith('/')) {      // the RN7 entry has no PON: '' -> no -pon
        def ponIdx = ponPath.endsWith('.gz') ? "${ponPath}.tbi" : "${ponPath}.idx"
        if (mutect2PonUser) {
            if (!file(ponPath).exists()) error "--mutect2_pon not found: ${ponPath}"
            if (!file(ponIdx).exists()) error "--mutect2_pon needs its index next to it: ${ponIdx}"
        }
        mutect2PonFiles = file(ponIdx).exists() ? [ponPath, ponIdx] : [ponPath]
    }
}
params.mutect2_pon_resolved = mutect2PonFiles
def mutect2PonMsg = params.mutect2_cohort_pon ?
    "Mutect2 PON      : cohort PON built in this run -> ${params.mutect2_cohort_pon_dir}/${params.mutect2_cohort_pon_name}.pon.vcf.gz" :
    "Mutect2 PON      : ${mutect2PonFiles ? mutect2PonFiles[0] : 'none'} (${mutect2PonUser ? '--mutect2_pon' : 'database'})"
log.info mutect2PonMsg
new File(params.log_file).append(mutect2PonMsg + "\n")
// <<< mutect2_pon

params.indel_vaf_script =
    params.indel_vaf_script ?: "${projectDir}/INDEL_VAF.py"

// Define which tools are available for the selected genome
def available_tools = params.genomes[params.genome].tools

// Create a function to write to the log file
def writeToLog(logFile, message) {
    new File(logFile).append(message + "\n")
}

// Generate comprehensive configuration log
def configInfo = """
Pipeline Configuration:
==============================================
Run ID         : ${workflow.runName}
Command Line   : ${workflow.commandLine}
Profile        : ${workflow.profile}
Container      : ${workflow.container}
Run as User    : ${workflow.userName}
Launch Dir     : ${workflow.launchDir}
Work Dir       : ${workflow.workDir}
Project Dir    : ${projectDir}
==============================================
Analysis Parameters:
==============================================
Analysis Type  : ${params.type}
Targets        : ${params.targets ?: 'none'} (padding ${params.target_padding} bp)
Contig rename  : ${params.contig_rename}
Results        : ${params.results_dir}
Genome         : ${params.genome}
==============================================
Reference Files:
==============================================
ref            : ${params.ref}
ref_fai        : ${params.ref_fai}
bed            : ${params.bed}
bed_tbi        : ${params.bed_tbi}
==============================================
Avaliable Tools:
==============================================
${available_tools.join(', ')}
==============================================
"""

params.bam_dir="${params.results_dir}/BAM"
params.report_dir="${params.results_dir}/REPORT"
params.manta_bed="${params.database_path}/GRCh38_ref/manta.bed.gz"
params.mkdup_temp_dir="$projectDir/mkdup_tmp"
params.FASTQC_dir="${params.results_dir}/FASTQC"
params.mkdup_dir="${params.results_dir}/MKDUP"
params.recal_dir="${params.results_dir}/RECALIBRATE"
params.SAGE_jar="${params.database_path}/SAGE/sage_v4.2.jar"
params.SAGE_dir="${params.results_dir}/SAGE"
params.strelka_dir="${params.results_dir}/STRELKA"
params.muse_dir="${params.results_dir}/MuSE"
params.MUTECT2_dir="${params.results_dir}/Mutect2"
params.mosdepth_dir="${params.results_dir}/mosdepth"
params.conpair_dir="${params.results_dir}/Conpair"
params.cnvkit_dir="${params.results_dir}/CNVkit"
params.Delly_dir="${params.results_dir}/Delly"
params.ascat_dir="${params.results_dir}/ASCAT"
params.conpair="${params.database_path}/Conpair-0.2"
params.manta_dir="${params.results_dir}/MANTA"
params.post_dir="${params.results_dir}/POST"
params.intervals_dir="${params.results_dir}/custom_intervals"
params.VAF_dir="${params.results_dir}/POST/VAF"
// >>> msi
params.msi_dir = "${params.results_dir}/MSI"
params.msi_env = "$projectDir/yml/msisensor_pro.yml"
// <<< msi

params.bwamem2_env = "$projectDir/yml/bwamem2.yml"
params.mkdup_env = "$projectDir/yml/picard.yml"
params.conpair_env = "$projectDir/yml/conpair.yml"
params.samtools_env = "$projectDir/yml/samtools.yml"
params.strelka_env = "$projectDir/yml/strelka_env.yml"
params.mosdepth_env = "$projectDir/yml/mosdepth_env.yml"
params.rename_bam_env = "$projectDir/yml/rename_bam.yml"
params.fastqc_env = "$projectDir/yml/fastqc_env.yml"
params.cnvkit_env = "$projectDir/yml/cnvkit.yml"
params.delly_env = "$projectDir/yml/delly.yml"
params.ascat_env = "$projectDir/yml/ascat.yml"
params.alleleCounter_env = "$projectDir/yml/alleleCounter.yml"
params.manta_env = "$projectDir/yml/manta.yml"
params.dkfz_env = "$projectDir/yml/dkfz.yml"
params.SAGE_java_env = "$projectDir/yml/sage_java.yml"
params.muse_env = "$projectDir/yml/muse.yml"
params.picard_merge_env = "$projectDir/yml/picard.yml"
params.gatk_env = "${projectDir}/yml/gatk_runtime.yml"
params.POST_env = "${projectDir}/yml/POST.yml"
params.POST_rename_env = "${projectDir}/yml/POST_rename.yml"
params.gatk = "${params.database_path}/gatk-4.6.2.0/gatk"
// >>> umi
params.umi_dir = "${params.results_dir}/UMI"
params.umi_env = "${projectDir}/yml/umi.yml"
// <<< umi

// Write configuration to log file
writeToLog(params.log_file, configInfo)

// Also display in terminal
log.info configInfo

if (params.cleanup) {
    def cleanupWarning = """
==============================================
WARNING: --cleanup is ENABLED
==============================================
Intermediate files will be permanently deleted during the run:
  - BWA_MEM work directories and raw BAMs
  - MKDUP BAMs and work directories
  - Recalibration chunk/merged BAMs and work directories

You will NOT be able to use -resume to reuse these steps.
Only final outputs in RESULTS/ will be kept.
==============================================
"""
    log.warn cleanupWarning
    writeToLog(params.log_file, cleanupWarning)
}

include { FASTQC } from './Modules/FASTQC'
include { BWA_MEM } from './Modules/BWA_MEM'
include { CHECK_BAM_BWA } from './Modules/CHECK_BAM_BWA'
include { CLEANUP_BWA } from './Modules/CLEANUP_BWA'
include { COMBINE_REPORTS_BWA } from './Modules/COMBINE_REPORTS_BWA'
include { RECALIBRATE_BaseRecal } from './Modules/RECALIBRATE_BaseRecal.nf'
include { RECALIBRATE_BaseRecal_exome } from './Modules/RECALIBRATE_BaseRecal_exome.nf'
include { RECALIBRATE_MergeReport } from './Modules/RECALIBRATE_MergeReport.nf'
include { RECALIBRATE_BQSR } from './Modules/RECALIBRATE_BQSR.nf'
include { RECALIBRATE_BQSR_exome } from './Modules/RECALIBRATE_BQSR_exome.nf'
include { RECALIBRATE_MergeBam } from './Modules/RECALIBRATE_MergeBam.nf'
include { RECALIBRATE_SortBam } from './Modules/RECALIBRATE_SortBam.nf'
include { CHECK_BAM_RECAL } from './Modules/CHECK_BAM_RECAL'
include { COMBINE_REPORTS_RECAL } from './Modules/COMBINE_REPORTS_RECAL'
include { MKDUP } from './Modules/MKDUP'
include { CLEANUP_MKDUP_RECAL } from './Modules/CLEANUP_MKDUP_RECAL'
include { CHECK_BAM_MKDUP } from './Modules/CHECK_BAM_MKDUP'
include { COMBINE_REPORTS_MKDUP } from './Modules/COMBINE_REPORTS_MKDUP'
include { MOSDEPTH } from './Modules/MOSDEPTH'
include { CONPAIR } from './Modules/CONPAIR'
include { SAGE } from './Modules/SAGE'
include { STRELKA } from './Modules/STRELKA'
include { MuSE } from './Modules/MuSE'

include { CNVkit_buildcnn } from './Modules/CNVkit_buildcnn'
include { CNVkit } from './Modules/CNVkit'
include { CNVkit_flatref } from './Modules/CNVkit_flatref'

include { Create_sample_tsv } from './Modules/Delly/Create_sample_tsv'
include { Delly_SVcalling } from './Modules/Delly/Delly_SVcalling'
include { Delly_Prefiltering } from './Modules/Delly/Delly_Prefiltering'
include { Delly_Filtering } from './Modules/Delly/Delly_Filtering'
include { Delly_Filtering_final } from './Modules/Delly/Delly_Filtering_final'

include { MANTA } from './Modules/MANTA'

include { ASCAT } from './Modules/ASCAT'
include { ASCAT_allelecount } from './Modules/ASCAT_exome/ASCAT_allelecount'
include { ASCAT_logrbaf } from './Modules/ASCAT_exome/ASCAT_logrbaf'
include { ASCAT_exome } from './Modules/ASCAT_exome/ASCAT_exome'

include { MUTECT2_CALLING } from './Modules/Mutect2/MUTECT2_CALLING'
include { GETpileUP } from './Modules/Mutect2/GETpileUP'
include { GETpileUP_Merge } from './Modules/Mutect2/GETpileUP_Merge'
include { LearnReadOrientationModel } from './Modules/Mutect2/LearnReadOrientationModel'
include { MergeMutectStats } from './Modules/Mutect2/MergeMutectStats'
include { MergeVcfs } from './Modules/Mutect2/MergeVcfs'
include { CalculateContamination } from './Modules/Mutect2/CalculateContamination'
include { FilterMutectCalls } from './Modules/Mutect2/FilterMutectCalls'
// >>> mutect2_pon
include { MUTECT2_COHORT_PON } from './Modules/MUTECT2_PON/MUTECT2_COHORT_PON'
// <<< mutect2_pon

include { RENAME_BAM_HEADER } from './Modules/RENAME_BAM_HEADER'
include { PREPARE_TARGETS } from './Modules/TARGETS/PREPARE_TARGETS'
include { GENERATE_TARGET_INTERVALS } from './Modules/TARGETS/GENERATE_TARGET_INTERVALS'
include { POST_TUMOR_ONLY } from './Modules/POST_TUMOR_ONLY'
include { GENERATE_INTERVALS } from './Modules/GENERATE_INTERVALS'

include { POST } from './Modules/POST.nf'
include { POST_rename } from './Modules/POST_rename.nf'
include { POST_scenario1 } from './Modules/POST_scenario1.nf'
include { POST_scenario2 } from './Modules/POST_scenario2.nf'
include { POST_scenario3 } from './Modules/POST_scenario3.nf'
include { POST_scenario4 } from './Modules/POST_scenario4.nf'
include { POST_scenario5 } from './Modules/POST_scenario5.nf'
include { POST_scenario6 } from './Modules/POST_scenario6.nf'
include { MERGE_HEADER } from './Modules/MERGE_HEADER.nf'
// >>> tmb
include { TMB } from './Modules/TMB/TMB.nf'
// <<< tmb
// >>> vep
include { VEP_ANNOTATE } from './Modules/VEP/VEP_ANNOTATE.nf'
include { VCF2MAF } from './Modules/VEP/VCF2MAF.nf'
include { VEP_COHORT_MAF } from './Modules/VEP/VEP_COHORT_MAF.nf'
// <<< vep
include { ALLELECOUNTER } from './Modules/ALLELECOUNTER.nf'

include { SAVE_CSV_FASTQC } from './Modules/SAVE_CSV/SAVE_CSV_FASTQC'
include { SAVE_CSV_BWA_MEM } from './Modules/SAVE_CSV/SAVE_CSV_BWA_MEM'
include { SAVE_CSV_MKDUP } from './Modules/SAVE_CSV/SAVE_CSV_MKDUP'
include { SAVE_CSV_RECAL } from './Modules/SAVE_CSV/SAVE_CSV_RECAL'
include { SAVE_CSV_SAGE } from './Modules/SAVE_CSV/SAVE_CSV_SAGE'
include { SAVE_CSV_STRELKA } from './Modules/SAVE_CSV/SAVE_CSV_STRELKA'
include { SAVE_CSV_MuSE } from './Modules/SAVE_CSV/SAVE_CSV_MuSE'
include { SAVE_CSV_Mutect2 } from './Modules/SAVE_CSV/SAVE_CSV_Mutect2'
include { SAVE_CSV_CONPAIR } from './Modules/SAVE_CSV/SAVE_CSV_CONPAIR'
include { SAVE_CSV_MOSDEPTH } from './Modules/SAVE_CSV/SAVE_CSV_MOSDEPTH'
// >>> msi
include { MSI_SCAN } from './Modules/MSI/MSI_SCAN'
include { MSI_PAIRED } from './Modules/MSI/MSI_PAIRED'
include { MSI_PRO_NORMAL } from './Modules/MSI/MSI_PRO_NORMAL'
include { MSI_BASELINE } from './Modules/MSI/MSI_BASELINE'
include { MSI_TUMOR_ONLY } from './Modules/MSI/MSI_TUMOR_ONLY'
include { MSI_SUMMARY } from './Modules/MSI/MSI_SUMMARY'
// <<< msi
// >>> umi
include { UMI_FASTQ_TO_BAM } from './Modules/UMI/UMI_FASTQ_TO_BAM'
include { UMI_ALIGN as UMI_ALIGN_RAW } from './Modules/UMI/UMI_ALIGN'
include { UMI_ALIGN as UMI_ALIGN_CONSENSUS } from './Modules/UMI/UMI_ALIGN'
include { UMI_GROUP } from './Modules/UMI/UMI_GROUP'
include { UMI_CONSENSUS } from './Modules/UMI/UMI_CONSENSUS'
include { CHECK_BAM_MKDUP as CHECK_BAM_UMI } from './Modules/CHECK_BAM_MKDUP'
include { COMBINE_REPORTS_MKDUP as COMBINE_REPORTS_UMI } from './Modules/COMBINE_REPORTS_MKDUP'
// <<< umi

// Sample-sheet checks and tumour-only routing (feature/clinical-bam).
// Reads only the sheet (IDs and paths), never a BAM. Tumours without a normal used to vanish silently at
// normal.cross(tumor); now they are called tumour-only from any entry point (GRCh38/GRCh37), or the run stops.
def sheetFile = file(params.sample)
def tumorOnlyPatients = [] as Set
def cnvkitRefMode = (params.cnvkit_reference == 'flat') ? 'flat' : 'pooled'
if (sheetFile.exists()) {
    def lines = sheetFile.readLines().findAll { it.trim() }
    def header = lines[0].split(',').collect { it.trim() }
    def iP = header.indexOf('patient')
    def iS = header.indexOf('sample')
    def iT = header.indexOf('status')
    if (iP < 0 || iS < 0 || iT < 0) {
        error "Sample sheet ${params.sample} needs patient, sample and status columns"
    }
    def problems = []
    def normals = [:].withDefault { 0 }
    def tumors = [:].withDefault { 0 }
    def seen = [] as Set
    lines.drop(1).eachWithIndex { line, k ->
        def f = line.split(',', -1) as List
        if (f.size() != header.size()) {
            problems << "line ${k + 2}: ${f.size()} fields but the header has ${header.size()} (comma inside a field?)"
            return
        }
        // splitCsv below keeps spaces and strips quotes, so IDs read here must already be clean or the
        // tumour-only routing (matched on patient ID) would miss them.
        if (f.any { it != it.trim() || it.contains('"') }) {
            problems << "line ${k + 2}: remove spaces, quotes or Windows line endings (\\r) around values"
            return
        }
        def pt = f[iP]
        def sm = f[iS]
        def st = f[iT]
        if (!(st in ['tumor', 'normal'])) {
            problems << "line ${k + 2}: status must be tumor or normal"
        }
        // ASCAT-exome splits file names (<patient>_<sample>_<status>_...) on _normal_ / _tumor_
        def ascatSelected = params.tool.toString().split(',').collect { it.trim() }.contains('ascat')
        if (ascatSelected && [pt, sm].any { "_${it}_".contains('_normal_') || "_${it}_".contains('_tumor_') }) {
            problems << "line ${k + 2}: 'tumor'/'normal' must not appear as an underscore-separated word in patient or sample IDs"
        }
        if (!seen.add([pt, sm, st])) problems << "line ${k + 2}: duplicate patient,sample,status"
        if (st == 'normal') normals[pt] += 1
        if (st == 'tumor') tumors[pt] += 1
    }
    // normal.cross(tumor) would emit one pair per normal x tumour
    normals.findAll { it.value > 1 }.each { problems << "a patient has ${it.value} normal rows; keep one" }
    if (problems) {
        error "Sample sheet ${params.sample}:\n  " + problems.join('\n  ')
    }
    tumorOnlyPatients = (tumors.keySet() - normals.keySet()) as Set
    // >>> mutect2_pon
    if (params.mutect2_cohort_pon == 'build') {
        if (normals.size() < 2) {
            error "--mutect2_cohort_pon build needs at least 2 normal rows in the sample sheet (found ${normals.size()})."
        }
        if (normals.size() < (params.mutect2_cohort_pon_warn_below as int)) {
            log.warn "Mutect2 cohort PON from only ${normals.size()} normals: a site enters the PON when >= 2 normals " +
                     "carry it, so few normals catch only the most recurrent artefacts. ~${params.mutect2_cohort_pon_warn_below}+ recommended."
        }
        log.info "Mutect2 cohort PON: building from ${normals.size()} normal BAM(s) (paired and normal-only rows)"
    }
    // <<< mutect2_pon
    def selectedTools = params.tool ? params.tool.toString().split(',').collect { it.trim() } : []
    // CNVkit reference: pooled from the run's normals, or flat (no normals). auto picks pooled when the sheet has
    // normals. A pooled reference with no normals would never be built (collect() on an empty channel emits
    // nothing) and CNVkit would silently not run, so that combination stops here.
    if (selectedTools.contains('cnvkit')) {
        if (params.cnvkit_reference == 'pooled' && normals.isEmpty()) {
            error "--cnvkit_reference pooled needs at least one normal row in the sample sheet. " +
                  "Use --cnvkit_reference flat (or auto) for a cohort without normals."
        }
        cnvkitRefMode = (params.cnvkit_reference == 'auto') ? (normals.isEmpty() ? 'flat' : 'pooled') : params.cnvkit_reference
        if (cnvkitRefMode == 'flat') {
            log.warn "CNVkit uses a FLAT reference (no normals): GC/repeat corrected but no per-bait capture correction. " +
                     "Large CNAs stay visible; small or low-amplitude calls are unreliable. Do not pool with pooled-reference results."
        }
    }
    def nPaired = tumors.keySet().intersect(normals.keySet()).size()
    log.info "Sample sheet: ${nPaired} patient(s) with a matched normal, ${tumorOnlyPatients.size()} tumour-only"
    // Tools that only run on tumour-normal pairs get no input when no patient is paired; say so instead of
    // finishing "successfully" without them.
    def pairOnly = selectedTools.findAll { it in ['manta', 'ascat', 'delly'] }
    if (nPaired == 0 && pairOnly) {
        log.warn "--tool ${pairOnly.join(',')} will not run: no patient has a matched normal, and these tools need one."
    }
    if (tumorOnlyPatients && !(params.genome in ['GRCh38', 'GRCh37'])) {
        error "${tumorOnlyPatients.size()} patient(s) have a tumour but no normal. Tumour-only calling needs " +
              "--genome GRCh38 or GRCh37 (panel of normals, germline resource, SAGE)."
    }
    if (tumorOnlyPatients) {
        def report = new File("${params.results_dir}/REPORT/tumor_only_patients.txt")
        report.parentFile.mkdirs()
        report.text = tumorOnlyPatients.sort().join('\n') + '\n'
        log.warn "${tumorOnlyPatients.size()} tumour-only patient(s): Mutect2 + SAGE tumour-only, consensus in POST/tumor_only " +
                 "(Strelka2, MuSE, Conpair, ASCAT, Delly and Manta need a normal). List: REPORT/tumor_only_patients.txt"
    }
}
params.tumor_only_patients = tumorOnlyPatients
params.cnvkit_ref_mode = cnvkitRefMode
// >>> msi
// MSI: tumour-only patients are scored against a baseline built from normals of the same assay (the run's
// normals, including normal-only pool rows) or a prebuilt --msi_baseline. Neither -> stop before any job runs.
def msiNormals = 0
if (params.tool && params.tool.toString().split(',').collect { it.trim() }.contains('msi')) {
    if (params.msi_baseline && !file(params.msi_baseline).exists()) {
        error "--msi_baseline file not found: ${params.msi_baseline}"
    }
    if (params.msi_sites && !file(params.msi_sites).exists()) {
        error "--msi_sites file not found: ${params.msi_sites}"
    }
    if (sheetFile.exists()) {
        def rows = sheetFile.readLines().findAll { it.trim() }
        def iStatus = rows[0].split(',').collect { it.trim() }.indexOf('status')
        msiNormals = rows.drop(1).count { (it.split(',', -1) as List)[iStatus] == 'normal' }
    }
    if (tumorOnlyPatients && !params.msi_baseline) {
        if (msiNormals == 0) {
            error "--tool msi: ${tumorOnlyPatients.size()} tumour-only patient(s) need an MSI baseline, but the sheet has no " +
                  "normals. Add normal rows from the same assay (normal-only rows are allowed) or pass --msi_baseline <file>."
        }
        if (msiNormals < 20) {
            log.warn "MSI tumour-only baseline is built from ${msiNormals} normal(s); MSIsensor-pro recommends >= 20 normals " +
                     "from the same assay and sequencer. Tumour-only MSI scores are less reliable with fewer."
        }
    }
}
params.msi_n_normals = msiNormals
// <<< msi

// Define a workflow completion handler
workflow.onComplete {
    def completionStatus = """
==============================================
Pipeline Execution Completed
==============================================
Status        : ${workflow.success ? 'SUCCESS' : 'FAILED'}
Completed at  : ${workflow.complete}
Duration      : ${workflow.duration}
Success       : ${workflow.success}
Exit status   : ${workflow.exitStatus}
Error report  : ${workflow.errorReport ?: 'None'}
==============================================
"""
    writeToLog(params.log_file, completionStatus)
    log.info completionStatus
}

// Helper function to check if a specific tool is selected
def isToolSelected(String tool) {
    if (params.tool) {
        def tools = params.tool.split(',').collect { it.trim() }
        return tools.contains(tool)
    }
    return false
}


workflow {
    chunk = Channel.of(1..20)

    // Mutect2/GetPileupSummaries shards: over --targets when given, else whole chromosomes of a custom ref,
    // else the database's prebuilt shards. Strelka call regions: padded targets, else the callable-regions --bed.
    if (params.targets) {
        PREPARE_TARGETS(file(params.targets_clean), file(params.ref_fai))
        GENERATE_TARGET_INTERVALS(PREPARE_TARGETS.out.padded_bed, file(params.ref), file(params.ref_fai), file(params.ref_dict))
        interval_dir_ch = GENERATE_TARGET_INTERVALS.out.interval_dir
        strelka_regions_ch = PREPARE_TARGETS.out.call_regions
    } else {
        if (params.custom_ref) {
            GENERATE_INTERVALS(
                file(params.ref),
                file(params.ref_fai),
                file(params.ref_dict)
            )
            interval_dir_ch = GENERATE_INTERVALS.out.interval_dir
        } else {
            interval_dir_ch = Channel.value(file(params.mutect2_interval_dir))
        }
        strelka_regions_ch = Channel.value([file(params.bed), file(params.bed_tbi)])
    }
    
    // Sanity check for ascat input format
    if (params.tool.toString().contains('ascat')) {
        def sampleFile = file(params.sample)
        def hasSexColumn = false
        
        if (sampleFile.exists()) {
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',')
            hasSexColumn = headers.contains("sex")
            
            if (!hasSexColumn) {
                log.error "ERROR: ASCAT tool requires 'sex' column in sample.csv file"
                exit 1
            }
        } else {
            log.error "ERROR: Cannot find input file ${params.sample}"
            exit 1
        }
    }

    // starts from scratch - mapping
    if (params.first_step in ['mapping']) {
        // Sanity check for sample sheet
        def requiredColumns = ['patient', 'sample', 'status', 'fastq_1', 'fastq_2']
        if (params.tool && isToolSelected('ascat')) {
            requiredColumns.add('sex')
        }
        
        // Check if sample sheet exists
        def sampleFile = file(params.sample)
        if (!sampleFile.exists()) {
            log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
            exit 1
        }
        
        // Check for required columns in the sample sheet
        def headerLine = sampleFile.readLines()[0]
        def headers = headerLine.split(',').collect { it.trim() }
        def missingColumns = requiredColumns.findAll { !headers.contains(it) }
        
        if (missingColumns) {
            def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
            log.error errorMsg
            exit 1
        }

        if (params.tool && params.tool.contains('ascat')) {
            Channel.fromPath(params.sample)
            | splitCsv(header:true)
            | map { row ->
                meta = row.subMap('patient', 'sample', 'status', 'fastq_1', 'fastq_2', 'sex')
                return meta
            } | set { sample_sheet }
        } else {
            Channel.fromPath(params.sample)
            | splitCsv(header:true)
            | map { row ->
                meta = row.subMap('patient', 'sample', 'status', 'fastq_1', 'fastq_2')
                return meta
            } | set { sample_sheet }
        }
        FASTQC(sample_sheet)
        // >>> umi
        // FASTQ -> uBAM (UMI in RX) -> bwa-mem2 + ZipperBams -> GroupReadsByUmi -> consensus -> re-align +
        // FilterConsensusReads. The consensus BAM stands in for MKDUP_out below; BWA_MEM, MKDUP and their checks
        // get an empty channel and do not run.
        if (params.umi) {
            def umiFastq = [:]
            if (params.umi_source == 'fastq') {
                def rows = file(params.sample).readLines().findAll { it.trim() }
                def h = rows[0].split(',').collect { it.trim() }
                rows.drop(1).each { line ->
                    def f = line.split(',', -1)
                    umiFastq["${f[h.indexOf('patient')]}|${f[h.indexOf('sample')]}|${f[h.indexOf('status')]}"] = f[h.indexOf('fastq_umi')]
                }
                def missingUmi = umiFastq.findAll { !it.value }.keySet()
                if (missingUmi) error "--umi_fastq_column: empty fastq_umi for ${missingUmi.join(', ')}"
            }
            UMI_FASTQ_TO_BAM(sample_sheet.map { meta -> [meta, umiFastq["${meta.patient}|${meta.sample}|${meta.status}"] ?: ''] })
            UMI_ALIGN_RAW(UMI_FASTQ_TO_BAM.out.ubam, 'raw')
            UMI_GROUP(UMI_ALIGN_RAW.out.raw)
            UMI_CONSENSUS(UMI_GROUP.out.grouped)
            UMI_ALIGN_CONSENSUS(UMI_CONSENSUS.out.ubam, 'consensus')
            CHECK_BAM_UMI(UMI_ALIGN_CONSENSUS.out.mkdup_bam)
            COMBINE_REPORTS_UMI(CHECK_BAM_UMI.out.individual_reports_mkdup.collect())
            sample_sheet = Channel.empty()
        }
        // <<< umi
        BWA_MEM(sample_sheet).set { BWA_MEM_out }
        CHECK_BAM_BWA(BWA_MEM_out.bam).set{ CHECK_BAM_BWA_out }
        COMBINE_REPORTS_BWA(CHECK_BAM_BWA_out.individual_reports_bwa.collect())
    }

    // starts from markduplicate
    if (params.first_step in ['mapping', 'markdup']) {
        if (params.first_step == "markdup"){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet }
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet }
            }
            RENAME_BAM_HEADER(sample_sheet)
            MKDUP(RENAME_BAM_HEADER.out.renamed_bam).set{ MKDUP_out }    
            if (params.cleanup) {
                CLEANUP_BWA(MKDUP.out.cleanup_trigger)
            }
        } else {
            MKDUP(BWA_MEM_out).set{ MKDUP_out }
            if (params.cleanup) {
                CLEANUP_BWA(MKDUP.out.cleanup_trigger)
            }            
            CHECK_BAM_MKDUP(MKDUP_out.mkdup_bam).set{ CHECK_BAM_MKDUP_out }
            COMBINE_REPORTS_MKDUP(CHECK_BAM_MKDUP_out.individual_reports_mkdup.collect())
        }
        // >>> umi
        if (params.umi) {
            MKDUP_out = UMI_ALIGN_CONSENSUS.out   // same emits as MKDUP: pair_mutect, mkdup_bam, cleanup_trigger
        }
        // <<< umi

        if (params.genome in ['mm39', 'RN7']){
            MKDUP_out.pair_mutect.filter{it[1].status == 'normal'}.set{normal}
            MKDUP_out.pair_mutect.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            RECALIBRATE_out = [pair_recal: MKDUP_out.pair_mutect]
        }
    }

    // starts from recalibration
    if (params.first_step in ['mapping', 'markdup', 'recalibration']) {
        if (params.first_step == "recalibration" && params.genome in ['GRCh38', 'GRCh37']){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map {row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet_raw }
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet_raw }
            }

            // Rename BAM headers to standardized format
            RENAME_BAM_HEADER(sample_sheet_raw)

            RENAME_BAM_HEADER.out.renamed_bam
            | map { meta, bam, bai ->
                [meta.patient, meta, bam, bai]
            } | set { sample_sheet }

            if (params.targeted && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal_exome(sample_sheet).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_out }
            } else if (params.type == "genome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal(sample_sheet, chunk, interval_dir_ch).set{ RECALIBRATE_BaseRecal_out }

                RECALIBRATE_BaseRecal_out.MergeReport_input.groupTuple(by:[0,1]).set { BaseRecal_out_pair }
                BaseRecal_out_pair.map{
                    [patient:it[0], status:it[1].status, meta:it[1], bam:it[2][0], table:it[3], bai:it[4][0]]
                }.set{ BaseRecal_out_MAP }

                RECALIBRATE_MergeReport(BaseRecal_out_MAP).set{ RECALIBRATE_MergeReport_out }

                RECALIBRATE_BQSR(RECALIBRATE_MergeReport_out.BQSR_input, chunk, interval_dir_ch).set { RECALIBRATE_BQSR_out }

                RECALIBRATE_BQSR_out.MergeBam_input.groupTuple(by:[0,1]).set { BaseRecal_BQSR_out_pair }
                BaseRecal_BQSR_out_pair.map{
                    [patient:it[0], meta:it[1], bam:it[2]]
                }.set{BaseRecal_BQSR_out_MAP}

                RECALIBRATE_MergeBam(BaseRecal_BQSR_out_MAP).set{ RECALIBRATE_MergeBam_out }

                RECALIBRATE_SortBam(RECALIBRATE_MergeBam_out.MergeBam_input).set{ RECALIBRATE_out }

                if (params.cleanup) {
                    CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
                }
            }
        } else {
            if (params.targeted && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal_exome(MKDUP_out.pair_mutect).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_out }
            } else if (params.type == "genome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal(MKDUP_out.pair_mutect, chunk, interval_dir_ch).set{ RECALIBRATE_BaseRecal_out }

                RECALIBRATE_BaseRecal_out.MergeReport_input.groupTuple(by:[0,1]).set { BaseRecal_out_pair }
                BaseRecal_out_pair.map{
                    [patient:it[0], status:it[1].status, meta:it[1], bam:it[2][0], table:it[3], bai:it[4][0]]
                }.set{ BaseRecal_out_MAP }

                RECALIBRATE_MergeReport(BaseRecal_out_MAP).set{ RECALIBRATE_MergeReport_out }

                RECALIBRATE_BQSR(RECALIBRATE_MergeReport_out.BQSR_input, chunk, interval_dir_ch).set { RECALIBRATE_BQSR_out }

                RECALIBRATE_BQSR_out.MergeBam_input.groupTuple(by:[0,1]).set { BaseRecal_BQSR_out_pair }
                BaseRecal_BQSR_out_pair.map{
                    [patient:it[0], meta:it[1], bam:it[2]]
                }.set{BaseRecal_BQSR_out_MAP}

                RECALIBRATE_MergeBam(BaseRecal_BQSR_out_MAP).set{ RECALIBRATE_MergeBam_out }

                RECALIBRATE_SortBam(RECALIBRATE_MergeBam_out.MergeBam_input).set{ RECALIBRATE_out }

                if (params.cleanup) {
                    CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
                }
            }
        }
        CHECK_BAM_RECAL(RECALIBRATE_out.pair_recal).set{ CHECK_BAM_RECAL_out }
        COMBINE_REPORTS_RECAL(CHECK_BAM_RECAL_out.individual_reports_recal.collect())
    }

    // starts from variant calling
    if (params.first_step in ['mapping', 'markdup', 'recalibration', 'variant_calling']) {
        if (params.first_step == "variant_calling"){
            // Sanity check for sample sheet
            def requiredColumns = ['patient', 'sample', 'status', 'bam', 'bai']
            if (params.tool && isToolSelected('ascat')) {
                requiredColumns.add('sex')
            }
            
            // Check if sample sheet exists
            def sampleFile = file(params.sample)
            if (!sampleFile.exists()) {
                log.error "ERROR: Sample sheet file '${params.sample}' does not exist"
                exit 1
            }
            
            // Check for required columns in the sample sheet
            def headerLine = sampleFile.readLines()[0]
            def headers = headerLine.split(',').collect { it.trim() }
            def missingColumns = requiredColumns.findAll { !headers.contains(it) }
            
            if (missingColumns) {
                def errorMsg = "ERROR: Sample sheet is missing required column(s): ${missingColumns.join(', ')}"
                log.error errorMsg
                exit 1
            }

            if (params.tool && params.tool.contains('ascat')) {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)  
                | map {row ->
                    meta = row.subMap('patient', 'sample', 'status', 'sex')
                    return [meta, file(row.bam), file(row.bai)]
                } | set {sample_sheet_raw}
            } else {
                Channel.fromPath(params.sample)
                | splitCsv(header:true)
                | map { row ->
                    meta = row.subMap('patient', 'sample', 'status')
                    return [meta, file(row.bam), file(row.bai)]
                } | set { sample_sheet_raw }
            }

            // Rename BAM headers to standardized format
            RENAME_BAM_HEADER(sample_sheet_raw)

            RENAME_BAM_HEADER.out.renamed_bam
            | map { meta, bam, bai ->
                [meta.patient, meta, bam, bai]
            } | set { sample_sheet }

            GETpileUP(sample_sheet, chunk, interval_dir_ch).set { GETpileUP_out }

            sample_sheet.filter{it[1].status == 'normal'}.set{normal}
            sample_sheet.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            // for CN analysis
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], gender:normal[1].sex, normal:normal[2],tumor:tumor[2], sample:tumor[1].sample]
            }.set{ RECALIBRATE_out_MAP_CN }

            // Collect all normal bams for CNVkit and Delly
            normal.map { patient, meta, bam, bai -> [bam, bai] }.flatten().collect().set { normal_bams }

            // Tumour-only patients (no normal row): same maps with normal = null. Mutect2, SAGE, mosdepth and
            // CNVkit (against the pooled normals) handle that; Strelka2, MuSE, Conpair, ASCAT, Delly and Manta get pairs only.
            tumor.filter { it[0] in params.tumor_only_patients }.map { patient, meta, bam, bai ->
                [patient:patient, normal:null, tumor:bam, tumor_meta:meta, normal_meta:null]
            }.set { TUMOR_ONLY_MAP }
            tumor.filter { it[0] in params.tumor_only_patients }.map { patient, meta, bam, bai ->
                [patient:patient, gender:meta.sex, normal:null, tumor:bam, sample:meta.sample]
            }.set { TUMOR_ONLY_MAP_CN }

            if (params.genome in ['GRCh38', 'GRCh37']) {
                SAGE(RECALIBRATE_out_MAP.mix(TUMOR_ONLY_MAP))
            }
            
            STRELKA(RECALIBRATE_out_MAP, strelka_regions_ch)
            MuSE(RECALIBRATE_out_MAP)
            
            if (params.genome in ['GRCh38', 'GRCh37']) {
                CONPAIR(RECALIBRATE_out_MAP)
                MOSDEPTH(RECALIBRATE_out_MAP.mix(TUMOR_ONLY_MAP))
            }

            // >>> mutect2_pon (was: MUTECT2_CALLING(..., chunk, interval_dir_ch))
            if (params.mutect2_cohort_pon == 'build') {
                mutect2_pon_ch = MUTECT2_COHORT_PON(normal, chunk, interval_dir_ch,
                    params.targets ? PREPARE_TARGETS.out.padded_bed : interval_dir_ch)
            } else {
                mutect2_pon_ch = Channel.value(params.mutect2_pon_resolved.collect { file(it) })
            }
            MUTECT2_CALLING(RECALIBRATE_out_MAP.mix(TUMOR_ONLY_MAP), chunk, interval_dir_ch, mutect2_pon_ch).set { MUTECT2_CALLING_out }
            // <<< mutect2_pon
            MUTECT2_CALLING_out.LearnReadOrientationModel_input.groupTuple(by:0).set { LearnReadOrientationModel_input_pair }
            MUTECT2_CALLING_out.MergeMutectStats_input.groupTuple(by:0).set { MergeMutectStats_input_pair }
            MUTECT2_CALLING_out.MergeVcfs_input.groupTuple(by:0).set { MergeVcfs_input_pair }
            LearnReadOrientationModel(LearnReadOrientationModel_input_pair).set { LearnReadOrientationModel_out }
            MergeMutectStats(MergeMutectStats_input_pair).set { MergeMutectStats_out }
            MergeVcfs(MergeVcfs_input_pair).set { MergeVcfs_out }

            GETpileUP_out.GETpileUP_Merge_input.map{patient, meta, table, chunk -> [patient, meta, [table, chunk]]}.groupTuple(by:[0,1],sort: {it[1]}).set { GETpileUP_out_pair }
            GETpileUP_out_pair.map{
            [patient:it[0], meta:it[1], mix:it[2]]
            }.set{ GETpileUP_out_pair_MAP }

            GETpileUP_Merge(GETpileUP_out_pair_MAP).set { GETpileUP_Merge_out }

            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
            GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
            }.set{ GETpileUP_out_MAP }

            GETpileUP_tumor.filter { it[0] in params.tumor_only_patients }.map { t ->
                [patient:t[0], normal:null, tumor:t[2], tumor_meta:t[1]]
            }.set { GETpileUP_tumor_only_MAP }

            CalculateContamination(GETpileUP_out_MAP.mix(GETpileUP_tumor_only_MAP)).set { CalculateContamination_out }
            def vcfChannel = MergeVcfs_out.MUTECT2_vcf
                .map{ map, vcf ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                }
            def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                .map{ map, table ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, table]]
                }
            def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                .map{ map, segments ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, segments]]
                }
            def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                .map{ map, orientation ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, orientation]]
                }
            def statsChannel = MergeMutectStats_out.MUTECT2_stats
                .map{ map, stats ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                }

            // Join by the composite key (patient_sample)
            vcfChannel
                .join(contaminationChannel)
                .join(segmentsChannel)
                .join(orientationChannel)
                .join(statsChannel)
                .map{ patient_sample, vcfData, contData, segData, oriData, statsData ->
                    [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                }
                .set{ filterInput }

            FilterMutectCalls(filterInput).set{ FILTER_OUT }
            // SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)

            // TOOL-SPECIFIC BLOCKS - MOVED INSIDE variant_calling STEP
            // ASCAT
            if (params.tool && params.tool.contains('ascat')) {
                chromosomes = Channel.of( *(1..22).collect { it.toString() } + ['X'] )
                if (params.targeted) {
                    ASCAT_allelecount(RECALIBRATE_out_MAP_CN, chromosomes)

                    ASCAT_allelecount.out.allelecount
                        .map { map, chr, file ->
                            def patient = map.patient
                            def sample = map.sample
                            def gender = map.gender
                            [patient, sample, gender, file]
                        }
                        .groupTuple(by: [0, 1, 2])  // Group by patient, sample, and gender
                        // Restructure from pairs to separate normal and tumor lists
                        .map { patient, sample, gender, files ->
                            // Extract normal and tumor files from the nested structure
                            def normalFiles = []
                            def tumorFiles = []

                            // Flatten any nested structure if present
                            def flatFiles = files.flatten()

                            // Separate normal and tumor files
                            flatFiles.each { file ->
                                if (file.toString().contains("_normal_")) {
                                    normalFiles << file
                                } else if (file.toString().contains("_tumor_")) {
                                    tumorFiles << file
                                }
                            }

                            [patient, sample, gender, normalFiles, tumorFiles]
                        }
                        .set { ascat_logrbaf_input }

                    ASCAT_logrbaf(ascat_logrbaf_input)
                    ASCAT_exome(ASCAT_logrbaf.out.ascat_input)

                } else if (params.type == "genome") {
                    ASCAT(RECALIBRATE_out_MAP_CN)
                }
            }

            // CNVkit
            if (params.tool && params.tool.contains('cnvkit')) {
                if (params.cnvkit_ref_mode == 'flat') {
                    cnvkit_ref_ch = CNVkit_flatref(Channel.value('flat')).CNVkit_ref_cnn
                } else {
                    cnvkit_ref_ch = CNVkit_buildcnn(normal_bams).CNVkit_ref_cnn
                }
                CNVkit(RECALIBRATE_out_MAP_CN.mix(TUMOR_ONLY_MAP_CN), cnvkit_ref_ch).set { CNVkit_out }
            }

            // Delly
            if (params.tool && params.tool.contains('delly')) {
                if (params.targeted) {
                    log.error "ERROR: Delly is not recommended for exome data. Please use a different SV caller (manta) for exome data"
                    exit 1  
                } else {
                    Create_sample_tsv(Channel.fromPath(params.sample, checkIfExists: true))
                    sample_tsv_collected = Create_sample_tsv.out.samples_tsv.collect()

                    Delly_SVcalling(RECALIBRATE_out_MAP_CN).set { SVcalling_out }
                    SVcalling_out.SVcalling_bcf
                    .combine(sample_tsv_collected)
                    .set { prefiltering_input }

                    Delly_Prefiltering(prefiltering_input).set { Prefiltering_out }
                    Delly_Filtering(Prefiltering_out.Delly_pre_bcf, normal_bams).set { Filtering_out }

                    Filtering_out.Delly_geno_bcf
                    .combine(sample_tsv_collected)
                    .set { filtering_final_input }

                    Delly_Filtering_final(filtering_final_input)
                }
            }

            // Manta
            if (params.tool && params.tool.contains('manta')) {
                MANTA(RECALIBRATE_out_MAP_CN)
            }

            // >>> msi
            // MSIsensor-pro: pairs with msi (tumour vs its normal), tumour-only with pro against a baseline of normals
            if (params.tool && isToolSelected('msi')) {
                msi_sites_ch = params.msi_sites ? Channel.value(file(params.msi_sites)) :
                    MSI_SCAN(params.targets ? PREPARE_TARGETS.out.padded_bed : Channel.value('genome')).sites
                MSI_PAIRED(RECALIBRATE_out_MAP, msi_sites_ch)
                if (params.msi_baseline) {
                    msi_baseline_ch = Channel.value(file(params.msi_baseline))
                } else {
                    MSI_PRO_NORMAL(normal.map { patient, meta, bam, bai -> ["${patient}_${meta.sample}", bam] }, msi_sites_ch)
                    msi_baseline_ch = MSI_BASELINE(MSI_PRO_NORMAL.out.all.collect(), msi_sites_ch).baseline
                }
                MSI_TUMOR_ONLY(TUMOR_ONLY_MAP, msi_baseline_ch)
                MSI_SUMMARY(MSI_PAIRED.out.status.map { it[1] }.mix(MSI_TUMOR_ONLY.out.status.map { it[1] }).collect())
            }
            // <<< msi

        } else {

            RECALIBRATE_out.pair_recal.filter{it[1].status == 'normal'}.set{normal}
            RECALIBRATE_out.pair_recal.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            // for CN analysis
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], gender:normal[1].sex, normal:normal[2],tumor:tumor[2], sample:tumor[1].sample]
            }.set{ RECALIBRATE_out_MAP_CN }

            // Collect all normal bams for CNVkit and Delly
            normal.map { patient, meta, bam, bai -> [bam, bai] }.flatten().collect().set { normal_bams }

            // Tumour-only patients (no normal row), as in the variant_calling entry: normal = null maps for
            // Mutect2, SAGE, mosdepth and CNVkit; Strelka2, MuSE, Conpair, ASCAT, Delly and Manta get pairs only.
            tumor.filter { it[0] in params.tumor_only_patients }.map { patient, meta, bam, bai ->
                [patient:patient, normal:null, tumor:bam, tumor_meta:meta, normal_meta:null]
            }.set { TUMOR_ONLY_MAP }
            tumor.filter { it[0] in params.tumor_only_patients }.map { patient, meta, bam, bai ->
                [patient:patient, gender:meta.sex, normal:null, tumor:bam, sample:meta.sample]
            }.set { TUMOR_ONLY_MAP_CN }

            if (params.genome in ['GRCh38', 'GRCh37']) {
                SAGE(RECALIBRATE_out_MAP.mix(TUMOR_ONLY_MAP))
            }

            STRELKA(RECALIBRATE_out_MAP, strelka_regions_ch)
            MuSE(RECALIBRATE_out_MAP)
            
            if (params.genome in ['GRCh38', 'GRCh37']) {
                CONPAIR(RECALIBRATE_out_MAP)
                MOSDEPTH(RECALIBRATE_out_MAP.mix(TUMOR_ONLY_MAP))
            }

            GETpileUP(RECALIBRATE_out.pair_recal, chunk, interval_dir_ch).set { GETpileUP_out }
            // >>> mutect2_pon (was: MUTECT2_CALLING(..., chunk, interval_dir_ch))
            if (params.mutect2_cohort_pon == 'build') {
                mutect2_pon_ch = MUTECT2_COHORT_PON(normal, chunk, interval_dir_ch,
                    params.targets ? PREPARE_TARGETS.out.padded_bed : interval_dir_ch)
            } else {
                mutect2_pon_ch = Channel.value(params.mutect2_pon_resolved.collect { file(it) })
            }
            MUTECT2_CALLING(RECALIBRATE_out_MAP.mix(TUMOR_ONLY_MAP), chunk, interval_dir_ch, mutect2_pon_ch).set { MUTECT2_CALLING_out }
            // <<< mutect2_pon
            MUTECT2_CALLING_out.LearnReadOrientationModel_input.groupTuple(by:0).set { LearnReadOrientationModel_input_pair }
            MUTECT2_CALLING_out.MergeMutectStats_input.groupTuple(by:0).set { MergeMutectStats_input_pair }
            MUTECT2_CALLING_out.MergeVcfs_input.groupTuple(by:0).set { MergeVcfs_input_pair }
            LearnReadOrientationModel(LearnReadOrientationModel_input_pair).set { LearnReadOrientationModel_out }
            MergeMutectStats(MergeMutectStats_input_pair).set { MergeMutectStats_out }
            MergeVcfs(MergeVcfs_input_pair).set { MergeVcfs_out }

            GETpileUP_out.GETpileUP_Merge_input.map{patient, meta, table, chunk -> [patient, meta, [table, chunk]]}.groupTuple(by:[0,1],sort: {it[1]}).set { GETpileUP_out_pair }
            GETpileUP_out_pair.map{
            [patient:it[0], meta:it[1], mix:it[2]]
            }.set{ GETpileUP_out_pair_MAP }

            GETpileUP_Merge(GETpileUP_out_pair_MAP).set { GETpileUP_Merge_out }

            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'normal'}.set{ GETpileUP_normal }
            GETpileUP_Merge_out.CalculateContamination_input.filter{it[1].status== 'tumor'}.set{ GETpileUP_tumor }
            GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2],tumor:tumor[2], tumor_meta:tumor[1]]
            }.set{ GETpileUP_out_MAP }

            GETpileUP_tumor.filter { it[0] in params.tumor_only_patients }.map { t ->
                [patient:t[0], normal:null, tumor:t[2], tumor_meta:t[1]]
            }.set { GETpileUP_tumor_only_MAP }

            CalculateContamination(GETpileUP_out_MAP.mix(GETpileUP_tumor_only_MAP)).set { CalculateContamination_out }
            def vcfChannel = MergeVcfs_out.MUTECT2_vcf
                .map{ map, vcf ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
                }
            def contaminationChannel = CalculateContamination_out.MUTECT2_contamination_table
                .map{ map, table ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, table]]
                }
            def segmentsChannel = CalculateContamination_out.MUTECT2_segments_table
                .map{ map, segments ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, segments]]
                }
            def orientationChannel = LearnReadOrientationModel_out.MUTECT2_read_orientation
                .map{ map, orientation ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, orientation]]
                }
            def statsChannel = MergeMutectStats_out.MUTECT2_stats
                .map{ map, stats ->
                    ["${map.patient}_${map.tumor_meta.sample}", [map, stats]]
                }

            // Join by the composite key (patient_sample)
            vcfChannel
                .join(contaminationChannel)
                .join(segmentsChannel)
                .join(orientationChannel)
                .join(statsChannel)
                .map{ patient_sample, vcfData, contData, segData, oriData, statsData ->
                    [vcfData[0], vcfData[1], contData[1], segData[1], oriData[1], statsData[1]]
                }
                .set{ filterInput }

            FilterMutectCalls(filterInput).set{ FILTER_OUT }
            // SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)


            // TOOL-SPECIFIC BLOCKS - ALSO FOR NON variant_calling START STEPS
	    // ASCAT
	    if (params.tool && params.tool.contains('ascat')) {
	        chromosomes = Channel.of( *(1..22).collect { it.toString() } + ['X'] )
	        if (params.targeted) {
	            ASCAT_allelecount(RECALIBRATE_out_MAP_CN, chromosomes)
	
	            ASCAT_allelecount.out.allelecount
	                .map { map, chr, file ->
	                    def patient = map.patient
	                    def sample = map.sample
	                    def gender = map.gender
	                    [patient, sample, gender, file]
	                }
	                .groupTuple(by: [0, 1, 2])  // Group by patient, sample, and gender
	                // Restructure from pairs to separate normal and tumor lists
	                .map { patient, sample, gender, files ->
	                    // Extract normal and tumor files from the nested structure
	                    def normalFiles = []
	                    def tumorFiles = []
	
	                    // Flatten any nested structure if present
	                    def flatFiles = files.flatten()
	
	                    // Separate normal and tumor files
	                    flatFiles.each { file ->
	                        if (file.toString().contains("_normal_")) {
	                            normalFiles << file
	                        } else if (file.toString().contains("_tumor_")) {
	                            tumorFiles << file
	                        }
	                    }
	
	                    [patient, sample, gender, normalFiles, tumorFiles]
	                }
	                .set { ascat_logrbaf_input }
	
	            ASCAT_logrbaf(ascat_logrbaf_input)
	            ASCAT_exome(ASCAT_logrbaf.out.ascat_input)
	
	        } else if (params.type == "genome") {
	            ASCAT(RECALIBRATE_out_MAP_CN)
	        }
	    }
	
	    // CNVkit
	    if (params.tool && params.tool.contains('cnvkit')) {
	        if (params.cnvkit_ref_mode == 'flat') {
	            cnvkit_ref_ch = CNVkit_flatref(Channel.value('flat')).CNVkit_ref_cnn
	        } else {
	            cnvkit_ref_ch = CNVkit_buildcnn(normal_bams).CNVkit_ref_cnn
	        }
	        CNVkit(RECALIBRATE_out_MAP_CN.mix(TUMOR_ONLY_MAP_CN), cnvkit_ref_ch).set { CNVkit_out }
	    }
	
	    // Delly
	    if (params.tool && params.tool.contains('delly')) {
	        if (params.targeted && params.tool.toString().contains('delly')) {
	            log.error "ERROR: Delly is not recommended for exome data. Please use a different SV caller (manta) for exome data"
	            exit 1  
	        } else {
	            Create_sample_tsv(Channel.fromPath(params.sample, checkIfExists: true))
	            sample_tsv_collected = Create_sample_tsv.out.samples_tsv.collect()
	
	            Delly_SVcalling(RECALIBRATE_out_MAP_CN).set { SVcalling_out }
	            SVcalling_out.SVcalling_bcf
	            .combine(sample_tsv_collected)
	            .set { prefiltering_input }
	
	            Delly_Prefiltering(prefiltering_input).set { Prefiltering_out }
	            Delly_Filtering(Prefiltering_out.Delly_pre_bcf, normal_bams).set { Filtering_out }
	
	            Filtering_out.Delly_geno_bcf
	            .combine(sample_tsv_collected)
	            .set { filtering_final_input }
	
	            Delly_Filtering_final(filtering_final_input)
	        }
	    }
	
	    // Manta
	    if (params.tool && params.tool.contains('manta')) {
	        MANTA(RECALIBRATE_out_MAP_CN)
	    }

            // >>> msi
            // MSIsensor-pro, as in the variant_calling entry
            if (params.tool && isToolSelected('msi')) {
                msi_sites_ch = params.msi_sites ? Channel.value(file(params.msi_sites)) :
                    MSI_SCAN(params.targets ? PREPARE_TARGETS.out.padded_bed : Channel.value('genome')).sites
                MSI_PAIRED(RECALIBRATE_out_MAP, msi_sites_ch)
                if (params.msi_baseline) {
                    msi_baseline_ch = Channel.value(file(params.msi_baseline))
                } else {
                    MSI_PRO_NORMAL(normal.map { patient, meta, bam, bai -> ["${patient}_${meta.sample}", bam] }, msi_sites_ch)
                    msi_baseline_ch = MSI_BASELINE(MSI_PRO_NORMAL.out.all.collect(), msi_sites_ch).baseline
                }
                MSI_TUMOR_ONLY(TUMOR_ONLY_MAP, msi_baseline_ch)
                MSI_SUMMARY(MSI_PAIRED.out.status.map { it[1] }.mix(MSI_TUMOR_ONLY.out.status.map { it[1] }).collect())
            }
            // <<< msi
        }
    }

// POST - filtering and consensus calling
    if (params.first_step in ['mapping', 'markdup', 'recalibration', 'variant_calling'] && 
        params.genome in ['GRCh38', 'GRCh37']) {
        
        // Use the outputs that are guaranteed to exist from variant calling steps
        mutect2Channel = FilterMutectCalls.out.Mutect2_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        museChannel = MuSE.out.MuSE_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        strelkaChannel = STRELKA.out.STRELKA_out
            .map { map, snv, indel ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, snv, indel]]
            }
        
        sageChannel = SAGE.out.SAGE_out
            .map { map, vcf ->
                ["${map.patient}_${map.tumor_meta.sample}", [map, vcf]]
            }
        
        // Tumour-only consensus (Mutect2 + SAGE). Paired samples go on to the 4-caller scenarios below;
        // tumour-only ones drop out of those joins because MuSE/Strelka2 never ran for them.
        FilterMutectCalls.out.Mutect2_out
            .filter { map, vcf -> !map.normal }
            .map { map, vcf -> ["${map.patient}_${map.tumor_meta.sample}", vcf] }
            .join(SAGE.out.SAGE_out
                .filter { map, vcf -> !map.normal }
                .map { map, vcf -> ["${map.patient}_${map.tumor_meta.sample}", vcf] })
            .set { tumorOnlyPostInput }
        POST_TUMOR_ONLY(tumorOnlyPostInput)

        // Recreate recalChannel from the available data
        // Choose the source based on first_step
        sourceChannel = (params.first_step == "variant_calling") ? 
            sample_sheet : 
            RECALIBRATE_out.pair_recal

        // Create recalChannel from the chosen source
        normal_post = sourceChannel.filter{it[1].status == 'normal'}
        tumor_post = sourceChannel.filter{it[1].status == 'tumor'}

        recalChannel = normal_post.cross(tumor_post){it[0]}.map{
            normal, tumor ->
            ["${tumor[1].patient}_${tumor[1].sample}", [tumor[1], tumor[2], tumor[3]]]
        }
        
        postInput = mutect2Channel
            .join(museChannel)
            .join(strelkaChannel)
            .join(sageChannel)
            .join(recalChannel)
            .map { patient_sample, mutect2_vcf, muse_vcf, strelka_vcf, sage_vcf, recal_bam ->
                [patient_sample, mutect2_vcf[1], muse_vcf[1], strelka_vcf[1], strelka_vcf[2], sage_vcf[1], recal_bam[1], recal_bam[2]]
            }
        
        // POST(postInput)
        POST_rename(postInput)

        POST_scenario1(POST_rename.out.renamed_postInput)
        POST_scenario2(POST_rename.out.renamed_postInput)
        POST_scenario3(POST_rename.out.renamed_postInput)
        POST_scenario4(POST_rename.out.renamed_postInput)
        POST_scenario5(POST_rename.out.renamed_postInput)
        POST_scenario6(POST_rename.out.renamed_postInput)

        merge_3 = POST_scenario3.out.merge_input.map { file -> tuple(file, 'scenario_3') }
        merge_4 = POST_scenario4.out.merge_input.map { file -> tuple(file, 'scenario_4') }
        merge_5 = POST_scenario5.out.merge_input.map { file -> tuple(file, 'scenario_5') }
        merge_6 = POST_scenario6.out.merge_input.map { file -> tuple(file, 'scenario_6') }

        merge_all = merge_3.mix(merge_4, merge_5, merge_6)

		MERGE_HEADER(merge_all)

        // >>> vep
        // VEP + vcf2maf: the paired standard consensus VCF of --vep_scenario and the tumour-only PASS consensus,
        // one MAF per sample, then MAF/cohort.maf with an Analysis_Mode column.
        if (isToolSelected('vep')) {
            // Sample IDs as the callers name their VCF columns (Mutect2/SAGE): <patient>_<sample>_tumor|_normal
            vep_pair_ids = normal_post.cross(tumor_post){it[0]}.map { normal, tumor ->
                ["${tumor[1].patient}_${tumor[1].sample}".toString(),
                 "${tumor[1].patient}_${tumor[1].sample}_tumor".toString(),
                 "${normal[1].patient}_${normal[1].sample}_normal".toString()]
            }
            def vep_suffix = "_${params.vep_scenario}_standard_PASSed.vcf".toString()
            vep_paired = MERGE_HEADER.out.final_standard_vcf
                .filter { vcf -> vcf.name.endsWith(vep_suffix) }
                .map { vcf -> [vcf.name.substring(0, vcf.name.length() - vep_suffix.length()), vcf] }
                .join(vep_pair_ids)
                .map { patient_sample, vcf, tumor_id, normal_id -> [patient_sample, 'paired', tumor_id, normal_id, vcf] }
            vep_tumor_only = POST_TUMOR_ONLY.out.consensus_pass
                .map { patient_sample, vcf, tbi -> [patient_sample, 'tumor_only', "${patient_sample}_tumor".toString(), 'NORMAL', vcf] }
            VEP_ANNOTATE(vep_paired.mix(vep_tumor_only))
            VCF2MAF(VEP_ANNOTATE.out.vcf)
            VEP_COHORT_MAF(VCF2MAF.out.maf.map { patient_sample, mode, maf -> maf }.collect())
        }
        // <<< vep

		

        // AlleleCounter for VAF calculation
        def vaf_source = (params.first_step == "variant_calling") ? 
            RENAME_BAM_HEADER.out.renamed_bam.map { meta, bam, bai -> [meta.patient, meta, bam, bai] } :
            RECALIBRATE_out.pair_recal

        vaf_source.filter{it[1].status == 'normal'}.set{normal_vaf}
        vaf_source.filter{it[1].status == 'tumor'}.set{tumor_vaf}

        normal_vaf.cross(tumor_vaf){it[0]}.map{
            normal, tumor ->
            [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1], tumor_bai:tumor[3], normal_bai:normal[3]]
        }.set{ RECALIBRATE_out_MAP_VAF }

        RECALIBRATE_out_MAP_VAF
            .map { map ->
                def patient_sample = "${map.patient}_${map.tumor_meta.sample}"
                [patient_sample, map]
            }
            .set { recal_map_keyed }

        POST_scenario4.out.final_vcf
            .join(recal_map_keyed)
            .map { patient_sample, vcf_file, recal_map -> 
                [
                    patient_sample: patient_sample,
                    final_vcf: vcf_file,
                    tumor: recal_map.tumor,
                    normal: recal_map.normal,
                    tumor_bai: recal_map.tumor_bai,
                    normal_bai: recal_map.normal_bai
                ]
            }
            .set { allelecounter_input_4 }


        POST_scenario5.out.final_vcf
            .join(recal_map_keyed)
            .map { patient_sample, vcf_file, recal_map -> 
                [
                    patient_sample: patient_sample,
                    final_vcf: vcf_file,
                    tumor: recal_map.tumor,
                    normal: recal_map.normal,
                    tumor_bai: recal_map.tumor_bai,
                    normal_bai: recal_map.normal_bai
                ]
            }
            .set { allelecounter_input_5 }

        POST_scenario6.out.final_vcf
            .join(recal_map_keyed)
            .map { patient_sample, vcf_file, recal_map -> 
                [
                    patient_sample: patient_sample,
                    final_vcf: vcf_file,
                    tumor: recal_map.tumor,
                    normal: recal_map.normal,
                    tumor_bai: recal_map.tumor_bai,
                    normal_bai: recal_map.normal_bai
                ]
            }
            .set { allelecounter_input_6 }


        ch_out1 = allelecounter_input_4.map { file -> tuple(file, 'scenario_4') }
        ch_out2 = allelecounter_input_5.map { file -> tuple(file, 'scenario_5') }
        ch_out3 = allelecounter_input_6.map { file -> tuple(file, 'scenario_6') }

        ch_all_outputs = ch_out1.mix(ch_out2, ch_out3)

        ALLELECOUNTER(ch_all_outputs)

        // >>> tmb
        // One block serves every entry point (mapping/markdup/recalibration via RECALIBRATE_out, variant_calling via
        // RENAME_BAM_HEADER): vaf_source / recal_map_keyed above already pick the right BAMs for each path.
        if (params.tool && isToolSelected('tmb')) {
            def tmbScenario = params.tmb_scenario.toString()
            def tmbPairedVcf = [
                '4': POST_scenario4.out.final_vcf,
                '5': POST_scenario5.out.final_vcf,
                '6': POST_scenario6.out.final_vcf
            ][tmbScenario]
            // HOOK (annotation branch): to count coding-only on an annotated VCF, replace tmbPairedVcf /
            // tmbTumorOnlyVcf with channels of [patient_sample, annotated_vcf] (CSQ or ANN in the header).
            def tmbTumorOnlyVcf = POST_TUMOR_ONLY.out.consensus_pass.map { ps, vcf, tbi -> [ps, vcf] }

            // HOOK (PureCN branch): replace this with the in-pipeline PureCN [patient_sample, variants.csv] channel.
            def tmbGermline = Channel.empty()
            if (params.tmb_germline_calls) {
                tmbGermline = tumor_vaf
                    .filter { it[0] in params.tumor_only_patients }
                    .map { p, meta, bam, bai -> "${p}_${meta.sample}" }
                    .map { ps ->
                        def f = file("${params.tmb_germline_calls}/${ps}_variants.csv")
                        if (!f.exists()) log.warn "TMB ${ps}: no ${f.name} in --tmb_germline_calls; no germline filter for this sample"
                        [ps, f.exists() ? f : []]
                    }
            }

            def tmbPaired = tmbPairedVcf
                .join(recal_map_keyed)
                .map { ps, vcf, m -> [ps, 'paired', "scenario_${tmbScenario}", vcf, m.tumor, m.tumor_bai, m.normal, m.normal_bai, []] }

            def tmbTumorOnlyBams = tumor_vaf
                .filter { it[0] in params.tumor_only_patients }
                .map { p, meta, bam, bai -> ["${p}_${meta.sample}", bam, bai] }
            def tmbTumorOnly = tmbTumorOnlyVcf
                .join(tmbTumorOnlyBams)
                .join(tmbGermline, remainder: true)
                .filter { it.size() == 5 && it[1] != null && it[2] != null }
                .map { ps, vcf, bam, bai, germ -> [ps, 'tumor_only', 'tumor_only_consensus', vcf, bam, bai, [], [], germ ?: []] }

            TMB(tmbPaired.mix(tmbTumorOnly), file(params.tmb_regions ?: params.targets_clean))

            TMB.out.tmb_tsv
                .map { ps, tsv -> tsv }
                .collectFile(name: 'tmb_cohort.tsv', keepHeader: true, skip: 1, sort: true, storeDir: params.tmb_dir)
            TMB.out.tmb_tsv.subscribe { ps, tsv ->
                def rows = tsv.readLines()
                def h = rows[0].split('\t') as List
                def v = rows[1].split('\t') as List
                def w = v[h.indexOf('warning')]
                if (w != 'none') log.warn "TMB ${ps}: ${w}"
            }
        }
        // <<< tmb
    }
}

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
// >>> qc
params.somalier_sites = "release"          // --tool somalier sites: release (public somalier sites) | targets (common gnomAD SNPs in the padded targets) | <sites.vcf.gz>
params.somalier_release_sites = null       // override the release sites path (default <database>/Databases/<genome>/somalier/sites.<build>.vcf.gz)
params.somalier_min_af = 0.2               // targets mode: gnomAD AF window
params.somalier_max_af = 0.8
params.somalier_swap_max_relatedness = 0.6 // same patient below this -> swap (loose: tumour CN/LOH lowers T-N relatedness)
params.somalier_dup_min_relatedness = 0.7  // different patients at or above this -> duplicate/mislabel
params.somalier_min_sites = 50             // fewer shared sites -> not judged (same-patient pairs flagged insufficient_sites)
params.somalier_min_x_sites = 10           // informative chrX sites needed to infer sex
params.somalier_relate_args = ""
params.multiqc = true                      // one MultiQC report at the end: --multiqc false to skip
params.multiqc_config = null
// <<< qc
// >>> purecn
params.purecn_normaldb = null           // dir with normalDB*.rds + intervals.txt (a previous run's PureCN/normaldb/); else built from this run's normals
params.purecn_mappability = null        // optional mappability bigWig/BED for IntervalFile.R (same contig names as --ref)
params.purecn_off_target = false        // add off-target bins (IntervalFile.R --off-target)
params.purecn_assay = ""                // assay name used in the normal DB file name
params.purecn_min_pop_af = 0.001        // gnomAD AF threshold for the germline SNP sites (and the DB flag)
params.purecn_min_mapq = 20             // allele counting at SNP sites (bcftools mpileup -q)
params.purecn_min_bq = 20               // allele counting at SNP sites (bcftools mpileup -Q)
params.purecn_min_purity = 0.15         // PureCN purity search range (always; a sheet purity does not narrow it)
params.purecn_max_purity = 0.95
params.purecn_purity_window = 0.15      // with a sheet `purity`: pick the most likely PureCN solution within purity +- this
params.purecn_model = "betabin"         // PureCN.R --model (beta | betabin)
params.purecn_sex = "?"                 // PureCN.R --sex (? detect | F | M | diploid)
params.purecn_matched_normal_coverage = false  // paired: normalise against the matched normal instead of the normal DB
params.purecn_args = ""                 // extra PureCN.R arguments, passed verbatim
params.purecn_gc_norm = true            // GC loess normalisation of coverage (Coverage.R); keep true for real data
// <<< purecn
// >>> pcgr
// PCGR 2.3.2 clinical interpretation of the consensus calls (--tool pcgr). Env/data: Modules/PCGR/install_pcgr.sh
params.pcgr_env = null                    // prebuilt conda env <dir>/pcgr (sibling <dir>/pcgrr must exist)
params.pcgr_refdata = null                // PCGR bundle dir holding data/grch38 (bundle 20260620 for PCGR 2.3.2)
params.pcgr_vep_dir = null                // VEP cache parent dir holding homo_sapiens/115_GRCh38
params.pcgr_scenario = 5                  // paired input: POST/standard_vcf scenario 3|4|5|6
params.pcgr_assay = null                  // WGS|WES|TARGETED; default from --type (genome|exome|panel)
params.pcgr_effective_target_size_mb = null // default: merged size of --targets (capped at 34), else 34
params.pcgr_tumor_site = 0                // PCGR primary-site code used when the sheet has no tumor_site value (0 = any)
params.pcgr_msi = true                    // only applied to WES/WGS tumour-normal pairs (PCGR's own restriction)
params.pcgr_signatures = true
params.pcgr_min_mutations_signatures = 200  // PCGR default; >= 100. Panels rarely reach it, so fitting is skipped
params.pcgr_tumor_dp_min = null
params.pcgr_tumor_af_min = null
params.pcgr_tmb_dp_min = null
params.pcgr_tmb_af_min = null
params.pcgr_gnomad_popmax_af = 0.001      // tumour-only germline filter (PCGR default)
params.pcgr_tonly_germline_filters = true // tumour-only: --exclude_likely_hom_germline --exclude_likely_het_germline
params.pcgr_pon_vcf = null                // HOOK (cohort PON): VCF with INFO flag PANEL_OF_NORMALS, tumour-only only
params.pcgr_cna_from_purecn = false       // HOOK (PureCN): feed PureCN loh.csv segments as --input_cna
params.pcgr_extra_args = ''               // passed verbatim to pcgr
// <<< pcgr
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
// >>> gene_status
// Per-gene somatic status (--tool gene_status): one row per covered panel gene and tumour, combining callable
// coverage (mosdepth), the consensus SNV/indel calls (VEP-annotated when --tool vep is on), CNVkit copy number
// and PureCN integer CN/LOH when those tools ran. Covered genes = genes BED ∩ padded --targets. The genes BED
// is --genes_bed, else the 4th column of --targets (bin/make_genes_bed.py builds one from an Ensembl/GENCODE
// GTF). Paired variants come from the POST consensus of --tmb_scenario (or --vep_scenario when vep is on).
params.genes_bed = null                       // 4-column BED: chrom, start, end, gene name
params.gene_status_min_depth = 50             // tumour depth for a targeted base of a gene to count as callable
params.gene_status_min_callable = 0.5         // callable fraction below this -> status NOT_CALLABLE
params.gene_status_min_mapq = 20              // reads counted by mosdepth for depth
params.gene_status_dir = "${params.results_dir}/REPORT/gene_status"
params.gene_status_env = "$projectDir/yml/gene_status.yml"
// <<< gene_status
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
        error "--vep_scenario must be scenario_3, scenario_4, scenario_5 or scenario_6 (the POST_SCENARIOS final VCFs)."
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

if (!(params.type in ['genome', 'exome', 'panel'])) {
    error "ERROR: --type must be 'genome', 'exome' or 'panel' (got: ${params.type == null ? 'not set' : "'${params.type}'"}).\n" +
          "Usage: nextflow run main.nf --type [genome|exome|panel] --genome [GRCh38|GRCh37|mm39|RN7] ..."
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
// >>> gene_status
// Resolve the genes BED once (gunzip, drop track/browser/# lines, keep 4 columns, check contigs against --ref)
// so the module and the -preview DAG read one clean file: --genes_bed, else the 4th column of --targets.
def geneStatusGenes = null
if (params.tool && params.tool.toString().split(',').collect { it.trim() }.contains('gene_status')) {
    if (!params.targets) {
        error "--tool gene_status needs --targets: covered genes are the genes BED intersected with the panel targets."
    }
    if (!(params.genome in ['GRCh38', 'GRCh37'])) {
        error "--tool gene_status reads the POST consensus, which runs only for --genome GRCh38 or GRCh37."
    }
    def gsSrc, gsWhat
    if (params.genes_bed) {
        gsSrc = file(params.genes_bed)
        if (!gsSrc.exists()) error "--genes_bed file not found: ${params.genes_bed}"
        gsWhat = "--genes_bed ${params.genes_bed}"
    } else {
        gsSrc = file(params.targets)
        gsWhat = "the 4th column of --targets ${params.targets}"
    }
    def gsRaw = gsSrc.name.endsWith('.gz') ?
        new java.util.zip.GZIPInputStream(gsSrc.newInputStream()).newReader().readLines() : gsSrc.readLines()
    def gsRows = gsRaw.findAll { it.trim() && !(it ==~ /^(track|browser|#).*/) }
        .collect { it.split('\t') }
        .findAll { it.size() >= 4 && it[3].trim() }
    if (!gsRows) {
        error "--tool gene_status found no gene names in ${gsWhat}. Give --genes_bed <BED: chrom,start,end,gene> " +
              "or a --targets BED whose 4th column is the gene name (bin/make_genes_bed.py builds a genes BED " +
              "from an Ensembl/GENCODE GTF)."
    }
    def gsFaiContigs = file(params.ref_fai).readLines().collect { it.split('\t')[0] } as Set
    def gsMissing = gsRows.collect { it[0] }.unique().findAll { !(it in gsFaiContigs) }
    if (gsMissing) {
        error "gene_status genes BED uses contigs that are not in --ref: ${gsMissing.take(5).join(' ')} " +
              "(chr prefix mismatch? bin/make_genes_bed.py --add-chr renames Ensembl contigs)"
    }
    def gsOut = new File("${params.results_dir}/custom_intervals/genes.clean.bed")
    gsOut.parentFile.mkdirs()
    def gsText = gsRows.collect { "${it[0]}\t${it[1]}\t${it[2]}\t${it[3].trim()}" }.join('\n') + '\n'
    if (!gsOut.exists() || gsOut.text != gsText) {
        gsOut.text = gsText
    }
    geneStatusGenes = gsOut.absolutePath
    log.info "gene_status: ${gsRows.collect { it[3].trim() }.unique().size()} gene(s) from ${gsWhat}; " +
             "min depth ${params.gene_status_min_depth}, min callable fraction ${params.gene_status_min_callable}"
}
params.gene_status_genes = geneStatusGenes
// <<< gene_status
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

// >>> qc
params.somalier_dir = "${params.results_dir}/somalier"
params.multiqc_dir = "${params.results_dir}/MultiQC"
params.somalier_env = "${projectDir}/yml/somalier.yml"
params.multiqc_env = "${projectDir}/yml/multiqc.yml"
def somalierSelected = params.tool && params.tool.toString().split(',').collect { it.trim() }.contains('somalier')
def somalierSitesFile = null
if (somalierSelected) {
    if (!(params.genome in ['GRCh38', 'GRCh37'])) {
        error "--tool somalier needs --genome GRCh38 or GRCh37 (human common-SNP sites)."
    }
    if (params.somalier_sites == 'targets') {
        if (!params.targets) error "--somalier_sites targets needs --targets."
        if (!file(params.mutect2_germline).exists() || !file("${params.mutect2_germline}.tbi").exists()) {
            error "--somalier_sites targets reads common SNPs from ${params.mutect2_germline} (+ .tbi), which is missing."
        }
    } else {
        def build = (params.genome == 'GRCh38') ? 'hg38' : 'GRCh37'
        somalierSitesFile = (params.somalier_sites == 'release') ?
            (params.somalier_release_sites ?: "${params.genome_database_dir}/somalier/sites.${build}.vcf.gz") : params.somalier_sites
        if (!file(somalierSitesFile).exists()) {
            error "somalier sites not found: ${somalierSitesFile}\n" +
                  "  Download the public sites (GRCh38: https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz, " +
                  "GRCh37: https://github.com/brentp/somalier/files/3412455/sites.GRCh37.vcf.gz) and pass " +
                  "--somalier_release_sites <file>, or use --somalier_sites targets with --targets."
        }
    }
}
params.somalier_sites_file = somalierSitesFile
// <<< qc

params.vaf_script = params.vaf_script ?: "${projectDir}/VAF.py"

// Custom reference: every packaged resource this run will use must share the reference's contig names and
// lengths. severity = what the consuming tool does on a mismatch: GATK aborts ('error'); SAGE, mosdepth, MuSE,
// Strelka/Manta and Delly silently ignore or mis-annotate ('warn'). A missing file is always an error.
if (params.custom_ref) {
    def human = params.genome in ['GRCh38', 'GRCh37']
    def tools = params.tool ? params.tool.toString().split(',')*.trim() : []
    def resources = [
        [label: '--ref_dict',       path: params.ref_dict,         kind: 'dict'],
        [label: '--bed',            path: params.bed,              kind: 'bed', severity: 'warn'],
        [label: 'Mutect2 germline', path: params.mutect2_germline, kind: 'vcf'],
        // cohort-build runs never read the database/user PON, so skip its contig check (empty path = skipped)
        [label: 'Mutect2 PON',      path: params.mutect2_cohort_pon ? '' : params.mutect2_pon, kind: 'vcf'],
        [label: 'MuSE dbSNP',       path: params.muse_vcf,         kind: 'vcf', severity: 'warn'],
    ]
    if (params.targeted)
        resources << [label: 'mosdepth exome BED', path: params.mosdepth_bed, kind: 'bed', severity: 'warn']
    if (human && params.first_step != 'variant_calling') {
        resources << [label: 'BQSR known sites 1', path: params.recal_knownsite1, kind: 'vcf']
        resources << [label: 'BQSR known sites 2', path: params.recal_knownsite2, kind: 'vcf']
        if (params.targeted)
            resources << [label: 'BQSR exome intervals', path: params.recal_interval_wes, kind: 'interval_list']
    }
    if (human) {
        resources << [label: 'SAGE hotspots',            path: params.SAGE_hotspots,            kind: 'vcf', severity: 'warn']
        resources << [label: 'SAGE panel BED',           path: params.SAGE_panel_bed,           kind: 'bed', severity: 'warn']
        resources << [label: 'SAGE high-confidence BED', path: params.SAGE_high_confidence_bed, kind: 'bed', severity: 'warn']
        resources << [label: 'Conpair markers',          path: params.conpair_marker,           kind: 'bed']
    }
    if ('delly' in tools)
        resources << [label: 'Delly exclude list', path: params.delly_excl, kind: 'bed', severity: 'warn', allowUnknown: true]

    def check = ReferenceCheck.run(params.ref_fai, resources)
    check.warnings.each { log.warn "Custom reference check: ${it}" }
    if (check.errors) {
        error "ERROR: custom reference ${params.ref} is not compatible with the packaged resources for genome '${params.genome}':\n  - " +
              check.errors.join('\n  - ') +
              "\nProvide resources built for this reference through the matching --<resource> options, or use the packaged reference."
    }
    log.info "Custom reference check: ${resources.count { it.path }} packaged resources checked against ${params.ref}, ${check.warnings.size()} warning(s)"
}

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
// >>> pcgr
params.pcgr_dir="${params.results_dir}/PCGR"
// <<< pcgr
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
params.bcftools_env = "$projectDir/yml/bcftools.yml"   // clinical: POST_TUMOR_ONLY (upstream removed POST_rename.yml)
params.strelka_env = "$projectDir/yml/strelka_env.yml"
params.mosdepth_env = "$projectDir/yml/mosdepth_env.yml"
params.rename_bam_env = "$projectDir/yml/rename_bam.yml"
params.fastqc_env = "$projectDir/yml/fastqc_env.yml"
params.cnvkit_env = "$projectDir/yml/cnvkit.yml"
params.delly_env = "$projectDir/yml/delly.yml"
params.ascat_env = "$projectDir/yml/ascat.yml"
params.alleleCounter_env = "$projectDir/yml/alleleCounter.yml"
params.manta_env = "$projectDir/yml/manta.yml"
params.SAGE_java_env = "$projectDir/yml/sage_java.yml"
params.muse_env = "$projectDir/yml/muse.yml"
params.picard_merge_env = "$projectDir/yml/picard.yml"
params.gatk_env = "${projectDir}/yml/gatk_runtime.yml"
params.POST_env = "${projectDir}/yml/POST.yml"
params.gatk = "${params.database_path}/gatk-4.6.2.0/gatk"
// >>> umi
params.umi_dir = "${params.results_dir}/UMI"
params.umi_env = "${projectDir}/yml/umi.yml"
// <<< umi
// >>> purecn
params.purecn_dir = "${params.results_dir}/PureCN"
params.purecn_env = "${projectDir}/yml/purecn.yml"
// PureCN (--tool purecn): targeted assays on GRCh38/GRCh37 only. Needs a normal DB: built from >= 2 normal rows of the
// sheet (paired or normal-only), or a prebuilt one (--purecn_normaldb). With exactly one normal and no DB, only that
// patient's paired tumour can run (against its matched normal). Optional sheet column `purity` (0-1, or 1-100 %)
// selects among that tumour's PureCN solutions (see Modules/PURECN/PURECN.nf). Reads only the sheet (IDs, status, purity), never a BAM.
// Set before the includes: modules only see params that exist when they are included.
params.purecn_genome = [GRCh38: 'hg38', GRCh37: 'hg19'][params.genome]
params.purecn_targets = params.targets_clean ?: params.genomes[params.genome].mosdepth_bed
def purecnPurity = [:]
def purecnBuildDb = false
if (params.tool && params.tool.toString().split(',').collect { it.trim() }.contains('purecn')) {
    if (!params.purecn_genome) {
        error "--tool purecn needs --genome GRCh38 or GRCh37 (PureCN genome hg38/hg19)."
    }
    if (!params.targeted) {
        error "--tool purecn is for targeted data: use --type panel (with --targets) or --type exome."
    }
    if (!file(params.purecn_targets).exists()) {
        error "--tool purecn: target BED not found: ${params.purecn_targets}"
    }
    if (!file(params.mutect2_germline).exists() || !file("${params.mutect2_germline}.tbi").exists()) {
        error "--tool purecn needs the gnomAD germline resource and its .tbi: ${params.mutect2_germline}"
    }
    if (params.purecn_normaldb) {
        def dbDir = file(params.purecn_normaldb)
        if (!dbDir.isDirectory() || !dbDir.listFiles().any { it.name ==~ /normalDB.*\.rds/ } ||
            !file("${params.purecn_normaldb}/intervals.txt").exists()) {
            error "--purecn_normaldb must be a directory with normalDB*.rds and intervals.txt " +
                  "(the PureCN/normaldb/ folder of an earlier run of the same assay): ${params.purecn_normaldb}"
        }
    }
    if (params.purecn_mappability && !file(params.purecn_mappability).exists()) {
        error "--purecn_mappability file not found: ${params.purecn_mappability}"
    }
    def pSheet = file(params.sample)
    def pLines = pSheet.exists() ? pSheet.readLines().findAll { it.trim() } : []
    def pHead = pLines ? pLines[0].split(',').collect { it.trim() } : []
    def iP = pHead.indexOf('patient'), iS = pHead.indexOf('sample'), iT = pHead.indexOf('status')
    def iPur = pHead.indexOf('purity')
    def nNormals = 0
    def pTumorPts = [] as Set, pNormalPts = [] as Set
    def purityProblems = []
    if (iP < 0 || iS < 0 || iT < 0) pLines = []   // the sheet checks below report the missing columns
    pLines.drop(1).each { line ->
        def f = line.split(',', -1) as List
        if (f.size() != pHead.size()) return
        if (f[iT] == 'normal') { nNormals++; pNormalPts << f[iP] }
        if (f[iT] == 'tumor') pTumorPts << f[iP]
        if (iPur >= 0 && f[iT] == 'tumor' && f[iPur]) {
            def v = f[iPur].isNumber() ? (f[iPur] as double) : -1d
            if (v > 1 && v <= 100) v = v / 100
            if (v <= 0 || v > 1) purityProblems << "${f[iP]},${f[iS]}: purity '${f[iPur]}' is not a fraction (0-1] or a percentage"
            else purecnPurity[(f[iP] + "_" + f[iS]).toString()] = v
        }
    }
    if (purityProblems) {
        error "Sample sheet ${params.sample}, column purity:\n  " + purityProblems.join('\n  ')
    }
    if (!params.purecn_normaldb) {
        if (nNormals == 0) {
            error "--tool purecn needs normals for its normal DB: add normal rows (normal-only rows are fine) or give " +
                  "--purecn_normaldb <dir> from an earlier run of the same assay."
        }
        purecnBuildDb = nNormals >= 2
        if (!purecnBuildDb && (pTumorPts - pNormalPts)) {
            error "--tool purecn: the sheet has one normal, and tumour-only patients need a normal DB of >= 2 normals. " +
                  "Add normal-only rows or give --purecn_normaldb <dir>."
        }
        if (!purecnBuildDb) {
            log.warn "PureCN: one normal only, so no normal DB; the paired tumour is normalised against its matched normal."
        }
    }
    log.info "PureCN: genome ${params.purecn_genome}, targets ${params.purecn_targets}, normal DB " +
             (params.purecn_normaldb ?: (purecnBuildDb ? "built from ${nNormals} normals" : 'none')) +
             (purecnPurity ? ", sheet purity for ${purecnPurity.size()} tumour(s)" : '')
}
params.purecn_purity_map = purecnPurity
params.purecn_build_normaldb = purecnBuildDb
// <<< purecn

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
// >>> qc
include { SOMALIER_SITES } from './Modules/SOMALIER/SOMALIER_SITES'
include { SOMALIER_EXTRACT } from './Modules/SOMALIER/SOMALIER_EXTRACT'
include { SOMALIER_RELATE } from './Modules/SOMALIER/SOMALIER_RELATE'
include { BCFTOOLS_STATS } from './Modules/MULTIQC/BCFTOOLS_STATS'
include { MULTIQC } from './Modules/MULTIQC/MULTIQC'
// <<< qc
// >>> purecn
include { PURECN_WF } from './Modules/PURECN/PURECN_WF'
// <<< purecn

// >>> tmb
include { TMB } from './Modules/TMB/TMB.nf'
// <<< tmb
// >>> gene_status
include { GENE_STATUS; GENE_STATUS_MATRIX } from './Modules/GENE_STATUS/GENE_STATUS.nf'
// <<< gene_status
// >>> vep
include { VEP_ANNOTATE } from './Modules/VEP/VEP_ANNOTATE.nf'
include { VCF2MAF } from './Modules/VEP/VCF2MAF.nf'
include { VEP_COHORT_MAF } from './Modules/VEP/VEP_COHORT_MAF.nf'
// <<< vep
include { POST_SCENARIOS } from './Modules/POST_scenarios.nf'
include { POST_FINAL_VCF } from './Modules/POST_clinical/POST_FINAL_VCF.nf'
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
// >>> pcgr
include { PCGR } from './Modules/PCGR/PCGR'
// <<< pcgr

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

// >>> pcgr
// PCGR settings resolved before any job starts: assay, TMB target size, per-tumour primary site (optional sheet
// column tumor_site, PCGR codes 0-30). Reads only the sheet header and that column, never a BAM/VCF.
def pcgrSites = [:]
def pcgrAssay = params.pcgr_assay ?: (params.type == 'panel' ? 'TARGETED' : (params.type == 'exome' ? 'WES' : 'WGS'))
def pcgrSizeMb = 34.0
if (params.tool && params.tool.toString().split(',').collect { it.trim() }.contains('pcgr')) {
    def problems = []
    if (!(params.genome in ['GRCh38', 'GRCh37'])) problems << "--tool pcgr needs --genome GRCh38 or GRCh37"
    if (!(pcgrAssay in ['WGS', 'WES', 'TARGETED'])) problems << "--pcgr_assay must be WGS, WES or TARGETED"
    if (!(params.pcgr_scenario.toString() in ['3', '4', '5', '6'])) problems << "--pcgr_scenario must be 3, 4, 5 or 6"
    def asm = params.genome == 'GRCh37' ? 'grch37' : 'grch38'
    if (!params.pcgr_env || !file("${params.pcgr_env}/bin/pcgr").exists()) {
        problems << "--pcgr_env must point to the PCGR conda env (<dir>/pcgr with bin/pcgr; see Modules/PCGR/install_pcgr.sh)"
    } else if (!file("${file(params.pcgr_env).parent}/pcgrr/bin/Rscript").exists()) {
        problems << "PCGR needs its R env next to --pcgr_env: ${file(params.pcgr_env).parent}/pcgrr not found"
    }
    if (!params.pcgr_refdata || !file("${params.pcgr_refdata}/data/${asm}").isDirectory()) {
        problems << "--pcgr_refdata must contain data/${asm} (PCGR reference bundle)"
    }
    if (!params.pcgr_vep_dir || !file("${params.pcgr_vep_dir}/homo_sapiens").isDirectory() ||
        !file("${params.pcgr_vep_dir}/homo_sapiens").list().any { it.endsWith("_${params.genome}") }) {
        problems << "--pcgr_vep_dir must contain homo_sapiens/<release>_${params.genome} (VEP cache)"
    }
    if (params.pcgr_effective_target_size_mb) {
        pcgrSizeMb = params.pcgr_effective_target_size_mb as double
    } else if (params.targets_clean) {
        // merged length of the clean (unpadded) targets
        def iv = file(params.targets_clean).readLines().collect { it.split('\t') }
            .collect { [it[0], it[1] as long, it[2] as long] }.sort { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] }
        long bp = 0; def cur = null
        iv.each { r ->
            if (cur && cur[0] == r[0] && r[1] <= cur[2]) { cur[2] = Math.max(cur[2], r[2]) }
            else { if (cur) bp += cur[2] - cur[1]; cur = [r[0], r[1], r[2]] }
        }
        if (cur) bp += cur[2] - cur[1]
        pcgrSizeMb = Math.round(bp / 1e4) / 100.0
        if (pcgrSizeMb > 34) {
            log.warn "PCGR: targets cover ${pcgrSizeMb} Mb; PCGR caps the TMB denominator at the 34 Mb coding exome."
            pcgrSizeMb = 34.0
        }
        if (pcgrSizeMb < 0.01) pcgrSizeMb = 0.01
    }
    if (pcgrAssay == 'TARGETED' && pcgrSizeMb < 1) {
        log.warn "PCGR: effective target size ${pcgrSizeMb} Mb (< 1 Mb). PCGR's TMB is not robust below 1 Mb; MSI is not " +
                 "predicted for TARGETED assays or tumour-only samples, and signature fitting needs >= ${params.pcgr_min_mutations_signatures} SNVs."
    }
    if (sheetFile.exists()) {
        def lines = sheetFile.readLines().findAll { it.trim() }
        def header = lines[0].split(',').collect { it.trim() }
        def iSite = header.indexOf('tumor_site')
        if (iSite >= 0) {
            lines.drop(1).eachWithIndex { line, k ->
                def f = line.split(',', -1) as List
                if (f[header.indexOf('status')] != 'tumor' || !f[iSite]) return
                if (!(f[iSite] ==~ /\d+/) || (f[iSite] as int) > 30) {
                    problems << "sample sheet line ${k + 2}: tumor_site must be a PCGR site code 0-30 (see pcgr --help)"
                } else {
                    pcgrSites["${f[header.indexOf('patient')]}_${f[header.indexOf('sample')]}"] = f[iSite] as int
                }
            }
        }
    }
    if (params.pcgr_pon_vcf && !file(params.pcgr_pon_vcf).exists()) problems << "--pcgr_pon_vcf not found"
    if (problems) error "PCGR:\n  " + problems.join('\n  ')
    log.info "PCGR: assay ${pcgrAssay}, effective target size ${pcgrSizeMb} Mb, paired input scenario ${params.pcgr_scenario}, " +
             "${pcgrSites.size()} tumour(s) with a tumor_site code. Reports in ${params.results_dir}/PCGR hold patient-level results: keep them on the cluster."
}
params.pcgr_assay_run = pcgrAssay
params.pcgr_target_size_mb_run = pcgrSizeMb
params.pcgr_tumor_sites = pcgrSites
// <<< pcgr

// Define a workflow completion handler
workflow.onComplete {
    def ignored = workflow.stats.ignoredCount
    def status  = !workflow.success ? 'FAILED'
                : ignored > 0       ? "COMPLETED WITH ${ignored} IGNORED TASK FAILURE(S) - outputs missing, see .nextflow.log"
                :                     'SUCCESS'
    def completionStatus = """
==============================================
Pipeline Execution Completed
==============================================
Status        : ${status}
Tasks         : ${workflow.stats.succeedCount} succeeded, ${workflow.stats.failedCount} failed, ${ignored} ignored
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

def validatePairing(sheet) {
    def sampleFile = file(sheet)
    if (!sampleFile.exists()) error "ERROR: Sample sheet file '${sheet}' does not exist"
    def header = sampleFile.readLines()[0].split(',')*.trim()
    if (!['patient', 'sample', 'status'].every { it in header }) return      // reported by the per-step column check
    def rows = sampleFile.splitCsv(header: true)
    def problems = []
    rows.groupBy { [it.patient?.trim(), it.sample?.trim(), it.status?.trim()] }.each { k, v ->
        if (v.size() > 1) problems << "patient ${k[0]}, sample ${k[1]}, status ${k[2]}: listed ${v.size()} times"
    }
    rows.groupBy { it.patient?.trim() }.each { patient, prows ->
        def status  = prows.collect { it.status?.trim() }
        def unknown = status.findAll { !(it in ['tumor', 'normal']) }.unique()
        if (unknown)                  problems << "patient ${patient}: status must be 'tumor' or 'normal' (got: ${unknown.join(', ')})"
        def nNormal = status.count('normal'), nTumor = status.count('tumor')
        // Clinical: 0 normals = tumour-only patient, 0 tumours = normal-only pool row (both handled downstream);
        // more than one normal stays ambiguous.
        if (nNormal > 1)              problems << "patient ${patient}: ${nNormal} normal sample(s), at most 1 is allowed"
    }
    if (problems) {
        error "ERROR: tumor/normal pairing of sample sheet ${sheet} is ambiguous:\n  - " + problems.join('\n  - ') +
              "\nEach patient has at most one row with status 'normal'; tumour-only and normal-only patients are allowed."
    }
    def byPatient = rows.groupBy { it.patient?.trim() }
    def nPairs = byPatient.count { p, r -> r*.status*.trim().contains('normal') && r*.status*.trim().contains('tumor') }
    log.info "Sample sheet ${sheet}: ${nPairs} paired patient(s), " +
             "${byPatient.count { p, r -> !r*.status*.trim().contains('normal') }} tumour-only, " +
             "${byPatient.count { p, r -> !r*.status*.trim().contains('tumor') }} normal-only"
}


workflow {
    validatePairing(params.sample)
    chunk = Channel.of(1..20)

    // BAMs removed by samtools quickcheck (CHECK_BAM_BWA / CHECK_BAM_MKDUP), as
    // tuple(patient, meta, stage). Used by the pair filter after CHECK_BAM_RECAL.
    qc_failed = Channel.empty()

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

        Channel.fromPath(params.sample)
        | splitCsv(header:true)
        | map { row ->
            def keys = ['patient', 'sample', 'status'] + (isToolSelected('ascat') ? ['sex'] : [])
            def meta = row.subMap(keys)
            [meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
        } | set { sample_sheet }
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

        // Number of BAMs per patient in the sample sheet (used by the pair filter)
        sample_sheet
            .map { meta, fq1, fq2 -> meta.patient }
            .toList()
            .map { patients -> patients.countBy { it } }
            .set { bams_per_patient }

        // Only BAMs that pass quickcheck go on to MarkDuplicates
        CHECK_BAM_BWA_out.checked
            .branch { meta, bam, qc ->
                pass: qc == 'PASS'
                fail: true
            }
            .set { bwa_qc }
        bwa_qc.pass
            .map { meta, bam, qc -> tuple(meta, bam) }
            .set { BWA_MEM_good }
        qc_failed = qc_failed.mix(
            bwa_qc.fail.map { meta, bam, qc ->
                log.warn "quickcheck FAIL after BWA: removing ${meta.patient}_${meta.sample}_${meta.status} (${bam.name})"
                tuple(meta.patient, meta, 'BWA')
            }
        )
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
            sample_sheet
                .map { meta, bam, bai -> meta.patient }
                .toList()
                .map { patients -> patients.countBy { it } }
                .set { bams_per_patient }

            RENAME_BAM_HEADER(sample_sheet)
            MKDUP(RENAME_BAM_HEADER.out.renamed_bam).set{ MKDUP_out }    
            if (params.cleanup) {
                CLEANUP_BWA(MKDUP.out.cleanup_trigger)
            }
        } else {
            MKDUP(BWA_MEM_good).set{ MKDUP_out }
            if (params.cleanup) {
                CLEANUP_BWA(MKDUP.out.cleanup_trigger)
            }            
        }
        // >>> umi
        if (params.umi) {
            MKDUP_out = UMI_ALIGN_CONSENSUS.out   // same emits as MKDUP: pair_mutect, mkdup_bam, cleanup_trigger
        }
        // <<< umi

        CHECK_BAM_MKDUP(MKDUP_out.mkdup_bam).set{ CHECK_BAM_MKDUP_out }
        COMBINE_REPORTS_MKDUP(CHECK_BAM_MKDUP_out.individual_reports_mkdup.collect())

        // Only BAMs that pass quickcheck go on to recalibration.
        // MKDUP_good has the same shape as MKDUP's pair_mutect output: tuple(patient, meta, bam, bai)
        CHECK_BAM_MKDUP_out.checked
            .branch { meta, bam, bai, qc ->
                pass: qc == 'PASS'
                fail: true
            }
            .set { mkdup_qc }
        mkdup_qc.pass
            .map { meta, bam, bai, qc -> tuple(meta.patient, meta, bam, bai) }
            .set { MKDUP_good }
        qc_failed = qc_failed.mix(
            mkdup_qc.fail.map { meta, bam, bai, qc ->
                log.warn "quickcheck FAIL after MarkDuplicates: removing ${meta.patient}_${meta.sample}_${meta.status} (${bam.name})"
                tuple(meta.patient, meta, 'MKDUP')
            }
        )

        if (params.genome in ['mm39', 'RN7']){
            MKDUP_good.filter{it[1].status == 'normal'}.set{normal}
            MKDUP_good.filter{it[1].status == 'tumor'}.set{tumor}
            normal.cross(tumor){it[0]}.map{
                normal, tumor ->
                [patient:normal[0], normal:normal[2], tumor:tumor[2], tumor_meta:tumor[1], normal_meta:normal[1]]
            }.set{ RECALIBRATE_out_MAP }

            RECALIBRATE_out = [pair_recal: MKDUP_good]
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

            sample_sheet_raw
                .map { meta, bam, bai -> meta.patient }
                .toList()
                .map { patients -> patients.countBy { it } }
                .set { bams_per_patient }

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

                RECALIBRATE_MergeBam_out.MergeBam_input
                    .map { patient, meta, merged -> [[patient, meta.sample, meta.status], patient, meta, merged] }
                    .join( RECALIBRATE_MergeReport_out.BQSR_input.map { patient, meta, bam, table, bai -> [[patient, meta.sample, meta.status], bam, bai] }, failOnMismatch: true )
                    .map { key, patient, meta, merged, mkdup_bam, mkdup_bai -> [patient, meta, merged, mkdup_bam, mkdup_bai] }
                    .set { SortBam_input }
                RECALIBRATE_SortBam(SortBam_input).set{ RECALIBRATE_out }

                if (params.cleanup) {
                    CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
                }
            }
        } else {
            if (params.targeted && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal_exome(MKDUP_good).set{ RECALIBRATE_BaseRecal_out }
                RECALIBRATE_BQSR_exome(RECALIBRATE_BaseRecal_out.BQSR_input).set { RECALIBRATE_out }
            } else if (params.type == "genome" && params.genome in ['GRCh38', 'GRCh37']) {
                RECALIBRATE_BaseRecal(MKDUP_good, chunk, interval_dir_ch).set{ RECALIBRATE_BaseRecal_out }

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

                RECALIBRATE_MergeBam_out.MergeBam_input
                    .map { patient, meta, merged -> [[patient, meta.sample, meta.status], patient, meta, merged] }
                    .join( RECALIBRATE_MergeReport_out.BQSR_input.map { patient, meta, bam, table, bai -> [[patient, meta.sample, meta.status], bam, bai] }, failOnMismatch: true )
                    .map { key, patient, meta, merged, mkdup_bam, mkdup_bai -> [patient, meta, merged, mkdup_bam, mkdup_bai] }
                    .set { SortBam_input }
                RECALIBRATE_SortBam(SortBam_input).set{ RECALIBRATE_out }

                if (params.cleanup) {
                    CLEANUP_MKDUP_RECAL(RECALIBRATE_out.cleanup_trigger)
                }
            }
        }
        CHECK_BAM_RECAL(RECALIBRATE_out.pair_recal).set{ CHECK_BAM_RECAL_out }
        COMBINE_REPORTS_RECAL(CHECK_BAM_RECAL_out.individual_reports_recal.collect())

        // Pair filter: keep only BAMs that passed every quickcheck, and only for
        // patients that still have at least one normal and one tumor BAM.
        // Every BAM in the sample sheet ends here exactly once, either as a
        // recalibrated BAM or as a quickcheck failure from an earlier step, so each
        // patient is released as soon as all of its own BAMs are accounted for.
        CHECK_BAM_RECAL_out.checked
            .map { patient, meta, bam, bai, qc ->
                if (qc != 'PASS') {
                    log.warn "quickcheck FAIL after recalibration: removing ${meta.patient}_${meta.sample}_${meta.status} (${bam.name})"
                }
                tuple(patient, meta, qc == 'PASS' ? [bam, bai] : null, qc == 'PASS' ? 'PASS' : 'RECAL')
            }
            .mix(qc_failed.map { patient, meta, stage -> tuple(patient, meta, null, stage) })
            .combine(bams_per_patient)
            .map { patient, meta, files, stage, counts -> tuple(groupKey(patient, counts[patient]), meta, files, stage) }
            .groupTuple()
            .flatMap { patient, metas, files, stages ->
                def rows = [metas, files, stages].transpose()
                def good = rows.findAll { m, f, s -> s == 'PASS' }
                def has_normal = good.any { m, f, s -> m.status == 'normal' }
                def has_tumor  = good.any { m, f, s -> m.status == 'tumor' }
                if (!(has_normal && has_tumor)) {
                    def missing = [has_normal ? null : 'normal', has_tumor ? null : 'tumor'].findAll()
                    log.warn "Skipping patient ${patient}: no passing ${missing.join(' or ')} BAM"
                    return []
                }
                good.collect { m, f, s -> tuple(m.patient, m, f[0], f[1]) }
            }
            .set { RECAL_good }
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
                    .map { map, chr, tumor_file, normal_file -> [map.patient, map.sample, map.gender, normal_file, tumor_file] }
                    .groupTuple(by: [0, 1, 2])   // -> [patient, sample, gender, normalFiles, tumorFiles]
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
            // >>> purecn
            // PureCN: every BAM (normals incl. normal-only rows, tumours) + this run's Mutect2 filtered calls
            if (isToolSelected('purecn')) {
                PURECN_WF(normal.mix(tumor), FILTER_OUT.Mutect2_out)
            }
            // <<< purecn

        } else {

            RECAL_good.filter{it[1].status == 'normal'}.set{normal}
            RECAL_good.filter{it[1].status == 'tumor'}.set{tumor}
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

            GETpileUP(RECAL_good, chunk, interval_dir_ch).set { GETpileUP_out }
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
                    .map { map, chr, tumor_file, normal_file -> [map.patient, map.sample, map.gender, normal_file, tumor_file] }
                    .groupTuple(by: [0, 1, 2])   // -> [patient, sample, gender, normalFiles, tumorFiles]
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
            // >>> purecn
            if (isToolSelected('purecn')) {
                PURECN_WF(normal.mix(tumor), FILTER_OUT.Mutect2_out)
            }
            // <<< purecn
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
            RECAL_good

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
		
        POST_SCENARIOS(postInput)

        // >>> clinical_post
        // Upstream POST_SCENARIOS (all six scenarios in one task) replaced the old per-scenario scripts and MERGE_HEADER.
        // Clinical additions on its outputs: sorted/indexed final VCFs with the exact commands in their headers,
        // the "no consensus mutations" report, and the VEP/PCGR wiring.
        def post_final = [
            'scenario_3': POST_SCENARIOS.out.final_vcf_s3, 'scenario_4': POST_SCENARIOS.out.final_vcf_s4,
            'scenario_5': POST_SCENARIOS.out.final_vcf_s5, 'scenario_6': POST_SCENARIOS.out.final_vcf_s6 ]
        post_final_all = post_final.collect { sc, ch -> ch.map { ps, vcf -> [ps, sc, vcf] } }
            .inject { a, b -> a.mix(b) }
        POST_FINAL_VCF(post_final_all.combine(POST_SCENARIOS.out.provenance, by: 0))
        POST_FINAL_VCF.out.status
            .collectFile(name: 'consensus_status.tsv', storeDir: "${params.results_dir}/REPORT", keepHeader: true, skip: 1, sort: true)
            .subscribe { f ->
                def none = f.readLines().drop(1).findAll { it.contains('NO_CONSENSUS_MUTATIONS') }
                if (none) {
                    log.warn "${none.size()} tumour/scenario row(s) have NO consensus mutations (run continued): REPORT/consensus_status.tsv"
                }
            }
        // <<< clinical_post

        // >>> vep
        // VEP + vcf2maf: the paired consensus VCF of --vep_scenario and the tumour-only PASS consensus,
        // one MAF per sample, then MAF/cohort.maf with an Analysis_Mode column.
        if (isToolSelected('vep')) {
            // Sample IDs as the callers name their VCF columns (Mutect2/SAGE): <patient>_<sample>_tumor|_normal
            vep_pair_ids = normal_post.cross(tumor_post){it[0]}.map { normal, tumor ->
                ["${tumor[1].patient}_${tumor[1].sample}".toString(),
                 "${tumor[1].patient}_${tumor[1].sample}_tumor".toString(),
                 "${normal[1].patient}_${normal[1].sample}_normal".toString()]
            }
            vep_paired = post_final[params.vep_scenario]
                .join(vep_pair_ids)
                .map { patient_sample, vcf, tumor_id, normal_id -> [patient_sample, 'paired', tumor_id, normal_id, vcf] }
            vep_tumor_only = POST_TUMOR_ONLY.out.consensus_pass
                .map { patient_sample, vcf, tbi -> [patient_sample, 'tumor_only', "${patient_sample}_tumor".toString(), 'NORMAL', vcf] }
            VEP_ANNOTATE(vep_paired.mix(vep_tumor_only))
            VCF2MAF(VEP_ANNOTATE.out.vcf)
            VEP_COHORT_MAF(VCF2MAF.out.maf.map { patient_sample, mode, maf -> maf }.collect())
        }
        // <<< vep
        // >>> pcgr
        // PCGR on every tumour's consensus: pairs from the chosen scenario's final VCF (+ the caller VCFs whose
        // FORMAT fields give tumour/normal depth and AF), tumour-only from POST_TUMOR_ONLY's PASS VCF.
        if (isToolSelected('pcgr')) {
            def pcgrNone = file("${projectDir}/Modules/PCGR/assets/NO_FILE")
            pcgr_paired = post_final["scenario_${params.pcgr_scenario}".toString()]
                .join(postInput.map { ps, m2, muse, st_snv, st_indel, sage, bam, bai -> [ps, m2, sage, st_snv, st_indel] })
                .map { ps, vcf, m2, sage, st_snv, st_indel -> [ps, 'paired', vcf, m2, sage, st_snv, st_indel] }
            pcgr_tonly = POST_TUMOR_ONLY.out.consensus_pass
                .map { ps, vcf, tbi -> [ps, 'tumor_only', vcf, pcgrNone, pcgrNone, pcgrNone, pcgrNone] }
            // HOOK (PureCN, --pcgr_cna_from_purecn): replace Channel.empty() with the PureCN allele-specific segments.
            pcgr_cna = Channel.empty()
            if (params.pcgr_cna_from_purecn) {
                log.warn "--pcgr_cna_from_purecn: PureCN is not wired on this branch yet; PCGR runs without CNA input."
            }
            pcgr_paired.mix(pcgr_tonly)
                .join(pcgr_cna, remainder: true)
                .filter { it.size() == 8 && it[1] in ['paired', 'tumor_only'] }
                .map { ps, mode, vcf, m2, sage, st_snv, st_indel, cna ->
                    def site = params.pcgr_tumor_sites.containsKey(ps) ? params.pcgr_tumor_sites[ps] : params.pcgr_tumor_site
                    [ps, mode, site, params.pcgr_assay_run, params.pcgr_target_size_mb_run, vcf, m2, sage, st_snv, st_indel, cna ?: pcgrNone]
                }
                .set { pcgr_input }
            PCGR(pcgr_input)
        }
        // <<< pcgr
        
        // AlleleCounter for VAF calculation
        def vaf_source = (params.first_step == "variant_calling") ? 
            RENAME_BAM_HEADER.out.renamed_bam.map { meta, bam, bai -> [meta.patient, meta, bam, bai] } :
            RECAL_good

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

        POST_SCENARIOS.out.final_vcf_s4
            .filter { ps, vcf -> vcf.readLines().any { !it.startsWith('#') } }  // clinical: skip AlleleCounter for an empty consensus (it segfaults); see REPORT/consensus_status.tsv
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


        POST_SCENARIOS.out.final_vcf_s5
            .filter { ps, vcf -> vcf.readLines().any { !it.startsWith('#') } }  // clinical: skip AlleleCounter for an empty consensus (it segfaults); see REPORT/consensus_status.tsv
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

        POST_SCENARIOS.out.final_vcf_s6
            .filter { ps, vcf -> vcf.readLines().any { !it.startsWith('#') } }  // clinical: skip AlleleCounter for an empty consensus (it segfaults); see REPORT/consensus_status.tsv
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
                '4': POST_SCENARIOS.out.final_vcf_s4,
                '5': POST_SCENARIOS.out.final_vcf_s5,
                '6': POST_SCENARIOS.out.final_vcf_s6
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
        // >>> gene_status
        // Per-gene somatic status for every tumour (paired and tumour-only) from every entry point, reusing the
        // channel sources of the TMB block above (recal_map_keyed / tumor_vaf pick the right BAMs per path).
        // Variants: the VEP-annotated consensus when --tool vep is on (SYMBOL + worst consequence; paired follows
        // --vep_scenario), else the plain consensus of --tmb_scenario (paired) / the tumour-only PASS consensus.
        // Copy number: CNVkit *_calls.cns and PureCN *_loh.csv when those tools ran (optional joins).
        if (params.tool && isToolSelected('gene_status')) {
            def gsScenario = params.tmb_scenario.toString()
            def gsVcfs
            if (isToolSelected('vep')) {
                gsVcfs = VEP_ANNOTATE.out.vcf.map { ps, mode, tid, nid, vcf ->
                    [ps.toString(), mode, (mode == 'paired' ? "${params.vep_scenario}_vep" : 'tumor_only_consensus_vep').toString(), vcf]
                }
            } else {
                def gsPairedVcf = [
                    '4': POST_SCENARIOS.out.final_vcf_s4,
                    '5': POST_SCENARIOS.out.final_vcf_s5,
                    '6': POST_SCENARIOS.out.final_vcf_s6
                ][gsScenario]
                gsVcfs = gsPairedVcf.map { ps, vcf -> [ps.toString(), 'paired', "scenario_${gsScenario}".toString(), vcf] }
                    .mix(POST_TUMOR_ONLY.out.consensus_pass.map { ps, vcf, tbi -> [ps.toString(), 'tumor_only', 'tumor_only_consensus', vcf] })
            }
            def gsBams = recal_map_keyed.map { ps, m -> [ps.toString(), m.tumor, m.tumor_bai] }
                .mix(tumor_vaf.filter { it[0] in params.tumor_only_patients }
                    .map { p, meta, bam, bai -> ["${p}_${meta.sample}".toString(), bam, bai] })
            def gsMain = gsVcfs.join(gsBams)    // [ps, mode, source, vcf, bam, bai]

            def gsCnvkit = isToolSelected('cnvkit') ?
                CNVkit_out.CNVkit_cn_files.flatten()
                    .filter { it.name.endsWith('_calls.cns') }
                    .map { f -> [(f.name - ~/_calls\.cns$/).toString(), f] } :
                Channel.empty()
            def gsPurecn = isToolSelected('purecn') ?
                PURECN_WF.out.loh.map { m, f -> ["${m.patient}_${m.tumor_meta.sample}".toString(), f] } :
                Channel.empty()

            def gsInput = gsMain.join(gsCnvkit, remainder: true)
                .filter { it[1] in ['paired', 'tumor_only'] }
                .map { row -> [row[0], row[1], row[2], row[3], row[4], row[5],
                               (row.size() > 6 && row[6]) ? row[6] : []] }
                .join(gsPurecn, remainder: true)
                .filter { it[1] in ['paired', 'tumor_only'] }
                .map { row -> [row[0], row[1], row[2], row[3], row[4], row[5], row[6],
                               (row.size() > 7 && row[7]) ? row[7] : []] }

            GENE_STATUS(gsInput, file(params.gene_status_genes), PREPARE_TARGETS.out.padded_bed.first())
            GENE_STATUS_MATRIX(GENE_STATUS.out.gene_status_tsv.map { ps, tsv -> tsv }.collect())
        }
        // <<< gene_status
    }

    // >>> qc
    // Sample identity (somalier) over every BAM of the run (tumours, normals, normal-only pool rows), from either
    // entry path, then one MultiQC report over all QC files.
    def humanGenome = params.genome in ['GRCh38', 'GRCh37']
    def qcBams = (params.first_step == 'variant_calling') ? sample_sheet : RECALIBRATE_out.pair_recal
    def mqc_files = Channel.empty()
    if (isToolSelected('somalier')) {
        def somalier_sites_ch = (params.somalier_sites == 'targets') ?
            SOMALIER_SITES(PREPARE_TARGETS.out.padded_bed,
                           [file(params.mutect2_germline), file("${params.mutect2_germline}.tbi")]).sites.first() :
            Channel.value(file(params.somalier_sites_file))
        SOMALIER_EXTRACT(qcBams, somalier_sites_ch)
        SOMALIER_RELATE(SOMALIER_EXTRACT.out.somalier.collect(), file(params.sample))
        SOMALIER_RELATE.out.flags.subscribe { f ->
            def n = f.readLines().findAll { it.trim() }.size() - 1
            def msg = n > 0 ?
                "Sample identity: ${n} flag(s) (swap / duplicate / sex / insufficient sites). See ${params.report_dir}/identity_flags.tsv" :
                "Sample identity: no flags (somalier, ${params.somalier_dir})"
            if (n > 0) log.warn msg else log.info msg
            writeToLog(params.log_file, msg)
        }
        mqc_files = mqc_files.mix(SOMALIER_RELATE.out.pairs, SOMALIER_RELATE.out.samples)
    }
    if (params.multiqc) {
        if (params.first_step == 'mapping') {
            mqc_files = mqc_files.mix(FASTQC.out.zip.flatten())
        }
        if (params.first_step in ['mapping', 'markdup']) {
            mqc_files = mqc_files.mix(MKDUP.out.metrics)
        }
        if (params.first_step in ['mapping', 'markdup', 'recalibration'] && humanGenome) {
            def bqsr_tables = params.targeted ?
                RECALIBRATE_BaseRecal_exome.out.BQSR_input.map { it[3] } :
                RECALIBRATE_MergeReport.out.BQSR_input.map { it[3] }
            mqc_files = mqc_files.mix(bqsr_tables)
        }
        if (humanGenome) {
            mqc_files = mqc_files.mix(MOSDEPTH.out.coverage.map { it[1] }.flatten())
            BCFTOOLS_STATS(POST_FINAL_VCF.out.vcf.map { it[2] }.mix(POST_TUMOR_ONLY.out.consensus_pass.map { it[1] }))
            mqc_files = mqc_files.mix(BCFTOOLS_STATS.out.stats)
        }
        MULTIQC(mqc_files.collect())
    }
    // <<< qc
}

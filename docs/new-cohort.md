# Running your own cohort

[&larr; back to the main README](../README.md)

The path from "a folder of clinical panel BAMs" to results, in six steps. Every helper lives in
[`scripts/cohort_prep/`](../scripts/cohort_prep). Steps 1-3 read your BAMs/BEDs, so run them yourself as the
data owner; their printed output is aggregate-only by design. Keep everything derived from the BAMs on your
protected filesystem.

## 0. What you need

- The BAMs (or FASTQs), with `.bai` next to each BAM.
- The SomaticBunny reference database for your genome build (see [installation.md](installation.md)).
- Ideally: the vendor's target BED and, for `--tool gene_status`, a genes BED
  (`bin/make_genes_bed.py` builds one from any public GENCODE/Ensembl GTF).
- Nextflow ≥ 25.04 and conda. Keep the pinned envs in `yml/` as they are ([tools.md](tools.md) explains why).

## 1. What are these BAMs?

```bash
# manifest.tsv: patient<TAB>sample<TAB>status<TAB>bam  (status: tumor|normal)
bash scripts/cohort_prep/check_bams_aggregate.sh manifest.tsv
```

This answers, without exposing any patient fields: the genome build (chr1 length), chr-prefix style, whether
the contigs match your `--ref`, and what processing the BAMs carry (dup-marking / BQSR / UMI consensus).
Decisions that follow: the entry point (`--first_step variant_calling` for processed BAMs, `mapping` for
FASTQs), and whether UMIs are already collapsed (vendor BAMs usually are; raw UMI FASTQs want `--umi`).

## 2. If the BAMs are b37-style (contigs `1, 2, MT`)

Make a chr-named copy of the BAMs' **own** FASTA (same sequence and contig order; GATK checks order):

```bash
SRC=/path/to/bams_reference.fasta OUT=/path/to/ref GATK=... bash scripts/cohort_prep/make_chr_reference.sh
```

Then run the pipeline with the printed `--ref/--ref_fai/--ref_dict` plus `--contig_rename auto` — the BAM
headers are renamed on the fly; the BAM files themselves are never modified. Any other resource you pass
(`--targets`, `--genes_bed`, a prebuilt PON) must use the same chr names; `bed_to_ref_names.sh` converts BEDs.

## 3. Targets

- **Vendor BED available:** convert its contig names to match `--ref` if needed:
  `bash scripts/cohort_prep/bed_to_ref_names.sh vendor.bed ref.fasta.fai targets.chr.bed`
- **No BED:** infer one from coverage (one assay at a time, 8-12 BAMs):
  `bash scripts/cohort_prep/infer_targets.sh bam_list.txt ref.fai outdir` — then sanity-check the printed
  footprint against the vendor's stated panel size.
- **Mixed assays in one cohort?** One run per assay, each with its own targets. Mixing captures in one run
  corrupts the per-assay resources (PON, MSI baseline, CNVkit reference).

## 4. The sample sheet

One CSV for the whole assay — paired, tumour-only and normal-only patients together
(see the main README for the format). Tumour-only patients just lack a normal row. Add every spare normal of
the assay as a normal-only row: they power the PON, the PureCN normal DB and the MSI baseline.

## 5. Pilot, then the full run

Run ONE pair + one tumour-only patient + all the normals first, with the resource-building flags on:

```bash
nextflow run main.nf -profile conda --sample pilot.csv --outdir PILOT \
    --type panel --targets targets.chr.bed --genome GRCh37 --first_step variant_calling \
    --database_path <db> --ref ... --ref_fai ... --ref_dict ... --contig_rename auto \
    --tool cnvkit,purecn,tmb,msi,vep,somalier,pcgr,gene_status --genes_bed genes.bed \
    --mutect2_cohort_pon build --somalier_sites targets
```

Review the pilot (`REPORT/`, MultiQC, [example-output.md](example-output.md) shows what to expect), then run
the full sheet reusing the pilot's resources:

```bash
    ... --mutect2_pon PILOT/Mutect2_PON/cohort.pon.vcf.gz \
        --purecn_normaldb PILOT/PureCN/normaldb/normalDB*.rds \
        --msi_baseline PILOT/MSI/baseline/msi_baseline.tsv --msi_sites PILOT/MSI/sites/msi_sites.tsv
```

(`--mutect2_cohort_pon build` and `--mutect2_pon` are mutually exclusive: build in the pilot, reuse in full.)

## 6. Cohorts with no normals at all

Everything still runs tumour-only: `--cnvkit_reference auto` (flat reference), the database PON, no Conpair.
PureCN needs a normal DB and MSI a baseline from the *same assay* — leave `purecn`/`msi` off unless you have
them from a previous cohort. Expect the documented tumour-only caveats (rare germline retained; see
[clinical-panels.md](clinical-panels.md)).

## Gotchas collected from real onboarding

- **Java:** GATK needs Java 17+; many cluster nodes default to Java 8. Export `JAVA_HOME` before dictionary
  creation and in scheduler wrappers.
- **Conda:** keep the pinned env files; don't "simplify" them ([tools.md](tools.md) has the war story). Put the
  conda cache on a filesystem with space - env builds fail on full home quotas.
- **MSI baseline needs ≥ 20 normals** from the same assay; below that the pipeline warns and the tumour-only
  cut-off should be calibrated on samples of known MSI status.
- **TMB on small panels** (< 1 Mb callable) is flagged unreliable - that's the warning working, not a bug.
- **Assay mismatches hide in capture tables:** verify any "WGS" classification against the actual covered
  footprint before believing it (a panel BAM classified as WGS usually means a mislabel).

# Example output folder

[&larr; back to the main README](../README.md)

The tree below is the real `--outdir` of the synthetic integration test (three patients: one pair, one
tumour-only, one zero-mutation pair, six pool normals; every tool on). The small tables in
[`example_output/`](example_output/) are actual files from that run - all simulated data, no patient
information.

```text
RESULTS/
├── REPORT/                      # run-level tables: read these first
│   ├── consensus_status.tsv     #   SNV/indel counts + NO_CONSENSUS_MUTATIONS per tumour and scenario
│   ├── identity_flags.tsv       #   somalier duplicates / swaps / sex flags
│   └── tumor_only_patients.txt  #   which patients ran without a normal
├── POST/                        # consensus calls
│   ├── <patient>_<sample>/scenario_{1..6}/   # paired: final.vcf.gz (+tbi, commands in header), snvs/, indels/
│   ├── tumor_only/              #   tumour-only consensus (PASS in Mutect2 ∩ SAGE)
│   └── VAF/                     #   AlleleCounter allele fractions at consensus sites
├── Mutect2/  STRELKA/  MuSE/  SAGE/          # per-caller VCFs
├── Mutect2_PON/                 # cohort panel of normals (--mutect2_cohort_pon build)
├── CNVkit/                      # log-ratio segments + calls (pooled or flat reference)
├── PureCN/
│   ├── <patient>_<sample>/      #   purity, ploidy, integer CN, LOH (csv per tumour)
│   ├── normaldb/                #   per-assay normal DB (reuse with --purecn_normaldb)
│   ├── coverage/  intervals/
├── TMB/                         # per-sample .tmb.tsv + tmb_cohort.tsv (+ callable regions, variant lists)
├── MSI/
│   ├── paired/  tumor_only/     #   <id>.msi_status.tsv per tumour
│   ├── baseline/                #   per-assay baseline (reuse with --msi_baseline)
│   └── sites/                   #   scanned microsatellite sites (reuse with --msi_sites)
├── VEP/  MAF/                   # annotated VCFs; per-sample MAFs + MAF/cohort.maf
├── somalier/                    # pairs/samples tables, identity_sex.tsv, interactive html
├── PCGR/<patient>_<sample>/     # <id>.pcgr.grch38.html + .xlsx + tiered TSVs per tumour
├── MultiQC/multiqc_report.html  # one QC report over the whole run
├── Conpair/  mosdepth/          # pair concordance; depth over targets
└── custom_intervals/            # the interval shards the run used
```

## Sample files from the synthetic run

| File | What it shows |
|---|---|
| [`consensus_status.tsv`](example_output/consensus_status.tsv) | the zero-mutation tumour reported as `NO_CONSENSUS_MUTATIONS` while the run continues |
| [`identity_flags.tsv`](example_output/identity_flags.tsv) | both planted duplicate samples caught at relatedness 1.000 |
| [`tmb_cohort.tsv`](example_output/tmb_cohort.tsv) | cohort TMB with callable Mb and the small-panel warning |
| [`PAIR1_A.tmb.tsv`](example_output/PAIR1_A.tmb.tsv) | one tumour's TMB record |
| [`PAIR1_A.msi_status.tsv`](example_output/PAIR1_A.msi_status.tsv) | paired MSI call (MSS) |
| [`TONLY1_A.msi_status.tsv`](example_output/TONLY1_A.msi_status.tsv) | tumour-only MSI call against the baseline |

On real data, most of these files are patient-level clinical results: keep `--outdir` on a protected
filesystem and never commit or share them. This example is publishable only because every read is simulated.

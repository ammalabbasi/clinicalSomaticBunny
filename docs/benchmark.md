# Validation and benchmarks

[&larr; back to the main README](../README.md)

All validation uses simulated reads on public GRCh38 with planted truth; no patient data.

## Merged-branch integration test (every tool on at once)

542 tasks, 0 failures. Paired consensus recovered the planted variants; a zero-mutation pair completed with
`NO_CONSENSUS_MUTATIONS` reported; every final VCF parses and carries its generating commands in the header;
input BAMs were byte-identical afterwards. somalier flagged both planted duplicate samples (relatedness 1.000)
with no false flags. PCGR reports held 20/20 planted somatic variants (paired; 0/1472 planted germline) and
18/20 (tumour-only; 8 rare germline retained, the documented tumour-only caveat).

## Per-tool synthetic validation

| Tool | Test | Result |
|---|---|---|
| Cohort PON | 5 recurrent artefacts planted in every normal | 5/5 removed, 17/18 true variants kept |
| PureCN | tumour at purity 0.6: gain, amp, del-LOH, CN-LOH, neutral, loss | purity 0.600 in all 3 modes, CN 6/6, all LOH states |
| MSI | MSI-high (~40% loci shifted) vs stable, 150 MS targets | 28.0% vs 0.9%, separated paired and tumour-only |
| TMB | expected value hand-computed from truth | exact match; small-panel warning fires |
| UMI (fgbio) | PCR families + single-read errors, 5 variants at 2-3% VAF | false positives 6 to 0; low-VAF recovery 0/5 to 4/5 |
| VEP + MAF | TP53-region consensus VCFs | MAF rows = PASS variants; TP53 annotated |
| somalier | planted swap + duplicate | both flagged, 0 false, sex 5/5 |
| Target inference | 40 planted targets, coverage only | 40/40 recovered, 41.4/42 kb |

## Comparison with nf-core/sarek

Same synthetic reads, references, PON and targets through sarek 3.x: paired consensus calls were identical;
in tumour-only mode this pipeline leaked 0/15 planted germline SNPs vs sarek's 8/15.

## Benchmark

Code to reproduce the benchmark against nf-core/sarek is in [`benchmark/`](benchmark).

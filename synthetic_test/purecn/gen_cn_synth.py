#!/usr/bin/env python3
"""Synthetic copy-number truth set for testing PureCN (and CNVkit) in SomaticBunny panel mode.

Everything comes from the PUBLIC GRCh38 reference and the PUBLIC gnomAD af-only VCF; no sample data is used.

Design
  * Six 2 Mb regions on chr1, chr3, chr7, chr8, chr12, chr17, each with 200 targets of 200 bp (one per 10 kb slot).
    A target is centred on a common gnomAD SNP (AF 0.2-0.8) when its slot has one, so most targets carry >= 1 SNP.
  * Germline: the tumour's individual (I0) is heterozygous at every such SNP inside a target, with random phase.
    The five pool-only normals (I1..I5) get Hardy-Weinberg genotypes from the gnomAD AF.
  * Capture efficiency per target: log-normal (sdlog 0.3), identical in every sample.
  * Tumour (purity 0.6), allele-specific copy number (hapA, hapB) per region:
      chr1 (2,1) gain CN3 | chr3 (1,0) loss CN1, LOH | chr7 (4,1) amp CN5 | chr8 (2,0) copy-neutral LOH
      chr12 (1,1) CN2     | chr17 (1,0) loss CN1, LOH
  * ~20 somatic SNVs on hapA; expected VAF = p*m / (2(1-p) + p*CN), m = mutated copies.
  * Reads: 150 bp paired-end, fragment ~N(300, 30), substitution errors 0.1 %. Per-target fragment counts are
    proportional to efficiency x copies (tumour: (1-p)*2 + p*CN), ~300x tumour / ~150x normals.

Usage: gen_cn_synth.py <ref.fa> <gnomad_af_only.vcf.gz> <bcftools> <samtools> <out_dir> [seed]
Writes <out_dir>/truth/* and <out_dir>/fastq/<sample>_{1,2}.fq
"""
import bisect
import math
import os
import random
import subprocess
import sys

REGIONS = [  # name, chrom, start (1-based), nA, nB
    ("chr1", "chr1", 30_000_001, 2, 1),
    ("chr3", "chr3", 50_000_001, 1, 0),
    ("chr7", "chr7", 100_000_001, 4, 1),
    ("chr8", "chr8", 120_000_001, 2, 0),
    ("chr12", "chr12", 50_000_001, 1, 1),
    ("chr17", "chr17", 40_000_001, 1, 0),
]
REGION_LEN = 2_000_000
N_TARGETS = 200
TARGET_LEN = 200
SLOT = REGION_LEN // N_TARGETS
PURITY = 0.6
EFF_SDLOG = 0.3
DEPTH_TUMOR = 300 * float(os.environ.get("PCN_SCALE", "1"))
DEPTH_NORMAL = 150 * float(os.environ.get("PCN_SCALE", "1"))
READ_LEN = 150
FRAG_MEAN, FRAG_SD = 300, 30
ERROR = 0.001
N_NORMALS = 6
SOMATIC_PER_REGION = [4, 3, 4, 3, 3, 3]   # 20 in total
COMP = str.maketrans("ACGTN", "TGCAN")


def sh(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def read_region(samtools, ref, chrom, start):
    out = sh([samtools, "faidx", ref, f"{chrom}:{start}-{start + REGION_LEN - 1}"])
    return "".join(out.split("\n")[1:]).upper()


def common_snps(bcftools, gnomad, chrom, start, seq):
    out = sh([bcftools, "query", "-r", f"{chrom}:{start}-{start + REGION_LEN - 1}", "-i",
              'TYPE="snp" && N_ALT=1 && INFO/AF>=0.2 && INFO/AF<=0.8', "-f", "%POS\t%REF\t%ALT\t%INFO/AF\n", gnomad])
    snps = []
    for line in out.splitlines():
        pos, ref, alt, af = line.split("\t")
        off = int(pos) - start                     # 0-based offset in region
        if 0 <= off < REGION_LEN and seq[off] == ref and len(alt) == 1 and alt in "ACGT":
            snps.append((off, ref, alt, float(af)))
    return snps


def poisson(rng, lam):
    if lam < 30:
        l, k, p = math.exp(-lam), 0, 1.0
        while True:
            p *= rng.random()
            if p <= l:
                return k
            k += 1
    return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))


def main(ref, gnomad, bcftools, samtools, out_dir, seed=7):
    rng = random.Random(seed)
    tdir, fdir = os.path.join(out_dir, "truth"), os.path.join(out_dir, "fastq")
    os.makedirs(tdir, exist_ok=True)
    os.makedirs(fdir, exist_ok=True)

    regions = []
    for (name, chrom, start, na, nb), n_som in zip(REGIONS, SOMATIC_PER_REGION):
        seq = read_region(samtools, ref, chrom, start)
        snps = common_snps(bcftools, gnomad, chrom, start, seq)
        snp_offs = [s[0] for s in snps]
        targets = []
        for k in range(N_TARGETS):
            lo, hi = k * SLOT + 500, (k + 1) * SLOT - 500 - TARGET_LEN
            i0, i1 = bisect.bisect_left(snp_offs, lo + 60), bisect.bisect_right(snp_offs, hi + TARGET_LEN - 60)
            cands = [o for o in snp_offs[i0:i1]]
            t0 = None
            for _ in range(20):
                t = (rng.choice(cands) - TARGET_LEN // 2) if cands else rng.randrange(lo, hi)
                t = min(max(t, lo), hi)
                if "N" not in seq[t - 600:t + TARGET_LEN + 600]:
                    t0 = t
                    break
            if t0 is not None:
                targets.append((t0, t0 + TARGET_LEN))
        eff = [rng.lognormvariate(0, EFF_SDLOG) for _ in targets]
        # germline SNPs inside targets (5 bp margin): I0 het with random phase; I1..I5 Hardy-Weinberg
        germ = []
        for t0, t1 in targets:
            for s in snps[bisect.bisect_left(snp_offs, t0 + 5):bisect.bisect_left(snp_offs, t1 - 5)]:
                off, r, a, af = s
                gts = [(1, 0) if rng.random() < 0.5 else (0, 1)]
                for _ in range(N_NORMALS - 1):
                    gts.append((int(rng.random() < af), int(rng.random() < af)))
                germ.append((off, r, a, af, gts))
        germ_offs = {g[0] for g in germ}
        # somatic SNVs on hapA, away from SNPs and from each other
        som = []
        while len(som) < n_som:
            t0, t1 = rng.choice(targets)
            off = rng.randrange(t0 + 20, t1 - 20)
            if any(abs(off - g) < 10 for g in germ_offs) or any(abs(off - s[0]) < 50 for s in som):
                continue
            r = seq[off]
            if r not in "ACGT":
                continue
            a = rng.choice([b for b in "ACGT" if b != r])
            m = 1
            if name == "chr7" and len(som) == 0:
                m = 4          # early event on every amplified copy
            if name == "chr8" and len(som) == 0:
                m = 2          # before the copy-neutral LOH duplication
            som.append((off, r, a, m))
        cn = na + nb
        regions.append(dict(name=name, chrom=chrom, start=start, seq=seq, na=na, nb=nb, cn=cn,
                            targets=targets, eff=eff, germ=germ, som=som))

    # ---- truth files
    with open(os.path.join(tdir, "targets.bed"), "w") as fh:
        for R in regions:
            for t0, t1 in R["targets"]:
                fh.write(f"{R['chrom']}\t{R['start'] - 1 + t0}\t{R['start'] - 1 + t1}\n")
    with open(os.path.join(tdir, "regions.tsv"), "w") as fh:
        fh.write("region\tchrom\tstart\tend\tnA\tnB\tCN\tLOH\tn_targets\tn_het_snps_I0\n")
        for R in regions:
            loh = int(R["nb"] == 0)
            fh.write(f"{R['name']}\t{R['chrom']}\t{R['start']}\t{R['start'] + REGION_LEN - 1}\t{R['na']}\t{R['nb']}\t"
                     f"{R['cn']}\t{loh}\t{len(R['targets'])}\t{len(R['germ'])}\n")
    with open(os.path.join(tdir, "efficiency.tsv"), "w") as fh:
        fh.write("chrom\tstart0\tend\tefficiency\n")
        for R in regions:
            for (t0, t1), e in zip(R["targets"], R["eff"]):
                fh.write(f"{R['chrom']}\t{R['start'] - 1 + t0}\t{R['start'] - 1 + t1}\t{e:.4f}\n")
    with open(os.path.join(tdir, "germline.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\tgnomad_af\t" + "\t".join(f"I{i}" for i in range(N_NORMALS)) + "\n")
        for R in regions:
            for off, r, a, af, gts in R["germ"]:
                fh.write(f"{R['chrom']}\t{R['start'] + off}\t{r}\t{a}\t{af}\t" +
                         "\t".join(f"{x}|{y}" for x, y in gts) + "\n")
    with open(os.path.join(tdir, "somatic.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\tregion\tCN\tmultiplicity\texpected_vaf\n")
        for R in regions:
            for off, r, a, m in R["som"]:
                vaf = PURITY * m / (2 * (1 - PURITY) + PURITY * R["cn"])
                fh.write(f"{R['chrom']}\t{R['start'] + off}\t{r}\t{a}\t{R['name']}\t{R['cn']}\t{m}\t{vaf:.4f}\n")

    # ---- reads
    samples = [("PCN01_T_tumor", 0, True)] + [(f"PCN0{i + 1}_N_normal", i, False) for i in range(N_NORMALS)]
    with open(os.path.join(tdir, "samples.tsv"), "w") as fh:
        fh.write("fastq_prefix\tindividual\tis_tumor\n")
        for s, ind, is_t in samples:
            fh.write(f"{s}\tI{ind}\t{int(is_t)}\n")
    for sname, ind, is_tumor in samples:
        r1 = open(os.path.join(fdir, f"{sname}_1.fq"), "w")
        r2 = open(os.path.join(fdir, f"{sname}_2.fq"), "w")
        nread = 0
        for R in regions:
            seq = R["seq"]
            # per-haplotype variant tables: offset -> base
            hap = [{}, {}]
            for off, r, a, af, gts in R["germ"]:
                for h in (0, 1):
                    if gts[ind][h]:
                        hap[h][off] = a
            som = {off: (a, m) for off, r, a, m in R["som"]} if is_tumor else {}
            hap_keys = [sorted(hap[0]), sorted(hap[1])]
            som_keys = sorted(som)
            if is_tumor:
                comps = [((1 - PURITY), 0, False), ((1 - PURITY), 1, False),
                         (PURITY * R["na"], 0, True), (PURITY * R["nb"], 1, True)]
                rel = ((1 - PURITY) * 2 + PURITY * R["cn"]) / 2
                depth = DEPTH_TUMOR
            else:
                comps = [(1.0, 0, False), (1.0, 1, False)]
                rel, depth = 1.0, DEPTH_NORMAL
            wsum = sum(c[0] for c in comps)
            for (t0, t1), e in zip(R["targets"], R["eff"]):
                lam = 1.47 * depth * e * rel * (t1 - t0 + FRAG_MEAN - 60) / (TARGET_LEN + FRAG_MEAN - 60)
                for _ in range(poisson(rng, lam)):
                    flen = int(min(max(rng.gauss(FRAG_MEAN, FRAG_SD), 200), 450))
                    s = rng.randrange(t0 - flen + 30, t1 - 30)
                    if s < 0 or s + flen > len(seq):
                        continue
                    x = rng.random() * wsum
                    for w, h, tum in comps:
                        x -= w
                        if x <= 0:
                            break
                    frag = list(seq[s:s + flen])
                    for off in _offs_in(hap_keys[h], s, s + flen):
                        frag[off - s] = hap[h][off]
                    if tum and h == 0:
                        for off in _offs_in(som_keys, s, s + flen):
                            a, m = som[off]
                            if rng.random() < m / R["na"]:
                                frag[off - s] = a
                    frag = "".join(frag)
                    a_read = frag[:READ_LEN]
                    b_read = frag[-READ_LEN:].translate(COMP)[::-1]
                    if rng.random() < 0.5:
                        a_read, b_read = b_read, a_read
                    nread += 1
                    name = f"{sname}:{R['name']}:{nread}"
                    for fh, rd in ((r1, a_read), (r2, b_read)):
                        rd = _errors(rng, rd)
                        fh.write(f"@{name}\n{rd}\n+\n{'F' * len(rd)}\n")
        r1.close()
        r2.close()
        print(f"{sname}: {nread} read pairs", flush=True)


def _offs_in(ks, lo, hi):
    return ks[bisect.bisect_left(ks, lo):bisect.bisect_left(ks, hi)]


def _errors(rng, rd):
    k = 0
    lam = ERROR * len(rd)
    l, p = math.exp(-lam), rng.random()
    while p > l:
        k += 1
        p *= rng.random()
    if not k:
        return rd
    rd = list(rd)
    for _ in range(k):
        i = rng.randrange(len(rd))
        rd[i] = rng.choice([b for b in "ACGT" if b != rd[i]])
    return "".join(rd)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2], a[3], a[4], int(a[5]) if len(a) > 5 else 7)

#!/usr/bin/env python3
"""Infer a panel target BED from coverage: keep regions well covered in most samples of ONE assay.

usage: consensus_targets.py --out OUT_DIR --fai REF.fai --min-frac 0.5 --min-len 50 --gap 100 \
           sample1.quantized.bed.gz sample2.quantized.bed.gz ...

Inputs are mosdepth --quantize outputs (see infer_targets.sbatch), where label HIGH marks bases at or above the
depth cutoff. Only primary chromosomes (1-22, X, Y, with or without chr) are kept.

Writes OUT_DIR/inferred_targets.bed and OUT_DIR/inferred_targets.summary.tsv.
Prints and writes ONLY aggregate numbers (bp, interval counts, per-sample overlap by index), never intervals,
so the summary can be shared; the BED itself is derived from the BAMs and stays on the cluster.
"""
import argparse
import gzip
import math
import os
import re

PRIMARY = re.compile(r"^(chr)?([0-9]+|X|Y)$")


def read_high(path):
    """Merged HIGH intervals per contig from a mosdepth quantized BED."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, label = line.rstrip("\n").split("\t")[:4]
            if label != "HIGH" or not PRIMARY.match(c):
                continue
            iv = out.setdefault(c, [])
            s, e = int(s), int(e)
            if iv and s <= iv[-1][1]:
                iv[-1][1] = max(iv[-1][1], e)
            else:
                iv.append([s, e])
    return out


def total(ivs):
    return sum(e - s for v in ivs.values() for s, e in v)


def overlap_bp(a, b):
    bp = 0
    for c in a:
        if c not in b:
            continue
        x, y = a[c], b[c]
        i = j = 0
        while i < len(x) and j < len(y):
            lo, hi = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
            if lo < hi:
                bp += hi - lo
            if x[i][1] < y[j][1]:
                i += 1
            else:
                j += 1
    return bp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--min-frac", type=float, default=0.5, help="fraction of samples a base must be HIGH in")
    ap.add_argument("--min-len", type=int, default=50, help="drop consensus intervals shorter than this")
    ap.add_argument("--gap", type=int, default=100, help="merge consensus intervals closer than this")
    ap.add_argument("beds", nargs="+")
    a = ap.parse_args()

    order = [l.split("\t")[0] for l in open(a.fai) if PRIMARY.match(l.split("\t")[0])]
    samples = [read_high(p) for p in a.beds]
    n = len(samples)
    need = max(1, math.ceil(a.min_frac * n))

    consensus = {}
    for c in order:
        events = []
        for smp in samples:
            for s, e in smp.get(c, []):
                events.append((s, 1))
                events.append((e, -1))
        if not events:
            continue
        events.sort(key=lambda t: (t[0], t[1]))
        depth, start, iv = 0, None, []
        for pos, d in events:
            before = depth
            depth += d
            if before < need <= depth:
                start = pos
            elif before >= need > depth and start is not None:
                if pos > start:
                    if iv and start - iv[-1][1] <= a.gap:
                        iv[-1][1] = pos
                    else:
                        iv.append([start, pos])
                start = None
        iv = [x for x in iv if x[1] - x[0] >= a.min_len]
        if iv:
            consensus[c] = iv

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "inferred_targets.bed"), "w") as fh:
        for c in order:
            for s, e in consensus.get(c, []):
                fh.write(f"{c}\t{s}\t{e}\n")

    cons_bp = total(consensus)
    lines = [
        ("samples", n),
        ("min_samples_per_base", need),
        ("consensus_intervals", sum(len(v) for v in consensus.values())),
        ("consensus_bp", cons_bp),
        ("consensus_Mb", round(cons_bp / 1e6, 3)),
        ("contigs_with_targets", len(consensus)),
    ]
    rows = []
    for i, smp in enumerate(samples, 1):
        bp = total(smp)
        ov = overlap_bp(smp, consensus)
        rows.append((i, bp, ov, round(ov / bp, 3) if bp else 0, round(ov / cons_bp, 3) if cons_bp else 0))
    with open(os.path.join(a.out, "inferred_targets.summary.tsv"), "w") as fh:
        for k, v in lines:
            fh.write(f"{k}\t{v}\n")
        fh.write("\nsample_index\thigh_bp\tin_consensus_bp\tfrac_of_sample_in_consensus\tfrac_of_consensus_covered\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(a.out, "inferred_targets.summary.tsv")) as fh:
        print(fh.read())


if __name__ == "__main__":
    main()

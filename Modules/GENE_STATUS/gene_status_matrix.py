#!/usr/bin/env python3
"""Cohort gene-status matrix from per-sample *.gene_status.tsv files (Modules/GENE_STATUS/gene_status.py).

Rows = genes ordered by chrom/start (coordinates from the per-sample TSVs), columns = samples (sorted).
A gene without a row for a sample (e.g. mixed assays) holds '.'. Stdlib only.
"""
import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("tsv", nargs="+")
    a = ap.parse_args()

    status = {}      # gene -> {sample: status}
    coords = {}      # gene -> (chrom, start)
    samples = set()
    for path in a.tsv:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                need = ("sample", "gene", "chrom", "start", "status")
                if any(k not in r or r[k] is None for k in need):
                    sys.exit(f"{path}: not a gene_status TSV (needs columns {', '.join(need)})")
                samples.add(r["sample"])
                status.setdefault(r["gene"], {})[r["sample"]] = r["status"]
                c = (r["chrom"], int(r["start"]))
                coords[r["gene"]] = min(coords.get(r["gene"], c), c)

    sample_order = sorted(samples)
    gene_order = sorted(status, key=lambda g: (coords[g][0], coords[g][1], g))
    with open(a.out, "w") as fh:
        fh.write("\t".join(["gene"] + sample_order) + "\n")
        for g in gene_order:
            fh.write("\t".join([g] + [status[g].get(s, ".") for s in sample_order]) + "\n")


if __name__ == "__main__":
    main()

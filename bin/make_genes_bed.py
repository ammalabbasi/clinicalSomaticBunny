#!/usr/bin/env python3
"""Build a 4-column genes BED (chrom, start, end, gene name) for --genes_bed from a public Ensembl or
GENCODE GTF (plain or .gz), e.g.:

    python3 bin/make_genes_bed.py --gtf gencode.v46.annotation.gtf.gz --out genes.bed
    python3 bin/make_genes_bed.py --gtf Homo_sapiens.GRCh38.112.gtf.gz --genes panel_genes.txt \\
        --add-chr --out genes.bed

One line per `gene` feature (GTF is 1-based inclusive; the BED is 0-based half-open). Options:
  --genes FILE        keep only these gene names (one per line, case-insensitive); unmatched names are
                      reported on stderr and the script exits non-zero so a panel typo cannot pass silently
  --protein-coding    keep only gene_type/gene_biotype == protein_coding
  --attr NAME         attribute holding the gene name (default gene_name, falling back to gene_id)
  --add-chr           prefix contigs with 'chr' (Ensembl GTFs name contigs 1..22,X; GENCODE uses chr1..)
The BED must use the same contig names as the pipeline's --ref (main.nf checks and stops otherwise).
"""
import argparse
import gzip
import re
import sys


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def attr(attrs, name):
    m = re.search(name + r'\s+"([^"]+)"', attrs) or re.search(name + r"\s+([^;]+);", attrs)
    return m.group(1).strip() if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", help="file with one gene name per line to keep")
    ap.add_argument("--attr", default="gene_name")
    ap.add_argument("--protein-coding", action="store_true")
    ap.add_argument("--add-chr", action="store_true")
    a = ap.parse_args()

    wanted = None
    if a.genes:
        with open(a.genes) as fh:
            wanted = {w.strip().upper() for w in fh if w.strip() and not w.startswith("#")}
        if not wanted:
            sys.exit(f"--genes {a.genes}: empty")

    rows, found = [], set()
    with open_text(a.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            if a.protein_coding:
                bt = attr(f[8], "gene_type") or attr(f[8], "gene_biotype")
                if bt != "protein_coding":
                    continue
            name = attr(f[8], a.attr) or attr(f[8], "gene_id")
            if not name:
                continue
            if wanted is not None:
                if name.upper() not in wanted:
                    continue
                found.add(name.upper())
            chrom = ("chr" + f[0]) if (a.add_chr and not f[0].startswith("chr")) else f[0]
            rows.append((chrom, int(f[3]) - 1, int(f[4]), name))

    if not rows:
        sys.exit(f"{a.gtf}: no matching 'gene' features (is this an Ensembl/GENCODE GTF?)")
    rows.sort()
    with open(a.out, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"{a.out}: {len(rows)} gene interval(s), {len({r[3] for r in rows})} gene name(s)", file=sys.stderr)

    if wanted is not None:
        missing = sorted(wanted - found)
        if missing:
            print(f"ERROR: {len(missing)} requested gene(s) not in the GTF: {', '.join(missing[:20])}"
                  + (" ..." if len(missing) > 20 else ""), file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()

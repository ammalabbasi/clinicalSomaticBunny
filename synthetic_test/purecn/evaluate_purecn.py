#!/usr/bin/env python3
"""Score PureCN (and CNVkit) results on the synthetic copy-number truth of gen_cn_synth.py.

  evaluate_purecn.py --truth <dir> --out <tsv> --sample <label>=<PureCN result prefix> ... --cnvkit <label>=<calls.cns> ...

<PureCN result prefix> is <results>/PureCN/<patient>_<sample> (the dir); files <dir>/<basename>.csv, _dnacopy.seg,
_loh.csv and _variants.csv are read. Checks:
  purity  within 0.10 of the truth (0.6)
  ploidy  within 0.35 of the truth (bp-weighted mean tumour CN of the six regions = 14/6 = 2.33)
  CN      per region: bp-weighted majority integer C of the PureCN segments == truth CN
  LOH     per region: >= 50 % of the region in segments with minor copy number M == 0 (and C > 0), vs truth
  somatic truth somatic SNVs with POSTERIOR.SOMATIC > 0.5 / truth het germline SNPs with POSTERIOR.SOMATIC < 0.5
  CNVkit  per region: bp-weighted mean segment log2 -> gain (> +0.15), loss (< -0.15) or neutral, vs the truth
          direction of the tumour/normal coverage ratio (CN > 2 gain, CN < 2 loss, CN 2 neutral)
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

PURITY = 0.6


def chrom_key(c):
    c = str(c)
    return c[3:] if c.startswith("chr") else c


def read_regions(truth):
    with open(os.path.join(truth, "regions.tsv")) as fh:
        return [dict(r, start=int(r["start"]), end=int(r["end"]), CN=int(r["CN"]), LOH=int(r["LOH"]),
                     nA=int(r["nA"]), nB=int(r["nB"])) for r in csv.DictReader(fh, delimiter="\t")]


def overlap(r, chrom, s, e):
    if chrom_key(chrom) != chrom_key(r["chrom"]):
        return 0
    return max(0, min(e, r["end"]) - max(s, r["start"]) + 1)


def read_table(path, delim):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def num(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def score_purecn(label, prefix, regions, truth_dir, rows):
    base = prefix if os.path.exists(prefix + ".csv") else os.path.join(prefix, os.path.basename(prefix.rstrip("/")))
    if not os.path.exists(base + ".csv"):
        rows.append([label, "run", "PureCN output", "missing", base + ".csv", "FAIL"])
        return
    cur = read_table(base + ".csv", ",")[0]
    purity, ploidy = num(cur.get("Purity")), num(cur.get("Ploidy"))
    true_ploidy = sum(r["CN"] for r in regions) / len(regions)
    rows.append([label, "purity", "all", f"{PURITY:.2f}", f"{purity:.3f}", "PASS" if abs(purity - PURITY) <= 0.10 else "FAIL"])
    rows.append([label, "ploidy", "all", f"{true_ploidy:.2f}", f"{ploidy:.3f}",
                 "PASS" if abs(ploidy - true_ploidy) <= 0.35 else "FAIL"])
    rows.append([label, "flagged", "all", "", f"{cur.get('Flagged')} {cur.get('Comment', '')}".strip(), "INFO"])

    seg = read_table(base + "_dnacopy.seg", "\t")
    for r in regions:
        w = defaultdict(int)
        for s in seg:
            ov = overlap(r, s["chrom"], int(float(s["loc.start"])), int(float(s["loc.end"])))
            if ov and s.get("C") not in (None, "", "NA"):
                w[int(round(float(s["C"])))] += ov
        call = max(w, key=w.get) if w else None
        detail = ",".join(f"C{k}:{v / sum(w.values()):.2f}" for k, v in sorted(w.items())) if w else "no segment"
        rows.append([label, "CN", r["region"], str(r["CN"]), f"{call} ({detail})", "PASS" if call == r["CN"] else "FAIL"])

    loh_path = base + "_loh.csv"
    if os.path.exists(loh_path):
        loh = read_table(loh_path, ",")
        for r in regions:
            tot = lohbp = 0
            mstates = defaultdict(int)
            for s in loh:
                ov = overlap(r, s["chr"], int(float(s["start"])), int(float(s["end"])))
                if not ov:
                    continue
                c, m = num(s.get("C")), num(s.get("M"))
                tot += ov
                mstates[f"{s.get('C')}/{s.get('M')}"] += ov
                if m is not None and c is not None and m == 0 and c > 0:
                    lohbp += ov
            call = int(tot > 0 and lohbp / tot >= 0.5)
            detail = ",".join(f"{k}:{v / tot:.2f}" for k, v in sorted(mstates.items())) if tot else "no segment"
            rows.append([label, "LOH", r["region"], str(r["LOH"]), f"{call} (C/M {detail})",
                         "PASS" if tot and call == r["LOH"] else "FAIL"])
            rows.append([label, "allele_CN", r["region"], f"{r['nA']}/{r['nB']}",
                         max(mstates, key=mstates.get) if mstates else "none",
                         "PASS" if mstates and max(mstates, key=mstates.get) == f"{r['CN']}/{min(r['nA'], r['nB'])}" else "INFO"])
    else:
        rows.append([label, "LOH", "all", "", "no _loh.csv", "FAIL"])

    var_path = base + "_variants.csv"
    if os.path.exists(var_path):
        var = {(chrom_key(v["chr"]), int(float(v["start"]))): num(v.get("POSTERIOR.SOMATIC"))
               for v in read_table(var_path, ",")}
        som = [(chrom_key(t["chrom"]), int(t["pos"])) for t in read_table(os.path.join(truth_dir, "somatic.tsv"), "\t")]
        germ = [(chrom_key(t["chrom"]), int(t["pos"])) for t in read_table(os.path.join(truth_dir, "germline.tsv"), "\t")]
        s_in = [var[k] for k in som if k in var and var[k] is not None]
        g_in = [var[k] for k in germ if k in var and var[k] is not None]
        s_ok = sum(p > 0.5 for p in s_in)
        g_ok = sum(p < 0.5 for p in g_in)
        rows.append([label, "somatic_posterior", "truth somatic", str(len(som)),
                     f"{s_ok}/{len(s_in)} scored > 0.5 ({len(som) - len(s_in)} not in PureCN output)",
                     "PASS" if s_in and s_ok >= 0.8 * len(s_in) else "FAIL"])
        rows.append([label, "germline_posterior", "truth het germline", str(len(germ)),
                     f"{g_ok}/{len(g_in)} scored < 0.5 ({len(germ) - len(g_in)} not in PureCN output)",
                     "PASS" if g_in and g_ok >= 0.95 * len(g_in) else "FAIL"])


def score_cnvkit(label, path, regions, rows):
    if not os.path.exists(path):
        rows.append([label, "cnvkit", "all", "", "missing " + path, "FAIL"])
        return
    seg = read_table(path, "\t")
    for r in regions:
        tot = acc = 0.0
        for s in seg:
            ov = overlap(r, s["chromosome"], int(s["start"]) + 1, int(s["end"]))
            if ov:
                tot += ov
                acc += ov * float(s["log2"])
        if not tot:
            rows.append([label, "CNVkit_log2", r["region"], "", "no segment", "FAIL"])
            continue
        l2 = acc / tot
        got = "gain" if l2 > 0.15 else ("loss" if l2 < -0.15 else "neutral")
        want = "gain" if r["CN"] > 2 else ("loss" if r["CN"] < 2 else "neutral")
        rows.append([label, "CNVkit_log2", r["region"], want, f"{got} ({l2:+.3f})", "PASS" if got == want else "FAIL"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", action="append", default=[])
    ap.add_argument("--cnvkit", action="append", default=[])
    a = ap.parse_args()
    regions = read_regions(a.truth)
    rows = []
    for s in a.sample:
        label, prefix = s.split("=", 1)
        score_purecn(label, prefix, regions, a.truth, rows)
    for s in a.cnvkit:
        label, path = s.split("=", 1)
        score_cnvkit(label, path, regions, rows)
    with open(a.out, "w") as fh:
        fh.write("run\tcheck\tregion\texpected\tobserved\tresult\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    w = [max(len(str(r[i])) for r in rows) for i in range(6)] if rows else []
    for r in rows:
        print("  ".join(str(x).ljust(w[i]) for i, x in enumerate(r)))
    n_fail = sum(r[5] == "FAIL" for r in rows)
    print(f"\n{sum(r[5] == 'PASS' for r in rows)} PASS, {n_fail} FAIL -> {a.out}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
